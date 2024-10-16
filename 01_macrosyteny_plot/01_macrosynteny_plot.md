# Macro-synteny analysis
Author: Ying Chen</br>

### Step 1. Download genome assembly from NCBI and run BUSCO

```bash
#!/bin/bash
LINEAGE_tetr=/project/lougheedlab/chorus/02_HiFi_HiC/04_contig_QC/BUSCO/busco_downloads/lineages/tetrapoda_odb10
for i in /project/lougheedlab/chorus/frog_genome/seq/add_more/*.fasta;
do
    OUT=results/$(basename ${i%.fasta})
    conda run -n env_busco --live-stream busco -c 88 -i $i -o $OUT -m genome -l $LINEAGE_tetr
done

# plot BUSCO results
mkdir BUSCO_summaries
cp results/*/short_summary.*.txt ./BUSCO_summaries/
conda activate env_busco
generate_plot.py -wd BUSCO_summaries
conda deactivate
```


### Step 2. Make synteny plots

[jcvi](https://github.com/tanghaibao/jcvi) macro-syteny plot input files:</br>
- `.bed` gene position information for both complete and fragmented BUSCOs for each species. Six columns: Chromosome, Gene Start, Gene End, Busco id, Busco Score, Strand</br>
- `.simple` synteny blocks (using a single gene will give you warning of coordinates not found). Six columns: start and stop gene of genome 1, start and stop gene of genome 2, score, orientation. tab-delineated</br>
- `seqids` chromosome/scaffold name<br>
- `layout` plot layout parameters

**Step 2.1**: create `.bed` file
```bash
# create a full list of all frog species
printf '%s\n' /project/lougheedlab/chorus/frog_genome/seq/*.fasta | sed 's!.*/!!' | sed 's/.fasta//g' > frog.list

# convert BUSCO tsv to bed file format
while read p; do
  IN=/project/lougheedlab/chorus/frog_genome/01_BUSCO/results/$p/run_tetrapoda_odb10/full_table.tsv
  OUT=$p.bed
  paste <(cut -f3-5 $IN) <(cut -f1-2 $IN) <(cut -f7 $IN)  <(cut -f6 $IN) | grep -e 'Complete' -e 'Fragmented' | cut -f1-4,6-7 | sort -k1,1 -k2,2n > tempfile
  # modify sequence name to keep chromosome name only, make sure start position is smaller than end position
  paste <(cut -d: -f1 tempfile) <(cut -f2-6 tempfile) | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' > $OUT
  rm tempfile
done < frog.list

# manually edit frog.list to add species name abbreviation and chromosome number
for i in *.bed; do 
    # grep species name
    PREF=`grep ${i%.bed} frog_chr.list | cut -f 2`
    CHR=`grep ${i%.bed} frog_chr.list | cut -f 3`
    echo "working on $PREF $CHR chromosomes"
    
    # extract the original chromosome name
    awk -v OFS='\t' -v var="${PREF}_chr" '{print $0,var NR}' <(cut -f 1 $i | sort | uniq | head -n $CHR) > ${PREF}.chr
done
# manually examine the chr names

# replace the chromosome names in bed file
for i in *.bed; do
    echo "working on $i"
    # grep species name
    PREF=`grep ${i%.bed} frog_chr.list | cut -f 2`
    # replace chromosome names
    awk -v OFS='\t' 'FNR==NR {id[$1]=$2; next} {if ($1 in id) $1=id[$1];print}' ${PREF}.chr $i > t; mv t $i
    mv $i ${PREF}.bed
done

# create bed file for JT22M_Mar12 (P. triseriata)
IN=/project/lougheedlab/chorus/02_HiFi_HiC/11_curation/04_QC/BUSCO/JT22M_Mar12/run_tetrapoda_odb10/full_table.tsv
paste <(cut -f3-5 $IN) <(cut -f1-2 $IN) <(cut -f7 $IN)  <(cut -f6 $IN) | grep -e 'Complete' -e 'Fragmented' | cut -f1-4,6-7 | sort -k1,1 -k2,2n > tempfile
paste <(cut -d: -f1 tempfile) <(cut -f2-6 tempfile) | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' > JT22M_Mar12.bed
rm tempfile

# prep `all.chr` file: column 1 original chromosome name, column 2 modified chromosome name, manually add P. triseriata
cat *.chr > all.chr
```

**Step 2.2**: Using R script `synteny_block_id.R` to create `.simple` and `.synteny` file. The codes extract synteny blocks spanning at least 3 BUSCO genes (both complete and fragmented). Note that the synteny is pairwise comparison and the order depends on the phylogeny so that synteny plots can be placed right beside the phylogeny for context. 

**Step 2.3**: Manually create `seqids` file (use chromosome names in the bed files). Note the order of frogs (each row) is the same as your phylogeny and the `layout` file; the order of chromosomes is based on your synteny so that they align well and look pretty.

**Step 2.4**: Manually create the `layout` file, make sure there is no blank line in the end. 

**Step 2.5**: run the code:
```bash
python -m jcvi.graphics.karyotype seqids layout --chrstyle=roundrect --nocircles
# run the command without arguments will print the help page
```
Note: the length of the chromosomes on the plot depends on number of genes rather than the actual length of the chromosomes [link](https://github.com/tanghaibao/jcvi/issues/427)

**Step 2.6**: based on the synteny, now you can add color to the synteny blocks to make the plot prettier. You need a `frog_synteny_color.csv` to list the color and coordinates of the synteny blocks (Note if one species have multiple blocks with the same color, separate the coordinates with "," no space). The codes to write new `.simple.color` are in `synteny_block_color.R`. 

Now run the final code with modified file `layout_color`:
```bash
python -m jcvi.graphics.karyotype seqids layout_color --chrstyle=roundrect --nocircles 
#--format=png --outfile=karyotype.png
```

