# 'Rebel' BUSCO gene in anurans
Author: Ying Chen

Network-based synteny analysis using BUSCO genes of 29 chromosome-level frog assemblies (28 NCBI frogs + *P. triseriata*) using [SynNet](https://almeidasilvaf.github.io/syntenet/articles/syntenet.html) and [PanSyn](https://github.com/yhw320/PanSyn). 

### Step 1. Data preparation

BUSCO analyses already conducted. Here to prepare the following files for PanSyn pipeline.

- [x] Protein sequence (*.pep) **Note: one protein sequence per gene**
- [x] Gene coordinates (*_simplified.gff): chromosomeID, geneID, start position, end position and strand (five columns without a header). **Note: one protein sequence per gene**
- [x] chromosome length (*.len) (no need to SynNet)
- [x] Gene coordinates for chromosomes only (*_simplified_chr.gff)

```bash
# create a list
printf '%s\n' /project/lougheedlab/chorus/frog_genome/01_BUSCO/results/* > frog_folder.list

# prep .pep and .gff file for each species
while read p; do 
    temp=$(basename $p)
    PREFIX=${temp%_GC*}
    echo "working on" $PREFIX
    # prep *.pep 
    cat $p/run_tetrapoda_odb10/busco_sequences/*/*.faa > $PREFIX.pep
    # prep *.gff
    grep -e "Complete" -e "Duplicated" -e "Fragmented" $p/run_tetrapoda_odb10/full_table.tsv | awk  -v OFS='\t' '{print $3,$3,$4,$5,$6}'> test.gff
    paste <(cut -f 1 test.gff | awk -F":" '/:/{print $1}') <(cut -f 2-5 test.gff) | sort -k1,1 -k3n,3 | awk '$3 > $4 { temp = $4; $4 = $3; $3 = temp } 1' OFS='\t' > $PREFIX.gff
    rm test*
done < frog_folder.list

# create a `frog_file.name`: first column abbreviated species name, second column the current species name

# replace file names
for i in *.gff; do 
mv $i $(grep ${i%.gff} frog_file.name | cut -f 1).gff
done
for i in *.pep; do 
mv $i $(grep ${i%.pep} frog_file.name | cut -f 1).pep
done

# prepare *.len file
for i in *.gff; do
cut -f 1 $i | uniq > temp
awk -v OFS='\t' -v pref="${i%.gff}" '{print pref"_chr"NR, $0}' temp > ${i%.gff}.chr 
rm temp
done

# only keep chromosomes in *.chr, manually check each file

# prep *_simplified.gff: modify chromosome names for each gff 
for i in *.gff; do 
    awk -v OFS='\t' 'FNR==NR {id[$2]=$1; next} {if ($1 in id) $1=id[$1];print}' ${i%.gff}.chr $i > ${i%.gff}_simplified.gff
done

# prepare .len: add length of each chromosomes (no need for SynNet)
for i in /project/lougheedlab/chorus/frog_genome/seq/*.fai; do
    NAME=$(basename $i)
    PREFIX=$(grep ${NAME%_GC*} frog_file.name | cut -f 1)
    awk -v OFS='\t' 'FNR==NR {len[$1]=$2; next} {if ($2 in len) print $0, len[$2]}' $i $PREFIX.chr > $PREFIX.len
done
```

### Step 2 (PanSyn step 2C). SynNet
- [x] Protein sequence (*.pep)
- [x] Gene coordinates (*_simplified.gff)
- [x] Species list (Species_list.txt)
- [x] provided script pheatmap.r

```bash
mkdir peps gffs 
cp ../../data_BUSCO_29frogs/*_simplified.gff ./gffs/
cp ../../data_BUSCO_29frogs/*.pep ./peps/
```

`SynetBuild-X` has three parameters: 
- n1 = Number of best non-self alignment hits (e.g., 5)
- n2 = Minimum of anchors for a synteny block; higher, stricter; (e.g., 5)
- n3 = Maximum of genes allowed as the GAP between Anchors; fewer, stricter; (e.g., 25)
- n4 = Number of threads (e.g., 12)

k5s5m15 following mammal genomes in [Zhao and Schranz 2019](https://www.pnas.org/doi/full/10.1073/pnas.1801757116).

```bash
mkdir k5s5m15_output
conda activate pansyn
prepare_for_synnet -g gffs/ -p peps/ -o k5s5m15_output/
# create Species_list.txt
for i in ./gffs/*.gff; do 
PREF=$(basename ${i%_simplified.gff})
echo $PREF
done > Species_list.txt

# run SynNet
cd k5s5m15_output/prepare_data/
SynetBuild-X 5 5 15 88 ../../Species_list.txt
```

Output:<br>
`output/prepare_data/SynNetBuild*-SynNet-k5s5m15/SynNet-k5s5m15` network database with four columns: Block_ID, Block_Score, Gene1, and Gene2 (Gene 1 and Gene 2 are a syntenic gene pair)
&nbsp; 

### Step 3 (PanSyn step 4B). SynNet phylogenomic profiling

Use infomap script in PanSyn pipeline to get `SynNet-k5s5m15_2cols` and `SynNet-k5s5m15_2cols_infoclusters` files

```bash
mkdir out
cp /project/lougheedlab/chorus/05_comparative_genomics/04_PanSyn/02C/BUSCO_29frogs/k5s5m15_output/prepare_data/SynNetBuild20240716_1749-SynNet-k5s5m15/SynNet-k5s5m15 ./
awk '{print $3"\t"$4}' SynNet-k5s5m15 > SynNet-k5s5m15_2cols
conda activate pansyn
infomap SynNet-k5s5m15_2cols SynNet-k5s5m15_2cols_infoclusters
sed -i '1s/^/Anchor1\tAnchor2\n/' SynNet-k5s5m15_2cols
sed -i '1s/.*/Gene\tCluster/' SynNet-k5s5m15_2cols_infoclusters
```
Use `SynNet-k5s5m15_2cols` and `SynNet-k5s5m15_2cols_infoclusters` to make phylogenomic profiling plots using script `01_profiling_plot/plot.r`. [syntenet tutorial](https://almeidasilvaf.github.io/syntenet/articles/syntenet.html)

### Step 4 (PanSyn step 4B). 'rebel' genes

To be able to check which BUSCO genes have two or more clusters in R, I need to know which BUSCO ID each sequence name corresponds to. Create a file BUSCO_ID_sequence.list: first column BUSCO ID, second column sequence name

```bash
echo -e "BUSCO_ID\tSequence" > BUSCOgeneID_SppSeqName.list
while read p; do 
    temp=$(basename $p)
    temp2=${temp%_GC*}
    PREF=$(grep $temp2 /project/lougheedlab/chorus/05_comparative_genomics/04_PanSyn/data_BUSCO_29frogs/frog_file.name | cut -f 1)
    echo "working on" $PREF
    # grep BUSCO ID and the sequence name for each species
    grep -e "Complete" -e "Duplicated" -e "Fragmented" $p/run_tetrapoda_odb10/full_table.tsv | awk  -v OFS='\t' -v var="${PREF}_" '{print $1,var $3}' | sed 's/:/b/g' >> BUSCOgeneID_SppSeqName.list
done < /project/lougheedlab/chorus/05_comparative_genomics/04_PanSyn/data_BUSCO_29frogs/frog_folder.list
```


