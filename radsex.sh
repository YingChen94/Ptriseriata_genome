# script to run radsex pipeline to find sex chromosome in P. triseriata
# July 2024

# RADSex 
# paper: https://doi.org/10.1111/1755-0998.13360
# github: https://github.com/SexGenomicsToolkit/radsex
# tutorial: https://sexgenomicstoolkit.github.io/html/index.html

# sample size: I think the program requires a specific number of samples (>30 per sex) to be able to have statistically significant markers
# too divergent samples will have diverged markers (i.e. no identical sequences across the same sex) and thus lead to false negative. 



# step 1. prep data
# - transfer ddRADseq data to Europa. Note RADsex only takes R1 reads, and the reads need to be the same length across samples. 
#	I chose samples of both maculata east and triseriata mitotype, including many US samples. 
# - create a file: column 1 sample file name, column 2 sex. e.g. `JT22F.1	F`

# step 2. Computing the marker depths table. Each thread per sample.
/home/ying_local/radsex/bin/radsex process --input-dir ./samples --output-file markers_table.tsv --threads 60 --min-depth 1


# step 3. Computing the distribution of markers between sexes
/home/ying_local/radsex/bin/radsex distrib --markers-table markers_table.tsv --output-file distribution.tsv --popmap popmap.txt --min-depth 10 --groups M,F

# Extracting markers significantly associated with sex
/home/ying_local/radsex/bin/radsex signif --markers-table markers_table.tsv --output-file signif_markers.tsv --popmap popmap.txt --min-depth 10 --groups M,F
/home/ying_local/radsex/bin/radsex signif --markers-table markers_table.tsv --output-file signif_markers.fasta --popmap popmap.txt --min-depth 10 --groups M,F --output-fasta


# step 4. Aligning markers to a genome
REF=/project/lougheedlab/chorus/02_HiFi_HiC/12_ref/JT22M.softmasked.fasta
/home/ying_local/radsex/bin/radsex map --markers-file markers_table.tsv --output-file alignment_results.tsv --popmap popmap.txt --genome-file $REF --min-quality 20 --min-depth 10 --groups M,F

# check which chromosome the significant markers are mapped to. 
tail -n +3 alignment_results.tsv | sort --version-sort | grep -B2 -A2 "True"

# significant marker position (note the end position is based on blast)
# Marker_id		position					Notes
# 11843592		chr1:393734341-393734482	part of LINE/N1
# 19742804		chr1:393746144-393746258	part of unknown repeat
# 18173673		chr1:394119305-394119446	part of DNA/hAT-Ac and unknown repeat
# 2575995		chr1:394119442-394119583	part of unknown repeat and DNA/TcMar-Tc2
# 12561139		chr1:394472712-394472853	part of DNA/TcMar-Tc2
# 18904198		chr1:394483143-394483284	contain LTR/Pao, part of DNA/TcMar-Tc2	
# 19053281		chr1:394506082-394506210	part of unknown repeats
# 21808671		chr1:394592296-394592437	contain simple repeats
# 20847340		chr1:394613540-394613681	part of LINE/Penelope
# 19330785		chr1:394721434-394721554 	end of the region is DNA/hAT-Tip100
# 9676207		chr1:394873959-394874100	no repeats
# note that 6 markers between these regions are not significantly associated with sex
# easier to see the repeats in JBrowse (i.e. add repeat and gene track to reference genome)
# sex region: 393714122-394874100

# step 5. manual check sex-linked region

# 5.1. check whether all females truly don't have the sequence or just low number of reads (i.e. below threshold 10)
tail -n +2 signif_markers.tsv | cut -f 4,6,8,9,20,22,29,32,45,50 > signif_female.tsv
# except for a few sequence depth of 1, the rest are 0
# -> females truly don't have those sequences 

# 5.2. check ddRAD bam file (i.e. mapping each sample to reference genome) to see whether females are homozygotes and males are heterozygotes in these regions
# Frontenac cluster
cd /global/project/hpcg1145/04_chorus_frog/05_sex_chr
salloc --time=3:00:00
module load StdEnv/2023 vcftools/0.1.16 bcftools/1.19 gcc/12.3 samtools/1.20
# create `regions.txt` to store loci region: chr1:start-end
# create `female_bam.list` to store the bam file path : /global/project/hpcg1145/04_chorus_frog/01_ddRAD_demultiplex/02_plate6789/04_bwa_mem/run_June/out_align
cat female_bam.list | while read bam; do
	samtools index -c $bam
	PREF=$(basename ${bam%.bam})
	xargs -a regions.txt -I {} samtools view $bam {} >> ${PREF}.sam
done
# create `male_bam.list` to store the bam file path : /global/project/hpcg1145/04_chorus_frog/01_ddRAD_demultiplex/02_plate6789/04_bwa_mem/run_June/out_align
cat male_bam.list | while read bam; do
	samtools index -c $bam
	PREF=$(basename ${bam%.bam})
	xargs -a regions.txt -I {} samtools view $bam {} >> ${PREF}.sam
done
# females have either low mapping quality alignments or no alignment in those regions
# males mostly have mapping quality of 60
# -> females have deletions in these regions; Or females have very divergent sequences that don't get mapped to reference genome


# 5.3. check hap2 (maternal haplotype) to see whether the sex-linked region are deletions in maternal haplotype

# 5.3.1 extract the region from hap2 and map to reference genome
REF=/project/lougheedlab/chorus/02_HiFi_HiC/12_ref/JT22M_ref.fa
HAP2=/project/lougheedlab/chorus/02_HiFi_HiC/11_curation/02_remove_mito/hap2/hap2_nomito.fa
samtools faidx $HAP2
samtools faidx $HAP2 "scaffold_1:402000000-407000000" > hap2_region.fasta
# align two fasta
minimap2 -ax asm5 -L $REF hap2_region.fasta -t 88 | samtools sort -O bam > hap2_scaf1_402M_407M.bam
samtools index -c hap2_scaf1_402M_407M.bam
samtools view hap2_scaf1_402M_407M.bam "chr1:393734000-393735000" | less
# -> the sex-linked regions have no mapping

# 5.3.2 do whole genome alignment and check paf file in JBrowse
# alignment is in folder: `/project/lougheedlab/chorus/02_HiFi_HiC/11_curation/05_synteny`
# download `hap2_ref.paf` to local computer
grep -w "chr1" hap2_ref.paf | sort -n -k8 > chr1.paf
grep -w "chr1" hap2_ref.paf | sort -n -k8 | awk '{if ($8>393700000 && $8 < 394900000) {print $0}}' > chr1_393.7M_394.9M.paf
# -> no synteny in the sex-linked region.
# add tracks of repeats and genes
# -> lots of repeats, overlap with `g1992.t1` gene


# 5.3.3 blast to check whether unplaced scaffolds have those sequences
#!/bin/bash
HAP2=/project/lougheedlab/chorus/02_HiFi_HiC/11_curation/02_remove_mito/hap2/hap2_nomito.fa
makeblastdb -in $HAP2 -dbtype nucl -out db_hap2
sed -i 's/:/_/g' signif_markers.fasta
SIGMARKER=../signif_markers.fasta
# blast against reference genome
db_Ptri=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/01_dmrt/02_blast/db_frog/db_Ptri
blastn -db $db_Ptri -query $SIGMARKER -out ref_signi_marker_blast.out -evalue 1e-5 -outfmt "7 std" -num_threads 44 -max_target_seqs 2 -max_hsps 2
# all on chromosome 1
# blast against HAP2
blastn -db db_hap2/db_hap2 -query $SIGMARKER -out hap2_signi_marker_blast.out -evalue 1e-5 -outfmt "7 std" -num_threads 44 -max_target_seqs 5 -max_hsps 5
# top hit is not on chromosome 1


# 5.3.4. check whether two assemblies have gaps in this region, especially hap2, if there are gaps, might be different region to sequence
# reference genome (check the rapid curation `hap1_nomito_v3.tpf`)
seqtk gap JT22M_ref.fa > JT22M_ref_gaps.bed
# gap : 393845946bp
# Hap2
cd /project/lougheedlab/chorus/02_HiFi_HiC/11_curation/02_remove_mito/hap2
seqtk gap hap2_nomito.fa > hap2_nomito_gaps.bed
mv hap2_nomito_gaps.bed /project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/03_radsex/07_assembly_gaps_check/
# gap : 399521129bp
# gap : 401234546bp (= 392.5M in reference)
# gap : 404978053bp (= 397.3M in reference)
# gap : 405908495bp 
# gap : 407694595bp 
# Note: `bed` file can be visualized in JBrowse2. 


# 5.4. check CF8 alignment bam file
BAM=/project/lougheedlab/chorus/04_phylogenomics/02_Bankeretal2020/02_map_WGS/CF8/CF8_sort.bam
samtools view -b $BAM "chr1:393000000-395000000" > CF8_chr1_393M_395M.bam
samtools index -c CF8_chr1_393M_395M.bam
# open in IGV, couldn't see any obvious difference in sequencing depth in those regions, lots of low quality mapping due to repeats


# 5.5. check JT22M alignment bam file
#!/bin/bash
BAM=/project/lougheedlab/chorus/06_conservation_genomics/01_PSMC/JT22M/JT22M_sort.bam
samtools view -b $BAM "chr1:393000000-395000000" > JT22M_chr1_393M_395M.bam
samtools index -c JT22M_chr1_393M_395M.bam
# open in IGV, couldn't see any obvious difference in those regions, lots of low quality mapping due to repeats



# step 6. check whether any testes RNA seq mapped to the sex-linked region
#!/bin/bash
BAM1=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_testis_map/GNRH.bam
BAM2=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_testis_map/SAL.bam
samtools view $BAM1 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > GNRH_393M_395M.sam
samtools view $BAM2 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > SAL_393M_395M.sam

# add header to the sam, otherwise it can't be indexed
samtools view -H ../05.4_CF8/CF8_chr1_393M_395M.bam > header.sam  # extract header only
cat header.sam GNRH_393M_395M.sam | samtools view -b - > GNRH_393M_395M.bam
cat header.sam SAL_393M_395M.sam | samtools view -b - > SAL_393M_395M.bam

# grep `g1992.t1` `g1991.t1` protein sequence
#!/bin/bash
nano sex_genes.txt
DAT=/project/lougheedlab/chorus/03_annotation/03_braker/BRAKER3_compleasm_out/better.aa
seqtk subseq $DAT sex_genes.txt > sex_gene_protein.fasta 
DAT=/project/lougheedlab/chorus/03_annotation/03_braker/BRAKER3_compleasm_out/better.codingseq
seqtk subseq $DAT sex_genes.txt > sex_gene_codingseq.fasta 
# blast to see whether X chromosome also have g1992.t1
blastn -db ../05.3_hap2/db_hap2/db_hap2 -query sex_gene_codingseq.fasta -out hap2_g1992.out -evalue 1e-5 -outfmt "7 std" -num_threads 44 -max_target_seqs 5 -max_hsps 5


# check other tissues
BAM1=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_RNA_map/JT22M-eye.bam
samtools view $BAM1 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > JT22M-eye_393M_395M.sam
cat header.sam JT22M-eye_393M_395M.sam | samtools view -b - > JT22M-eye_393M_395M.bam
samtools index JT22M-eye_393M_395M.bam

BAM2=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_RNA_map/JT22M-kidney.bam
samtools view $BAM2 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > JT22M-kidney_393M_395M.sam
cat header.sam JT22M-kidney_393M_395M.sam | samtools view -b - > JT22M-kidney_393M_395M.bam
samtools index JT22M-kidney_393M_395M.bam

BAM3=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_RNA_map/JT22M-liver.bam
samtools view $BAM3 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > JT22M-liver_393M_395M.sam
cat header.sam JT22M-liver_393M_395M.sam | samtools view -b - > JT22M-liver_393M_395M.bam
samtools index JT22M-liver_393M_395M.bam

BAM4=/project/lougheedlab/chorus/05_comparative_genomics/03_sex_chr/02_RNA_map/YC22Cec1-liver.bam
samtools view $BAM4 | awk -v OFS='\t' '{if ($3=="chr1" && $4>393000000 && $4<395000000) {print $0}}' > YC22Cec1-liver_393M_395M.sam
cat header.sam YC22Cec1-liver_393M_395M.sam | samtools view -b - > YC22Cec1-liver_393M_395M.bam
samtools index YC22Cec1-liver_393M_395M.bam





