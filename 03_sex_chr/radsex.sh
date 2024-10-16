# script to run radsex pipeline to find sex chromosome in P. triseriata
# Ying chen

# RADSex 
# paper: https://doi.org/10.1111/1755-0998.13360
# github: https://github.com/SexGenomicsToolkit/radsex
# tutorial: https://sexgenomicstoolkit.github.io/html/index.html


# step 1. prep data
# - ddRADseq data. Note RADsex only takes R1 reads, and the reads need to be the same length across samples. 
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



