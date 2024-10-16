# R script to look for `rebel` genes that have lineage-specific synteny
# tutorial: https://almeidasilvaf.github.io/syntenet/articles/syntenet.html

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(pheatmap)
library(cluster)
library(vegan)
# BiocManager::install("syntenet")
library(syntenet)
library(dplyr)
library(ggplot2)
library(tidyr)

# Input file is a two-column output from infomap clustering
data <- read.table("../../SynNet-k5s5m15_2cols_infoclusters",header=T)
# add a column of species
data$species <- substr(data$Gene,1,4) # extract first three letters of the gene name as species id
colnames(data) <-c("Gene","Cluster","Species")
data <- data[,c(1,3,2)]

# Create a named vector of custom species order to plot
# this is needed for species_annotation metadata
species_order <- setNames(
  # vector elements
  c("Ptri","Hsar","Debr","Rimi","Ecoq","Bbuf","Bgar","Bvir","Lfus","Epus",
    "Rkuk","Rtem","Rmus","Lsyl","Pads","Gcar","Pcor","Mfle","Sbom","Pcul",
    "Llei","Lail","Xbor","Xlae","Xtro","Hboe","Dpic","Bbom","Atru"),
  # vector names
  c("P. triseriata", "H. sarda", "D. ebraccatus","R. imitator","E. coqui",
    "B. bufo", "B. gargarizans", "B. viridis","L. fuscus","E. pustulosus",
    "R. kukunoris","R. temporaria","R. muscosa", "L. sylvaticus","P. adspersus",
    "G. carolinensis", "P. corroboree", "M. fleayi","S. bombifrons","P. cultripes",
    "L. leishanense","L. ailaonicum","X. borealis","X. laevis", "X. tropicalis", 
    "H. Boettgeri", "D. pictus", "B. bombina", "A. truei")
)


#---------------------------
# phylogenomic profiling, needed to find group-specific clusters

# cluster by rows, Species by columns
profiles <- as.data.frame.matrix(table(data$Cluster,data$Species))
head(profiles)

#---------------------------
# read the file with BUSCO_ID and corresponding chromosome for each species
busco_pos_chr_element <- read.csv("data/busco_pos_chr_element.csv")


#---------------------------
# find genes with lineage-specific network clusters
# you will get the cluster name here, you can get the actual BUSCO gene ID in next session

#--
# group species into family
species_annotation <- data.frame(
  Species = species_order,
  Family = c(
    rep("Hylidae", 3), "Dendrobatidae","Eleutherodactylidae",
    rep("Bufonidae", 3), rep("Leptodactylidae",2),
    rep("Ranidae",4),"Pyxicephalidae",
    "Microhylidae",rep("Myobatrachidae",2),rep("Pelobatidae",2),
    rep("Megophryidae",2),rep("Pipidae", 4),
    "Alytidae","Bombinatoridae", "Ascaphidae")
)
head(species_annotation)
colnames(species_annotation)[2] <- " "

# Find family-specific clusters
gs_clusters_family <- find_GS_clusters(profiles, species_annotation)
head(gs_clusters_family)
# How many family-specific clusters are there?
nrow(gs_clusters_family)
# 185


#--
# group species into three suborders
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Mesobatrachia", 8),  # this is no longer used but helpful to see split between Neo+Mesobatrachia vs Archarobatrachia
    rep("Archaeobatrachia", 3))
)
head(species_annotation)
colnames(species_annotation)[2] <- " "

# Find group-specific clusters
gs_clusters_suborder <- find_GS_clusters(profiles, species_annotation)
head(gs_clusters_suborder)
# How many suborder-specific clusters are there?
nrow(gs_clusters_suborder)
# 75

#--
# group species into two suborders
# Archaeobatrachia has too few species
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Archaeobatrachia", 11))
)
head(species_annotation)
colnames(species_annotation)[2] <- " "

# Find group-specific clusters
gs_clusters_suborder <- find_GS_clusters(profiles, species_annotation)
head(gs_clusters_suborder)
# How many suborder-specific clusters are there?
nrow(gs_clusters_suborder)
# 69

#--
# group species into suborder but split Hyloidea and Raniodea
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Hyloidea", 10),
    rep("Ranoidea", 6),
    rep("Neobatrachia_split", 2),
    rep("Archaeobatrachia", 11))
)
head(species_annotation)
colnames(species_annotation)[2] <- " "

# Find group-specific clusters
gs_clusters_group <- find_GS_clusters(profiles, species_annotation)
head(gs_clusters_group)
# How many family-specific clusters are there?
nrow(gs_clusters_group)
# 89

#--
# combine three into one dataset
#gs_clusters <- rbind(gs_clusters_family,gs_clusters_suborder,gs_clusters_group)
gs_clusters <- rbind(gs_clusters_suborder,gs_clusters_group)

# keep only clusters with >90% percent
gs_clusters <- gs_clusters %>% filter(Percentage > 90) %>% filter(Group != "Neobatrachia_split")

# add cluster information
busco_pos_chr_element2 <- separate_rows(busco_pos_chr_element, Cluster)
gs_clusters <- left_join(gs_clusters,busco_pos_chr_element2,by="Cluster")
gs_clusters <- distinct(gs_clusters)

write.csv(gs_clusters, "rebel_gene_90per.csv",row.names = F)



