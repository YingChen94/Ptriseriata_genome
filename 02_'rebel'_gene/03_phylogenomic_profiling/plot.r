### Rscripts Phylogenomic_Profiling.r modified
# toturial: https://almeidasilvaf.github.io/syntenet/articles/syntenet.html

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(pheatmap)
library(cluster)
library(vegan)
# BiocManager::install("syntenet")
library(syntenet)
library(dplyr)
library(ggplot2)

# read network file
network <- read.table("../SynNet-k5s5m15_2cols",header=T)

# Input file is a two-column output from infomap clustering
data <- read.table("../SynNet-k5s5m15_2cols_infoclusters",header=T)
clusters <- data
data$species <- substr(clusters$Gene,1,4) # extract first three letters of the gene name as species id
colnames(data) <-c("Gene","Cluster","Species")
data <- data[,c(1,3,2)]

# Create a named vector of custom species order to plot
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


# Create a metadata to group species into family/suborder
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
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Mesobatrachia", 8),
    rep("Archaeobatrachia", 3))
)
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Archaeobatrachia", 11))
)
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Hyloidea", 10),
    rep("Ranoidea", 6),
    rep("Neobatrachia", 2),
    rep("Mesobatrachia", 8),
    rep("Archaeobatrachia", 3))
)
head(species_annotation)
colnames(species_annotation)[2] <- " "


#---------------------------
# phylogenomic profiling plot

# cluster by rows, Species by columns
profiles <- as.data.frame.matrix(table(data$Cluster,data$Species))
head(profiles)

# plot
png("profile_superorder_2.png", units = "in", width = 8, height = 5, res=400)
plot_profiles(
  profiles, 
  species_annotation, 
  cluster_species = species_order, 
  dist_function = labdsv::dsvdis,
  dist_params = list(index = "ruzicka")
)
# take a bit of time
dev.off()


# Find group-specific clusters
gs_clusters <- find_GS_clusters(profiles, species_annotation)
head(gs_clusters)
# How many family-specific clusters are there?
nrow(gs_clusters)

# Filter profiles matrix to only include group-specific clusters
idx <- rownames(profiles) %in% gs_clusters$Cluster
p_gs <- profiles[idx, ]

# Plot heatmap
png("profile_superfamily_gs.png", units = "in", width = 8, height = 5, res=400)
plot_profiles(
  p_gs, species_annotation, 
  cluster_species = species_order, 
  cluster_columns = TRUE,
  #border_color=NA,
)
dev.off()











