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
busco_pos_chr_element <- read.csv("../01_prep_BUSCO_dataset/busco_pos_chr_element.csv")


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



#---------------------------
# OLD codes
#---------------------------
#---------------------------
#---------------------------
# plot network for particular genes

# read network file
network <- read.table("../SynNet-k5s5m15_2cols",header=T)
# read clusters file
clusters <- read.table("../SynNet-k5s5m15_2cols_infoclusters",header=T)

# color by species
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)

# color by family
#gene_df <- merge(gene_df, species_annotation)[, c("Gene", " ")] # column name was replaced as empty
#head(gene_df)


plot_network(network, clusters, 
             cluster_id = c("5161","5258","3466"), 
             color_by = gene_df)



# --- 
# busco genes with 2 clusters
busco_2 <- busco_count %>% filter(n==2) 
cluster_2 <- cluster_busco[which(cluster_busco$BUSCO_ID %in% busco_2$BUSCO_ID),]

# check any gene with 2 clusters overlap with the lineage-specific list
overlap <- cluster_2$Cluster[which(cluster_2$Cluster %in% gs_clusters$Cluster)]
i=3
id <- cluster_2$Cluster[which(cluster_2$BUSCO_ID==cluster_2$BUSCO_ID[which(cluster_2$Cluster==overlap[i])])]
#png("plot_cluster_122359at32523.png", units = "in", width = 5.5, height = 4.5, res=400)
plot_network(network, clusters, 
             cluster_id = id, 
             color_by = gene_df)
#dev.off()
id %in% gs_clusters$Cluster


b <- "122359at32523"
b <- "301105at32523"
(id <- cluster_busco$Cluster[which(cluster_busco$BUSCO_ID == b)])
plot_network(network, clusters, 
             cluster_id = id, 
             color_by = gene_df)

# grep the BUSCO gene ID
id # the culster ID
(temp_busco <- cluster_2$BUSCO_ID[which(cluster_2$Cluster==overlap[i])])
(temp_cluster_id <- cluster_2$Cluster[which(cluster_2$BUSCO_ID==temp_busco)])
busco$Gene[which(busco$BUSCO_ID=="122359at32523")]
# split between neobatrachia vs archeao-meso
# BUSCO gene, cluster number
# 122359at32523, 1298 4890, UTP6, small subunit processome component
# 125111at32523, 1974 5050, UDP-GalNAc:beta-1, 3-N-acetylgalactosaminyltransferase 2
# 211392at32523, 1698 4951
# 211410at32523, 4770 4907
# 249397at32523, 2775 5003
# 268229at32523, 3253 4881
# 68460at32523, 4926 5163
# 83450at32523, 109 5030
# 92890at32523, 4782 4993

# split between archaeobatrachia vs meso-neo
# 301105at32523, 4823 4983 - spermatogenesis



#---
# busco genes with 3 clusters
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Mesobatrachia", 8),
    rep("Archaeobatrachia", 3))
)
colnames(species_annotation)[2] <- " "
gs_clusters <- find_GS_clusters(profiles, species_annotation)
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)
gene_df <- merge(gene_df, species_annotation)[, c("Gene", " ")] # column name was replaced as empty



busco_3 <- busco_count %>% filter(n==3) %>% select(BUSCO_ID)
cluster_3 <- cluster_busco[which(cluster_busco$BUSCO_ID %in% busco_3$BUSCO_ID),]

overlap <- cluster_3$Cluster[which(cluster_3$Cluster %in% gs_clusters$Cluster)]
i=20
id <- cluster_3$Cluster[which(cluster_3$BUSCO_ID==cluster_3$BUSCO_ID[which(cluster_3$Cluster==overlap[i])])]
#png("plots/cluster_174944at32523.png", units = "in", width = 5.5, height = 4.5, res=400)
plot_network(network, clusters, cluster_id=c(4855,4899,5219), color_by=gene_df)
#dev.off()
id %in% gs_clusters$Cluster

# check the BUSCO ID
id
cluster_3$BUSCO_ID[which(cluster_3$Cluster==overlap[i])]
busco$Gene[which(busco$BUSCO_ID=="256818at32523")]
cluster_3$Cluster[which(cluster_3$BUSCO_ID=="181890at32523")]

# split between 3 suborder
# 174944at32523, 4815 4896 5215

# Hyloidea and Ranoidea together
# 181890at32523, 3285 4406 4903
# 256818at32523, 1417 4901 1313

# neobatrachia together
# 50206at32523, 4766 5017 5223
# 161343at32523, 4855 4899 5219

# Hyloidea specific
# 223932at32523, 4851 5011 5203
# 275471at32523, 2396 4934 5303 - mitochondrial large ribosomal subunit
png("plots/cluster_223932at32523.png", units = "in", width = 5.5, height = 4.5, res=400)
plot_network(network, clusters, 
             cluster_id = c(2396,4934,5303), 
             color_by = gene_df)
dev.off()

#---
# busco genes with 4 clusters
busco_4 <- busco_count %>% filter(n==4) %>% select(BUSCO_ID)
cluster_4 <- cluster_busco[which(cluster_busco$BUSCO_ID %in% busco_4$BUSCO_ID),]

overlap <- cluster_4$Cluster[which(cluster_4$Cluster %in% gs_clusters$Cluster)]
i=3
id <- cluster_4$Cluster[which(cluster_4$BUSCO_ID==cluster_4$BUSCO_ID[which(cluster_4$Cluster==overlap[i])])]
plot_network(network, clusters, cluster_id = id, color_by = gene_df)
id %in% gs_clusters$Cluster

# check the BUSCO ID
id
cluster_4$BUSCO_ID[which(cluster_4$Cluster==overlap[i])]
busco$Gene[which(busco$BUSCO_ID=="106876at32523")]

# neobatrachia together
# 262758at32523, 4803 5014 5265 5418

# Hyloidea and Ranoidea partially together
# 266561at32523, 2695 5012 5228 5325


#---
# family-specific network
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
colnames(species_annotation)[2] <- " "
gs_clusters <- find_GS_clusters(profiles, species_annotation)
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)
gene_df <- merge(gene_df, species_annotation)[, c("Gene", " ")] # column name was replaced as empty


# check which busco genes
cluster_busco[which(cluster_busco$Cluster==5239),]
id <- cluster_busco$Cluster[which(cluster_busco$BUSCO_ID=="60551at32523")]
plot_network(network, clusters, cluster_id = id, color_by = gene_df)



