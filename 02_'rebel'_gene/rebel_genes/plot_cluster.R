#---------------------------
# plot network for particular genes

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(syntenet)

# read network file
network <- read.table("../../SynNet-k5s5m15_2cols",header=T)
# read clusters file
clusters <- read.table("../../SynNet-k5s5m15_2cols_infoclusters",header=T)


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

#----------
# color by species
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)

#----------
# color by family
# group species into family
species_annotation <- data.frame(
  Species = species_order,
  Family = c(
    rep("Hylidae", 3), "Dendrobatidae","Eleutherodactylidae",
    rep("Bufonidae", 3), rep("Leptodactylidae",2),
    rep("Ranidae",4),"Pyxicephalidae",
    "Microhylidae",rep("Myobatrachidae",2),rep("Pelobatidae",2),
    rep("Megophryidae",2),rep("Pipidae", 4),
    "Alytidae","Bombinatoridae", "Ascaphidae"))
head(species_annotation)
colnames(species_annotation)[2] <- " "
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)
gene_df_fam <- merge(gene_df, species_annotation)[, c("Gene", " ")] # column name was replaced as empty
head(gene_df_fam)

#----------
# color by suborder
# group species into suborder
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neobatrachia", 18),
    rep("Archaeobatrachia", 11)))
species_annotation <- data.frame(
  Species = species_order,
  Suborder = c(
    rep("Neo_Hyloidea", 10),
    rep("Neobatrachia", 8),
    rep("Archaeobatrachia", 11)))
head(species_annotation)
colnames(species_annotation)[2] <- " "
genes <- unique(c(network$Anchor1, network$Anchor2))
gene_df <- data.frame(
  Gene = genes,
  Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)
gene_df_suborder <- merge(gene_df, species_annotation)[, c("Gene", " ")] # column name was replaced as empty
head(gene_df_suborder)
gene_df_suborder[gene_df_suborder=="Neobatrachia"] <- "A"
#----------

# read BUSCO_ID and corresponding chromosome for each species
busco_pos_chr_element <- read.csv("data/busco_pos_chr_element.csv")
# get cluster number from BUSCO ID
busco_id <- "275471at32523"
(id <- busco_pos_chr_element$Cluster[which(busco_pos_chr_element$BUSCO_ID==busco_id)])
id <- unlist(strsplit(id,,split=", "))
# plot
png("gene_275471at32523.png", units = "in", width = 6, height =5, res=400)
plot_network(network, clusters, 
             cluster_id = id, 
             color_by = gene_df_suborder)
dev.off()



