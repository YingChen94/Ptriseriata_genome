# scripts to identify BUSCO synteny blocks
# Ying Chen
# Mar 2024

library(dplyr)

rm(list=ls())
# set working directory to the same as the R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Note that those frog bed files chromosome names are original
# the genome order is according to the phylogeny in order to make nice plots later
genome1 <- read.delim("frog_bed_files/JT22M_Mar12.bed",header=F)
genome2 <- read.delim("frog_bed_files/Hsar.bed",header=F)
genome3 <- read.delim("frog_bed_files/Debr.bed",header=F)
genome4 <- read.delim("frog_bed_files/Rimi.bed",header=F)
genome5 <- read.delim("frog_bed_files/Bbuf.bed",header=F)
genome6 <- read.delim("frog_bed_files/Lfus.bed",header=F)
genome7 <- read.delim("frog_bed_files/Epus.bed",header=F)
genome8 <- read.delim("frog_bed_files/Ecoq.bed",header=F)
genome9 <- read.delim("frog_bed_files/Pcor.bed",header=F)
genome10 <- read.delim("frog_bed_files/Mfle.bed",header=F)
genome11 <- read.delim("frog_bed_files/Rkuk.bed",header=F)
genome12 <- read.delim("frog_bed_files/Rtem.bed",header=F)
genome13 <- read.delim("frog_bed_files/Lsyl.bed",header=F)
genome14 <- read.delim("frog_bed_files/Pads.bed",header=F)
genome15 <- read.delim("frog_bed_files/Gcar.bed",header=F)
genome16 <- read.delim("frog_bed_files/Pcul.bed",header=F)
genome17 <- read.delim("frog_bed_files/Llei.bed",header=F)
genome18 <- read.delim("frog_bed_files/Sbom.bed",header=F)
genome19 <- read.delim("frog_bed_files/Xtro.bed",header=F)
genome20 <- read.delim("frog_bed_files/Hboe.bed",header=F)
genome21 <- read.delim("frog_bed_files/Dpic.bed",header=F)
genome22 <- read.delim("frog_bed_files/Bbom.bed",header=F)
genome23 <- read.delim("frog_bed_files/Atru.bed",header=F)


genomeLIST <- list(genome1,genome2,genome3,genome4,genome5,genome6,genome7,genome8,genome9,genome10,
                   genome11,genome12,genome13,genome14,genome15,genome16,genome17,genome18,genome19,genome20,
                   genome21,genome22,genome23)
genomeLIST <- lapply(genomeLIST, setNames, nm=c("chromosome","gene_start","gene_end","busco_id","score","orientation"))

speciesLIST <- c("JT22M","Hsar","Debr","Rimi",
                 "Bbuf","Lfus","Epus","Ecoq","Pcor","Mfle",
                 "Rkuk","Rtem","Lsyl","Pads","Gcar",
                 "Pcul","Llei","Sbom","Xtro",
                 "Hboe","Dpic","Bbom","Atru")

# read the chromosome name list
chr_list <- read.delim("frog_bed_files/all.chr",header=F)

for (k in 1:(length(genomeLIST)-1)) {
  cat("now working on",speciesLIST[k],speciesLIST[k+1],"\n")
  genome1 <- genomeLIST[[k]]
  genome2 <- genomeLIST[[k+1]]
  
  # only keep genes on main chromosomes
  genome1 <- genome1 %>% filter(chromosome %in% chr_list$V2)
  genome2 <- genome2 %>% filter(chromosome %in% chr_list$V2)
  
  # only keep shared busco genes 
  shared_busco_list <- intersect(genome1$busco_id,genome2$busco_id)
  genome1 <- genome1 %>% filter(busco_id %in% shared_busco_list)
  genome2 <- genome2 %>% filter(busco_id %in% shared_busco_list)
  
  # sort the bed file by chromosome and gene positions
  genome1 <- genome1[order(genome1[,1],genome1[,2]),]
  genome2 <- genome2[order(genome2[,1],genome2[,2]),]
  
  # make .simple dataset 
  dat_synteny <- data.frame(matrix(ncol=12, nrow = 10000))
  colnames(dat_synteny) <- c("chromosome","start","end","lower","upper","chromosome2","start2","end2","lower2","upper2","score","strand")
  x=1 # x tracks the dat_synteny lines
  
  # identify synteny blocks
  for (i in 1:length(genome1$chromosome)){ # i tracks the lines in genome 1
    if (i==1){
      dat_synteny$chromosome[x] <- genome1$chromosome[i]
      dat_synteny$start[x] <- genome1$busco_id[i]
      dat_synteny$chromosome2[x] <- genome2$chromosome[which(genome2$busco_id==genome1$busco_id[i])]
      next
    } else if (which(genome2$busco_id==genome1$busco_id[i])==(which(genome2$busco_id==genome1$busco_id[i-1])+1) & # in synteny
               genome1$chromosome[i]==dat_synteny$chromosome[x] & # check genes are on the same chromosome
               genome2$chromosome[which(genome2$busco_id==genome1$busco_id[i])]==dat_synteny$chromosome2[x]) { 
      dat_synteny$strand[x] <- "+"
      next
    } else if (which(genome2$busco_id==genome1$busco_id[i])==(which(genome2$busco_id==genome1$busco_id[i-1])-1) & # in synteny
               genome1$chromosome[i]==dat_synteny$chromosome[x] & # check genes are on the same chromosome
               genome2$chromosome[which(genome2$busco_id==genome1$busco_id[i])]==dat_synteny$chromosome2[x]) { 
      dat_synteny$strand[x] <- "-"
      next
    } else {
      dat_synteny$end[x] <- genome1$busco_id[i-1]
      x <- x+1
      dat_synteny$chromosome[x] <- genome1$chromosome[i]
      dat_synteny$start[x] <- genome1$busco_id[i]
      dat_synteny$chromosome2[x] <- genome2$chromosome[which(genome2$busco_id==genome1$busco_id[i])]
    }
  }
  
  # fill in empty cell in the synteny block datasheet 
  dat_synteny <- dat_synteny %>% filter(!is.na(start))
  dat_synteny$strand[is.na(dat_synteny$strand)] <- "+"
  dat_synteny$score <- "99"
  for (i in 1:length(dat_synteny$start)){
    if (dat_synteny$strand[i]=="-"){
      dat_synteny$start2[i] <- dat_synteny$end[i]
      dat_synteny$end2[i] <- dat_synteny$start[i]
    } else {
      dat_synteny$start2[i] <- dat_synteny$start[i]
      dat_synteny$end2[i] <- dat_synteny$end[i]
    }
  }
  dat_synteny <- dat_synteny %>% filter(!is.na(end))
  
  # fill in the coordinates
  for (m in 1:length(dat_synteny$chromosome)){
    # get the lower boundary of the synteny block for genome 1
    dat_synteny$lower[m] <- min(c(genome1$gene_start[which(genome1$busco_id==dat_synteny$start[m])],
                                  genome1$gene_end[which(genome1$busco_id==dat_synteny$start[m])],
                                  genome1$gene_start[which(genome1$busco_id==dat_synteny$end[m])],
                                  genome1$gene_end[which(genome1$busco_id==dat_synteny$end[m])]))
    # get the upper boundary of the synteny block for genome 1
    dat_synteny$upper[m] <- max(c(genome1$gene_start[which(genome1$busco_id==dat_synteny$start[m])],
                                  genome1$gene_end[which(genome1$busco_id==dat_synteny$start[m])],
                                  genome1$gene_start[which(genome1$busco_id==dat_synteny$end[m])],
                                  genome1$gene_end[which(genome1$busco_id==dat_synteny$end[m])]))
    # get the lower boundary of the synteny block for genome 2
    dat_synteny$lower2[m] <- min(c(genome2$gene_start[which(genome2$busco_id==dat_synteny$start[m])],
                                   genome2$gene_end[which(genome2$busco_id==dat_synteny$start[m])],
                                   genome2$gene_start[which(genome2$busco_id==dat_synteny$end[m])],
                                   genome2$gene_end[which(genome2$busco_id==dat_synteny$end[m])]))
    # get the upper boundary of the synteny block for genome 2
    dat_synteny$upper2[m] <- max(c(genome2$gene_start[which(genome2$busco_id==dat_synteny$start[m])],
                                   genome2$gene_end[which(genome2$busco_id==dat_synteny$start[m])],
                                   genome2$gene_start[which(genome2$busco_id==dat_synteny$end[m])],
                                   genome2$gene_end[which(genome2$busco_id==dat_synteny$end[m])]))
  }

  
  # now check how many genes each synteny block include
  dat_synteny$n_gene_in_synteny <- NA
  for (j in 1:length(dat_synteny$chromosome)){
    dat_synteny$n_gene_in_synteny[j] <- abs(which(genome1$busco_id==dat_synteny$start[j]) - which(genome1$busco_id==dat_synteny$end[j]))+1
  }
  cat("total synteny blocks",length(dat_synteny$chromosome),"\n")
  
  # discard blocks with only 1 gene/2 genes
  dat_synteny <- dat_synteny %>% filter(n_gene_in_synteny > 2)
  cat("total synteny blocks spanning at least 2 genes",length(dat_synteny$chromosome),"\n")
  
  # write out this full table 
  write.table(dat_synteny,paste0("synteny_files_v2_timetree/",speciesLIST[k],"_",speciesLIST[k+1],".synteny",sep=""),sep = "\t", row.names=FALSE, col.names=T, quote=FALSE)
  
  # make .simple file
  dat_simple <- dat_synteny[,c(-1,-4,-5,-6,-9,-10,-13)] # delete chromosome column and number of genes columns
  write.table(dat_simple,paste0("synteny_files_v2_timetree/",speciesLIST[k],"_",speciesLIST[k+1],".simple",sep=""),sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}





