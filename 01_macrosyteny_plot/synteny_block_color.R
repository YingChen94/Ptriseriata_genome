# scripts to identify BUSCO synteny blocks
# Ying Chen

library(dplyr)

rm(list=ls())
# set working directory to the same as the R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read the full synteny table with chromosome and coordinates
synteny1 <- read.delim("synteny_files_v2_timetree/JT22M_Hsar.synteny",header=T)
synteny2 <- read.delim("synteny_files_v2_timetree/Hsar_Debr.synteny",header=T)
synteny3 <- read.delim("synteny_files_v2_timetree/Debr_Rimi.synteny",header=T)
synteny4 <- read.delim("synteny_files_v2_timetree/Rimi_Bbuf.synteny",header=T)
synteny5 <- read.delim("synteny_files_v2_timetree/Bbuf_Lfus.synteny",header=T)
synteny6 <- read.delim("synteny_files_v2_timetree/Lfus_Epus.synteny",header=T)
synteny7 <- read.delim("synteny_files_v2_timetree/Epus_Ecoq.synteny",header=T)
synteny8 <- read.delim("synteny_files_v2_timetree/Ecoq_Pcor.synteny",header=T)
synteny9 <- read.delim("synteny_files_v2_timetree/Pcor_Mfle.synteny",header=T)
synteny10 <- read.delim("synteny_files_v2_timetree/Mfle_Rkuk.synteny",header=T)
synteny11 <- read.delim("synteny_files_v2_timetree/Rkuk_Rtem.synteny",header=T)
synteny12 <- read.delim("synteny_files_v2_timetree/Rtem_Lsyl.synteny",header=T)
synteny13 <- read.delim("synteny_files_v2_timetree/Lsyl_Pads.synteny",header=T)
synteny14 <- read.delim("synteny_files_v2_timetree/Pads_Gcar.synteny",header=T)
synteny15 <- read.delim("synteny_files_v2_timetree/Gcar_Pcul.synteny",header=T)
synteny16 <- read.delim("synteny_files_v2_timetree/Pcul_Llei.synteny",header=T)
synteny17 <- read.delim("synteny_files_v2_timetree/Llei_Sbom.synteny",header=T)
synteny18 <- read.delim("synteny_files_v2_timetree/Sbom_Xtro.synteny",header=T)
synteny19 <- read.delim("synteny_files_v2_timetree/Xtro_Hboe.synteny",header=T)
synteny20 <- read.delim("synteny_files_v2_timetree/Hboe_Dpic.synteny",header=T)
synteny21 <- read.delim("synteny_files_v2_timetree/Dpic_Bbom.synteny",header=T)
synteny22 <- read.delim("synteny_files_v2_timetree/Bbom_Atru.synteny",header=T)

syntenyLIST <- list(synteny1,synteny2,synteny3,synteny4,synteny5,synteny6,synteny7,synteny8,synteny9,synteny10,
                    synteny11,synteny12,synteny13,synteny14,synteny15,synteny16,synteny17,synteny18,synteny19,synteny20,
                    synteny21,synteny22)

speciesLIST <- c("JT22M","Hsar","Debr","Rimi",
                 "Bbuf","Lfus","Epus","Ecoq","Pcor","Mfle",
                 "Rkuk","Rtem","Lsyl","Pads","Gcar",
                 "Pcul","Llei","Sbom","Xtro",
                 "Hboe","Dpic","Bbom","Atru")

# read the csv file with color info
col_coord_list <- read.csv("synteny_files_v2_timetree/frog_synteny_color.csv",header=T)

# add color to the synteny file
for (k in 2:(length(syntenyLIST)+1)) { # go through each pairwise synteny
  synteny_to_color <- syntenyLIST[[k-1]]
  species <- speciesLIST[k-1]
  
  for (i in 1:length(col_coord_list$color)){ # go through each synteny color
    if (col_coord_list[i,k]!=""){ # !is.na(col_coord_list[i,k])
      for (sec in 1:length(strsplit(col_coord_list[i,k],",")[[1]])){ # there might be multiple segments, go through each segment
        x <- strsplit(col_coord_list[i,k],",")[[1]][sec]
        chr_to_color <- paste(strsplit(x,"_")[[1]][1],"_",strsplit(x,"_")[[1]][2],sep="")
        left_coord <- as.numeric(strsplit(x,"_")[[1]][3])-1
        right_coord <- as.numeric(strsplit(x,"_")[[1]][4])+1
        rows_2_color <- which(synteny_to_color$chromosome==chr_to_color &
                                synteny_to_color$lower > left_coord &
                                synteny_to_color$upper < right_coord)
        for (m in 1:length(rows_2_color)){
          synteny_to_color$start[rows_2_color[m]] <- paste0(col_coord_list$color[i],"*",synteny_to_color$start[rows_2_color[m]])
        }
      }
    }
  }
  
  dat_simple_col <- synteny_to_color[,c(-1,-4,-5,-6,-9,-10,-13)] # delete chromosome column and number of genes columns
  cat("write out",paste0(speciesLIST[k-1],"_",speciesLIST[k],".simple.color",sep=""),"\n")
  write.table(dat_simple_col,paste0("synteny_files_v2_timetree/",speciesLIST[k-1],"_",speciesLIST[k],".simple.color",sep=""),sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}















