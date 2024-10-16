# scripts to count percent of BUSCO genes with 1,2,3... clusters

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# count number of clusters for each BUSCO gene

# read a file: BUSCO_ID and its position info on each species
busco <- read.table("data/BUSCOgeneID_SppSeqName.list",header=T)
colnames(busco) <- c("BUSCO_ID","Gene")
head(busco)
length(unique(busco$BUSCO_ID))
# total 5233 BUSCO genes (77 BUSCO genes are missing in all species)

# create a dataset to store BUSCO ID and its corresponding cluster number
clusters <- read.table("SynNet-k5s5m15_2cols_infoclusters",header=T)
cluster_busco <- left_join(clusters,busco,by="Gene")[,c(3,2)]
cluster_busco <- left_join(busco,clusters,by="Gene")[,c(1,3)]
cluster_busco <- cluster_busco[order(cluster_busco$BUSCO_ID),]
cluster_busco <- distinct(cluster_busco)
cluster_busco <- drop_na(cluster_busco)

# count how many cluster each BUSCO gene has
busco_count <- cluster_busco %>% count(BUSCO_ID) 

library(janitor)
freq <- tabyl(busco_count$n, sort = TRUE)
freq[length(freq$`busco_count$n`)+1,] <- c(">5",sum(freq$n[6:14]),sum(freq$percent[6:14]))
freq <- freq[c(-6:-14),]
freq$percent <- as.numeric(freq$percent)

#write.csv(freq,"freq.csv",row.names = F)

# plot
#png("plot_cluster_per.png", units = "in", width = 3, height = 2, res=400)
ggplot(data = freq, aes(x = `busco_count$n`,y=percent)) +
  geom_bar(stat = "identity",fill = "skyblue")+
  #scale_x_discrete(limits = c("1","2","3","4","5",">5"))+
  scale_x_discrete(limits = c(">5","5","4","3","2","1"),position = "bottom")+
  xlab("")+
  #ylab("Percentage")+
  ylab("")+
  ggtitle("# of synteny clusters")+
  coord_flip()+
  #scale_y_reverse(expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5))
#dev.off()
