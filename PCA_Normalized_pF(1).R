#####R version: 3.6.2#####
library(ggplot2)
library(vegan)
library(grid)
library(gridExtra)

geneName = "GAPDH"
#geneName = "ACTB"
#geneName = "EDF1"
#geneName = "GSTP1"
#geneName = "OAZ1"
#geneName = "UBB"
#geneName = "YBX1"


##### total pF 9vs9
pF.9vs9 <- read.csv("../../../1_Paper_Data/Total_All_pF(our)/18_pFs_TPM_SVA_All_genes.csv", header=T, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
sample.table1 <- read.csv("../../../RNA_Total_RNASeq/2017-12-All_data_v19/EC-DC-2542_metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
sample.table1 <- sample.table1[order(sample.table1$group, decreasing = TRUE),]
sample.table1$Batch <- "Old"
sample.table2 <- read.csv("../../../RNA_Total_RNASeq/2017-12-All_data_v19/EC-DC-4731_metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
sample.table2 <- sample.table2[order(sample.table2$group, decreasing = TRUE),]
sample.table2 <- sample.table2[rownames(sample.table2)!="pF-01"&rownames(sample.table2)!="pF-08",]
sample.table2$Batch <- "New"
## Combind two corhorts
sample.table <- rbind(sample.table1[,c(3,6)], sample.table2[,c(2,3)])
sample.table <- sample.table[order(sample.table$group, decreasing = TRUE),]
pF.9vs9 <- pF.9vs9[,match(rownames(sample.table), colnames(pF.9vs9))]
pF.9vs9 <- apply(pF.9vs9, 2, function(x) x/x[names(x)==geneName])
sample.table[sample.table$group=="Y",]$group <- "Young"
sample.table[sample.table$group=="O",]$group <- "Old"
sample.table$cohort <- "Total pF 9vs9"
sample.pF9vs9 <- sample.table[,-2]
remove(sample.table1, sample.table2, sample.table)

##### polyA+ pF 4vs4
sample.table <- read.csv("../../../1_Paper_Data/RawCount/EC-JM-1753_metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
PolyA_pF.4vs4 <- read.csv("../../../1_Paper_Data/TPM/EC-JM-1753_Gencode_TPM.csv", header=T, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
sample.table <- sample.table[(sample.table$group=="Y"|sample.table$group=="O"),]
sample.table <- sample.table[sample.table$Type=="PF",]
sample.table <- sample.table[order(rownames(sample.table)),]
sample.table <- sample.table[order(sample.table$group, decreasing = TRUE),]
PolyA_pF.4vs4 <- PolyA_pF.4vs4[,match(rownames(sample.table), colnames(PolyA_pF.4vs4))]
PolyA_pF.4vs4 <- apply(PolyA_pF.4vs4, 2, function(x) x/x[names(x)==geneName])
sample.table[sample.table$group=="Y",]$group <- "Young"
sample.table[sample.table$group=="O",]$group <- "Old"
sample.table$cohort <- "polyA+ pF 4vs4"
sample.PolyA_pF4vs4 <- sample.table[, c(3,6)]
remove(sample.table)

##### polyA+ Gage pF
sample.table <- read.csv("../../../1_Paper_Data/RawCount/Gage_pF_metadata.csv", header=T, stringsAsFactors=FALSE,row.names=1)
Gage_pF.15vs70 <- read.csv("../../../1_Paper_Data/TPM/Gage-pF_Gencode_TPM.csv", header=T, stringsAsFactors=FALSE, row.names=1, check.names=FALSE)
sample.table[sample.table$group=="Y",]$group<-"M"
sample.table[sample.table$Age<=15,]$group<-"Y"
sample.table <- sample.table[(sample.table$group=="Y"|sample.table$group=="O"),]
sample.table <- sample.table[order(sample.table$group, decreasing = TRUE),]
Gage_pF.15vs70 <- Gage_pF.15vs70[,match(rownames(sample.table), colnames(Gage_pF.15vs70))]
Gage_pF.15vs70 <- apply(Gage_pF.15vs70, 2, function(x) x/x[names(x)==geneName])
sample.table[sample.table$group=="Y",]$group <- "Young"
sample.table[sample.table$group=="O",]$group <- "Old"
sample.table$cohort <- "Gage's pF 6vs4"
sample.Gage_pF15vs70 <- sample.table[, c(1,5)]
remove(sample.table)

combineData <- function(count.matrix, sample.table){
  grid.arrange(plotPCA("L1_common_all", count.matrix, sample.table), 
               plotPCA("L2_common_all_pF", count.matrix, sample.table), 
               plotPCA("L6_pan_Brain", count.matrix, sample.table), 
               plotPCA("L7_FC", count.matrix, sample.table),
               nrow = 2)
}

plotPCA <- function(markersName, count.matrix, sample.table){
  Block.markers <- read.csv(paste("../../Figure1/Figure1C/", markersName,"_pvals.csv",sep=""), header=T, stringsAsFactors=FALSE, row.names=1)
  Block.markers <- Block.markers[rownames(Block.markers) %in% rownames(pF.9vs9),]
  Block.markers <- Block.markers[rownames(Block.markers) %in% rownames(PolyA_pF.4vs4),]
  Block.markers <- Block.markers[rownames(Block.markers) %in% rownames(Gage_pF.15vs70),]
  Block.markers <- Block.markers[order(Block.markers$edgington),]
  Block.markers <- Block.markers[1:100,]
  count.matrix <- cbind(pF.9vs9[match(rownames(Block.markers), rownames(pF.9vs9)),], 
                        PolyA_pF.4vs4[match(rownames(Block.markers), rownames(PolyA_pF.4vs4)),], 
                        Gage_pF.15vs70[match(rownames(Block.markers), rownames(Gage_pF.15vs70)),])
  sample.table <- rbind(sample.pF9vs9, sample.PolyA_pF4vs4, sample.Gage_pF15vs70)
  
  pc <- prcomp(t(count.matrix))
  percentVar <- pc$sdev^2 / sum( pc$sdev^2)
  
  PCi<-data.frame(pc$x,Group=sample.table$group, Cohort=sample.table$cohort)
  
  g1 <- ggplot(PCi,aes(x=PC2,y=PC3,col=Group, shape=Cohort))+
    xlab(paste0("PC2 (", round(percentVar[2]*100), "% variance)")) + 
    ylab(paste0("PC3 (", round(percentVar[3]*100), "% variance)")) + 
    geom_point(size=4,alpha=0.7)+ #Size and alpha just for fun
    scale_color_manual(values = c("#DC143C", "#1E90FF"))+ 
    theme_classic() + ggtitle(markersName)
  return(g1)
}

pdf(paste(gsub("-","",Sys.Date()),"_", geneName, "_pF_PC_2_3.pdf",sep=""), width=10.5, height=8, useDingbats = FALSE)
  combineData(count.matrix, sample.table)
dev.off()


