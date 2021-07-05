# Title: PCA plots
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date last modified: 7/5/21

library("DESeq2")
library("tidyverse")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("ggplot2")
library("gplots")
library("limma")
library("spdep") 
library("adegenet") 
library("goseq")
library("gridExtra")
library("clusterProfiler")
library(stringr)
library("heatmaply")

# Acerv PCA

acerv <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv/acerv_sub_unique.sig.list_20210219.csv")


test <- estimateSizeFactors(acerv)


SFtest <- estimateSizeFactors(acerv)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(acerv, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# PCA plot of diff-expressed genes 
acerv_sub_DEG_PCA <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_acerv_sub <- round(100*attr(acerv_sub_DEG_PCA, "percentVar")) #plot PCA of samples with all data
acerv_sub_DEG_PCA_plot <- ggplot(acerv_sub_DEG_PCA, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=8) +
  #geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_acerv_sub[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv_sub[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  #scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  scale_color_manual(values = c(control="gray", Treatment1="darkslategray2", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  #ggtitle("A. cervicornis") +
  theme_bw() + #Set background color
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size=25),
        #title = element_text(size=30),
        legend.position = "none",
        panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_sub_DEG_PCA_plot
# PCA plot is of differentially expressed genes only
#PC.info <- mcav_DEGPCAplot$data
ggsave("~/Desktop/acerv_sub_DEGs_PCA_20210219.jpeg", acerv_sub_DEG_PCA_plot, width = 30, height = 20,, units = "cm")
