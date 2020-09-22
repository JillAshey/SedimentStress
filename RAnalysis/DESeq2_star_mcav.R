# Title: DESeq2 with Mcav samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/16/20

# Code for Francois sedimentation data. M. cav only samples analyzed here aligned against M. cav. STAR was read aligner with gff annotation file from Matz lab (pers. comm.)

# Load packages
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

# Load gene count matrix
mcav_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_mcav_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(mcav_counts) # 25142 x 15
head(mcav_counts)
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(mcav_counts)[col])
}
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  gsub("X", "", colnames(mcav_counts)[col])
}
# In Hollie's code, she read a functional annotation file in here 
annot <- read.csv("~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3",header = FALSE, sep="\t", skip=2)
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
# annot$gene <- annot$attr
annot <- annot[!grepl("##", annot$scaffold),]
annot$gene <-gsub(";.*", "", annot$attr)
annot$gene <-gsub("ID=", "", annot$gene)














































# Load gene count matrix
mcav_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_mcav_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(mcav_counts) # 25142 x 15
head(mcav_counts)
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(mcav_counts)[col])
}
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  gsub("X", "", colnames(mcav_counts)[col])
}

# Load metadata
metadata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata.csv", header = TRUE)
dim(metadata) # 45 by 12
head(metadata)
# Selecting only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Acerv species only 
mcav_metadata <- subset(metadata, Species=="Montastraea cavernosa")

# Renaming treatments
mcav_metadata$Treatment <- gsub("Ctl", "control", mcav_metadata$Treatment)
mcav_metadata$Treatment <- gsub("T1", "Treatment1", mcav_metadata$Treatment)
mcav_metadata$Treatment <- gsub("T2", "Treatment2", mcav_metadata$Treatment)
mcav_metadata$Treatment <- gsub("T3", "Treatment3", mcav_metadata$Treatment)
mcav_metadata$Treatment <- gsub("T4", "Treatment4", mcav_metadata$Treatment)
# Removing unwanted text from SampleID
mcav_metadata$SampleID <- gsub(".txt.gz", "", mcav_metadata$SampleID)
mcav_metadata$SampleID <- gsub(";.*", "", mcav_metadata$SampleID)
mcav_metadata$SampleID <- gsub(".fastq.gz", "", mcav_metadata$SampleID)
mcav_metadata$SampleID <- sub("\\.", "", mcav_metadata$SampleID)

# Making sampleID as rownames in metadata 
rownames(mcav_metadata) <- mcav_metadata$SampleID

# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(mcav_metadata) %in% colnames(mcav_counts)) # must come out TRUE



## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(mcav_counts, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- mcav_counts[gfilt,]  
dim(gkeep) # 12873 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
mcav_gcount_filt <- as.data.frame(mcav_counts[which(rownames(mcav_counts) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(mcav_gcount_filt)
dim(mcav_gcount_filt) # 12873 x 15 -- only 12873 genes kept after filtering 

# Write acerv metadata and counts tables with corrected column and row names and filtered gene counts
write.csv(mcav_metadata, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_mcav_filtered.csv")
write.csv(mcav_gcount_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_mcav_star_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(mcav_metadata) %in% colnames(mcav_gcount_filt))




## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
mcav_metadata$Treatment <- factor(mcav_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
head(mcav_metadata)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_mcav <- DESeqDataSetFromMatrix(countData = mcav_gcount_filt,
                                     colData = mcav_metadata,
                                     design = ~Treatment)
gdds_mcav





## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_mcav <- estimateSizeFactors(gdds_mcav) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_mcav
print(sizeFactors(SF.gdds_mcav)) #view size factors

# size factors all less than 4, can use VST
# vst not working, using rlog 
gvst_mcav <- vst(gdds_mcav, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_mcav))
dim(gvst_mcav)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_mcav))) # calculate distance matrix, t returns transpose of assay(gvst_ofav)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_mcav) # assign row names 
colnames(gsampleDistsMatrix) <- NULL # assign col names 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
mcav_heatmap <- pheatmap(gsampleDistsMatrix, # plot matrix
                          clustering_distance_rows = gsampleDists, # cluster rows
                          clustering_distance_cols = gsampleDists, # cluster cols
                          col = colors) # set colors 
## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_mcav, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
mcav_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  #geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() + 
  ggtitle("M. cav (star)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
mcav_heatmap_PCA <- grid.arrange(mcav_PCAplot, mcav_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/mcav_star_heatmap_PCA.pdf", mcav_heatmap_PCA, width = 8, height = 8, units = c("in"))






## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_mcav <- DESeq(gdds_mcav) #run differential expression test by group using the Wald model
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
res_DEG_mcav <- results(DEG_mcav)
res_DEG_mcav_Ordered <- res_DEG_mcav[order(res_DEG_mcav$pvalue),]
DEG_mcav$Treatment
resultsNames(DEG_mcav)
# [1] "Intercept"                       "Treatment_Treatment1_vs_control"
# [3] "Treatment_Treatment2_vs_control" "Treatment_Treatment3_vs_control"
# [5] "Treatment_Treatment4_vs_control"

# Explore significant p-values for treatments 
# Control vs treatment1
DEG_mcav_results_control_vs_T1 <- results(DEG_mcav, name = "Treatment_Treatment1_vs_control") # results only for control vs treatment1 treatments 
results_ordered_DEG_mcav_results_control_vs_T1 <- DEG_mcav_results_control_vs_T1[order(DEG_mcav_results_control_vs_T1$pvalue),] # order from smallest pvalue 
summary(DEG_mcav_results_control_vs_T1) # view summary of results with adj p < 0.1
mcav_sig.num.control_vs_T1 <- sum(DEG_mcav_results_control_vs_T1$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 19 significantly differentially expressed genes between control and T1 that are less than 0.05
mcav_DEGs.control_vs_T1 <- subset(DEG_mcav_results_control_vs_T1, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.control_vs_T1 <- as.data.frame(mcav_DEGs.control_vs_T1) # make df
mcav_DEGs.control_vs_T1$contrast <- as_factor(c("Treatment_Treatment1_vs_control")) # set contrast as a factor 
mcav_DEGs.control_vs_T1 <- cbind(gene_id = rownames(mcav_DEGs.control_vs_T1), mcav_DEGs.control_vs_T1) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.control_vs_T1) <- NULL # remove row names 
mcav_DEGs.control_vs_T1
dim(mcav_DEGs.control_vs_T1) # 19 by 8
write.csv(mcav_DEGs.control_vs_T1, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.control_vs_T1.csv")
mcav_sig.num.control_vs_T1

# Control vs treatment2
DEG_mcav_results_control_vs_T2 <- results(DEG_mcav, name = "Treatment_Treatment2_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_mcav_results_control_vs_T2 <- DEG_mcav_results_control_vs_T2[order(DEG_mcav_results_control_vs_T2$pvalue),] # order from smallest pvalue 
summary(DEG_mcav_results_control_vs_T2) # view summary of results with adj p < 0.1
mcav_sig.num.control_vs_T2 <- sum(DEG_mcav_results_control_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 43 significantly differentially expressed genes between control and T2 that are less than 0.05
mcav_DEGs.control_vs_T2 <- subset(DEG_mcav_results_control_vs_T2, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.control_vs_T2 <- as.data.frame(mcav_DEGs.control_vs_T2) # make df
mcav_DEGs.control_vs_T2$contrast <- as_factor(c("Treatment_Treatment2_vs_control")) # set contrast as a factor 
mcav_DEGs.control_vs_T2 <- cbind(gene_id = rownames(mcav_DEGs.control_vs_T2), mcav_DEGs.control_vs_T2) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.control_vs_T2) <- NULL # remove row names 
mcav_DEGs.control_vs_T2
dim(mcav_DEGs.control_vs_T2) # 43 by 8
write.csv(mcav_DEGs.control_vs_T2, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.control_vs_T2.csv")
mcav_sig.num.control_vs_T2

# Control vs treatment3
DEG_mcav_results_control_vs_T3 <- results(DEG_mcav, name = "Treatment_Treatment3_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_mcav_results_control_vs_T3 <- DEG_mcav_results_control_vs_T3[order(DEG_mcav_results_control_vs_T3$pvalue),] # order from smallest pvalue 
summary(DEG_mcav_results_control_vs_T3) # view summary of results with adj p < 0.1
mcav_sig.num.control_vs_T3 <- sum(DEG_mcav_results_control_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 26 significantly differentially expressed genes between control and T3 that are less than 0.05
mcav_DEGs.control_vs_T3 <- subset(DEG_mcav_results_control_vs_T3, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.control_vs_T3 <- as.data.frame(mcav_DEGs.control_vs_T3) # make df
mcav_DEGs.control_vs_T3$contrast <- as_factor(c("Treatment_Treatment3_vs_control")) # set contrast as a factor 
mcav_DEGs.control_vs_T3 <- cbind(gene_id = rownames(mcav_DEGs.control_vs_T3), mcav_DEGs.control_vs_T3) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.control_vs_T3) <- NULL # remove row names 
mcav_DEGs.control_vs_T3
dim(mcav_DEGs.control_vs_T3) # 26 by 8
write.csv(mcav_DEGs.control_vs_T3, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.control_vs_T3.csv")
mcav_sig.num.control_vs_T3

# Control vs treatment4
DEG_mcav_results_control_vs_T4 <- results(DEG_mcav, name = "Treatment_Treatment4_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_mcav_results_control_vs_T4 <- DEG_mcav_results_control_vs_T4[order(DEG_mcav_results_control_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_mcav_results_control_vs_T4) # view summary of results with adj p < 0.1
mcav_sig.num.control_vs_T4 <- sum(DEG_mcav_results_control_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 18 significantly differentially expressed genes between control and T4 that are less than 0.05
mcav_DEGs.control_vs_T4 <- subset(DEG_mcav_results_control_vs_T4, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.control_vs_T4 <- as.data.frame(mcav_DEGs.control_vs_T4) # make df
mcav_DEGs.control_vs_T4$contrast <- as_factor(c("Treatment_Treatment4_vs_control")) # set contrast as a factor 
mcav_DEGs.control_vs_T4 <- cbind(gene_id = rownames(mcav_DEGs.control_vs_T4), mcav_DEGs.control_vs_T4) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.control_vs_T4) <- NULL # remove row names 
mcav_DEGs.control_vs_T4
dim(mcav_DEGs.control_vs_T4) # 26 by 8
write.csv(mcav_DEGs.control_vs_T4, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.control_vs_T4.csv")
mcav_sig.num.control_vs_T4

# Treatment1 vs treatment2
DEG_mcav_results_T1_vs_T2 <- results(DEG_mcav, contrast = c("Treatment", "Treatment1", "Treatment2")) # results only for control vs treatment2 treatments 
results_ordered_DEG_mcav_results_T1_vs_T2 <- DEG_mcav_results_T1_vs_T2[order(DEG_mcav_results_T1_vs_T2$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T1_vs_T2) # view summary of results with adj p < 0.1
mcav_sig.num.T1_vs_T2 <- sum(DEG_mcav_results_T1_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between control and T2 that are less than 0.05

# Treatment1 vs treatment3
DEG_mcav_results_T1_vs_T3 <- results(DEG_mcav, contrast = c("Treatment", "Treatment1", "Treatment3")) # results only for T1 vs T3 treatments 
results_ordered_DEG_mcav_results_T1_vs_T3 <- DEG_mcav_results_T1_vs_T3[order(DEG_mcav_results_T1_vs_T3$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T1_vs_T3) # view summary of results with adj p < 0.1
mcav_sig.num.T1_vs_T3 <- sum(DEG_mcav_results_T1_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 1 significantly differentially expressed genes between T1 and T3 that are less than 0.05
mcav_DEGs.T1_vs_T3 <- subset(DEG_mcav_results_T1_vs_T3, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.T1_vs_T3 <- as.data.frame(mcav_DEGs.T1_vs_T3) # make df
mcav_DEGs.T1_vs_T3$contrast <- as_factor(c("Treatment_Treatment1_vs_Treatment3")) # set contrast as a factor 
mcav_DEGs.T1_vs_T3 <- cbind(gene_id = rownames(mcav_DEGs.T1_vs_T3), mcav_DEGs.T1_vs_T3) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.T1_vs_T3) <- NULL # remove row names 
mcav_DEGs.T1_vs_T3
dim(mcav_DEGs.T1_vs_T3) # 1 by 8
write.csv(mcav_DEGs.T1_vs_T3, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.T1_vs_T3.csv")
mcav_sig.num.T1_vs_T3

# Treatment1 vs treatment4
DEG_mcav_results_T1_vs_T4 <- results(DEG_mcav, contrast = c("Treatment", "Treatment1", "Treatment4")) # results only for T1 vs T4 treatments 
results_ordered_DEG_mcav_results_T1_vs_T4 <- DEG_mcav_results_T1_vs_T4[order(DEG_mcav_results_T1_vs_T4$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T1_vs_T4) # view summary of results with adj p < 0.1
mcav_sig.num.T1_vs_T4 <- sum(DEG_mcav_results_T1_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between T1 and T4 that are less than 0.05

# Treatment2 vs treatment3
DEG_mcav_results_T2_vs_T3 <- results(DEG_mcav, contrast = c("Treatment", "Treatment2", "Treatment3")) # results only for T1 vs T3 treatments 
results_ordered_DEG_mcav_results_T2_vs_T3 <- DEG_mcav_results_T2_vs_T3[order(DEG_mcav_results_T2_vs_T3$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T2_vs_T3) # view summary of results with adj p < 0.1
mcav_sig.num.T2_vs_T3 <- sum(DEG_mcav_results_T2_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 1 significantly differentially expressed genes between T2 and T3 that are less than 0.05
mcav_DEGs.T2_vs_T3 <- subset(DEG_mcav_results_T2_vs_T3, padj<0.05) # subset only <0.05 padj values
mcav_DEGs.T2_vs_T3 <- as.data.frame(mcav_DEGs.T2_vs_T3) # make df
mcav_DEGs.T2_vs_T3$contrast <- as_factor(c("Treatment_Treatment2_vs_Treatment3")) # set contrast as a factor 
mcav_DEGs.T2_vs_T3 <- cbind(gene_id = rownames(mcav_DEGs.T2_vs_T3), mcav_DEGs.T2_vs_T3) # make gene id a row and bind it to the rest of the df
rownames(mcav_DEGs.T2_vs_T3) <- NULL # remove row names 
mcav_DEGs.T2_vs_T3
dim(mcav_DEGs.T2_vs_T3) # 1 by 8
write.csv(mcav_DEGs.T2_vs_T3, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEGs.T2_vs_T3.csv")
mcav_sig.num.T2_vs_T3

# Treatment2 vs treatment4
DEG_mcav_results_T2_vs_T4 <- results(DEG_mcav, contrast = c("Treatment", "Treatment2", "Treatment4")) # results only for T2 vs T4 treatments 
results_ordered_DEG_mcav_results_T2_vs_T4 <- DEG_mcav_results_T2_vs_T4[order(DEG_mcav_results_T2_vs_T4$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T2_vs_T4) # view summary of results with adj p < 0.1
mcav_sig.num.T2_vs_T4 <- sum(DEG_mcav_results_T2_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between T2 and T3 that are less than 0.05

# Treatment3 vs treatment4
DEG_mcav_results_T3_vs_T4 <- results(DEG_mcav, contrast = c("Treatment", "Treatment3", "Treatment4")) # results only for T2 vs T4 treatments 
results_ordered_DEG_mcav_results_T3_vs_T4 <- DEG_mcav_results_T3_vs_T4[order(DEG_mcav_results_T3_vs_T4$padj),] # order from smallest pvalue 
summary(DEG_mcav_results_T3_vs_T4) # view summary of results with adj p < 0.1
mcav_sig.num.T3_vs_T4 <- sum(DEG_mcav_results_T3_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between T3 and T3 that are less than 0.05

# Combining all DEGs among all treatment comparison
DEGs_mcav_all_treatments <- bind_rows(mcav_DEGs.control_vs_T1, mcav_DEGs.control_vs_T2, mcav_DEGs.control_vs_T3, mcav_DEGs.control_vs_T4, mcav_DEGs.T1_vs_T3, mcav_DEGs.T2_vs_T3) # bind CvsT1, CvsT2, CvsT3, CvsT4, T1vsT3, T2vsT3 results together by row - comparisons where there were DEGs
dim(DEGs_mcav_all_treatments) # 108 by 8
DEG_mcav.sg.num <- sum(DEGs_mcav_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_mcav.sg.num # 89 adjusted p-values 
summary(DEGs_mcav_all_treatments)
write.csv(DEGs_mcav_all_treatments, file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/DEGs_mcav_all_treatments.csv")





## Visualize diff-expressed genes

# Subset and log-transform count data 
# Subset list of genes by those which padj>0.
dim(DEGs_mcav_all_treatments)
DEGs_mcav <- DEGs_mcav_all_treatments$gene_id # list all gene names 
DEGs_mcav <- unique(DEGs_mcav) # select only unique gene names 
DEG_mcav_list <- gdds_mcav[which(rownames(gdds_mcav) %in% DEGs_mcav)] # filter gdds_acerv DESeq2 object by unique gene names
dim(DEG_mcav_list) # 62 x 15
print(counts(DEG_mcav_list))

# As determined above, size factors all less than 4, so proceed with VST
#apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_mcav_list, blind = FALSE, nsub = nrow(counts(DEG_mcav_list)))
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
dim(DEGvst) # 62 x 15
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
# Testing if the first two command in the heatmap is necessary given that we are clustering the columns anyways.
mcav_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_mcav.sg.num) #sort by decreasing sig ?
mat_mcav <- assay(DEGvst)[mcav_topVarGenes, ] #make an expression object
mat_mcav <- mat_mcav - rowMeans(mat_mcav) #diff_mcav in expression compared to average across all samples
dim(mat_mcav)
ann_colors <- list(Treatment= c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange"))
df_DEG_mcav <- as.data.frame(colData(DEGvst)[c("Treatment")]) #make dataframe for column naming and associated treatment
mcav_DEGheatmap <- pheatmap(mat_mcav, scale= "row", legend=TRUE, annotation_legend=TRUE, annotation_col=df_DEG_mcav, annotation_colors = ann_colors,
                             clustering_distance_rows="euclidean", clustering_method = "average",
                             show_rownames =FALSE,
                             show_colnames =TRUE,
                             cluster_cols = TRUE)
mcav_DEGheatmap 

# PCA plot of diff-expressed genes 
mcav_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_mcav <- round(100*attr(mcav_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
mcav_DEGPCAplot <- ggplot(mcav_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_mcav[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_mcav[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
mcav_DEGPCAplot
# PCA plot is of differentially expressed genes only

# Save results
write.csv(counts(DEG_mcav_list), file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/mcav_DEG_list_unique.csv")
mcav_DEGs_heatmap_PCA <- grid.arrange(mcav_DEGPCAplot, mcav_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/mcav_DEGs_heatmap_PCA.pdf", mcav_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))

## Can I make DEG PCAs with individual contrasts? eg control vs T1, T1 vs T4, etc 




