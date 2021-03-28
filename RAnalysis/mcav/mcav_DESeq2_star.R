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
mcav_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/mcav/mcav_gene_count_matrix.csv", header = TRUE, row.names = "gene_id")
dim(mcav_counts) # 25142 x 15
head(mcav_counts)
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(mcav_counts)[col])
}
for ( col in 1:ncol(mcav_counts)){
  colnames(mcav_counts)[col] <-  gsub("X", "", colnames(mcav_counts)[col])
}

# functional annotation gff
# annot <- read.csv("~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3",header = FALSE, sep="\t")
# colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
# # annot$gene <- annot$attr
# annot <- annot[!grepl("##", annot$scaffold),]
# annot$gene <-gsub(";.*", "", annot$attr)
# annot$gene <-gsub("ID=", "", annot$gene)
# annot$gene <- gsub("-.*", "", annot$gene)

# Load metadata
metadata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/FL_sediment_metadata.csv", header = TRUE)
dim(metadata) # 45 by 15
head(metadata)
# Selecting only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Mcav species only 
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

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(mcav_counts, filt) # create filter for counts data 
keep <- mcav_counts[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
mcav_counts_filt <- as.matrix(mcav_counts[which(rownames(mcav_counts) %in% gn.keep),]) 
write.csv(mcav_counts_filt, "~/Desktop/mcav_counts_filtered.csv")
storage.mode(mcav_counts_filt) <- "integer" # stores count data as integer 
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(mcav_metadata) %in% colnames(mcav_counts_filt)) # must come out TRUE
# Set Treatment as a factor
mcav_metadata$Treatment <- factor(mcav_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = mcav_counts_filt, colData = mcav_metadata, design = ~ Treatment)

# Expression visualization
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.data <- estimateSizeFactors(data) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.data
print(sizeFactors(SF.data)) #view size factors

vst <- vst(data, blind = FALSE) 

head(assay(vst), 3) # view data
sampleDists <- dist(t(assay(vst))) # calculate distance matrix
sampleDistMatrix <- as.matrix(sampleDists) # create distance matrix
rownames(sampleDistMatrix) <- colnames(vst) # assign row names
colnames(sampleDistMatrix) <- NULL # assign col names 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) # assign colors 
pheatmap(sampleDistMatrix, # plot matrix
         clustering_distance_rows = sampleDists, # cluster rows
         clustering_distance_cols = sampleDists, # cluster cols
         col=colors) # set colors
plotPCA(vst, intgroup = c("Treatment")) # plot PCA of samples with all data 
# Treatment2 has a weird outlier, may remove?

# Differential gene expression analysis 
DEG.int <- DESeq(data) # run differential expression test by treatment (?) using wald test 
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
DEG.int.res <- results(DEG.int) # save DE results ; why does it say 'Wald test p-value: Treatment Treatment4 vs control' for DEG.int.res? Is it only looking at treatment 4 and control? In DESeq object created above, it says that design is Treatment
resultsNames(DEG.int) # view DE results 
#[1] "Intercept"                       "Treatment_Treatment1_vs_control" "Treatment_Treatment2_vs_control"
#[4] "Treatment_Treatment3_vs_control" "Treatment_Treatment4_vs_control"
# NAs in padj: https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na

# Compare C vs T1
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1")) # results of DESeq2 comparing C and T1
DEG_control_vs_T1 
DEG_control_vs_T1 <- as.data.frame(DEG_control_vs_T1) # make results into a df
DEG_control_vs_T1["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
write.csv(DEG_control_vs_T1, file = "~/Desktop/mcav_control_vs_T1_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 19 DEGs
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T1.sig.list_full <- cbind(DEG_control_vs_T1.sig, DEG_control_vs_T1.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T1.sig.list_full, file = "~/Desktop/mcav_control_vs_T1_DEG_full.csv") # write out csv
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2")) # results of DESeq2 comparing C and T1
DEG_control_vs_T2 
DEG_control_vs_T2 <- as.data.frame(DEG_control_vs_T2) # make results into a df
DEG_control_vs_T2["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
write.csv(DEG_control_vs_T2, file = "~/Desktop/mcav_control_vs_T2_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T2.sig.num
# 43 DEGs
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T2.sig.list <- as.data.frame(counts(DEG_control_vs_T2.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T2.sig.list_full <- cbind(DEG_control_vs_T2.sig, DEG_control_vs_T2.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T2.sig.list_full, file = "~/Desktop/mcav_control_vs_T2_DEG_full.csv") # write out csv
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3")) # results of DESeq2 comparing C and T1
DEG_control_vs_T3
DEG_control_vs_T3 <- as.data.frame(DEG_control_vs_T3) # make results into a df
DEG_control_vs_T3["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
write.csv(DEG_control_vs_T3, file = "~/Desktop/mcav_control_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T3.sig.num
# 26 DEGs
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T3.sig.list <- as.data.frame(counts(DEG_control_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T3.sig.list_full <- cbind(DEG_control_vs_T3.sig, DEG_control_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T3.sig.list_full, file = "~/Desktop/mcav_control_vs_T3_DEG_full.csv") # write out csv
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4")) # results of DESeq2 comparing C and T1
DEG_control_vs_T4
DEG_control_vs_T4 <- as.data.frame(DEG_control_vs_T4) # make results into a df
DEG_control_vs_T4["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
write.csv(DEG_control_vs_T4, file = "~/Desktop/mcav_control_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T4.sig.num
# 18 DEGs
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T4.sig.list_full <- cbind(DEG_control_vs_T4.sig, DEG_control_vs_T4.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T4.sig.list_full, file = "~/Desktop/mcav_control_vs_T4_DEG_full.csv") # write out csv
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2")) # results of DESeq2 comparing C and T1
DEG_T1_vs_T2
DEG_T1_vs_T2 <- as.data.frame(DEG_T1_vs_T2) # make results into a df
DEG_T1_vs_T2["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
write.csv(DEG_T1_vs_T2, file = "~/Desktop/mcav_T1_vs_T2_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T2.sig.num
# 0 DEGs

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3")) # results of DESeq2 comparing C and T1
DEG_T1_vs_T3
DEG_T1_vs_T3 <- as.data.frame(DEG_T1_vs_T3) # make results into a df
DEG_T1_vs_T3["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
write.csv(DEG_T1_vs_T3, file = "~/Desktop/mcav_T1_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T3.sig.num
# 1 DEGs
DEG_T1_vs_T3.sig <- subset(DEG_T1_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T3.sig["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
DEG_T1_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_T1_vs_T3.sig.list <- as.data.frame(counts(DEG_T1_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_T1_vs_T3.sig.list_full <- cbind(DEG_T1_vs_T3.sig, DEG_T1_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T1_vs_T3.sig.list_full, file = "~/Desktop/mcav_T1_vs_T3_DEG_full.csv") # write out csv
DEG_T1_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4")) # results of DESeq2 comparing C and T1
DEG_T1_vs_T4
DEG_T1_vs_T4 <- as.data.frame(DEG_T1_vs_T4) # make results into a df
DEG_T1_vs_T4["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
write.csv(DEG_T1_vs_T4, file = "~/Desktop/mcav_T1_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T4.sig.num
# 0 DEGs

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3")) # results of DESeq2 comparing C and T1
DEG_T2_vs_T3
DEG_T2_vs_T3 <- as.data.frame(DEG_T2_vs_T3) # make results into a df
DEG_T2_vs_T3["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
write.csv(DEG_T2_vs_T3, file = "~/Desktop/mcav_T2_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T2_vs_T3.sig.num
# 1 DEG
DEG_T2_vs_T3.sig <- subset(DEG_T2_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T2_vs_T3.sig["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
DEG_T2_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T2_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_T2_vs_T3.sig.list <- as.data.frame(counts(DEG_T2_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_T2_vs_T3.sig.list_full <- cbind(DEG_T2_vs_T3.sig, DEG_T2_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T2_vs_T3.sig.list_full, file = "~/Desktop/mcav_T2_vs_T3_DEG_full.csv") # write out csv
DEG_T2_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_T3_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4")) # results of DESeq2 comparing C and T1
DEG_T2_vs_T4
DEG_T2_vs_T4 <- as.data.frame(DEG_T2_vs_T4) # make results into a df
DEG_T2_vs_T4["Treatment_Compare"] <- "T2vsT4" # adding treatment comparison column
write.csv(DEG_T2_vs_T4, file = "~/Desktop/mcav_T2_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T2_vs_T4.sig.num
# 0 DEG

# Compare T2 vs T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4")) # results of DESeq2 comparing C and T1
DEG_T3_vs_T4
DEG_T3_vs_T4 <- as.data.frame(DEG_T3_vs_T4) # make results into a df
DEG_T3_vs_T4["Treatment_Compare"] <- "T3vsT4" # adding treatment comparison column
write.csv(DEG_T3_vs_T4, file = "~/Desktop/mcav_T3_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T3_vs_T4.sig.num
# 0 DEG


# Make full list of genes and treatments: CvsT1, CvsT2, CvsT3, CvsT4, T1vsT3, T2vsT3
DEG_control_vs_T1.sig.list_full$gene_id <- rownames(DEG_control_vs_T1.sig.list_full)
rownames(DEG_control_vs_T1.sig.list_full) <- NULL
DEG_control_vs_T2.sig.list_full$gene_id <- rownames(DEG_control_vs_T2.sig.list_full)
rownames(DEG_control_vs_T2.sig.list_full) <- NULL
DEG_control_vs_T3.sig.list_full$gene_id <- rownames(DEG_control_vs_T3.sig.list_full)
rownames(DEG_control_vs_T3.sig.list_full) <- NULL
DEG_control_vs_T4.sig.list_full$gene_id <- rownames(DEG_control_vs_T4.sig.list_full)
rownames(DEG_control_vs_T4.sig.list_full) <- NULL
DEG_T1_vs_T3.sig.list_full$gene_id <- rownames(DEG_T1_vs_T3.sig.list_full)
rownames(DEG_T1_vs_T3.sig.list_full) <- NULL
DEG_T2_vs_T3.sig.list_full$gene_id <- rownames(DEG_T2_vs_T3.sig.list_full)
rownames(DEG_T2_vs_T3.sig.list_full) <- NULL

DEGs.all <- rbind(DEG_control_vs_T1.sig.list_full, 
                  DEG_control_vs_T2.sig.list_full,
                  DEG_control_vs_T3.sig.list_full, 
                  DEG_control_vs_T4.sig.list_full,
                  DEG_T1_vs_T3.sig.list_full,
                  DEG_T2_vs_T3.sig.list_full)

dim(DEGs.all) # 108 x 23
length(unique(DEGs.all$gene_id)) # 62 unique genes between all treatments 
write.csv(DEGs.all, file = "~/Desktop/mcav_DEGs.all_treatment_20210208.csv")

## Find intersections and unique results between treatments 
# interactions
int1 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T2.sig.list_full$gene_id)
length(unique(int1)) # 15 DEGs shared between CvT1 and CvT2
int2 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T3.sig.list_full$gene_id)
length(unique(int2)) # 16 DEGs shared between CvT1 and CvT3
int3 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int3)) # 8 DEGs shared between CvT1 and CvT2
int4 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_T1_vs_T3.sig.list_full$gene_id)
length(unique(int4)) # 0 DEGs shared between CvT1 and T1vT3
int5 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int5)) # 0 DEGs shared between CvT1 and T2vT3
int6 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_control_vs_T3.sig.list_full$gene_id)
length(unique(int6)) # 15 DEGs shared between CvT2 and CvT3
int7 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int7)) # 10 DEGs shared between CvT2 and CvT4
int8 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_T1_vs_T3.sig.list_full$gene_id)
length(unique(int8)) # 0 DEGs shared between CvT2 and T1vT3
int9 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int9)) # 1 DEGs shared between CvT2 and T2vT3
int10 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int10)) # 11 DEGs shared between CvT3 and T2vT4
int11 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_T1_vs_T3.sig.list_full$gene_id)
length(unique(int11)) # 0 DEGs shared between CvT3 and T1vT3
int12 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int12)) # 0 DEGs shared between CvT3 and T2vT3
int13 <- intersect(DEG_control_vs_T4.sig.list_full$gene_id, DEG_T1_vs_T3.sig.list_full$gene_id)
length(unique(int13)) # 0 DEGs shared between CvT4 and T1vT3
int14 <- intersect(DEG_control_vs_T4.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int14)) # 0 DEGs shared between CvT4 and T2vT3
int15 <- intersect(DEG_T1_vs_T3.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int15)) # 0 DEGs shared between T1vT3 and T2vT3


##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT3, T2vsT3
DEGs.all_mcav <- DEGs.all$gene_id
DEGs.all_mcav <- unique(DEGs.all_mcav)
DEGs.all_mcav <- as.data.frame(DEGs.all_mcav) 
dim(DEGs.all_mcav) # 62 unique DEGs among treatment comparisons 

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_mcav$DEGs), ] # subset list of sig transcripts from original count data
write.csv(counts(unique.sig.list), file = "~/Desktop/mcav_unique.sig.list_20210208.csv")
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# PCA.plot <- plotPCA(unique.rsig, intgroup = "Treatment") # plot PCA of all samples for DEG only 
# PCA.plot
# PC.info <- PCA.plot$data

# PCA plot of diff-expressed genes 
mcav_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_mcav <- round(100*attr(mcav_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
mcav_DEGPCAplot <- ggplot(mcav_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=8) +
  #geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_mcav[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_mcav[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  #scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  scale_color_manual(values = c(control="gray", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  #ggtitle("M. cavernosa") +
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
mcav_DEGPCAplot
# PCA plot is of differentially expressed genes only
#PC.info <- mcav_DEGPCAplot$data
ggsave("~/Desktop/mcav_DEGs_PCA.png", mcav_DEGPCAplot, width = 30, height = 20,, units = "cm")




## Heatmap of DEGs
# This heatmap is going to be wild. Grouping columns by treatment, putting gene id on left hand side and GO term on right hand side
# Also taking out legend 

# Need a file with counts, gene names, GO IDs, term, and ontology 
# This file is acerv genes with GO terms associated with them. One GO term per line, so multiple acerv gene names sometimes if genes have multiple GO terms
# This file includes all gene names and GO IDs
mcav_sig <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/mcav_GOterms.unique.csv", header = TRUE)
mcav_sig <- select(mcav_sig, -X)
colnames(mcav_sig)[1] <-"gene"
head(mcav_sig)

# Getting col order to order unique counts data
list(mcav_metadata)
list <- mcav_metadata[order(mcav_metadata$Treatment),] # need to order them so it will group by treatment in plot
list(list$SampleID) # look at sample IDs and use that list to make col.order
col.order <- c("22_ctl2_Mc_TWF_1",
               "28_ctl1_Mc_GBM_1",
               "42_ctl3_Mc_MGR_1",
               "20_T12_Mc_PWC",
               "39_T13_Mc_FJE",
               "61_T11_Mc_RAP",
               "29_T23_Mc_PND",
               "34_T22_Mc_SVS",
               "58_T21_Mc_EAH",
               "21_T33_Mc_EOU",
               "49_T31_Mc_SWQ",
               "55_T32_Mc_TWP",
               "33_T43_Mc_RFV",
               "46_T41_Mc_QYH_1",
               "56_T42_Mc_JAW") 
  


























# Now I will order the counts data so the samples will group by treatment
unique.DEG.annot <- as.data.frame(counts(unique.sig.list)) # make df of sig genes with counts and sample IDs
list(colnames(unique.DEG.annot))
unique.DEG.annot2 <- unique.DEG.annot[, col.order]
unique.DEG.annot2$gene <- rownames(unique.DEG.annot2)

# Now we will take the unique.DEG.annot2 and merge it with acerv_sig
# The unique.DEG.annot2 file includes gene names for DEGs and counts data
test_merge <- merge(unique.DEG.annot2, mcav_sig, by = "gene", all.x = TRUE)
# test_merge now holds gene names for DEGs, counts data, and GO.IDs

# Now we need info about term and ontology 
GO_all <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/mcav_GO_ALL.csv", header = TRUE) 
GO_all <- select(GO_all, -X)
colnames(GO_all)[1] <-"GO.ID"
GO_merge <- merge(test_merge, GO_all, by = "GO.ID", all.x = T)
# Great! GO_merge now contains GO IDs, gene names for DEGs, counts data, over and under represented pvalue, numCat, term, and ontology.
# All of the genes are in there (sometimes duplicated because multple GO/term/ontology info per gene) and some don't have any GO/term/ontology info. 
# *** I could also just plot GO_merge, but have duplicate genes for some that have multiple multple GO/term/ontology info per gene. Probably not the best idea tho

# Trying to aggregate based on GO terms. Hopefully this works because I also want the term and ontology to also aggregate but i think they may just go.
# Well maybe I could aggregate multiple times and then bind them? Lets see
agg_GO <- aggregate(GO_merge$GO.ID, list(GO_merge$gene), paste, collapse = ",") # aggregate GO terms 
colnames(agg_GO) <- c("gene", "GO.ID")
agg_term <- aggregate(GO_merge$term, list(GO_merge$gene), paste, collapse = ",") # aggregate term
colnames(agg_term) <- c("gene", "term")
agg_ont <- aggregate(GO_merge$ontology, list(GO_merge$gene), paste, collapse = ",") # aggregate ontology
colnames(agg_ont) <- c("gene", "ontology")
agg_over <- aggregate(GO_merge$over_represented_pvalue, list(GO_merge$gene), paste, collapse = ",")
colnames(agg_over) <- c("gene", "over_represented_pvalue")

# Now I'll merge them all together!
merge_all <- merge(agg_GO, agg_term, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, agg_ont, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, agg_over, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, unique.DEG.annot2, by = "gene", all.x = TRUE)
write.csv(merge_all, file = "~/Desktop/mcav_GO_DEG.csv") # maybe include gene counts too?



## Hooray! Now I have a lovely file with counts, gene names, GO IDs, term, and ontology 
# Now I must put it in the heatmap...........

# First, lets make a matrix of gene counts 
rownames(merge_all) <- merge_all$gene
mat <- select(merge_all, -c("gene", "GO.ID", "term", "ontology", "over_represented_pvalue", "term2"))
mat <- as.matrix(mat)

# Now lets make df of only treatment and sample ID
#df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
#colnames(df) <- "Treatment"
#df <- df[order(df$Treatment),]
#df <- as.data.frame(df)
#colnames(df) <- "Treatment"
df <- select(mcav_metadata, c("Treatment"))
#df <- df[order(df$Treatment),]
# probably just easier to take treatment info straight from metadata file

# Now lets make a df of only gene names 
df_gene <- as.data.frame(merge_all$gene)
colnames(df_gene) <- "DEG"
rownames(df_gene) <- df_gene$DEG

# Some genes have multiple terms, so I am going to select the first term for every gene 
merge_all$term2 <- merge_all$term
merge_all$term <- gsub(",.*", "", merge_all$term)

# Some genes have NAs, so subbing blank for NA to see the actual terms
merge_all[is.na(merge_all$term)] <- " "

#Set colors for treatment
ann_colors <- list(Treatment = c(control="gray", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray"))

## Plot heatmap
mcav_sub_heatmap <- pheatmap(mat, 
                              annotation_col = df,
                              #annotation_row = df_gene,
                              annotation_colors = ann_colors,
                              annotation_legend = F,
                              cluster_rows = F,
                              show_rownames = T,
                              cluster_cols = F,
                              show_colnames = T,
                              scale = "row",
                              fontsize_row = 8,
                              labels_row = merge_all$term)
mcav_sub_heatmap
ggsave("~/Desktop/mcav_heatmap.png", mcav_sub_heatmap, width = 30, height = 20,, units = "cm")












## Going to remove the control from the list and see how DEG PCA looks without it--can look at spread of mid and high 
# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
metadata_mcav_treatment <- subset(mcav_metadata, Treatment=="Treatment1" | Treatment=="Treatment2" | Treatment=="Treatment3" | Treatment=="Treatment4")
mcav_ID_treatment <- metadata_mcav_treatment$SampleID
count_mcav_treatment <- select(mcav_counts, all_of(mcav_ID_treatment))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_mcav_treatment, filt) # create filter for counts data 
keep <- count_mcav_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
mcav_treatment_counts_filt <- as.matrix(count_mcav_treatment[which(rownames(count_mcav_treatment) %in% gn.keep),]) 
storage.mode(mcav_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(mcav_treatment_counts_filt, "~/Desktop/plob_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_mcav_treatment) %in% colnames(mcav_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_mcav_treatment$Treatment <- factor(metadata_mcav_treatment$Treatment, levels = c("Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = mcav_treatment_counts_filt, colData = metadata_mcav_treatment, design = ~ Treatment)

# Expression visualization
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.data <- estimateSizeFactors(data) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.data
print(sizeFactors(SF.data)) #view size factors
# size factors all less than 4, can use VST

vst <- vst(data, blind = FALSE) 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

head(assay(vst), 3) # view data
sampleDists <- dist(t(assay(vst))) # calculate distance matrix
sampleDistMatrix <- as.matrix(sampleDists) # create distance matrix
rownames(sampleDistMatrix) <- colnames(vst) # assign row names
colnames(sampleDistMatrix) <- NULL # assign col names 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) # assign colors 
pheatmap(sampleDistMatrix, # plot matrix
         clustering_distance_rows = sampleDists, # cluster rows
         clustering_distance_cols = sampleDists, # cluster cols
         col=colors) # set colors
plotPCA(vst, intgroup = c("Treatment")) # plot PCA of samples with all data 

# Differential gene expression analysis 
DEG.int <- DESeq(data) # run differential expression test by treatment (?) using wald test 
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
DEG.int.res <- results(DEG.int) # save DE results ; why does it say 'Wald test p-value: Treatment Treatment4 vs control' for DEG.int.res? Is it only looking at treatment 4 and control? In DESeq object created above, it says that design is Treatment
resultsNames(DEG.int) # view DE results 
#[1] "Intercept"          "Treatment_Treatment2_vs_Treatment1" "Treatment_Treatment3_vs_Treatment1" "Treatment_Treatment4_vs_Treatment1"

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2"))
DEG_T1_vs_T2
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T2.sig.num
# 0 DEGs

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3"))
DEG_T1_vs_T3
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T3.sig.num
# 2 DEGs
DEG_T1_vs_T3.sig <- subset(DEG_T1_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T3.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T1_vs_T3.sig.list)
print(sizeFactors(SFtest))
DEG_T1_vs_T3.rsig <- varianceStabilizingTransformation(DEG_T1_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning messages:
#  In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#              Estimated rdf < 1.0; not estimating variance
DEG_T1_vs_T3.sig.list <- as.data.frame(counts(DEG_T1_vs_T3.sig.list))
DEG_T1_vs_T3.sig.list["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
#write.csv(DEG_T1_vs_T3.sig.list, file = "~/Desktop/mcav_T1_vs_T3_DEG.csv")

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num
# 2 DEGs
DEG_T1_vs_T4.sig <- subset(DEG_T1_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T4.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T1_vs_T4.sig.list)
print(sizeFactors(SFtest))
DEG_T1_vs_T4.rsig <- varianceStabilizingTransformation(DEG_T1_vs_T4.sig.list, blind = FALSE) 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning messages:
#  In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#              Estimated rdf < 1.0; not estimating variance
DEG_T1_vs_T4.sig.list <- as.data.frame(counts(DEG_T1_vs_T4.sig.list))
DEG_T1_vs_T4.sig.list["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
#write.csv(DEG_T1_vs_T4.sig.list, file = "~/Desktop/mcav_T1_vs_T4_DEG.csv")

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3"))
DEG_T2_vs_T3
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T3.sig.num
# 1 DEGs
DEG_T2_vs_T3.sig <- subset(DEG_T2_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T2_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T2_vs_T3.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T2_vs_T3.sig.list)
print(sizeFactors(SFtest))
# Kinda close to 4 but not quite...going to stick with varianceStabilizingTransformation
DEG_T2_vs_T3.rsig <- varianceStabilizingTransformation(DEG_T2_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# Error in estimateDispersionsFit(object, fitType, quiet = TRUE) : 
# all gene-wise dispersion estimates are within 2 orders of magnitude
# from the minimum value, and so the standard curve fitting techniques will not work.
# One can instead use the gene-wise estimates as final estimates:
#   dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
# ...then continue with testing using nbinom WaldTest or nbinomLRT
DEG_T2_vs_T3.sig.list <- as.data.frame(counts(DEG_T2_vs_T3.sig.list))
DEG_T2_vs_T3.sig.list["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
#write.csv(DEG_T2_vs_T3.sig.list, file = "~/Desktop/mcav_T2_vs_T3_DEG.csv")

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4"))
DEG_T2_vs_T4
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Comapre T3 vs T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num
# 0 DEGs

##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT2, T1vsT3, T1vsT4, T2vsT3, T2vsT4, T3vsT4
DEGs_T1vsT3 <- as.data.frame(rownames(DEG_T1_vs_T3.sig.list), DEG_T1_vs_T3.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT3) <- "DEGs"
DEGs_T1vsT3$Treatment_Compare <- rownames(DEGs_T1vsT3)
DEGs_T1vsT4 <- as.data.frame(rownames(DEG_T1_vs_T4.sig.list), DEG_T1_vs_T4.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT4) <- "DEGs"
DEGs_T1vsT4$Treatment_Compare <- rownames(DEGs_T1vsT4)
DEGs_T2vsT3 <- as.data.frame(rownames(DEG_T2_vs_T3.sig.list), DEG_T2_vs_T3.sig.list$Treatment_Compare)
colnames(DEGs_T2vsT3) <- "DEGs"
DEGs_T2vsT3$Treatment_Compare <- rownames(DEGs_T2vsT3)

DEGs.all <- rbind(DEGs_T1vsT3,
                  DEGs_T1vsT4,
                  DEGs_T2vsT3) 
#write.csv(DEGs.all, file = "~/Desktop/mcav_DEGs.all_treatment.csv")
DEGs.all_mcav <- DEGs.all$DEGs
DEGs.all_mcav <- unique(DEGs.all_mcav)
DEGs.all_mcav <- as.data.frame(DEGs.all_mcav)

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_mcav$DEGs), ] # subset list of sig transcripts from original count data
#write.csv(counts(unique.sig.list), file = "~/Desktop/mcav_unique.sig.list.csv")
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
# some larger than 4, so use rlog
unique.vst.sig <- rlog(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning message:
#   In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#               Estimated rdf < 1.0; not estimating variance

# PCA plot of diff-expressed genes 
mcav_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_mcav <- round(100*attr(mcav_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
mcav_DEGPCAplot <- ggplot(mcav_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_mcav[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_mcav[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  scale_color_manual(values = c(Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  #scale_color_manual(values = c(Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  ggtitle("Mcav w/o control") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
mcav_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- mcav_DEGPCAplot$data
ggsave("~/Desktop/mcav_treatment_DEGs_PCA.pdf", mcav_DEGPCAplot)





## Try removing weird outliers in T1 and T2 and rerun
# Subset count data for only mcav samples based on SampleID and make sure rows of metadata = cols of count data
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("24_T12_Ac_FM", acerv_metadata$SampleID),]
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("45_T41_Ac_SC_1", acerv_metadata$SampleID),]
# For some reason, not getting the outliers to remove. Just going to subset them by row number 
mcav_metadata_sub <- acerv_metadata[-c(2,10),]
rownames(acerv_metadata_sub) <- acerv_metadata_sub$SampleID # rename row names to reflect subsetted samples 
acerv_ID_sub <- acerv_metadata_sub$SampleID
count_acerv_treatment <- select(acerv_counts, all_of(acerv_ID_sub))
# continue maybe?








## Connecting IPS annotations to full annotations
mcav_IPS <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/InterProScan/mcav.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(mcav_IPS$V1)) # 716241
colnames(mcav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
mcav_IPS_GO <- filter(mcav_IPS, grepl("GO:", attr)) # select only rows with GO terms

# Rename annotation gff cols to merge annot and interproscan file 
annot$gene <- paste0(annot$gene, "-RA")
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "prot")

# merge annot and interproscan file by protein
mcav_merge <- merge(annot, mcav_IPS_GO, by = "prot", all.x = TRUE)
mcav_merge <- na.omit(mcav_merge)

# subset bu Pfam predictor
mcav_pfam_merge <- subset(mcav_merge, Predict == "Pfam")

# Use prot id to merge with DEGs file with full annot file -- will give final annotation of DEGs with GO terms, etc
colnames(DEGs.all) <- "prot"
DEGs.all$prot <- paste0(DEGs.all$prot, "-RA")
mcav_full_annot <- merge(DEGs.all, mcav_pfam_merge, by = "prot", all.x = TRUE)
write.csv(mcav_full_annot, file = "~/Desktop/mcav_full_annot.csv")






