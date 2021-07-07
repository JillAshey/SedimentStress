# Title: DESeq2 with acerv samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date last modified: 2/8/21

# Data analyzed here is subsetted - removed 2 outliers

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
library("heatmaply")


# Load gene count matrix
acerv_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv/acerv_gene_count_matrix.csv", header = TRUE, row.names = "gene_id")
dim(acerv_counts) # 33715 x 15
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(acerv_counts)[col]) # removing excess from col names
}
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  gsub("X", "", colnames(acerv_counts)[col]) # remvoing X from col names
}
# the data has some gene names with 'TU' and some gene names with 'model'. I thought the TU and model were the same genes, just with one word switched 
# out for the other in certain files. Shoot i need go back to all acerv data and redo. 
# okay so the names with model in them get removed downstream because their counts are all 0 because they are not genes? need to check w/ star and stringtie
# so i'll remove the rows with model now so it doesn't confound anything downstream. 
acerv_counts <- acerv_counts[!grepl("model", rownames(acerv_counts)),]
dim(acerv_counts) # 18958 x 15

# Load annotation file
annot <- read.csv("~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", header = FALSE, sep = "\t") # gff annotation file 
# Looks like annotations with TU are the "Parent" or overall gene name and annotations with model are the gene parts (ie exon, CDS, etc)
# colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "blah") # dont look at blah col
# # okay to isolate gene names, I need to take the parent=xx from all ids except for gene and make it the gene name. 
# annot$gene <-sub("^([^;]*.[^;]*).*", "\\1", annot$attr) #remove everything after the third . in the attr column for all rows
# unique(annot$id)
# # [1] "gene"        "mRNA"        "exon"        "CDS"         "start_codon" "stop_codon"  "tRNA"       
# annot <- annot %>% separate(gene, c("ID", "Parent"), sep = ";") # separate ID and Parent in gene col
# annot <- na.omit(annot) # remove NAs
# # If id == anything BUT gene, remove Parent=; else, do nothing 
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "mRNA" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "exon" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "CDS" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "start_codon" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "stop_codon" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "tRNA" ,  
#                          gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# # If id == gene, remove Name=; else, do nothing
# annot <- annot %>% 
#   mutate(Parent = ifelse(id == "gene" ,  
#                          gsub("Name=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# # Make Parent names the same for each part in gene
# annot$Parent <- gsub("model", "TU", annot$Parent)
# annot$Parent <- gsub("gene", "evm", annot$Parent)
# annot <- select(annot, c(scaffold:attr, Parent)) # select cols 
# colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "gene") # rename cols 

# Load metadata
metadata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/FL_sediment_metadata.csv", header = TRUE)
dim(metadata) # 45 by 12
head(metadata)
# Select only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Rename cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Acerv species only 
acerv_metadata <- subset(metadata, Species=="Acropora cervicornis")
# Rename treatments
acerv_metadata$Treatment <- gsub("Ctl", "control", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T1", "Treatment1", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T2", "Treatment2", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T3", "Treatment3", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T4", "Treatment4", acerv_metadata$Treatment)
# Remove unwanted text from SampleID
acerv_metadata$SampleID <- gsub(".txt.gz", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- gsub(";.*", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- gsub(".fastq.gz", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- sub("\\.", "", acerv_metadata$SampleID)
# Make row names to be SampleID 
rownames(acerv_metadata) <- acerv_metadata$SampleID

## Try removing weird outliers in T1 and T2 and rerun
# Identified 24_T12_Ac_FM and 45_T41_Ac_SC_1 as outliers in previous script

# Subset acerv counts without the two outliers 
acerv_counts_sub <- select(acerv_counts, -c("24_T12_Ac_FM", "45_T41_Ac_SC_1"))

# Subset metadata to remove the two outliers 
rownames.remove <- c("24_T12_Ac_FM", "45_T41_Ac_SC_1")
acerv_metadata_sub <- acerv_metadata[!(rownames(acerv_metadata) %in% rownames.remove),]

# Check to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata_sub) %in% colnames(acerv_counts_sub)) # must come out TRUE


# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A = 0.85 and 5 (ie keep if 85% of samples have a count >5)
tfil <- genefilter(acerv_counts_sub, filt) # create filter for counts data 
keep <- acerv_counts_sub[tfil,] # identify genes to keep based on filter
# the gene names with model in them would get filtered out here because they do not have any counts 
gn.keep <- rownames(keep)
# Based on filt info, keep only the genes that pass in acerv_counts_filt
acerv_counts_sub_filt <- as.matrix(acerv_counts_sub[which(rownames(acerv_counts_sub) %in% gn.keep),]) 
write.csv(acerv_counts_sub_filt, "~/Desktop/acerv_counts_sub_filt_20210327.csv")
storage.mode(acerv_counts_sub_filt) <- "integer" # stores count data as integer 
# Check to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata_sub) %in% colnames(acerv_counts_sub_filt)) # must come out TRUE
# Set Treatment as a factor
acerv_metadata_sub$Treatment <- factor(acerv_metadata_sub$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
# create matrix that can be read by DESeq
data <- DESeqDataSetFromMatrix(countData = acerv_counts_sub_filt, colData = acerv_metadata_sub, design = ~ Treatment)

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

# do variance stabilizing transforamtion
vst <- vst(data, blind = FALSE) 

# Visualize filtered data with heatmaps and PCA
head(assay(vst), 3) # view data
sampleDists <- dist(t(assay(vst))) # calculate distance matrix
sampleDistMatrix <- as.matrix(sampleDists) # create distance matrix
rownames(sampleDistMatrix) <- colnames(vst) # assign row names
colnames(sampleDistMatrix) <- NULL # assign col names 
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) # assign colors 
pheatmap(sampleDistMatrix, # plot distance matrix
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
#[1] "Intercept"                       "Treatment_Treatment1_vs_control" "Treatment_Treatment2_vs_control"
#[4] "Treatment_Treatment3_vs_control" "Treatment_Treatment4_vs_control"
# NAs in padj: https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html 
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na


# Compare C vs T1
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1")) # results of DESeq2 comparing C and T1
DEG_control_vs_T1 
DEG_control_vs_T1 <- as.data.frame(DEG_control_vs_T1) # make results into a df
DEG_control_vs_T1["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
write.csv(DEG_control_vs_T1, file = "~/Desktop/acerv_sub_control_vs_T1_all_genes_20210219.csv") # maybe include gene counts too?
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 70 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=19). Trying to figure out why this is...Removed a sample from the pile? idk
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T1.sig.list_full <- cbind(DEG_control_vs_T1.sig, DEG_control_vs_T1.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T1.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T1_DEG_full_20210219.csv") # write out csv
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
## before running vst above, turn DEG_control_vs_T1.sig.list back into formal class DESeqTransform by running line DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),]
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2")) # results of DESeq2 comparing C and T2
DEG_control_vs_T2 
DEG_control_vs_T2 <- as.data.frame(DEG_control_vs_T2) # make results into a df
DEG_control_vs_T2["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
write.csv(DEG_control_vs_T2, file = "~/Desktop/acerv_sub_control_vs_T2_all_genes_20210219.csv") # maybe include gene counts too?
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T2.sig.num
# 141 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=35). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T2.sig.list <- as.data.frame(counts(DEG_control_vs_T2.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T2.sig.list_full <- cbind(DEG_control_vs_T2.sig, DEG_control_vs_T2.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T2.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T2_DEG_full_20210219.csv") # write out csv
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_control_vs_T3 
DEG_control_vs_T3 <- as.data.frame(DEG_control_vs_T3) # make results into a df
DEG_control_vs_T3["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
write.csv(DEG_control_vs_T3, file = "~/Desktop/acerv_sub_control_vs_T3_all_genes_20210219.csv") # maybe include gene counts too?
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T3.sig.num
# 56 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=20). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T3.sig.list <- as.data.frame(counts(DEG_control_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T3.sig.list_full <- cbind(DEG_control_vs_T3.sig, DEG_control_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T3.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T3_DEG_full_20210219.csv") # write out csv
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_control_vs_T4
DEG_control_vs_T4 <- as.data.frame(DEG_control_vs_T4) # make results into a df
DEG_control_vs_T4["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
write.csv(DEG_control_vs_T4, file = "~/Desktop/acerv_sub_control_vs_T4_all_genes_20210219.csv") # maybe include gene counts too?
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T4.sig.num
# 74 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=3). Trying to figure out why this is... Removed a sample from the pile? idk
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T4.sig.list_full <- cbind(DEG_control_vs_T4.sig, DEG_control_vs_T4.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T4.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T4_DEG_full_20210219.csv") # write out csv
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T2 
DEG_T1_vs_T2 <- as.data.frame(DEG_T1_vs_T2) # make results into a df
DEG_T1_vs_T2["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
write.csv(DEG_T1_vs_T2, file = "~/Desktop/acerv_sub_control_vs_T2_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T2.sig.num
# 51 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=0). Trying to figure out why this is...Removed a sample from the pile? idk
DEG_T1_vs_T2.sig <- subset(DEG_T1_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T2.sig["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # subset list of significant genes from original count data 
DEG_T1_vs_T2.sig.list <- as.data.frame(counts(DEG_T1_vs_T2.sig.list)) # make list of sig gene counts into a df
DEG_T1_vs_T2.sig.list_full <- cbind(DEG_T1_vs_T2.sig, DEG_T1_vs_T2.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T1_vs_T2.sig.list_full, file = "~/Desktop/acerv_sub_T1_vs_T2_DEG_full_20210219.csv") # write out csv
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_T1_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T3
DEG_T1_vs_T3 <- as.data.frame(DEG_T1_vs_T3) # make results into a df
DEG_T1_vs_T3["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
write.csv(DEG_T1_vs_T3, file = "~/Desktop/acerv_sub_T1_vs_T3_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T3.sig.num
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T4
DEG_T1_vs_T4 <- as.data.frame(DEG_T1_vs_T4) # make results into a df
DEG_T1_vs_T4["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
write.csv(DEG_T1_vs_T4, file = "~/Desktop/acerv_sub_T1_vs_T4_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T4.sig.num
# 0 DEGs

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_T2_vs_T3
DEG_T2_vs_T3 <- as.data.frame(DEG_T2_vs_T3) # make results into a df
DEG_T2_vs_T3["Treatment_Compare"] <- "T12sT3" # adding treatment comparison column
write.csv(DEG_T2_vs_T3, file = "~/Desktop/acerv_sub_T2_vs_T3_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T2_vs_T3.sig.num
# 3 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=0). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_T2_vs_T3.sig <- subset(DEG_T2_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T2_vs_T3.sig["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
DEG_T2_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T2_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_T2_vs_T3.sig.list <- as.data.frame(counts(DEG_T2_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_T2_vs_T3.sig.list_full <- cbind(DEG_T2_vs_T3.sig, DEG_T2_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T2_vs_T3.sig.list_full, file = "~/Desktop/acerv_sub_T2_vs_T3_DEG_full_20210219.csv") # write out csv
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # turn back into formal class DESeqTransform or else vst will not run
DEG_T1_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T2_vs_T4
DEG_T2_vs_T4 <- as.data.frame(DEG_T2_vs_T4) # make results into a df
DEG_T2_vs_T4["Treatment_Compare"] <- "T2vsT4" # adding treatment comparison column
write.csv(DEG_T2_vs_T4, file = "~/Desktop/acerv_sub_T2_vs_T4_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Compare T3 vs T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T3_vs_T4
DEG_T3_vs_T4 <- as.data.frame(DEG_T3_vs_T4) # make results into a df
DEG_T3_vs_T4["Treatment_Compare"] <- "T3vsT4" # adding treatment comparison column
write.csv(DEG_T3_vs_T4, file = "~/Desktop/acerv_sub_T3_vs_T4_all_genes_20210219.csv") # maybe include gene counts too?
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T3_vs_T4.sig.num
# 0 DEGs


# Make full list of genes and treatments
DEG_control_vs_T1.sig.list_full$gene_id <- rownames(DEG_control_vs_T1.sig.list_full)
rownames(DEG_control_vs_T1.sig.list_full) <- NULL
DEG_control_vs_T2.sig.list_full$gene_id <- rownames(DEG_control_vs_T2.sig.list_full)
rownames(DEG_control_vs_T2.sig.list_full) <- NULL
DEG_control_vs_T3.sig.list_full$gene_id <- rownames(DEG_control_vs_T3.sig.list_full)
rownames(DEG_control_vs_T3.sig.list_full) <- NULL
DEG_control_vs_T4.sig.list_full$gene_id <- rownames(DEG_control_vs_T4.sig.list_full)
rownames(DEG_control_vs_T4.sig.list_full) <- NULL
DEG_T1_vs_T2.sig.list_full$gene_id <- rownames(DEG_T1_vs_T2.sig.list_full)
rownames(DEG_T1_vs_T2.sig.list_full) <- NULL
DEG_T2_vs_T3.sig.list_full$gene_id <- rownames(DEG_T2_vs_T3.sig.list_full)
rownames(DEG_T2_vs_T3.sig.list_full) <- NULL

DEGs.all <- rbind(DEG_control_vs_T1.sig.list_full, 
                  DEG_control_vs_T2.sig.list_full,
                  DEG_control_vs_T3.sig.list_full, 
                  DEG_control_vs_T4.sig.list_full,
                  DEG_T1_vs_T2.sig.list_full,
                  DEG_T2_vs_T3.sig.list_full)
dim(DEGs.all) # 395 x 21
length(unique(DEGs.all$gene_id)) # 215 unique genes between all treatments 
write.csv(DEGs.all, file = "~/Desktop/acerv_sub_DEGs.all_treatment_20210219.csv")

## Find intersections and unique results between treatments 
# interactions
int1 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T2.sig.list_full$gene_id)
length(unique(int1)) # 45 DEGs shared between CvT1 and CvT2
int2 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T3.sig.list_full$gene_id)
length(unique(int2)) # 40 DEGs shared between CvT1 and CvT3
int3 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int3)) # 37 DEGs shared between CvT1 and CvT2
int4 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_T1_vs_T2.sig.list_full$gene_id)
length(unique(int4)) # 3 DEGs shared between CvT1 and T1vT2
int5 <- intersect(DEG_control_vs_T1.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int5)) # 0 DEGs shared between CvT1 and T2vT3
int6 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_control_vs_T3.sig.list_full$gene_id)
length(unique(int6)) # 50 DEGs shared between CvT2 and CvT3
int7 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int7)) # 67 DEGs shared between CvT2 and CvT4
int8 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_T1_vs_T2.sig.list_full$gene_id)
length(unique(int8)) # 7 DEGs shared between CvT2 and T1vT2
int9 <- intersect(DEG_control_vs_T2.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int9)) # 3 DEGs shared between CvT2 and T2vT3
int10 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_control_vs_T4.sig.list_full$gene_id)
length(unique(int10)) # 43 DEGs shared between CvT3 and T2vT4
int11 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_T1_vs_T2.sig.list_full$gene_id)
length(unique(int11)) # 0 DEGs shared between CvT3 and T1vT2
int12 <- intersect(DEG_control_vs_T3.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int12)) # 0 DEGs shared between CvT3 and T2vT3
int13 <- intersect(DEG_control_vs_T4.sig.list_full$gene_id, DEG_T1_vs_T2.sig.list_full$gene_id)
length(unique(int13)) # 0 DEGs shared between CvT4 and T1vT2
int14 <- intersect(DEG_control_vs_T4.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int14)) # 0 DEGs shared between CvT4 and T1vT2
int15 <- intersect(DEG_T1_vs_T2.sig.list_full$gene_id, DEG_T2_vs_T3.sig.list_full$gene_id)
length(unique(int15)) # 3 DEGs shared between T1vT2 and T3vT4


##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT2, T2vsT3
DEGs.all_acerv_sub <- DEGs.all$gene_id
DEGs.all_acerv_sub <- unique(DEGs.all_acerv_sub)
DEGs.all_acerv_sub <- as.data.frame(DEGs.all_acerv_sub) 
dim(DEGs.all_acerv_sub) # 215 unique DEGs among treatment comparisons 

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv_sub$DEGs), ] # subset list of sig transcripts from original count data
write.csv(counts(unique.sig.list), file = "~/Desktop/acerv_sub_unique.sig.list_20210219.csv")
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# PCA plot of diff-expressed genes 
acerv_sub_DEG_PCA <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_acerv_sub <- round(100*attr(acerv_sub_DEG_PCA, "percentVar")) #plot PCA of samples with all data
acerv_sub_DEG_PCA_plot <- ggplot(acerv_sub_DEG_PCA, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=8) +
  #geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_acerv_sub[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv_sub[2],"% variance")) +
  scale_color_manual(values = c(control="gray", Treatment1="darkslategray2", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  #ggtitle("A. cervicornis") +
  theme_bw() + #Set background color
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size=20),
        #title = element_text(size=30),
        legend.position = "right",
        panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_sub_DEG_PCA_plot
# PCA plot is of differentially expressed genes only
#PC.info <- mcav_DEGPCAplot$data
ggsave("~/Desktop/acerv_sub_DEGs_PCA_20210705.jpeg", acerv_sub_DEG_PCA_plot, width = 25, height = 25, units = "cm")
ggsave("~/Desktop/acerv_sub_DEGs_PCA_20210705.pdf", acerv_sub_DEG_PCA_plot, width = 25, height = 25, units = "cm")


## Heatmap of DEGs
# This heatmap is going to be wild. Grouping columns by treatment, putting gene id on left hand side and GO term on right hand side
# Also taking out legend 

# Need a file with counts, gene names, GO IDs, term, and ontology 
# This file is acerv genes with GO terms associated with them. One GO term per line, so multiple acerv gene names sometimes if genes have multiple GO terms
# This file includes all gene names and GO IDs
acerv_sig <- read.csv("~/Desktop/acerv_sub_GOterms.unique.csv", header = TRUE)
acerv_sig <- select(acerv_sig, -X)
colnames(acerv_sig)[1] <-"gene"
head(acerv_sig)


##### Volcano plots 
## Here, the log transformed adjusted p-values are plotted on the y-axis and log2 fold change values on the x-axis (https://hbctraining.github.io/Intro-to-R-with-DGE/lessons/B1_DGE_visualizing_results.html)
# Read in data 
acerv.DEG <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv/acerv_sub_DEGs.all_treatment_20210219.csv")
acerv.DEG <- select(acerv.DEG, -X)
acerv.DEG$diffexpressed <- "NA"
acerv.DEG$diffexpressed[acerv.DEG$log2FoldChange > 0] <- "Up"
acerv.DEG$diffexpressed[acerv.DEG$log2FoldChange < 0] <- "Down"

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

threshold <- acerv.DEG$padj < padj.cutoff & abs(acerv.DEG$log2FoldChange) > lfc.cutoff
length(which(threshold)) # this did not reduce anything, as the df only has DEGs in it?

# Add vector to df
acerv.DEG$threshold <- threshold   

# Volcano plot w/ DEGs
acerv.volcano <- ggplot(acerv.DEG) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), shape=Treatment_Compare, colour=diffexpressed), size = 2) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
acerv.volcano
ggsave("~/Desktop/acerv_volcano_20210705.pdf", acerv.volcano, width = 25, height = 25)
ggsave("~/Desktop/acerv_volcano_20210705.jpeg", acerv.volcano, width = 25, height = 25)


## trying volcano plot with GO term and GO slim data 
acerv_ByTreatment <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv/acerv_sub_ByTreatment_GO.terms_20210327.csv")
names(acerv_ByTreatment)[names(acerv_ByTreatment) == "category"] <- "GO.IDs"
acerv_ByTreatment$diffexpressed <- "NA"
acerv_ByTreatment$diffexpressed[acerv_ByTreatment$log2FoldChange > 0] <- "Up"
acerv_ByTreatment$diffexpressed[acerv_ByTreatment$log2FoldChange < 0] <- "Down"

# Read in GO slim info
go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
acerv_ByTreatment <- merge(acerv_ByTreatment, go.slim, by="GO.IDs", all = TRUE) # merge pdam info and GOslim
acerv_ByTreatment <- na.omit(acerv_ByTreatment)

# Volcano plot
acerv_go.volcano <- ggplot(acerv_ByTreatment) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, shape=GO.Slim.Term), size = 3) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
acerv_go.volcano
ggsave("~/Desktop/acerv_volcano.GOterms_20210705.pdf", acerv_go.volcano, width = 25, height = 25)
ggsave("~/Desktop/acerv_volcano.GOterms_20210705.jpeg", acerv_go.volcano, width = 25, height = 25)







