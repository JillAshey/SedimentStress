# Title: DESeq2 with Acerv samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/16/20

# Code for Francois sedimentation data. A. cerv only samples analyzed here aligned against A. cerv. STAR was read aligner with gff annotation file from Baums lab (pers. comm.)

# Code for Francois sedimentation data. A. cerv only samples analyzed here aligned against A. cerv. STAR was read aligner with gff annotation file from Baums lab (pers. comm.)
# Using Hollie code from pdam tawainn experiment as a guide. my code is being weird and i want to see if different code produces the same results 

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
acerv_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/gene_count_acerv_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(acerv_counts) # 33715 x 15
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(acerv_counts)[col]) # removing excess from col names
}
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  gsub("X", "", colnames(acerv_counts)[col]) # remvoing X from col names
}

# Load annotation file
annot <- read.csv("~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", header = FALSE, sep = "\t") # gff annotation file 
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "blah") # dont look at blah col
# okay to isolate gene names, I need to take the parent=xx from all ids except for gene and make it the gene name. 
annot$gene <-sub("^([^;]*.[^;]*).*", "\\1", annot$attr) #remove everything after the third . in the attr column for all rows
unique(annot$id)
# [1] "gene"        "mRNA"        "exon"        "CDS"         "start_codon" "stop_codon"  "tRNA"       
annot <- annot %>% separate(gene, c("ID", "Parent"), sep = ";") # separate ID and Parent in gene col
annot <- na.omit(annot) # remove NAs
# If id == anything BUT gene, remove Parent=; else, do nothing 
annot <- annot %>% 
  mutate(Parent = ifelse(id == "mRNA" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
annot <- annot %>% 
  mutate(Parent = ifelse(id == "exon" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
annot <- annot %>% 
  mutate(Parent = ifelse(id == "CDS" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
annot <- annot %>% 
  mutate(Parent = ifelse(id == "start_codon" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
annot <- annot %>% 
  mutate(Parent = ifelse(id == "stop_codon" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
annot <- annot %>% 
  mutate(Parent = ifelse(id == "tRNA" ,  
                         gsub("Parent=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# If id == gene, remove Name=; else, do nothing
annot <- annot %>% 
  mutate(Parent = ifelse(id == "gene" ,  
                         gsub("Name=", "", Parent, fixed = TRUE), gsub("", "", Parent)))
# Make Parent names the same for each part in gene
annot$Parent <- gsub("model", "TU", annot$Parent)
annot$Parent <- gsub("gene", "evm", annot$Parent)
annot <- select(annot, c(scaffold:attr, Parent)) # select cols 
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "gene") # rename cols 

# Load metadata
metadata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata.csv", header = TRUE)
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
# Make sampleID as rownames in metadata 
rownames(acerv_metadata) <- acerv_metadata$SampleID

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A = 0.85 and 5 (ie 85% of samples )
tfil <- genefilter(acerv_counts, filt) # create filter for counts data 
keep <- acerv_counts[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
# Based on filt info, keep only the genes that pass in acerv_counts_filt
acerv_counts_filt <- as.matrix(acerv_counts[which(rownames(acerv_counts) %in% gn.keep),]) 
write.csv(acerv_counts_filt, "~/Desktop/acerv_counts_filtered.csv")
storage.mode(acerv_counts_filt) <- "integer" # stores count data as integer 
# Check to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata) %in% colnames(acerv_counts_filt)) # must come out TRUE
# Set Treatment as a factor
acerv_metadata$Treatment <- factor(acerv_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
# create matrix that can be read in DESeq
data <- DESeqDataSetFromMatrix(countData = acerv_counts_filt, colData = acerv_metadata, design = ~ Treatment)

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

vst <- vst(data, blind = FALSE) # apply regularized log transformation to minimize effects of small counts and normalize wrt library
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
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing
# Save DE results
DEG.int.res <- results(DEG.int) # why does it say 'Wald test p-value: Treatment Treatment4 vs control' for DEG.int.res? Is it only looking at treatment 4 and control? In DESeq object created above, it says that design is Treatment
resultsNames(DEG.int) # view DE results 

# Compare C vs T1
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1")) # results of DESeq2 comparing C and T1
DEG_control_vs_T1 
DEG_control_vs_T1 <- as.data.frame(DEG_control_vs_T1) # make results into a df
DEG_control_vs_T1["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
write.csv(DEG_control_vs_T1, file = "~/Desktop/acerv_control_vs_T1_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 30 DEGs
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T1.sig.list_full <- cbind(DEG_control_vs_T1.sig, DEG_control_vs_T1.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T1.sig.list_full, file = "~/Desktop/acerv_control_vs_T1_DEG_full.csv") # write out csv
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2")) # results of DESeq2 comparing C and T1
DEG_control_vs_T2 
DEG_control_vs_T2 <- as.data.frame(DEG_control_vs_T2) # make results into a df
DEG_control_vs_T2["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
write.csv(DEG_control_vs_T2, file = "~/Desktop/acerv_control_vs_T2_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T2.sig.num
# 35 DEGs
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T2.sig.list <- as.data.frame(counts(DEG_control_vs_T2.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T2.sig.list_full <- cbind(DEG_control_vs_T2.sig, DEG_control_vs_T2.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T2.sig.list_full, file = "~/Desktop/acerv_control_vs_T2_DEG_full.csv") # write out csv
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3")) # results of DESeq2 comparing C and T1
DEG_control_vs_T3
DEG_control_vs_T3 <- as.data.frame(DEG_control_vs_T3) # make results into a df
DEG_control_vs_T3["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
write.csv(DEG_control_vs_T3, file = "~/Desktop/acerv_control_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T3.sig.num
# 20 DEGs
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T3.sig.list <- as.data.frame(counts(DEG_control_vs_T3.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T3.sig.list_full <- cbind(DEG_control_vs_T3.sig, DEG_control_vs_T3.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T3.sig.list_full, file = "~/Desktop/acerv_control_vs_T3_DEG_full.csv") # write out csv
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4"))
DEG_control_vs_T4
DEG_control_vs_T4 <- as.data.frame(DEG_control_vs_T4) # make results into a df
DEG_control_vs_T4["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
write.csv(DEG_control_vs_T4, file = "~/Desktop/acerv_control_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T4.sig.num
# 3 DEGs
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list))
DEG_control_vs_T4.sig.list_full <- cbind(DEG_control_vs_T4.sig, DEG_control_vs_T4.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T4.sig.list_full, file = "~/Desktop/acerv_control_vs_T4_DEG_full.csv") # write out csv
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2"))
DEG_T1_vs_T2 # why are there NAs? How are those generated?
DEG_T1_vs_T2 <- as.data.frame(DEG_T1_vs_T2) # make results into a df
DEG_T1_vs_T2["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
write.csv(DEG_T1_vs_T2, file = "~/Desktop/acerv_T1_vs_T2_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T2.sig.num
# 0 DEGs

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3"))
DEG_T1_vs_T3
DEG_T1_vs_T3 <- as.data.frame(DEG_T1_vs_T3) # make results into a df
DEG_T1_vs_T3["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
write.csv(DEG_T1_vs_T3, file = "~/Desktop/acerv_T1_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T3.sig.num
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4 <- as.data.frame(DEG_T1_vs_T4) # make results into a df
DEG_T1_vs_T4["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
write.csv(DEG_T1_vs_T4, file = "~/Desktop/acerv_T1_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num
# 1 DEGs
DEG_T1_vs_T4.sig <- subset(DEG_T1_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T4.sig["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
DEG_T1_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T4.sig)),] # subsey list of significant genes from original count data 
DEG_T1_vs_T4.sig.list <- as.data.frame(counts(DEG_T1_vs_T4.sig.list))
DEG_T1_vs_T4.sig.list_full <- cbind(DEG_T1_vs_T4.sig, DEG_T1_vs_T4.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T1_vs_T4.sig.list_full, file = "~/Desktop/acerv_T1_vs_T4_DEG_full.csv") # write out csv
DEG_T1_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3"))
DEG_T2_vs_T3
DEG_T2_vs_T3 <- as.data.frame(DEG_T2_vs_T3) # make results into a df
DEG_T2_vs_T3["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
write.csv(DEG_T2_vs_T3, file = "~/Desktop/acerv_T2_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T3.sig.num
# 0 DEGs

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4"))
DEG_T2_vs_T4
DEG_T2_vs_T4 <- as.data.frame(DEG_T2_vs_T4) # make results into a df
DEG_T2_vs_T4["Treatment_Compare"] <- "T2vsT4" # adding treatment comparison column
write.csv(DEG_T2_vs_T4, file = "~/Desktop/acerv_T2_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Compare T3 and T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4 <- as.data.frame(DEG_T3_vs_T4) # make results into a df
DEG_T3_vs_T4["Treatment_Compare"] <- "T3vsT4" # adding treatment comparison column
write.csv(DEG_T3_vs_T4, file = "~/Desktop/acerv_T3_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num
# 0 DEGs

##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT4
DEGs.all <- rbind(DEG_control_vs_T1.sig.list_full, 
              DEG_control_vs_T2.sig.list_full,
              DEG_control_vs_T3.sig.list_full, 
              DEG_control_vs_T4.sig.list_full,
              DEG_T1_vs_T4.sig.list_full
              )
write.csv(DEGs.all, file = "~/Desktop/acerv_DEGs.all_treatment.csv")
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all_acerv <- DEGs.all$DEGs
DEGs.all_acerv <- unique(DEGs.all_acerv)
DEGs.all_acerv <- as.data.frame(DEGs.all_acerv) # 89 unique DEGs among treatment comparisons

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv$DEGs), ] # subset list of sig transcripts from original count data
write.csv(counts(unique.sig.list), file = "~/Desktop/acerv_unique.sig.list.csv")
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# PCA.plot <- plotPCA(unique.rsig, intgroup = "Treatment") # plot PCA of all samples for DEG only 
# PCA.plot
# PC.info <- PCA.plot$data
# dev.off()
# jpeg(file="Output/Unique_PCA.DEG.jpg")
# plot(PC.info$PC1, PC.info$PC2, xlab="PC1 94%", ylab="PC2 2%", pch = c(15, 16)[as.numeric(sample.info$CO2)], col=c("gray", "black")[sample.info$Temperature], cex=1.3)
# legend(x="bottomright", 
#        bty="n",
#        legend = c("ATAC", "ATHC", "HTAC", "HTHC"),
#        pch = c(15, 16),
#        col=c("gray", "gray", "black", "black"),
#        cex=1)
# dev.off()

# PCA plot of diff-expressed genes 
acerv_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE) # create PCA object
percentVar_pca_acerv <- round(100*attr(acerv_DEGPCAdata, "percentVar")) # PC variance
acerv_DEGPCAplot <- ggplot(acerv_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_acerv[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  scale_color_manual(values = c(control="lightpink", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- acerv_DEGPCAplot$data
ggsave("~/Desktop/acerv_DEGs_PCA.pdf", acerv_DEGPCAplot)

# Heatmap of DEGs
df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")]) # make df of treatment by sample
colnames(df) <- "Treatment"
ann_colors <- list(Treatment = c(control="lightpink", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray"))
#ann_colors <- list(Treatment = c(control="lightpink", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4"))
col.order <- c("19_T33_Ac_WK",
               "24_T12_Ac_FM",
               "25_ctl1_Ac_GF_1",
               "27_ctl2_Ac_YG_1",
               "31_T22_Ac_UV",
               "35_T43_Ac_MT",
               "37_T13_Ac_ML",
               "38_T23_Ac_IN",
               "41_ctl3_Ac_RN_1",
               "45_T41_Ac_SC_1",
               "47_T31_Ac_JB",
               "52_T11_Ac_II",
               "53_T21_Ac_NH",
               "54_T42_Ac_JQ",
               "57_T32_Ac_NM") # specify column order for heatmap 

unique.DEG.annot <- as.data.frame(counts(unique.sig.list)) # make df of sig genes with counts and sample IDs
unique.DEG.annot$gene <- rownames(unique.DEG.annot) # make column with gene names
unique.DEG.annot <- merge(unique.DEG.annot, annot, by = "gene") # merge annotation file with df of sig genes by gene id
unique.DEG.annot <- unique.DEG.annot[!duplicated(unique.DEG.annot$gene),] # remove duplicate rows
rownames(unique.DEG.annot) <- unique.DEG.annot$gene
write.csv(unique.DEG.annot, file = "~/Desktop/acerv_unique_DEG_annotated.csv")

unique.DEG.annot <- unique.DEG.annot[ ,c(2:16)] # select only gene counts
rownames(df) <- colnames(unique.DEG.annot) # set df row names 
# unique.DEG.annot <- unique.DEG.annot[,-16]
mat <- as.matrix(unique.DEG.annot) # create matrix 

mat <- mat[,col.order] # order matrix by col order
#dev.off()
#pdf(file = "~/Desktop/Unique_Heatmap.DEG_Annotated.pdf")
acerv_heatmap <- pheatmap(mat, 
                          annotation_col = df,
                          annotation_colors = ann_colors,
                          scale = "row",
                          show_rownames = T,
                          fontsize_row = 4,
                          cluster_cols = T,
                          cluster_rows = T,
                          show_colnames = F
                          )
#dev.off()
# plot has all treatment comparisons 
ggsave("~/Desktop/acerv_DEGs_heatmap.pdf", acerv_heatmap)








## Going to remove the control from the list and see how DEG PCA looks without it--can look at spread of mid and high 
# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
metadata_acerv_treatment <- subset(acerv_metadata, Treatment=="Treatment1" | Treatment=="Treatment2" | Treatment=="Treatment3" | Treatment=="Treatment4")
rownames(metadata_acerv_treatment) <- metadata_acerv_treatment$SampleID # rename row names to reflect subsetted samples 
acerv_ID_treatment <- metadata_acerv_treatment$SampleID
count_acerv_treatment <- select(acerv_counts, all_of(acerv_ID_treatment))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_acerv_treatment, filt) # create filter for counts data 
keep <- count_acerv_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
acerv_treatment_counts_filt <- as.matrix(count_acerv_treatment[which(rownames(count_acerv_treatment) %in% gn.keep),]) 
storage.mode(acerv_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(acerv_treatment_counts_filt, "~/Desktop/acerv_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_acerv_treatment) %in% colnames(acerv_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_acerv_treatment$Treatment <- factor(metadata_acerv_treatment$Treatment, levels = c("Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = acerv_treatment_counts_filt, colData = metadata_acerv_treatment, design = ~ Treatment)

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
#[1] "Intercept"                          "Treatment_Treatment2_vs_Treatment1" "Treatment_Treatment3_vs_Treatment1"
#[4] "Treatment_Treatment4_vs_Treatment1"

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
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num
# 5 DEGs
DEG_T1_vs_T4.sig <- subset(DEG_T1_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T4.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T1_vs_T4.sig.list)
print(sizeFactors(SFtest))
DEG_T1_vs_T4.rsig <- varianceStabilizingTransformation(DEG_T1_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
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
# 0 DEGs

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4"))
DEG_T2_vs_T4
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Compare T3 vs T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num
# 0 DEGs

# Only one group with DEGs - T1vsT4
# PCA plot of diff-expressed genes 
acerv_DEGPCAdata <- plotPCA(DEG_T1_vs_T4.rsig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_acerv <- round(100*attr(acerv_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
acerv_DEGPCAplot <- ggplot(acerv_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_acerv[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  scale_color_manual(values = c(Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  #scale_color_manual(values = c(Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  ggtitle("Acerv no control treatment") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- acerv_DEGPCAplot$data
ggsave("~/Desktop/acerv_treatment_DEGs_PCA.pdf", acerv_DEGPCAplot)












## Try removing weird outliers in T1 and T2 and rerun
# Subset count data for only acerv samples based on SampleID and make sure rows of metadata = cols of count data
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("24_T12_Ac_FM", acerv_metadata$SampleID),]
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("45_T41_Ac_SC_1", acerv_metadata$SampleID),]
# For some reason, not getting the outliers to remove. Just going to subset them by row number 
acerv_metadata_sub <- acerv_metadata[-c(2,10),]
rownames(acerv_metadata_sub) <- acerv_metadata_sub$SampleID # rename row names to reflect subsetted samples 
acerv_ID_sub <- acerv_metadata_sub$SampleID
count_acerv_treatment <- select(acerv_counts, all_of(acerv_ID_sub))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_acerv_treatment, filt) # create filter for counts data 
keep <- count_acerv_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
acerv_treatment_counts_filt <- as.matrix(count_acerv_treatment[which(rownames(count_acerv_treatment) %in% gn.keep),]) 
storage.mode(acerv_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(acerv_treatment_counts_filt, "~/Desktop/acerv_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata_sub) %in% colnames(acerv_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
acerv_metadata_sub$Treatment <- factor(acerv_metadata_sub$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = acerv_treatment_counts_filt, colData = acerv_metadata_sub, design = ~ Treatment)

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
#[1] "Intercept"                       "Treatment_Treatment1_vs_control" "Treatment_Treatment2_vs_control"
#[4] "Treatment_Treatment3_vs_control" "Treatment_Treatment4_vs_control"

# Compare C vs T1
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1"))
DEG_control_vs_T1
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 70 DEGs
# why is it different than when outliers are included?
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list))
DEG_control_vs_T1.sig.list["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
#write.csv(DEG_control_vs_T1.sig.list, file = "~/Desktop/acerv_control_vs_T1_DEG.csv")

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2"))
DEG_control_vs_T2
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj< 0.05
DEG_control_vs_T2.sig.num
# 141 DEGs
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T2.sig.list <- as.data.frame(counts(DEG_control_vs_T2.sig.list))
DEG_control_vs_T2.sig.list["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
#write.csv(DEG_control_vs_T2.sig.list, file = "~/Desktop/acerv_control_vs_T2_DEG.csv")

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3"))
DEG_control_vs_T3
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T3.sig.num
# 56 DEGs
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T3.sig.list <- as.data.frame(counts(DEG_control_vs_T3.sig.list))
DEG_control_vs_T3.sig.list["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
#write.csv(DEG_control_vs_T3.sig.list, file = "~/Desktop/acerv_control_vs_T3_DEG.csv")

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4"))
DEG_control_vs_T4
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T4.sig.num
# 74 DEGs
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list))
DEG_control_vs_T4.sig.list["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
#write.csv(DEG_control_vs_T4.sig.list, file = "~/Desktop/acerv_control_vs_T4_DEG.csv")

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2"))
DEG_T1_vs_T2
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T2.sig.num
# 51 DEGs
DEG_T1_vs_T2.sig <- subset(DEG_T1_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # subsey list of significant genes from original count data 
DEG_T1_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# Error in estimateDispersionsFit(object, quiet = TRUE, fitType) : 
#   all gene-wise dispersion estimates are within 2 orders of magnitude
# from the minimum value, and so the standard curve fitting techniques will not work.
# One can instead use the gene-wise estimates as final estimates:
#   dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
# ...then continue with testing using nbinomWaldTest or nbinomLRT
DEG_T1_vs_T2.sig.list <- as.data.frame(counts(DEG_T1_vs_T2.sig.list))
DEG_T1_vs_T2.sig.list["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
#write.csv(DEG_T1_vs_T2.sig.list, file = "~/Desktop/acerv_T1_vs_T2_DEG.csv")

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3"))
DEG_T1_vs_T3
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T3.sig.num
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3"))
DEG_T2_vs_T3
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T3.sig.num
# 3 DEGs
DEG_T2_vs_T3.sig <- subset(DEG_T2_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T2_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T2_vs_T3.sig)),] # subsey list of significant genes from original count data 
DEG_T2_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_T2_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning message:
#   In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#               Estimated rdf < 1.0; not estimating variance
DEG_T2_vs_T3.sig.list <- as.data.frame(counts(DEG_T2_vs_T3.sig.list))
DEG_T2_vs_T3.sig.list["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
#write.csv(DEG_T2_vs_T3.sig.list, file = "~/Desktop/acerv_T2_vs_T3_DEG.csv")

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4"))
DEG_T2_vs_T4
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Compare T3 and T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num
# 0 DEGs

##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT4
DEGs_CvsT1 <- as.data.frame(rownames(DEG_control_vs_T1.sig.list), DEG_control_vs_T1.sig.list$Treatment_Compare)
colnames(DEGs_CvsT1) <- "DEGs"
DEGs_CvsT1$Treatment_Compare <- rownames(DEGs_CvsT1)
DEGs_CvsT2 <- as.data.frame(rownames(DEG_control_vs_T2.sig.list), DEG_control_vs_T2.sig.list$Treatment_Compare)
colnames(DEGs_CvsT2) <- "DEGs"
DEGs_CvsT2$Treatment_Compare <- rownames(DEGs_CvsT2)
DEGs_CvsT3 <- as.data.frame(rownames(DEG_control_vs_T3.sig.list), DEG_control_vs_T3.sig.list$Treatment_Compare)
colnames(DEGs_CvsT3) <- "DEGs"
DEGs_CvsT3$Treatment_Compare <- rownames(DEGs_CvsT3)
DEGs_CvsT4 <- as.data.frame(rownames(DEG_control_vs_T4.sig.list), DEG_control_vs_T4.sig.list$Treatment_Compare)
colnames(DEGs_CvsT4) <- "DEGs"
DEGs_CvsT4$Treatment_Compare <- rownames(DEGs_CvsT4)
DEGs_T1vsT2 <- as.data.frame(rownames(DEG_T1_vs_T2.sig.list), DEG_T1_vs_T2.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT2) <- "DEGs"
DEGs_T1vsT2$Treatment_Compare <- rownames(DEGs_T1vsT2)
DEGs_T2vsT3 <- as.data.frame(rownames(DEG_T2_vs_T3.sig.list), DEG_T2_vs_T3.sig.list$Treatment_Compare)
colnames(DEGs_T2vsT3) <- "DEGs"
DEGs_T2vsT3$Treatment_Compare <- rownames(DEGs_T2vsT3)

DEGs.all <- rbind(DEGs_CvsT1,
                  DEGs_CvsT2, 
                  DEGs_CvsT3, 
                  DEGs_CvsT4, 
                  DEGs_T1vsT2,
                  DEGs_T2vsT3) 
#write.csv(DEGs.all, file = "~/Desktop/acerv_DEGs.all_treatment.csv")
DEGs.all_acerv <- DEGs.all$DEGs
DEGs.all_acerv <- unique(DEGs.all_acerv)
DEGs.all_acerv <- as.data.frame(DEGs.all_acerv)

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv$DEGs), ] # subset list of sig transcripts from original count data
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
#write.csv(counts(unique.sig.list), file = "~/Desktop/acerv_unique.sig.list.csv")

# PCA plot of diff-expressed genes w/ outliers removed
acerv_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE) # create PCA object
percentVar_pca_acerv <- round(100*attr(acerv_DEGPCAdata, "percentVar")) # PC variance
acerv_DEGPCAplot <- ggplot(acerv_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  #geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_acerv[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  scale_color_manual(values = c(control="lightpink", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  ggtitle("Acerv w/o outliers") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- acerv_DEGPCAplot$data
ggsave("~/Desktop/acerv_outlierRemove_DEGs_PCA.pdf", acerv_DEGPCAplot)









## Now going to remove outliers AND remove control to see how DEG PCA looks 
# Subset count data for only acerv samples based on SampleID and make sure rows of metadata = cols of count data
metadata_acerv_treatment <- subset(acerv_metadata, Treatment=="Treatment1" | Treatment=="Treatment2" | Treatment=="Treatment3" | Treatment=="Treatment4")
metadata_acerv_treatment_sub <- metadata_acerv_treatment[-c(2,7),]
rownames(metadata_acerv_treatment_sub) <- metadata_acerv_treatment_sub$SampleID # rename row names to reflect subsetted samples 
acerv_ID_treatment_sub <- metadata_acerv_treatment_sub$SampleID
count_acerv_treatment <- select(acerv_counts, all_of(acerv_ID_treatment_sub))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_acerv_treatment, filt) # create filter for counts data 
keep <- count_acerv_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
acerv_treatment_counts_filt <- as.matrix(count_acerv_treatment[which(rownames(count_acerv_treatment) %in% gn.keep),]) 
storage.mode(acerv_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(acerv_treatment_counts_filt, "~/Desktop/acerv_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_acerv_treatment_sub) %in% colnames(acerv_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_acerv_treatment_sub$Treatment <- factor(metadata_acerv_treatment_sub$Treatment, levels = c("Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = acerv_treatment_counts_filt, colData = metadata_acerv_treatment_sub, design = ~ Treatment)


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
#[1] "Intercept"                          "Treatment_Treatment2_vs_Treatment1" "Treatment_Treatment3_vs_Treatment1"
#[4] "Treatment_Treatment4_vs_Treatment1"

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2"))
DEG_T1_vs_T2
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T2.sig.num
# 40 DEGs
DEG_T1_vs_T2.sig <- subset(DEG_T1_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T1_vs_T2.sig.list)
print(sizeFactors(SFtest))
DEG_T1_vs_T2.rsig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_T1_vs_T2.sig.list <- as.data.frame(counts(DEG_T1_vs_T2.sig.list))
DEG_T1_vs_T2.sig.list["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
#write.csv(DEG_T1_vs_T2.sig.list, file = "~/Desktop/acerv_T1_vs_T2_DEG.csv")

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3"))
DEG_T1_vs_T3
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T3.sig.num
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num
# 3 DEGs
DEG_T1_vs_T4.sig <- subset(DEG_T1_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T4.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T1_vs_T4.sig.list)
print(sizeFactors(SFtest))
DEG_T1_vs_T4.rsig <- varianceStabilizingTransformation(DEG_T1_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
#Warning message:
#  In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#              Estimated rdf < 1.0; not estimating variance
DEG_T1_vs_T4.sig.list <- as.data.frame(counts(DEG_T1_vs_T4.sig.list))
DEG_T1_vs_T4.sig.list["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
#write.csv(DEG_T1_vs_T4.sig.list, file = "~/Desktop/acerv_T1_vs_T4_DEG.csv")

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3"))
DEG_T2_vs_T3
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T3.sig.num
# 6 DEGs
DEG_T2_vs_T3.sig <- subset(DEG_T2_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_T2_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_T2_vs_T3.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_T2_vs_T3.sig.list)
print(sizeFactors(SFtest))
DEG_T2_vs_T3.rsig <- varianceStabilizingTransformation(DEG_T2_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
#Warning message:
#  In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#              Estimated rdf < 1.0; not estimating variance
DEG_T2_vs_T3.sig.list <- as.data.frame(counts(DEG_T2_vs_T3.sig.list))
DEG_T2_vs_T3.sig.list["Treatment_Compare"] <- "T2vsT3" # adding treatment comparison column
#write.csv(DEG_T2_vs_T3.sig.list, file = "~/Desktop/acerv_T2_vs_T3_DEG.csv")

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
DEGs_T1vsT2 <- as.data.frame(rownames(DEG_T1_vs_T2.sig.list), DEG_T1_vs_T2.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT2) <- "DEGs"
DEGs_T1vsT2$Treatment_Compare <- rownames(DEGs_T1vsT2)
DEGs_T1vsT4 <- as.data.frame(rownames(DEG_T1_vs_T4.sig.list), DEG_T1_vs_T4.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT4) <- "DEGs"
DEGs_T1vsT4$Treatment_Compare <- rownames(DEGs_T1vsT4)
DEGs_T2vsT3 <- as.data.frame(rownames(DEG_T2_vs_T3.sig.list), DEG_T2_vs_T3.sig.list$Treatment_Compare)
colnames(DEGs_T2vsT3) <- "DEGs"
DEGs_T2vsT3$Treatment_Compare <- rownames(DEGs_T2vsT3)

DEGs.all <- rbind(DEGs_T1vsT2,
                  DEGs_T1vsT4,
                  DEGs_T2vsT3) 
#write.csv(DEGs.all, file = "~/Desktop/acerv_DEGs.all_treatment.csv")
DEGs.all_acerv <- DEGs.all$DEGs
DEGs.all_acerv <- unique(DEGs.all_acerv)
DEGs.all_acerv <- as.data.frame(DEGs.all_acerv)

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv$DEGs), ] # subset list of sig transcripts from original count data
#write.csv(counts(unique.sig.list), file = "~/Desktop/mcav_unique.sig.list.csv")
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
# some larger than 4, so use rlog
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# PCA plot of diff-expressed genes 
acerv_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_acerv <- round(100*attr(acerv_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
acerv_DEGPCAplot <- ggplot(acerv_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_acerv[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv[2],"% variance")) +
  #scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  scale_color_manual(values = c(Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  #scale_color_manual(values = c(Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray")) +
  coord_fixed() +
  ggtitle("Acerv w/o control and outliers") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- acerv_DEGPCAplot$data
ggsave("~/Desktop/acerv_treatment_outlierRemove_DEGs_PCA.pdf", acerv_DEGPCAplot)










## Connecting IPS annotations to full annotations
acerv_IPS <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/InterProScan/acerv.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(acerv_IPS$V1)) # 981372
colnames(acerv_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
acerv_IPS_GO <- filter(acerv_IPS, grepl("GO:", attr)) # select only rows with GO terms

# Rename annotation gff cols to merge annot and interproscan file 
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr", "prot")
annot$prot <- gsub("TU", "model", annot$prot)   
            
# merge annot and interproscan file by protein
test <- merge(annot, acerv_IPS_GO, by = "prot")
test <- na.omit(test)

# subset bu Pfam predictor
pfam_test <- subset(test, Predict == "Pfam")

# modify DEG names to match those in the annotation file 
DEGs.all_acerv$DEGs <- gsub("TU", "model", DEGs.all_acerv$DEGs)
colnames(DEGs.all_acerv) <- "prot"

# Use prot id to merge with DEGs file with full annot file -- will give final annotation of DEGs with GO terms, etc
acerv_full_annot <- merge(DEGs.all_acerv, pfam_test, by = "prot", all.x = TRUE)
write.csv(acerv_full_annot, file = "~/Desktop/acerv_full_annot.csv")



























































































