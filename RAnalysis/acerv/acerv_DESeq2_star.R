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
acerv_counts <- read.csv("Output/DESeq2/acerv/acerv_gene_count_matrix.csv", header = TRUE, row.names = "gene_id")
dim(acerv_counts) # 33715 x 15
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(acerv_counts)[col]) # removing excess from col names
}
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  gsub("X", "", colnames(acerv_counts)[col]) # remvoing X from col names
}
# remove 24 and 45 samples 
acerv_counts <- acerv_counts[,-c(2,10)]

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
metadata <- read.csv("Data/FL_sediment_metadata.csv", header = TRUE)
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
# remove 24 and 45 samples 
acerv_metadata <- acerv_metadata[-c(2,10),]


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






##### Volcano plots 
## Here, the log transformed adjusted p-values are plotted on the y-axis and log2 fold change values on the x-axis (https://hbctraining.github.io/Intro-to-R-with-DGE/lessons/B1_DGE_visualizing_results.html)
# Read in data 
acerv.DEG <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv/acerv_sub_DEGs.all_treatment_20210219.csv")
acerv.DEG <- select(acerv.DEG, -X)

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

threshold <- acerv.DEG$padj < padj.cutoff & abs(acerv.DEG$log2FoldChange) > lfc.cutoff
length(which(threshold)) # this did not reduce anything, as the df only has DEGs in it?

# Add vector to df
acerv.DEG$threshold <- threshold   

# Volcano plot
acerv.volcano <- ggplot(acerv.DEG) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=Treatment_Compare)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
acerv.volcano
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/acerv/acerv_volcano.pdf", acerv.volcano, width = 28, height = 28, units = "cm")


## trying volcano plot with expanded data 
acerv_ByTreatment <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv/acerv_sub_ByTreatment_GO.terms_20210327.csv")
View(acerv_ByTreatment)

# Volcano plot
acerv.volcano <- ggplot(acerv_ByTreatment) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=term)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
acerv.volcano
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/acerv/acerv_volcano.GOterms.pdf", acerv.volcano, width = 28, height = 28, units = "cm")





