

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
acerv_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv_gene_count_matrix.csv", header = TRUE, row.names = "gene_id")
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
# Make sampleID as rownames in metadata 
rownames(acerv_metadata) <- acerv_metadata$SampleID


## Try removing weird outliers in T1 and T2 and rerun
# Identified 24_T12_Ac_FM and 45_T41_Ac_SC_1 as outliers

# Subset acerv counts without the two outliers 
acerv_counts_sub <- select(acerv_counts, -c("24_T12_Ac_FM", "45_T41_Ac_SC_1"))

# Subset metadata
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("24_T12_Ac_FM", acerv_metadata$SampleID),]
#acerv_metadata_sub$SampleID <- acerv_metadata[!grepl("45_T41_Ac_SC_1", acerv_metadata$SampleID),]
# For some reason, not getting the outliers to remove. Just going to subset them by row number 
acerv_metadata_sub <- acerv_metadata[-c(2,10),]
rownames(acerv_metadata_sub) <- acerv_metadata_sub$SampleID # rename row names to reflect subsetted samples 
#acerv_ID_sub <- acerv_metadata_sub$SampleID

# Check to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata_sub) %in% colnames(acerv_counts_sub)) # must come out TRUE


# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A = 0.85 and 5 (ie 85% of samples )
tfil <- genefilter(acerv_counts_sub, filt) # create filter for counts data 
keep <- acerv_counts_sub[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
# Based on filt info, keep only the genes that pass in acerv_counts_filt
acerv_counts_sub_filt <- as.matrix(acerv_counts_sub[which(rownames(acerv_counts_sub) %in% gn.keep),]) 
write.csv(acerv_counts_sub_filt, "~/Desktop/acerv_counts_sub_filt.csv")
storage.mode(acerv_counts_sub_filt) <- "integer" # stores count data as integer 
# Check to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata_sub) %in% colnames(acerv_counts_sub_filt)) # must come out TRUE
# Set Treatment as a factor
acerv_metadata_sub$Treatment <- factor(acerv_metadata_sub$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
# create matrix that can be read in DESeq
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
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1")) # results of DESeq2 comparing C and T1
DEG_control_vs_T1 
DEG_control_vs_T1 <- as.data.frame(DEG_control_vs_T1) # make results into a df
DEG_control_vs_T1["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
write.csv(DEG_control_vs_T1, file = "~/Desktop/acerv_sub_control_vs_T1_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 70 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=19). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T1.sig.list_full <- cbind(DEG_control_vs_T1.sig, DEG_control_vs_T1.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T1.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T1_DEG_full.csv") # write out csv
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2")) # results of DESeq2 comparing C and T2
DEG_control_vs_T2 
DEG_control_vs_T2 <- as.data.frame(DEG_control_vs_T2) # make results into a df
DEG_control_vs_T2["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
write.csv(DEG_control_vs_T2, file = "~/Desktop/acerv_sub_control_vs_T2_all_genes.csv") # maybe include gene counts too?
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
write.csv(DEG_control_vs_T2.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T2_DEG_full.csv") # write out csv
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_control_vs_T3 
DEG_control_vs_T3 <- as.data.frame(DEG_control_vs_T3) # make results into a df
DEG_control_vs_T3["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
write.csv(DEG_control_vs_T3, file = "~/Desktop/acerv_sub_control_vs_T3_all_genes.csv") # maybe include gene counts too?
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
write.csv(DEG_control_vs_T3.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T3_DEG_full.csv") # write out csv
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_control_vs_T4
DEG_control_vs_T4 <- as.data.frame(DEG_control_vs_T4) # make results into a df
DEG_control_vs_T4["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
write.csv(DEG_control_vs_T4, file = "~/Desktop/acerv_sub_control_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T4.sig.num
# 74 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=3). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_T4.sig.list_full <- cbind(DEG_control_vs_T4.sig, DEG_control_vs_T4.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_T4.sig.list_full, file = "~/Desktop/acerv_sub_control_vs_T4_DEG_full.csv") # write out csv
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T1 vs T2
DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T2 
DEG_T1_vs_T2 <- as.data.frame(DEG_T1_vs_T2) # make results into a df
DEG_T1_vs_T2["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
write.csv(DEG_T1_vs_T2, file = "~/Desktop/acerv_sub_control_vs_T2_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T2.sig.num
# 51 DEGs
# not sure why this # of DEGs differs from the one that includes the outliers (DEG=0). Trying to figure out why this is...
# Removed a sample from the pile? idk
DEG_T1_vs_T2.sig <- subset(DEG_T1_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T2.sig["Treatment_Compare"] <- "T1vsT2" # adding treatment comparison column
DEG_T1_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T2.sig)),] # subset list of significant genes from original count data 
DEG_T1_vs_T2.sig.list <- as.data.frame(counts(DEG_T1_vs_T2.sig.list)) # make list of sig gene counts into a df
DEG_T1_vs_T2.sig.list_full <- cbind(DEG_T1_vs_T2.sig, DEG_T1_vs_T2.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_T1_vs_T2.sig.list_full, file = "~/Desktop/acerv_sub_T1_vs_T2_DEG_full.csv") # write out csv
DEG_T1_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T1 vs T3
DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T3
DEG_T1_vs_T3 <- as.data.frame(DEG_T1_vs_T3) # make results into a df
DEG_T1_vs_T3["Treatment_Compare"] <- "T1vsT3" # adding treatment comparison column
write.csv(DEG_T1_vs_T3, file = "~/Desktop/acerv_sub_T1_vs_T3_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T3.sig.num
# 0 DEGs

# Compare T1 vs T4
DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T1_vs_T4
DEG_T1_vs_T4 <- as.data.frame(DEG_T1_vs_T4) # make results into a df
DEG_T1_vs_T4["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
write.csv(DEG_T1_vs_T4, file = "~/Desktop/acerv_sub_T1_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T1_vs_T4.sig.num
# 0 DEGs

# Compare T2 vs T3
DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3")) # results of DESeq2 comparing C and T2
DEG_T2_vs_T3
DEG_T2_vs_T3 <- as.data.frame(DEG_T2_vs_T3) # make results into a df
DEG_T2_vs_T3["Treatment_Compare"] <- "T12sT3" # adding treatment comparison column
write.csv(DEG_T2_vs_T3, file = "~/Desktop/acerv_sub_T2_vs_T3_all_genes.csv") # maybe include gene counts too?
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
write.csv(DEG_T2_vs_T3.sig.list_full, file = "~/Desktop/acerv_sub_T2_vs_T3_DEG_full.csv") # write out csv
DEG_T1_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.

# Compare T2 vs T4
DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T2_vs_T4
DEG_T2_vs_T4 <- as.data.frame(DEG_T2_vs_T4) # make results into a df
DEG_T2_vs_T4["Treatment_Compare"] <- "T2vsT4" # adding treatment comparison column
write.csv(DEG_T2_vs_T4, file = "~/Desktop/acerv_sub_T2_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T2_vs_T4.sig.num
# 0 DEGs

# Compare T3 vs T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4")) # results of DESeq2 comparing C and T2
DEG_T3_vs_T4
DEG_T3_vs_T4 <- as.data.frame(DEG_T3_vs_T4) # make results into a df
DEG_T3_vs_T4["Treatment_Compare"] <- "T3vsT4" # adding treatment comparison column
write.csv(DEG_T3_vs_T4, file = "~/Desktop/acerv_sub_T3_vs_T4_all_genes.csv") # maybe include gene counts too?
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_T3_vs_T4.sig.num
# 0 DEGs


##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT2, T2vsT3
DEGs.all <- rbind(DEG_control_vs_T1.sig.list_full, 
                  DEG_control_vs_T2.sig.list_full,
                  DEG_control_vs_T3.sig.list_full, 
                  DEG_control_vs_T4.sig.list_full,
                  DEG_T1_vs_T2.sig.list_full,
                  DEG_T2_vs_T3.sig.list_full
)
write.csv(DEGs.all, file = "~/Desktop/acerv_sub_DEGs.all_treatment.csv")
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all_acerv_sub <- DEGs.all$DEGs
DEGs.all_acerv_sub <- unique(DEGs.all_acerv_sub)
DEGs.all_acerv_sub <- as.data.frame(DEGs.all_acerv_sub) # 395 unique DEGs among treatment comparisons
# Not really sure what the heck is going on her with DEG counts

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv_sub$DEGs), ] # subset list of sig transcripts from original count data
write.csv(counts(unique.sig.list), file = "~/Desktop/acerv_sub_unique.sig.list.csv")
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
ggsave("~/Desktop/acerv_sub_DEGs_PCA.png", acerv_sub_DEG_PCA_plot, width = 30, height = 20,, units = "cm")


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

# Now we will take the unique.DEG.annot2 and merge it with acerv_sig
# The unique.DEG.annot2 file includes gene names for DEGs and counts data
test_merge <- merge(unique.DEG.annot2, acerv_sig, by = "gene", all.x = TRUE)
# test_merge now holds gene names for DEGs, counts data, and GO.IDs

# Now we need info about term and ontology 
GO_all <- read.csv("~/Desktop/acerv_sub_GO_ALL.csv", header = TRUE) 
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
write.csv(merge_all, file = "~/Desktop/acerv_sub_GO_DEG.csv") # maybe include gene counts too?



# Hooray! Now I have a lovely file with counts, gene names, GO IDs, term, and ontology 
# Now I must put it in the heatmap...........

# First, lets make a matrix of gene counts 
rownames(merge_all) <- merge_all$gene
mat <- select(merge_all, -c("gene", "GO.ID", "term", "ontology", "over_represented_pvalue"))
mat <- as.matrix(mat)

# Now lets make df of only treatment and sample ID
df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
colnames(df) <- "Treatment"

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
acerv_sub_heatmap <- pheatmap(mat, 
                              annotation_col = df,
                             # annotation_row = df_gene,
                              annotation_colors = ann_colors,
                              annotation_legend = F,
                              cluster_rows = F,
                              show_rownames = T,
                              cluster_cols = F,
                              show_colnames = T,
                              scale = "row",
                              fontsize_row = 8,
                              labels_row = merge_all$term,
                              naprint=F
                              )
acerv_sub_heatmap
ggsave("~/Desktop/acerv_sub_heatmap.png", acerv_sub_heatmap, width = 30, height = 20,, units = "cm")















merge_all <- merge(merge_all, agg_ont, by = "gene", all.x = TRUE)
write.csv(GO_all, file = "~/Desktop/acerv_sub_GO_all_gene_ID.csv") # maybe include gene counts too?












GO_test_merge <- merge(test_merge, acerv_sig, by = "gene", all.x = TRUE)
# So now we got GO.x and GO.y in this file, which are GO terms that came from the two files. some match and some dont. Why?? Okay I see
# In the acerv_sig file, there were repeats of genes because some genes had multiple GO terms associated. So Go.y is coming from that, so I'll keep that one
# Go.x is repeating
GO_test_merge <- select(GO_test_merge, -GO.ID.x)
colnames(GO_test_merge)[15] <-"category" # renaming GO col so that it can merge with the info generated from GOseq























# Creating df of only treatment
df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
colnames(df) <- "Treatment"
#df <- df[order(df$Treatment),]
#df <- as.data.frame(df)
#colnames(df) <- "Treatment"

# Getting col order 
list(acerv_metadata_sub)
test <- acerv_metadata_sub[order(acerv_metadata_sub$Treatment),] # need to order them so it will group by treatment in plot
list(test$SampleID) # look at sample IDs and use that list to make col.order
col.order <- c("25_ctl1_Ac_GF_1",
               "27_ctl2_Ac_YG_1",
               "41_ctl3_Ac_RN_1",
               "37_T13_Ac_ML",
               "52_T11_Ac_II",
               "31_T22_Ac_UV",
               "38_T23_Ac_IN",
               "53_T21_Ac_NH",
               "19_T33_Ac_WK",
               "47_T31_Ac_JB",
               "57_T32_Ac_NM",
               "35_T43_Ac_MT",
               "54_T42_Ac_JQ") 

# Getting unique DEGs
unique.DEG.annot <- as.data.frame(counts(unique.sig.list)) # make df of sig genes with counts and sample IDs
#unique.DEG.annot$gene <- rownames(unique.DEG.annot) # make column with gene names
#unique.DEG.annot <- merge(unique.DEG.annot, annot, by = "gene") # merge annotation file with df of sig genes by gene id
#unique.DEG.annot <- unique.DEG.annot[!duplicated(unique.DEG.annot$gene),] # remove duplicate rows
#rownames(unique.DEG.annot) <- unique.DEG.annot$gene
#write.csv(unique.DEG.annot, file = "~/Desktop/acerv_sub_unique_DEG_annotated.csv")
list(colnames(unique.DEG.annot))
unique.DEG.annot2 <- unique.DEG.annot[, col.order] # similarly to test above, need to order them so it will correctly group by treatment in plot
#unique.DEG.annot <- unique.DEG.annot[,2:14]
#rownames(df) <- colnames(unique.DEG.annot)

# Creating matrix to use as input for heatmap
mat <- as.matrix(unique.DEG.annot2)
mat <- mat[,col.order] # put in col.order
mat <- mat[order(rownames(mat)),]

#Set colors for treatment
ann_colors <- list(Treatment = c(control="lightpink", Treatment1="darkslategray1", Treatment2="darkslategray3", Treatment3="darkslategray4", Treatment4="darkslategray"))

# Now we need to add GO terms. 
acerv_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv_GOterms.unique.csv", header = TRUE)
acerv_GO <- select(acerv_GO, -X)
colnames(acerv_GO)[1] <-"gene"
# lets try merging acerv_GO terms with unique.DEG.annot2
unique.DEG.annot2$gene <- rownames(unique.DEG.annot2)
test_merge <- merge(unique.DEG.annot2, acerv_GO, by = "gene", all.x = TRUE)

# Now lets merge test_merge by acerv_sub significantly enriched GO terms 
acerv_sig <- read.csv("~/Desktop/acerv_sub_GOterms.unique.csv", header = TRUE)
acerv_sig <- select(acerv_sig, -X)
colnames(acerv_sig)[1] <-"gene"
GO_test_merge <- merge(test_merge, acerv_sig, by = "gene", all.x = TRUE)
# So now we got GO.x and GO.y in this file, which are GO terms that came from the two files. some match and some dont. Why?? Okay I see
# In the acerv_sig file, there were repeats of genes because some genes had multiple GO terms associated. So Go.y is coming from that, so I'll keep that one
# Go.x is repeating
GO_test_merge <- select(GO_test_merge, -GO.ID.x)
colnames(GO_test_merge)[15] <-"category" # renaming GO col so that it can merge with the info generated from GOseq

# Now we shall merge GO_test_merge with the associated GO terms as generated by GOSeq
GO_all <- read.csv("~/Desktop/acerv_sub_GO_ALL.csv", header = TRUE) 
GO_all <- select(GO_all, -X)
GO_all <- merge(GO_test_merge, GO_all, by = "category", all.x = TRUE)
write.csv(GO_all, file = "~/Desktop/acerv_sub_GO_all_gene_ID.csv") # maybe include gene counts too?

# Alright now we got a file that has the GO terms, terms, ontologies, counts, and gene ids
GO_all <- select(GO_all, c("category", "gene", "over_represented_pvalue", "under_represented_pvalue", "numDEInCat", "numInCat", "term", "ontology")) # removing counts to make it easier to look at 
# Now lets merge GO_all with unique.DEG.annot2 by gene
GO_DEG_all <- merge(unique.DEG.annot2, GO_all, by = "gene", all.x = TRUE)
GO_DEG_all <- unique(GO_DEG_all)
GO_DEG_all <- na.omit(GO_DEG_all)







# Filter based on the DEGs in unique.DEG.annot2
keep2 <- rownames(mat)

test <- as.data.frame(GO_DEG_all$gene[GO_DEG_all$gene %in% keep,]) 

# Based on filt info, keep only the genes that pass in acerv_counts_filt
acerv_counts_sub_filt <- as.matrix(acerv_counts_sub[which(rownames(acerv_counts_sub) %in% gn.keep),]) 








GO_DEG_all <- GO_DEG_all[order(GO_DEG_all$over_represented_pvalue),] # need to order them so it will group by treatment in plot
GO_DEG_all <- na.omit(GO_DEG_all)
# 218 rows, but only 217 unique DEGs...
head(GO_DEG_all)
# The 218 is a result of Acerv_evm.TU.Segkk2909_pilon.1 having 2 GO terms that are very similar: GO:0006427 (histidyl-tRNA aminoacylation, BP) and GO:0004821 (histidine-tRNA ligase activity, MF)
# These are quite similar...so I think I'll manually remove the first row. Need to find way to automate this
GO_DEG_all <- GO_DEG_all[-1,]
GO_DEG_all <- GO_DEG_all[order(GO_DEG_all$gene),]


## So now I have all my info I need to make the heatmap! now how do i put it all together 

list <- data.frame(GO_DEG_all$gene, rownames(mat))


# Trying to add terms to left of map
term_df <- data.frame(GO_DEG_all$term)




rownames(term_df) = rownames(mat) # match names 

# Plot
acerv_sub_heatmap <- pheatmap(mat, 
                              annotation_col = df,
                              annotation_colors = ann_colors,
                              scale = "row",
                              show_rownames = T,
                              fontsize_row = 4,
                              cluster_cols = F,
                              show_colnames = T,
                              cluster_rows = F)
acerv_sub_heatmap
#dev.off()
# plot has all treatment comparisons 
ggsave("~/Desktop/mcav_DEGs_heatmap.pdf", mcav_heatmap)
















library(heatmaply)
p <- heatmaply(mat, 
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               scale = "column",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = TRUE,
               branches_lwd = 0.1,
               label_names = c("Country", "Feature:", "Value"),
               fontsize_row = 5, fontsize_col = 5,
               labCol = colnames(mat),
               labRow = rownames(mat),
               heatmap_layers = theme(axis.line=element_blank())
)






