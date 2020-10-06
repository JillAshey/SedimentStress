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
DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1"))
DEG_control_vs_T1
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 5%FDR (padj<0.05) 
DEG_control_vs_T1.sig.num
# 30 DEGs
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T1.sig.list <- as.data.frame(counts(DEG_control_vs_T1.sig.list))
DEG_control_vs_T1.sig.list["Treatment_Compare"] <- "CvsT1" # adding treatment comparison column
write.csv(DEG_control_vs_T1.sig.list, file = "~/Desktop/acerv_control_vs_T1_DEG.csv")

# Compare C vs T2
DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2"))
DEG_control_vs_T2
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj< 0.05
DEG_control_vs_T2.sig.num
# 35 DEGs
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T2.sig.list <- as.data.frame(counts(DEG_control_vs_T2.sig.list))
DEG_control_vs_T2.sig.list["Treatment_Compare"] <- "CvsT2" # adding treatment comparison column
write.csv(DEG_control_vs_T2.sig.list, file = "~/Desktop/acerv_control_vs_T2_DEG.csv")

# Compare C vs T3
DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3"))
DEG_control_vs_T3
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T3.sig.num
# 20 DEGs
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_T3.sig.list <- as.data.frame(counts(DEG_control_vs_T3.sig.list))
DEG_control_vs_T3.sig.list["Treatment_Compare"] <- "CvsT3" # adding treatment comparison column
write.csv(DEG_control_vs_T3.sig.list, file = "~/Desktop/acerv_control_vs_T3_DEG.csv")

# Compare C vs T4
DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4"))
DEG_control_vs_T4
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T4.sig.num
# 3 DEGs
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
DEG_control_vs_T4.sig.list <- as.data.frame(counts(DEG_control_vs_T4.sig.list))
DEG_control_vs_T4.sig.list["Treatment_Compare"] <- "CvsT4" # adding treatment comparison column
write.csv(DEG_control_vs_T4.sig.list, file = "~/Desktop/acerv_control_vs_T4_DEG.csv")

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
# 1 DEGs
DEG_T1_vs_T4.sig <- subset(DEG_T1_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_T1_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_T1_vs_T4.sig)),] # subsey list of significant genes from original count data 
DEG_T1_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_T1_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# Error in estimateDispersionsFit(object, quiet = TRUE, fitType) :   ?????
#   all gene-wise dispersion estimates are within 2 orders of magnitude
# from the minimum value, and so the standard curve fitting techniques will not work.
# One can instead use the gene-wise estimates as final estimates:
#   dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
# ...then continue with testing using nbinomWaldTest or nbinomLRT
DEG_T1_vs_T4.sig.list <- as.data.frame(counts(DEG_T1_vs_T4.sig.list))
DEG_T1_vs_T4.sig.list["Treatment_Compare"] <- "T1vsT4" # adding treatment comparison column
write.csv(DEG_T1_vs_T4.sig.list, file = "~/Desktop/acerv_T1_vs_T4_DEG.csv")

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

# Compare T3 and T4
DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num

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
DEGs_T1vsT4 <- as.data.frame(rownames(DEG_T1_vs_T4.sig.list), DEG_T1_vs_T4.sig.list$Treatment_Compare)
colnames(DEGs_T1vsT4) <- "DEGs"
DEGs_T1vsT4$Treatment_Compare <- rownames(DEGs_T1vsT4)

DEGs.all <- rbind(DEGs_CvsT1,
                        DEGs_CvsT2, 
                        DEGs_CvsT3, 
                        DEGs_CvsT4, 
                        DEGs_T1vsT4) 
write.csv(DEGs.all, file = "~/Desktop/acerv_DEGs.all_treatment.csv")
DEGs.all_acerv <- DEGs.all$DEGs
DEGs.all_acerv <- unique(DEGs.all_acerv)
DEGs.all_acerv <- as.data.frame(DEGs.all_acerv)

unique.sig.list <- data[which(rownames(data) %in% DEGs.all_acerv$DEGs), ] # subset list of sig transcripts from original count data
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(unique.sig.list), file = "~/Desktop/acerv_unique.sig.list.csv")

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
acerv_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_acerv_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(acerv_counts) # 33715 x 15
head(acerv_counts)
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(acerv_counts)[col])
}
for ( col in 1:ncol(acerv_counts)){
  colnames(acerv_counts)[col] <-  gsub("X", "", colnames(acerv_counts)[col])
}
# Remove sample 24, as it may be an outlier to see how analysis goes without it 
acerv_counts <- acerv_counts[, -2] 

# Load metadata
metadata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata.csv", header = TRUE)
dim(metadata) # 45 by 12
head(metadata)
# Selecting only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Acerv species only 
acerv_metadata <- subset(metadata, Species=="Acropora cervicornis")
# Renaming treatments
acerv_metadata$Treatment <- gsub("Ctl", "control", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T1", "Treatment1", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T2", "Treatment2", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T3", "Treatment3", acerv_metadata$Treatment)
acerv_metadata$Treatment <- gsub("T4", "Treatment4", acerv_metadata$Treatment)
# Removing unwanted text from SampleID
acerv_metadata$SampleID <- gsub(".txt.gz", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- gsub(";.*", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- gsub(".fastq.gz", "", acerv_metadata$SampleID)
acerv_metadata$SampleID <- sub("\\.", "", acerv_metadata$SampleID)
# Making sampleID as rownames in metadata 
rownames(acerv_metadata) <- acerv_metadata$SampleID

# Removing sample 23 to see if outlier affects results 
acerv_metadata <- acerv_metadata[!grepl("24_T12_Ac_FM", acerv_metadata$SampleID),]

# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(acerv_metadata) %in% colnames(acerv_counts)) # must come out TRUE


## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(acerv_counts, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- acerv_counts[gfilt,]  
dim(gkeep) # 9055 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
acerv_gcount_filt <- as.matrix(acerv_counts[which(rownames(acerv_counts) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(acerv_gcount_filt)
dim(acerv_gcount_filt) # 9055 x 12 -- only 9055 genes kept after filtering 

# With sample 24 gone
dim(acerv_gcount_filt) # 9447 x 12 -- only 9447 genes kept after filtering 

# Write acerv metadata and counts tables with corrected column and row names and filtered gene counts
write.csv(acerv_metadata, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_acerv_filtered.csv")
write.csv(acerv_gcount_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_acerv_star_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(acerv_metadata) %in% colnames(acerv_gcount_filt))




## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
acerv_metadata$Treatment <- factor(acerv_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
head(acerv_metadata)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_acerv <- DESeqDataSetFromMatrix(countData = acerv_gcount_filt,
                                    colData = acerv_metadata,
                                    design = ~Treatment)
gdds_acerv





## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_acerv <- estimateSizeFactors(gdds_acerv) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_acerv
print(sizeFactors(SF.gdds_acerv)) #view size factors

# size factors all less than 4, can use VST
# vst not working, using rlog 
gvst_acerv <- vst(gdds_acerv, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
# when I run vst(): 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
head(assay(gvst_acerv))
dim(gvst_acerv)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_acerv))) # calculate distance matrix, t returns transpose of assay(gvst_ofav)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_acerv) # assign row names 
colnames(gsampleDistsMatrix) <- NULL # assign col names 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
acerv_heatmap <- pheatmap(gsampleDistsMatrix, # plot matrix
                         clustering_distance_rows = gsampleDists, # cluster rows
                         clustering_distance_cols = gsampleDists, # cluster cols
                         col = colors) # set colors 


## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_acerv, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
acerv_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  #geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() + 
  ggtitle("A. cerv (star)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
acerv_heatmap_PCA <- grid.arrange(mcav_PCAplot, mcav_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/acerv_star_heatmap_PCA.pdf", acerv_heatmap_PCA, width = 8, height = 8, units = c("in"))





## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_acerv <- DESeq(gdds_acerv) #run differential expression test by group using the Wald model
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing
res_DEG_acerv <- results(DEG_acerv)
res_DEG_acerv_Ordered <- res_DEG_acerv[order(res_DEG_acerv$pvalue),]
DEG_acerv$Treatment
resultsNames(DEG_acerv)
# [1] "Intercept"                       "Treatment_Treatment1_vs_control"
# [3] "Treatment_Treatment2_vs_control" "Treatment_Treatment3_vs_control"
# [5] "Treatment_Treatment4_vs_control"

# Explore significant p-values for treatments 
# Control vs treatment1
DEG_acerv_results_control_vs_T1 <- results(DEG_acerv, name = "Treatment_Treatment1_vs_control") # results only for control vs treatment1 treatments 
results_ordered_DEG_acerv_results_control_vs_T1 <- DEG_acerv_results_control_vs_T1[order(DEG_acerv_results_control_vs_T1$pvalue),] # order from smallest pvalue 
summary(DEG_acerv_results_control_vs_T1) # view summary of results with adj p < 0.1
acerv_sig.num.control_vs_T1 <- sum(DEG_acerv_results_control_vs_T1$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 30 significantly differentially expressed genes between control and mid that are less than 0.05
acerv_DEGs.control_vs_T1 <- subset(DEG_acerv_results_control_vs_T1, padj<0.05) # subset only <0.05 padj values
acerv_DEGs.control_vs_T1 <- as.data.frame(acerv_DEGs.control_vs_T1) # make df
acerv_DEGs.control_vs_T1$contrast <- as_factor(c("Treatment_Treatment1_vs_control")) # set contrast as a factor 
acerv_DEGs.control_vs_T1 <- cbind(gene_id = rownames(acerv_DEGs.control_vs_T1), acerv_DEGs.control_vs_T1) # make gene id a row and bind it to the rest of the df
rownames(acerv_DEGs.control_vs_T1) <- NULL # remove row names 
acerv_DEGs.control_vs_T1
dim(acerv_DEGs.control_vs_T1) # 30 by 8
write.csv(acerv_DEGs.control_vs_T1, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEGs.control_vs_T1.csv")
acerv_sig.num.control_vs_T1

# Control vs treatment2
DEG_acerv_results_control_vs_T2 <- results(DEG_acerv, name = "Treatment_Treatment2_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_acerv_results_control_vs_T2 <- DEG_acerv_results_control_vs_T2[order(DEG_acerv_results_control_vs_T2$pvalue),] # order from smallest pvalue 
summary(DEG_acerv_results_control_vs_T2) # view summary of results with adj p < 0.1
acerv_sig.num.control_vs_T2 <- sum(DEG_acerv_results_control_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 35 significantly differentially expressed genes between control and mid that are less than 0.05
acerv_DEGs.control_vs_T2 <- subset(DEG_acerv_results_control_vs_T2, padj<0.05) # subset only <0.05 padj values
acerv_DEGs.control_vs_T2 <- as.data.frame(acerv_DEGs.control_vs_T2) # make df
acerv_DEGs.control_vs_T2$contrast <- as_factor(c("Treatment_Treatment2_vs_control")) # set contrast as a factor 
acerv_DEGs.control_vs_T2 <- cbind(gene_id = rownames(acerv_DEGs.control_vs_T2), acerv_DEGs.control_vs_T2) # make gene id a row and bind it to the rest of the df
rownames(acerv_DEGs.control_vs_T2) <- NULL # remove row names 
acerv_DEGs.control_vs_T2
dim(acerv_DEGs.control_vs_T2) # 35 by 8
write.csv(acerv_DEGs.control_vs_T2, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEGs.control_vs_T2.csv")
acerv_sig.num.control_vs_T2

# Control vs treatment3
DEG_acerv_results_control_vs_T3 <- results(DEG_acerv, name = "Treatment_Treatment3_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_acerv_results_control_vs_T3 <- DEG_acerv_results_control_vs_T3[order(DEG_acerv_results_control_vs_T3$pvalue),] # order from smallest pvalue 
summary(DEG_acerv_results_control_vs_T3) # view summary of results with adj p < 0.1
acerv_sig.num.control_vs_T3 <- sum(DEG_acerv_results_control_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 20 significantly differentially expressed genes between control and mid that are less than 0.05
acerv_DEGs.control_vs_T3 <- subset(DEG_acerv_results_control_vs_T3, padj<0.05) # subset only <0.05 padj values
acerv_DEGs.control_vs_T3 <- as.data.frame(acerv_DEGs.control_vs_T3) # make df
acerv_DEGs.control_vs_T3$contrast <- as_factor(c("Treatment_Treatment3_vs_control")) # set contrast as a factor 
acerv_DEGs.control_vs_T3 <- cbind(gene_id = rownames(acerv_DEGs.control_vs_T3), acerv_DEGs.control_vs_T3) # make gene id a row and bind it to the rest of the df
rownames(acerv_DEGs.control_vs_T3) <- NULL # remove row names 
acerv_DEGs.control_vs_T3
dim(acerv_DEGs.control_vs_T3) # 20 by 8
write.csv(acerv_DEGs.control_vs_T3, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEGs.control_vs_T3.csv")
acerv_sig.num.control_vs_T3

# Control vs treatment4
DEG_acerv_results_control_vs_T4 <- results(DEG_acerv, name = "Treatment_Treatment4_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_acerv_results_control_vs_T4 <- DEG_acerv_results_control_vs_T4[order(DEG_acerv_results_control_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_acerv_results_control_vs_T4) # view summary of results with adj p < 0.1
acerv_sig.num.control_vs_T4 <- sum(DEG_acerv_results_control_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 3 significantly differentially expressed genes between control and mid that are less than 0.05
acerv_DEGs.control_vs_T4 <- subset(DEG_acerv_results_control_vs_T4, padj<0.05) # subset only <0.05 padj values
acerv_DEGs.control_vs_T4 <- as.data.frame(acerv_DEGs.control_vs_T4) # make df
acerv_DEGs.control_vs_T4$contrast <- as_factor(c("Treatment_Treatment4_vs_control")) # set contrast as a factor 
acerv_DEGs.control_vs_T4 <- cbind(gene_id = rownames(acerv_DEGs.control_vs_T4), acerv_DEGs.control_vs_T4) # make gene id a row and bind it to the rest of the df
rownames(acerv_DEGs.control_vs_T4) <- NULL # remove row names 
acerv_DEGs.control_vs_T4
dim(acerv_DEGs.control_vs_T4) # 3 by 8
write.csv(acerv_DEGs.control_vs_T4, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEGs.control_vs_T4.csv")
acerv_sig.num.control_vs_T4

# T1 vs T2
DEG_acerv_results_T1_vs_T2 <- results(DEG_acerv, contrast = c("Treatment", "Treatment1", "Treatment2")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T1_vs_T2 <- order(DEG_acerv_results_T1_vs_T2$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T1_vs_T2) # view summary of results with adj p < 0.1
acerv_sig.num.T1_vs_T2 <- sum(DEG_acerv_results_T1_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T1_vs_T2 # 0 significantly differentially expressed genes between control and mid 
# so do not have to move forward with this, as there are no sig DEGs

# T1 vs T3
DEG_acerv_results_T1_vs_T3 <- results(DEG_acerv, contrast = c("Treatment", "Treatment1", "Treatment3")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T1_vs_T3 <- order(DEG_acerv_results_T1_vs_T3$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T1_vs_T3) # view summary of results with adj p < 0.1
acerv_sig.num.T1_vs_T3 <- sum(DEG_acerv_results_T1_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T1_vs_T3 # 0 significantly differentially expressed genes between control and mid 
# so do not have to move forward with this, as there are no sig DEGs

# T1 vs T4
DEG_acerv_results_T1_vs_T4 <- results(DEG_acerv, contrast = c("Treatment", "Treatment1", "Treatment4")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T1_vs_T4 <- order(DEG_acerv_results_T1_vs_T4$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T1_vs_T4) # view summary of results with adj p < 0.1
acerv_sig.num.T1_vs_T4 <- sum(DEG_acerv_results_T1_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T1_vs_T4 # 1 significantly differentially expressed genes between control and mid 
acerv_DEGs.T1_vs_T4 <- subset(DEG_acerv_results_T1_vs_T4, padj<0.05) # subset only <0.05 padj values
acerv_DEGs.T1_vs_T4 <- as.data.frame(acerv_DEGs.T1_vs_T4) # make df
acerv_DEGs.T1_vs_T4$contrast <- as_factor(c("Treatment_Treatment1_vs_Treatment4")) # set contrast as a factor
acerv_DEGs.T1_vs_T4 <- cbind(gene_id = rownames(acerv_DEGs.T1_vs_T4), acerv_DEGs.T1_vs_T4) # make gene id a row and bind it to the rest of the df
rownames(acerv_DEGs.T1_vs_T4) <- NULL # remove row names 
acerv_DEGs.T1_vs_T4
dim(acerv_DEGs.T1_vs_T4) # 1 by 8
write.csv(acerv_DEGs.T1_vs_T4, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEGs.T1_vs_T4.csv")
acerv_sig.num.T1_vs_T4

# T2 vs T3
DEG_acerv_results_T2_vs_T3 <- results(DEG_acerv, contrast = c("Treatment", "Treatment2", "Treatment3")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T2_vs_T3 <- order(DEG_acerv_results_T2_vs_T3$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T2_vs_T3) # view summary of results with adj p < 0.1
acerv_sig.num.T2_vs_T3 <- sum(DEG_acerv_results_T2_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T2_vs_T3 # 0 significantly differentially expressed genes between T2 and T3 
# so do not have to move forward with this, as there are no sig DEGs

# T2 vs T4
DEG_acerv_results_T2_vs_T4 <- results(DEG_acerv, contrast = c("Treatment", "Treatment2", "Treatment4")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T2_vs_T4 <- order(DEG_acerv_results_T2_vs_T4$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T2_vs_T4) # view summary of results with adj p < 0.1
acerv_sig.num.T2_vs_T4 <- sum(DEG_acerv_results_T2_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T2_vs_T4 # 0 significantly differentially expressed genes between T2 and T3 
# so do not have to move forward with this, as there are no sig DEGs

# T3 vs T4
DEG_acerv_results_T3_vs_T4 <- results(DEG_acerv, contrast = c("Treatment", "Treatment3", "Treatment4")) # results only for control vs mid treatments 
results_ordered_DEG_acerv_results_T3_vs_T4 <- order(DEG_acerv_results_T3_vs_T4$pvalue) #Order p-values by smallest value first
summary(DEG_acerv_results_T3_vs_T4) # view summary of results with adj p < 0.1
acerv_sig.num.T3_vs_T4 <- sum(DEG_acerv_results_T3_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
acerv_sig.num.T3_vs_T4 # 0 significantly differentially expressed genes between T2 and T3 
# so do not have to move forward with this, as there are no sig DEGs

# Combining all DEGs among all treatment comparison
DEGs_acerv_all_treatments <- bind_rows(acerv_DEGs.control_vs_T1, acerv_DEGs.control_vs_T2, acerv_DEGs.control_vs_T3, acerv_DEGs.control_vs_T4, acerv_DEGs.T1_vs_T4) # bind CvsT1, CvsT2, CvsT3, CvsT4, T1vsT4 results together by row - comparisons where there were DEGs
dim(DEGs_acerv_all_treatments) # 89 by 8
DEG_acerv.sg.num <- sum(DEGs_acerv_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_acerv.sg.num # 89 adjusted p-values 
summary(DEGs_acerv_all_treatments)
write.csv(DEGs_acerv_all_treatments, file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/DEGs_acerv_all_treatments.csv")
         


 
## Visualize diff-expressed genes

# Subset and log-transform count data 
# Subset list of genes by those which padj>0.
dim(DEGs_acerv_all_treatments)
DEGs_acerv <- DEGs_acerv_all_treatments$gene_id # list all gene names 
DEGs_acerv <- unique(DEGs_acerv) # select only unique gene names 
DEG_acerv_list <- gdds_acerv[which(rownames(gdds_acerv) %in% DEGs_acerv)] # filter gdds_acerv DESeq2 object by unique gene names
dim(DEG_acerv_list) # 50 x 15
print(counts(DEG_acerv_list))

# As determined above, size factors all less than 4, so proceed with VST
#apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_acerv_list, blind = FALSE, nsub = nrow(counts(DEG_acerv_list)))
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
dim(DEGvst) # 50 x 15
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
# Testing if the first two command in the heatmap is necessary given that we are clustering the columns anyways.
acerv_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_acerv.sg.num) #sort by decreasing sig ?
mat_acerv <- assay(DEGvst)[acerv_topVarGenes, ] #make an expression object
mat_acerv <- mat_acerv - rowMeans(mat_acerv) #diff_pdam in expression compared to average across all samples
dim(mat_acerv)
ann_colors <- list(Treatment= c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange"))
df_DEG_acerv <- as.data.frame(colData(DEGvst)[c("Treatment")]) #make dataframe for column naming and associated treatment

acerv_DEGheatmap <- pheatmap(mat_acerv, scale= "row", legend=TRUE, annotation_legend=TRUE, annotation_col=df_DEG_acerv, annotation_colors = ann_colors,
                            clustering_distance_rows="euclidean", clustering_method = "average",
                            show_rownames =FALSE,
                            show_colnames =TRUE,
                            cluster_cols = TRUE)
acerv_DEGheatmap # this is a crazy map haha

# PCA plot of diff-expressed genes 
acerv_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_acerv <- round(100*attr(acerv_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
acerv_DEGPCAplot <- ggplot(acerv_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_acerv[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_acerv[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
acerv_DEGPCAplot
# PCA plot is of differentially expressed genes only

# Save results
write.csv(counts(DEG_acerv_list), file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/acerv_DEG_list_unique.csv")
acerv_DEGs_heatmap_PCA <- grid.arrange(acerv_DEGPCAplot, acerv_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/acerv_DEGs_heatmap_PCA.pdf", acerv_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))

## Can I make DEG PCAs with individual contrasts? eg control vs T1, T1 vs T4, etc 



