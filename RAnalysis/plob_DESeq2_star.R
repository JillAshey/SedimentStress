# Title: DESeq2 with plob samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. lob only samples analyzed here aligned against P. lutea. STAR was read aligner 

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
library("DataCombine")
library("VennDiagram")


# Load gene count matrix
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/gene_count_plob_only_matrix.csv")
dim(countdata) # 31126 x 17
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub("X", "", colnames(countdata)[col])
}
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub(".fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(countdata)[col])
}
rownames(countdata) <- countdata$gene_id

# Load functional annotation file 
annot <- read.csv("~/Desktop/GFFs/Plut.GFFannotation.fixed_transcript.gff",header = FALSE, sep="\t", skip=6)
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
annot$gene_id <-gsub(";.*", "", annot$attr)
annot$gene_id <-gsub("ID=", "", annot$gene_id)
# subset by gene ? not sure 

# Load plob metadata
metadata_plob <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_plob_raw_filtered.csv")
metadata_plob <- na.omit(metadata_plob)
metadata_plob$SampleID <- gsub("X", "", metadata_plob$SampleID)
rownames(metadata_plob) <- metadata_plob$SampleID

# Subset count data for only plob samples with metadata based on SampleID and make sure rows of metadata = cols of count data
plob_ID <- metadata_plob$SampleID
count_plob <- select(countdata, all_of(plob_ID))
all(rownames(metadata_plob) %in% colnames(count_plob)) # must come out TRUE

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_plob, filt) # create filter for counts data 
keep <- count_plob[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
plob_counts_filt <- as.matrix(count_plob[which(rownames(count_plob) %in% gn.keep),]) 
dim(plob_counts_filt) # 15369 x 12
storage.mode(plob_counts_filt) <- "integer" # stores count data as integer 
write.csv(plob_counts_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_plob) %in% colnames(plob_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_plob$Treatment <- factor(metadata_plob$Treatment, levels = c("control", "mid", "high"))
data <- DESeqDataSetFromMatrix(countData = plob_counts_filt, colData = metadata_plob, design = ~ Treatment)

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
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"

# Compare C vs mid
DEG_control_vs_mid <- results(DEG.int, contrast = c("Treatment", "control", "mid"))
DEG_control_vs_mid
DEG_control_vs_mid.sig.num <- sum(DEG_control_vs_mid$padj <0.05, na.rm = T) # identify # of significant pvalues with p< 0.05
DEG_control_vs_mid.sig.num
# 109 DEGs
DEG_control_vs_mid.sig <- subset(DEG_control_vs_mid, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_mid.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_mid.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_control_vs_mid.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_mid.rsig <- varianceStabilizingTransformation(DEG_control_vs_mid.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_mid.sig.list <- as.data.frame(counts(DEG_control_vs_mid.sig.list))
DEG_control_vs_mid.sig.list["Treatment_Compare"] <- "CvsMid" # adding treatment comparison column
write.csv(DEG_control_vs_mid.sig.list, file = "~/Desktop/plob_control_vs_mid_DEG.csv")

# Compare C vs high
DEG_control_vs_high <- results(DEG.int, contrast = c("Treatment", "control", "high"))
DEG_control_vs_high
DEG_control_vs_high.sig.num <- sum(DEG_control_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with p< 0.05
DEG_control_vs_high.sig.num
# 92 DEGs
DEG_control_vs_high.sig <- subset(DEG_control_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_high.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_control_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_high.rsig <- varianceStabilizingTransformation(DEG_control_vs_high.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_high.sig.list <- as.data.frame(counts(DEG_control_vs_high.sig.list))
DEG_control_vs_high.sig.list["Treatment_Compare"] <- "CvsHigh" # adding treatment comparison column
write.csv(DEG_control_vs_high.sig.list, file = "~/Desktop/plob_control_vs_high_DEG.csv")

# Compare mid vs high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with p< 0.05
DEG_mid_vs_high.sig.num
# 18 DEGs
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_mid_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_mid_vs_high.rsig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list))
DEG_mid_vs_high.sig.list["Treatment_Compare"] <- "MidvsHigh" # adding treatment comparison column
write.csv(DEG_mid_vs_high.sig.list, file = "~/Desktop/plob_mid_vs_high_DEG.csv")

##### Unique genes from intersections of DEG 
DEGs_CvsMid <- as.data.frame(DEG_control_vs_mid.sig.list)
#colnames(DEGs_CvsMid) <- "DEGs"
DEGs_CvsHigh <- as.data.frame(DEG_control_vs_high.sig.list)
#colnames(DEGs_CvsHigh) <- "DEGs"
DEGs_MidvsHigh <- as.data.frame(DEG_mid_vs_high.sig.list)
#colnames(DEGs_MidvsHigh) <- "DEGs"

DEGs.all <- rbind(DEGs_CvsMid, DEGs_CvsHigh, DEGs_MidvsHigh)
write.csv(DEGs.all, file = "~/Desktop/plob_DEGs.all_treatment.csv")
DEGs.all <- unique(DEGs.all) # 153 unique DEGs
unique.sig.num <- length(t(unique(DEGs.all)))

unique.sig.list <- data[which(rownames(data) %in% DEGs.all$DEGs), ] # subset list of sig transcripts from original count data
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(unique.sig.list), file = "~/Desktop/plob_unique.sig.list.csv")

# PCA plot of diff-expressed genes 
plob_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_plob <- round(100*attr(plob_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
plob_DEGPCAplot <- ggplot(plob_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_plob[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_plob[2],"% variance")) +
  scale_color_manual(values = c(control="black", mid = "pink", high = "darkgreen")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
plob_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- plob_DEGPCAplot$data
ggsave("~/Desktop/plob_DEGs_PCA.pdf", plob_DEGPCAplot)

df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
colnames(df) <- "Treatment"
colnames(count_plob)
col.order <- c("6_1",
               "7_1",
               "8_1",
               "9_1",
               "21_1",
               "22_1",
               "23_1",
               "25_1",
               "26_1",
               "27_1",
               "29_1",
               "34_1")
ann_colors <- list(Treatment = c(control="black", mid = "pink", high = "darkgreen"))

# Removing excess and isolating gene name              
unique.DEG.annot <- as.data.frame(counts(unique.sig.list))
unique.DEG.annot$gene_id <- rownames(unique.DEG.annot)
unique.DEG.annot <- merge(unique.DEG.annot, annot, by = "gene_id")
#rownames(unique.DEG.annot) <- unique.DEG.annot$gene_id
write.csv(unique.DEG.annot, file = "~/Desktop/plob_unique_DEG_annotated.csv")

unique.DEG.annot <- unique.DEG.annot[,1:13]
unique.DEG.annot <- unique(unique.DEG.annot)
rownames(unique.DEG.annot)<- unique.DEG.annot$gene_id
unique.DEG.annot <- select(unique.DEG.annot, -gene_id)
rownames(df) <- colnames(unique.DEG.annot)
mat <- as.matrix(unique.DEG.annot)

mat <- mat[,col.order]
#dev.off()
#pdf(file = "~/Desktop/Unique_Heatmap.DEG_Annotated.pdf")
plob_heatmap <- pheatmap(mat, 
                         annotation_col = df,
                         annotation_colors = ann_colors,
                         scale = "row",
                         show_rownames = T,
                         fontsize_row = 4,
                         cluster_cols = T,
                         show_colnames = T)
#dev.off()
# plot has all treatment comparisons 
ggsave("~/Desktop/plob_DEGs_heatmap.pdf", plob_heatmap)














## Going to remove the control from the list and see how DEG PCA looks without it--can look at spread of mid and high 
# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
metadata_plob_treatment <- subset(metadata_plob, Treatment=="mid" | Treatment=="high")
plob_ID_treatment <- metadata_plob_treatment$SampleID
count_plob_treatment <- select(countdata, all_of(plob_ID_treatment))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_plob_treatment, filt) # create filter for counts data 
keep <- count_plob_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
plob_treatment_counts_filt <- as.matrix(count_plob_treatment[which(rownames(count_plob_treatment) %in% gn.keep),]) 
storage.mode(plob_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(plob_treatment_counts_filt, "~/Desktop/plob_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_plob_treatment) %in% colnames(plob_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_plob_treatment$Treatment <- factor(metadata_plob_treatment$Treatment, levels = c("mid", "high"))
data <- DESeqDataSetFromMatrix(countData = plob_treatment_counts_filt, colData = metadata_plob_treatment, design = ~ Treatment)

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
#[1] "Intercept"             "Treatment_high_vs_mid"

# Compare mid and high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_mid_vs_high.sig.num
# 14 DEGs 
# why is it different than when control is included?
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subsey list of significant genes from original count data 
DEG_mid_vs_high.sig.list$contrast <- as_factor(c("Treatment_high_vs_mid")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_mid_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list))
DEG_mid_vs_high.sig.list["Treatment_Compare"] <- "MidsHigh" # adding treatment comparison column
#write.csv(DEG_mid_vs_high.sig.list, file = "~/Desktop/pdam_mid_vs_high_DEG.csv")

#write.csv(DEGs.all, file = "~/Desktop/pdam_DEGs.all_treatment.csv")
DEGs.all <- select(DEG_mid_vs_high.sig.list, -Treatment_Compare)
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all <- unique(DEGs.all) # 4 unique DEGs
#unique.sig.num <- length(t(unique(DEGs.all)))

# PCA plot of diff-expressed genes 
plob_DEGPCAdata <- plotPCA(DEG_mid_vs_high.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_plob <- round(100*attr(plob_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
plob_DEGPCAplot <- ggplot(plob_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_plob[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_plob[2],"% variance")) +
  scale_color_manual(values = c(mid = "pink", high = "darkgreen")) +
  coord_fixed() +
  ggtitle("Plob") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
plob_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- plob_DEGPCAplot$data
ggsave("~/Desktop/plob_treatment_DEGs_PCA.pdf", plob_DEGPCAplot)




























