# Title: DESeq2 with pdam samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. dam only samples analyzed here aligned against P. dam. STAR was read aligner with gff annotation file from NCBI.

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
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/gene_count_pdam_NCBI_matrix.csv")
dim(countdata) # 37630
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub("X", "", colnames(countdata)[col])
}
# Removing all gene_ids that are MSTRG - not helpful as they are novel spliced genes found by STAR. There may be some way to ID the specific genes associated with them, but not sure
countdata <- countdata[grep("LOC", countdata$gene_id), ]
countdata$gene_id <- gsub(".*L","", countdata$gene_id) # the | symbol being really annoying in subsetting, so removing everything up to L, then will add L back to gene_id
countdata$gene_id <- paste0("L", countdata$gene_id)
rownames(countdata) <- countdata$gene_id

# functional annotation file 
annot <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
# annot$gene <- annot$attr
annot <- annot[!grepl("##", annot$scaffold),]
annot$gene_id <- regmatches(annot$attr, gregexpr("(?<=gene=).*", annot$attr, perl = TRUE)) #removing everything in Symbol col up to LOC
annot$gene_id <-gsub(";.*", "", annot$gene_id)
annot <- annot[!grepl("character", annot$gene_id),]

metadata_pdam<- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_pdam_raw_filtered.csv")
metadata_pdam <- na.omit(metadata_pdam)
metadata_pdam$SampleID <- gsub("X", "", metadata_pdam$SampleID)
rownames(metadata_pdam) <- metadata_pdam$SampleID

# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
pdam_ID <- metadata_pdam$SampleID
count_pdam <- select(countdata, all_of(pdam_ID))
all(rownames(metadata_pdam) %in% colnames(count_pdam)) # must come out TRUE

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_pdam, filt) # create filter for counts data 
keep <- count_pdam[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
pdam_counts_filt <- as.matrix(count_pdam[which(rownames(count_pdam) %in% gn.keep),]) 
storage.mode(pdam_counts_filt) <- "integer" # stores count data as integer 
write.csv(pdam_counts_filt, "~/Desktop/pdam_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_pdam) %in% colnames(pdam_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_pdam$Treatment <- factor(metadata_pdam$Treatment, levels = c("control", "mid", "high"))
data <- DESeqDataSetFromMatrix(countData = pdam_counts_filt, colData = metadata_pdam, design = ~ Treatment)

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
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"

# Compare C and mid 
DEG_control_vs_mid <- results(DEG.int, contrast = c("Treatment", "control", "mid"))
DEG_control_vs_mid
DEG_control_vs_mid <- as.data.frame(DEG_control_vs_mid) # make full results into a df
DEG_control_vs_mid["Treatment_Compare"] <- "CvsMid" # add treatment comparison col
write.csv(DEG_control_vs_mid, file = "~/Desktop/pdam_control_vs_mid_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_mid.sig.num <- sum(DEG_control_vs_mid$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_mid.sig.num
# 317 DEGs
DEG_control_vs_mid.sig <- subset(DEG_control_vs_mid, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_mid.sig["Treatment_Compare"] <- "CvsMid" # adding treatment comparison column
DEG_control_vs_mid.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_mid.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_mid.sig.list <- as.data.frame(counts(DEG_control_vs_mid.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_mid.sig.list_full <- cbind(DEG_control_vs_mid.sig, DEG_control_vs_mid.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_mid.sig.list_full, file = "~/Desktop/pdam_control_vs_mid_DEG_full.csv") # write out csv
DEG_control_vs_mid.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_mid.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 





# Compare C and high
DEG_control_vs_high <- results(DEG.int, contrast = c("Treatment", "control", "high"))
DEG_control_vs_high
DEG_control_vs_high <- as.data.frame(DEG_control_vs_high) # make full results into a df
DEG_control_vs_high["Treatment_Compare"] <- "CvsHigh" # add treatment comparison col
write.csv(DEG_control_vs_high, file = "~/Desktop/pdam_control_vs_high_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_high.sig.num <- sum(DEG_control_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_high.sig.num
# 462 DEGs
DEG_control_vs_high.sig <- subset(DEG_control_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_high.sig["Treatment_Compare"] <- "CvsHigh" # adding treatment comparison column
DEG_control_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_high.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_high.sig.list <- as.data.frame(counts(DEG_control_vs_high.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_high.sig.list_full <- cbind(DEG_control_vs_high.sig, DEG_control_vs_high.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_high.sig.list_full, file = "~/Desktop/pdam_control_vs_high_DEG_full.csv") # write out csv
DEG_control_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 






# Compare mid and high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high <- as.data.frame(DEG_mid_vs_high) # make full results into a df
DEG_mid_vs_high["Treatment_Compare"] <- "MidvsHigh" # add treatment comparison col
write.csv(DEG_mid_vs_high, file = "~/Desktop/pdam_mid_vs_high_all_genes.csv") # maybe include gene counts too?
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_mid_vs_high.sig.num
# 8 DEGs
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig["Treatment_Compare"] <- "MidvsHigh" # adding treatment comparison column
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subset list of significant genes from original count data 
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list)) # make list of sig gene counts into a df
DEG_mid_vs_high.sig.list_full <- cbind(DEG_mid_vs_high.sig, DEG_mid_vs_high.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_mid_vs_high.sig.list_full, file = "~/Desktop/pdam_mid_vs_high_DEG_full.csv") # write out csv
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 



##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT3, T2vsT3
DEGs.all <- rbind(DEG_control_vs_mid.sig.list_full, 
                  DEG_control_vs_high.sig.list_full,
                  DEG_mid_vs_high.sig.list_full)
write.csv(DEGs.all, file = "~/Desktop/pdam_DEGs.all_treatment.csv")
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all_pdam <- DEGs.all$DEGs
DEGs.all_pdam <- unique(DEGs.all_pdam)
DEGs.all_pdam <- as.data.frame(DEGs.all_pdam) # 89 unique DEGs among treatment comparisons

unique.sig.list <- data[which(rownames(data) %in% DEGs.all$DEGs), ] # subset list of sig transcripts from original count data
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(unique.sig.list), file = "~/Desktop/pdam_unique.sig.list.csv")


# PCA plot of diff-expressed genes 
pdam_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_pdam <- round(100*attr(pdam_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
pdam_DEGPCAplot <- ggplot(pdam_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_pdam[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_pdam[2],"% variance")) +
  scale_color_manual(values = c(control="black", mid = "pink", high = "darkgreen")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
pdam_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- pdam_DEGPCAplot$data
ggsave("~/Desktop/pdam_DEGs_PCA.pdf", pdam_DEGPCAplot)

df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
colnames(df) <- "Treatment"
colnames(count_pdam)
col.order <- c( "1_2",
                "2_2",
                "4_2",
                "11_2",
                "28_2",
                "35_2",
                "36_2",
                "38_2",
                "39_2",
                "41_2",
                "42_2",
                "47_2")
ann_colors <- list(Treatment = c(control="black", mid = "pink", high = "darkgreen"))

# Removing excess and isolating gene name              
unique.DEG.annot <- as.data.frame(counts(unique.sig.list))
unique.DEG.annot$gene_id <- rownames(unique.DEG.annot)

unique.DEG.annot <- merge(unique.DEG.annot, annot, by = "gene_id")
#rownames(unique.DEG.annot) <- unique.DEG.annot$gene_id
write.csv(unique.DEG.annot, file = "~/Desktop/pdam_unique_DEG_annotated.csv")

unique.DEG.annot <- unique.DEG.annot[,1:13]
unique.DEG.annot <- unique(unique.DEG.annot)
rownames(unique.DEG.annot)<- unique.DEG.annot$gene_id
unique.DEG.annot <- select(unique.DEG.annot, -gene_id)
rownames(df) <- colnames(unique.DEG.annot)
mat <- as.matrix(unique.DEG.annot)

mat <- mat[,col.order]
#dev.off()
#pdf(file = "~/Desktop/Unique_Heatmap.DEG_Annotated.pdf")
pdam_heatmap <- pheatmap(mat, 
                         annotation_col = df,
                         annotation_colors = ann_colors,
                         scale = "row",
                         show_rownames = T,
                         fontsize_row = 4,
                         cluster_cols = T,
                         show_colnames = T)
#dev.off()
# plot has all treatment comparisons 
ggsave("~/Desktop/pdam_DEGs_heatmap.pdf", pdam_heatmap)









## Going to remove the control from the list and see how DEG PCA looks without it--can look at spread of mid and high 
# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
metadata_pdam_treatment <- subset(metadata_pdam, Treatment=="mid" | Treatment=="high")
pdam_ID_treatment <- metadata_pdam_treatment$SampleID
count_pdam_treatment <- select(countdata, all_of(pdam_ID_treatment))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_pdam_treatment, filt) # create filter for counts data 
keep <- count_pdam_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
pdam_treatment_counts_filt <- as.matrix(count_pdam_treatment[which(rownames(count_pdam_treatment) %in% gn.keep),]) 
storage.mode(pdam_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(pdam_treatment_counts_filt, "~/Desktop/pdam_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_pdam_treatment) %in% colnames(pdam_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_pdam_treatment$Treatment <- factor(metadata_pdam_treatment$Treatment, levels = c("mid", "high"))
data <- DESeqDataSetFromMatrix(countData = pdam_treatment_counts_filt, colData = metadata_pdam_treatment, design = ~ Treatment)

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
# 4 DEGs 
# why is it different than when control is included?
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subsey list of significant genes from original count data 
DEG_mid_vs_high.sig.list$contrast <- as_factor(c("Treatment_high_vs_mid")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_mid_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list))
DEG_mid_vs_high.sig.list["Treatment_Compare"] <- "MidsHigh" # adding treatment comparison column
#write.csv(DEG_mid_vs_high.sig.list, file = "~/Desktop/pdam_mid_vs_high_DEG.csv")

#write.csv(DEGs.all, file = "~/Desktop/pdam_DEGs.all_treatment.csv")
DEGs.all <- select(DEG_mid_vs_high.sig.list, -Treatment_Compare)
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all <- unique(DEGs.all) # 4 unique DEGs
#unique.sig.num <- length(t(unique(DEGs.all)))

# PCA plot of diff-expressed genes 
pdam_DEGPCAdata <- plotPCA(DEG_mid_vs_high.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_pdam <- round(100*attr(pdam_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
pdam_DEGPCAplot <- ggplot(pdam_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_pdam[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_pdam[2],"% variance")) +
  scale_color_manual(values = c(mid = "pink", high = "darkgreen")) +
  coord_fixed() +
  ggtitle("Pdam") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
pdam_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- pdam_DEGPCAplot$data
ggsave("~/Desktop/pdam_treatment_DEGs_PCA.pdf", pdam_DEGPCAplot)

















































































