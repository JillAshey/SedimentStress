# Title: DESeq2 with only known plob samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. lobata only samples analyzed here aligned against P. lutea. HISAT2 was read aligner.

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

## Obtaining and tidying data 

# Read in count data
countdata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/hisat2/gene_count_plob_hisat2_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 31126 x 64
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata)[col])
}

# Load metadata 
metadata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Data/sediment_metadata_allSamples_raw.csv", header = TRUE)
head(metadata)
# Replacing certain words/charaacters in columns
names(metadata)[names(metadata) == "File.Name.fastq"] <- "SampleID"
names(metadata)[names(metadata) == "Treatment.in.mg.L.of.sediment"] <- "Treatment"
names(metadata)[names(metadata) == "Time.point.in.days"] <- "Days"
metadata$Treatment <- gsub("400", "high", metadata$Treatment)
metadata$Treatment <- gsub("40", "mid", metadata$Treatment)
metadata$Treatment <- gsub("0", "control", metadata$Treatment)
metadata$SampleID <- gsub(".fastq.gz", "", metadata$SampleID)
metadata$SampleID <- paste0("X", metadata$SampleID) # add X to front so it matches countdata and isnt treated like a numerical variable 
metadata$Days <- gsub("7", "seven" ,metadata$Days)
metadata$Days <- gsub("4", "four" ,metadata$Days)
rownames(metadata) <- metadata$SampleID # make sample ids the row names
metadata <- metadata[-65,] # remove weird blank space
# Select plob species only 
metadata_plob <- subset(metadata, Species=="Porites lobata")

# Subset count data for only pcomp samples based on SampleID and make sure rows of metadata = cols of count data
plob_ID <- metadata_plob$SampleID
count_plob <- select(countdata, all_of(plob_ID))
all(rownames(metadata_plob) %in% colnames(count_plob)) 

## done before all plob samples were identified
# countdata_plob <- read.csv("/Users/jillashey/Desktop/gene_count_plob.csv", header=TRUE,row.names = "gene_id") 
# dim(countdata_plob) # 31126 x 12
# head(countdata_plob)
# # Take out extra characters in column headers 
# for ( col in 1:ncol(countdata_plob)){
#   colnames(countdata_plob)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata_plob)[col])
# }
# 
# # Load metadata for pdam
# metadata_plob <- read.csv("/Users/jillashey/Desktop/metadata_plob.csv", header=TRUE) 
# head(metadata_plob)
# 
# # Rename treatment column
# names(metadata_plob)[names(metadata_plob) == "Treatment.in.mg.L.of.sediment"] <- "Treatment"
# 
# # In treatment col, rename 0=control, 40=mid, 400=high
# metadata_plob$Treatment[metadata_plob$Treatment == "0"] <- "control"
# metadata_plob$Treatment[metadata_plob$Treatment == "40"] <- "mid"
# metadata_plob$Treatment[metadata_plob$Treatment == "400"] <- "high"
# head(metadata_plob)
# 
# # Remove fastq.gz from filename column
# metadata_plob$File.Name.fastq <- gsub(".fastq.gz", "", metadata_plob$File.Name.fastq)
# head(metadata_plob)
# 
# # Add X in front of lines in filename column 
# metadata_plob$File.Name.fastq <- paste("X", metadata_plob$File.Name.fastq, sep = "")
# head(metadata_plob)
# 
# # Make filename column the row names in metadata file
# rownames(metadata_plob) <- metadata_plob$File.Name.fastq
# 
# # Check to make sure rownames in metadata file match column names in countdata file - must return TRUE
# all(rownames(metadata_plob) %in% colnames(countdata_plob)) 




## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(count_plob, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_plob[gfilt,]  
dim(gkeep) # 12504

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(count_plob[which(rownames(count_plob) %in% gn.keep),])
head(gcount_filt)
dim(gcount_filt) # 12504 x 16 -- only 12054 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_plob, "~/Desktop/metadata_plob_raw_filtered.csv")
write.csv(gcount_filt, "~/Desktop/gene_count_plob_hisat2_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(metadata_plob) %in% colnames(gcount_filt))




## Construct DESeq2 dataset 

# Set group as a factor and give levels 
metadata_plob$Treatment <- factor(metadata_plob$Treatment, levels = c("control", "mid", "high", "unknown")) # unknown is treatment bc not sure what metadata for some samples are 
head(metadata_plob)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_plob <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                    colData = metadata_plob,
                                    design = ~Treatment)
gdds_plob




## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_plob <- estimateSizeFactors(gdds_plob) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_plob
print(sizeFactors(SF.gdds_plob)) #view size factors

# size factors all less than 4, can use VST
gvst_plob <- vst(gdds_plob, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_plob))
dim(gvst_plob)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_plob))) # calculate distance matrix, t returns transpose of assay(gvst_pdam)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_plob) # assign row names 
colnames(gsampleDistsMatrix) <- NULL # assign col names 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
plob_heatmap <- pheatmap(gsampleDistsMatrix, # plot matrix
         clustering_distance_rows = gsampleDists, # cluster rows
         clustering_distance_cols = gsampleDists, # cluster cols
         col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_plob, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
plob_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape =Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() +
  ggtitle("P. lobata (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
plob_heatmap_PCA <- grid.arrange(plob_PCAplot, plob_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/plob_heatmap_hisat2_PCA.pdf", plob_heatmap_PCA, width = 8, height = 8, units = c("in"))


















## Run DGE analysis

# After looking at plots to make sure all genes are looking okay, run DE analysis 
# Use Wald model 
DEG_plob <- DESeq(gdds_plob) #run differential expression test by group using the Wald model
res_DEG_plob <- results(DEG_plob)
res_DEG_plob_Ordered <- res_DEG_plob[order(res_DEG_plob$pvalue),]
head(res_DEG_plob_Ordered)

# # Summary all results - all DGE, regardless of treatment at the moment 
# #summary is just printing a table for you, you need to tell it what threshold you want
# help("summary",package="DESeq2")
# alpha <- 0.05 #set alpha to 0.05, this will control FDR
# summary(res_DEG_plob) #default FDR is still 0.1
# summary(res_DEG_plob, alpha) #no showing all genes with FRD < 0.05
# # To get significant genes, indepdent filtering in results() has alpha argument, usede to optimize a cutoff on mean normalized count
# res_DEG_plob_05 <- results(DEG_plob, alpha= alpha) #set FDR to 0.05 now
# res_DEG_plob_05_Sig <- res_DEG_plob[which(res_DEG_plob$padj < alpha),]
# summary(res_DEG_plob_05_Sig) #this is the significant ones!
# sum(res_DEG_plob_05_Sig$padj < 0.05, na.rm=TRUE) #singificant genes
# sig="significant"
# res_DEG_plob_05_Sig$Significance <- sig
# resS4_05_nonSig <- res_DEG_plob[which(res_DEG_plob$padj > alpha),] #create list of nonsig
# nonsig <- "non-significant"
# # Order results tables by adj pvalues
# head( res_DEG_plob_05_Sig[ order( res_DEG_plob_05_Sig$log2FoldChange ), ] ) #head for strongest downregulation
# tail( res_DEG_plob_05_Sig[ order( res_DEG_plob_05_Sig$log2FoldChange ), ] ) #tail for strongest up regulation


## Obtain results names in order to compare by treatment 
DEG_plob$Treatment
resultsNames(DEG_plob)
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"
# Why no mid v high? - bc its a contrast


## Explore significant p-values for treatments 
# Control vs mid
DEG_plob_results_control_vs_mid <- results(DEG_plob, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_plob_results_control_vs_mid <- DEG_plob_results_control_vs_mid[order(DEG_plob_results_control_vs_mid$pvalue),] # order from smallest pvalue 
summary(DEG_plob_results_control_vs_mid) # view summary of results with adj p < 0.1
plob_sig.num.control_vs_mid <- sum(DEG_plob_results_control_vs_mid$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 281 significantly differentially expressed genes between control and mid that are less than 0.05
plob_DEGs.control_vs_mid <- subset(DEG_plob_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
plob_DEGs.control_vs_mid <- as.data.frame(plob_DEGs.control_vs_mid) # make df
plob_DEGs.control_vs_mid$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
plob_DEGs.control_vs_mid <- cbind(gene_id = rownames(plob_DEGs.control_vs_mid), plob_DEGs.control_vs_mid) # make gene id a row and bind it to the rest of the df
rownames(plob_DEGs.control_vs_mid) <- NULL # remove row names 
plob_DEGs.control_vs_mid
dim(plob_DEGs.control_vs_mid) # 281 by 8
write.csv(plob_DEGs.control_vs_mid, "~/Desktop/plob_DEGs.control_vs_mid.csv")
plob_sig.num.control_vs_mid

# Control vs high
DEG_plob_results_control_vs_high <- results(DEG_plob, name = "Treatment_high_vs_control") # results only for control vs high treatments 
results_ordered_DEG_plob_results_control_vs_high <- DEG_plob_results_control_vs_high[order(DEG_plob_results_control_vs_high$pvalue),] # order from smallest pvalue 
summary(DEG_plob_results_control_vs_high) # view summary of results with adj p < 0.1
plob_sig.num.control_vs_high <- sum(DEG_plob_results_control_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 226 significantly differentially expressed genes between control and mid that are less than 0.05
plob_DEGs.control_vs_high <- subset(DEG_plob_results_control_vs_high, padj<0.05) # subset only <0.05 padj values
plob_DEGs.control_vs_high <- as.data.frame(plob_DEGs.control_vs_high) # make df
plob_DEGs.control_vs_high$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
plob_DEGs.control_vs_high <- cbind(gene_id = rownames(plob_DEGs.control_vs_high), plob_DEGs.control_vs_high) # make gene id a row and bind it to the rest of the df
rownames(plob_DEGs.control_vs_high) <- NULL # remove row names 
plob_DEGs.control_vs_high
dim(plob_DEGs.control_vs_high) # 226 by 8
write.csv(plob_DEGs.control_vs_high, "~/Desktop/plob_DEGs.control_vs_high.csv")
plob_sig.num.control_vs_high

# Mid vs high
DEG_plob_results_mid_vs_high <- results(DEG_plob, contrast = c("Treatment", "mid", "high")) # results only for control vs mid treatments 
results_ordered_DEG_plob_results_mid_vs_high <- order(DEG_plob_results_mid_vs_high$pvalue) #Order p-values by smallest value first
summary(DEG_plob_results_mid_vs_high) # view summary of results with adj p < 0.1
plob_sig.num.mid_vs_high <- sum(DEG_plob_results_mid_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
plob_sig.num.mid_vs_high # 6 significantly differentially expressed genes between control and mid 
plob_DEGs.mid_vs_high <- subset(DEG_plob_results_mid_vs_high, padj<0.05) # subset only <0.05 padj values
plob_DEGs.mid_vs_high <- as.data.frame(plob_DEGs.mid_vs_high) # make df
plob_DEGs.mid_vs_high$contrast <- as_factor(c("mid_vs_high")) # set contrast as a factor
plob_DEGs.mid_vs_high <- cbind(gene_id = rownames(plob_DEGs.mid_vs_high), plob_DEGs.mid_vs_high) # make gene id a row and bind it to the rest of the df
rownames(plob_DEGs.mid_vs_high) <- NULL # remove row names 
plob_DEGs.mid_vs_high
dim(plob_DEGs.mid_vs_high) # 6 by 8
write.csv(plob_DEGs.mid_vs_high, "~/Desktop/plob_DEGs.mid_vs_high.csv")
plob_sig.num.mid_vs_high

# Combine treatment comparisons together 
DEGs_plob_contrast_all_treatments <- bind_rows(plob_DEGs.control_vs_mid, plob_DEGs.control_vs_high, plob_DEGs.mid_vs_high) # bind mid_vs_control results, high_vs_control results, and mid_vs_high results together by row  
dim(DEGs_plob_contrast_all_treatments) # 513 by 8
DEG_plob.sg.num <- sum(DEGs_plob_contrast_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_plob.sg.num # 513 adjusted p-values 
summary(DEGs_plob_contrast_all_treatments)
write.csv(DEGs_plob_contrast_all_treatments, file="~/Desktop/DEGs_plob_contrast_all_treatments.csv")




## Visualize diff-expressed genes

# get rid of redundant genes 
dim(DEGs_plob_contrast_all_treatments)
DEGs_plob <- DEGs_plob_contrast_all_treatments$gene_id # list all gene names 
DEGs_plob <- unique(DEGs_plob) # select only unique gene names 
DEG_plob_list <- gdds_plob[which(rownames(gdds_plob) %in% DEGs_plob)] # filter gdds_plob DESeq2 object by unique gene names
dim(DEG_plob_list) # 369 x 12
DEG_plob_list_save <- print(counts(DEG_plob_list))
write.csv(DEG_plob_list_save, file="~/Desktop/DEG_plob_list_alltreatments_save.csv")


# As determined above, size factors all less than 4, so proceed with VST
# apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_plob_list, blind = FALSE, nsub = nrow(counts(DEG_plob_list)))
dim(DEGvst) # 369 by 12 
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
plob_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_plob.sg.num) #sort counts by decreasing sig ?
mat_plob <- assay(DEGvst)[plob_topVarGenes, ] #make an expression object
mat_plob <- mat_plob - rowMeans(mat_plob) #diff_plob in expression compared to average across all samples
dim(mat_plob)
ann_colors <- list(Treatment= c(control="cadetblue3", mid="palevioletred", high="darkgreen")) 
df_DEG_plob <- as.data.frame(colData(DEGvst)[c("Treatment")]) #make dataframe for column naming and associated treatment

plob_DEGheatmap <- pheatmap(mat_plob, scale= "row", legend=TRUE, annotation_legend=TRUE, annotation_col=df_DEG_plob, annotation_colors = ann_colors,
                            clustering_distance_rows="euclidean", clustering_method = "average",
                            show_rownames =FALSE,
                            show_colnames =TRUE,
                            cluster_cols = TRUE)
plob_DEGheatmap # this is a crazy map haha

# PCA plot of diff-expressed genes 
plob_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_plob <- round(100*attr(plob_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
plob_DEGPCAplot <- ggplot(plob_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_plob[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_plob[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
plob_DEGPCAplot
# PCA plot is of differentially expressed genes only 

# Save results
plob_DEGs_heatmap_PCA <- grid.arrange(plob_DEGPCAplot, plob_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/plob_DEGs_heatmap_PCA.pdf", plob_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))









