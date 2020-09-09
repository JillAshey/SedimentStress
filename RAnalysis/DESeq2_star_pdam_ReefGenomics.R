# Title: DESeq2 with pdam samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. dam only samples analyzed here aligned against P. dam. STAR was read aligner. ReefGenomics fasta and gff files were used.

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
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_pdam_rgGFF_star_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 1675 x 64
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(countdata)[col])
}

# Load metadata 
metadata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_HI_metadata_raw.csv", header = TRUE)
head(metadata)
# Renaming specific columns
names(metadata)[names(metadata) == "File.Name.fastq"] <- "SampleID" 
names(metadata)[names(metadata) == "Treatment.in.mg.L.of.sediment"] <- "Treatment"
names(metadata)[names(metadata) == "Time.point.in.days"] <- "Days"
# Replacing certain words/charaacters in columns
metadata$Treatment <- gsub("400", "high", metadata$Treatment)
metadata$Treatment <- gsub("40", "mid", metadata$Treatment)
metadata$Treatment <- gsub("0", "control", metadata$Treatment)
metadata$SampleID <- gsub(".fastq.gz", "", metadata$SampleID)
metadata$SampleID <- paste0("X", metadata$SampleID) # add X to front so it matches countdata and isnt treated like a numerical variable 
metadata$Days <- gsub("7", "seven" ,metadata$Days)
metadata$Days <- gsub("4", "four" ,metadata$Days)
# metadata$Treatment <- gsub("<NA>", "unknown", metadata$Treatment) - use if looking at unknown samples
metadata <- metadata %>% drop_na() # remove all rows with NAs - those are the unknown samples
rownames(metadata) <- metadata$SampleID # make sampleID the row names 
# metadata <- metadata[-65,] # remove random blank space at the end
# Select mcap species only 
metadata_pdam <- subset(metadata, Species=="Pocillopora damicornis")

# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
pdam_ID <- metadata_pdam$SampleID
count_pdam <- select(countdata, all_of(pdam_ID))
all(rownames(metadata_pdam) %in% colnames(count_pdam)) # must come out TRUE



## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(count_pdam, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pdam[gfilt,]  
dim(gkeep) # 969 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(count_pdam[which(rownames(count_pdam) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(gcount_filt)
dim(gcount_filt) # 969 x 17 -- only 969 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_pdam, "~/Desktop/metadata_pdam_filtered.csv")
write.csv(gcount_filt, "~/Desktop/genecount_pdam_rgGFF_star_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(metadata_pdam) %in% colnames(gcount_filt))




## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
metadata_pdam$Treatment <- factor(metadata_pdam$Treatment, levels = c("control", "mid", "high"))
head(metadata_pdam)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_pdam <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                    colData = metadata_pdam,
                                    design = ~Treatment)
gdds_pdam




## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_pdam <- estimateSizeFactors(gdds_pdam) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_pdam
print(sizeFactors(SF.gdds_pdam)) #view size factors

# size factors all less than 4, can use VST
# vst not working, using rlog 
gvst_pdam <- rlog(gdds_pdam, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_pdam))
dim(gvst_pdam)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_pdam))) # calculate distance matrix, t returns transpose of assay(gvst_pdam)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_pdam) # assign row names 
colnames(gsampleDistsMatrix) <- NULL # assign col names 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pdam_heatmap <- pheatmap(gsampleDistsMatrix, # plot matrix
                         clustering_distance_rows = gsampleDists, # cluster rows
                         clustering_distance_cols = gsampleDists, # cluster cols
                         col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_pdam, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
pdam_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  # geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen")) +
  coord_fixed() + 
  ggtitle("P. dam (star)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
pdam_heatmap_PCA <- grid.arrange(pdam_PCAplot, pdam_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/pdam_heatmap_star_PCA.pdf", pdam_heatmap_PCA, width = 8, height = 8, units = c("in"))



## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_pdam <- DESeq(gdds_pdam) #run differential expression test by group using the Wald model
res_DEG_pdam <- results(DEG_pdam)
res_DEG_pdam_Ordered <- res_DEG_pdam[order(res_DEG_pdam$pvalue),]

DEG_pdam$Treatment
resultsNames(DEG_pdam)
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"
# Why no mid v high? - bc its a contrast

# Control vs mid
DEG_pdam_results_control_vs_mid <- results(DEG_pdam, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_mid <- DEG_pdam_results_control_vs_mid[order(DEG_pdam_results_control_vs_mid$pvalue),] # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_mid) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_mid <- sum(DEG_pdam_results_control_vs_mid$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 1114 significantly differentially expressed genes between control and mid that are less than 0.05
pdam_DEGs.control_vs_mid <- subset(DEG_pdam_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_mid <- as.data.frame(pdam_DEGs.control_vs_mid) # make df
pdam_DEGs.control_vs_mid$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
pdam_DEGs.control_vs_mid <- cbind(gene_id = rownames(pdam_DEGs.control_vs_mid), pdam_DEGs.control_vs_mid) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_mid) <- NULL # remove row names 
pdam_DEGs.control_vs_mid
dim(pdam_DEGs.control_vs_mid) # 4 by 8
# write.csv(pdam_DEGs.control_vs_mid, "~/Desktop/pdam_DEGs.control_vs_mid.csv")
pdam_sig.num.control_vs_mid

# Control vs high
DEG_pdam_results_control_vs_high <- results(DEG_pdam, name = "Treatment_high_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_high <- order(DEG_pdam_results_control_vs_high$pvalue) # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_high <- sum(DEG_pdam_results_control_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
pdam_sig.num.control_vs_high # 1454 significantly differentially expressed genes between control and mid 
pdam_DEGs.control_vs_high <- subset(DEG_pdam_results_control_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_high <- as.data.frame(pdam_DEGs.control_vs_high) # make df
pdam_DEGs.control_vs_high$contrast <- as_factor(c("control_vs_high")) # set contrast as a factor
pdam_DEGs.control_vs_high <- cbind(gene_id = rownames(pdam_DEGs.control_vs_high), pdam_DEGs.control_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_high) <- NULL # remove row names 
pdam_DEGs.control_vs_high
dim(pdam_DEGs.control_vs_high) # 10 by 8
# write.csv(pdam_DEGs.control_vs_high, "~/Desktop/pdam_DEGs.control_vs_high.csv")
pdam_sig.num.control_vs_high

# Mid vs high
DEG_pdam_results_mid_vs_high <- results(DEG_pdam, contrast = c("Treatment", "mid", "high")) # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_mid_vs_high <- order(DEG_pdam_results_mid_vs_high$pvalue) #Order p-values by smallest value first
summary(DEG_pdam_results_mid_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.mid_vs_high <- sum(DEG_pdam_results_mid_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
pdam_sig.num.mid_vs_high # 33 significantly differentially expressed genes between control and mid 
pdam_DEGs.mid_vs_high <- subset(DEG_pdam_results_mid_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.mid_vs_high <- as.data.frame(pdam_DEGs.mid_vs_high) # make df
pdam_DEGs.mid_vs_high$contrast <- as_factor(c("mid_vs_high")) # set contrast as a factor
pdam_DEGs.mid_vs_high <- cbind(gene_id = rownames(pdam_DEGs.mid_vs_high), pdam_DEGs.mid_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.mid_vs_high) <- NULL # remove row names 
pdam_DEGs.mid_vs_high
dim(pdam_DEGs.mid_vs_high) # 2 by 8
# write.csv(pdam_DEGs.mid_vs_high, "~/Desktop/pdam_DEGs.mid_vs_high.csv")
pdam_sig.num.mid_vs_high

DEGs_pdam_contrast_all_treatments <- bind_rows(pdam_DEGs.control_vs_mid, pdam_DEGs.control_vs_high, pdam_DEGs.mid_vs_high) # bind mid_vs_control results, high_vs_control results, and mid_vs_high results together by row  
dim(DEGs_pdam_contrast_all_treatments) # 16 by 8
DEG_pdam.sg.num <- sum(DEGs_pdam_contrast_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_pdam.sg.num # 16 adjusted p-values 
summary(DEGs_pdam_contrast_all_treatments)
# write.csv(DEGs_pdam_contrast_all_treatments, file="~/Desktop/DEGs_pdam_contrast_all_treatments.csv")


## Visualize diff-expressed genes

# Subset and log-transform count data 
# Subset list of genes by those which padj>0.
dim(DEGs_pdam_contrast_all_treatments)
DEGs_pdam <- DEGs_pdam_contrast_all_treatments$gene_id # list all gene names 
DEGs_pdam <- unique(DEGs_pdam) # select only unique gene names 
DEG_pdam_list <- gdds_pdam[which(rownames(gdds_pdam) %in% DEGs_pdam)] # filter gdds_pdam DESeq2 object by unique gene names
dim(DEG_pdam_list) # 1766 x 12
print(counts(DEG_pdam_list))

# As determined above, size factors all less than 4, so proceed with VST
#apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_pdam_list, blind = FALSE, nsub = nrow(counts(DEG_pdam_list)))
dim(DEGvst) # 11 by 12 
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
# Testing if the first two command in the heatmap is necessary given that we are clustering the columns anyways.
pdam_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_pdam.sg.num) #sort by decreasing sig ?
mat_pdam <- assay(DEGvst)[pdam_topVarGenes, ] #make an expression object
mat_pdam <- mat_pdam - rowMeans(mat_pdam) #diff_pdam in expression compared to average across all samples
dim(mat_pdam)
ann_colors <- list(Treatment= c(control="cadetblue3", mid="palevioletred", high="darkgreen")) 
df_DEG_pdam <- as.data.frame(colData(DEGvst)[c("Treatment")]) #make dataframe for column naming and associated treatment

pdam_DEGheatmap <- pheatmap(mat_pdam, scale= "row", legend=TRUE, annotation_legend=TRUE, annotation_col=df_DEG_pdam, annotation_colors = ann_colors,
                            clustering_distance_rows="euclidean", clustering_method = "average",
                            show_rownames =FALSE,
                            show_colnames =TRUE,
                            cluster_cols = TRUE)
pdam_DEGheatmap # this is a crazy map haha

# PCA plot of diff-expressed genes 
pdam_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_pdam <- round(100*attr(pdam_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
pdam_DEGPCAplot <- ggplot(pdam_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_pdam[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_pdam[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
pdam_DEGPCAplot
# PCA plot is of differentially expressed genes only 


























