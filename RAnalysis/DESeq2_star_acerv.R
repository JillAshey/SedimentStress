# Title: DESeq2 with Acerv samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/16/20

# Code for Francois sedimentation data. A. cerv only samples analyzed here aligned against A. cerv. STAR was read aligner with gff annotation file from Baums lab (pers. comm.)

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



