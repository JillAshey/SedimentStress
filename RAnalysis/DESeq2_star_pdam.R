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
DEG_control_vs_mid.sig.num <- sum(DEG_control_vs_mid$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_mid.sig.num
# 317 DEGs
DEG_control_vs_mid.sig <- subset(DEG_control_vs_mid, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_mid.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_mid.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_mid.sig.list$contrast <- as_factor(c("Treatment_mid_vs_control")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_control_vs_mid.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_mid.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_mid.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_mid.sig.list <- as.data.frame(counts(DEG_control_vs_mid.sig.list))
DEG_control_vs_mid.sig.list["Treatment_Compare"] <- "CvsMid" # adding treatment comparison column
write.csv(DEG_control_vs_mid.sig.list, file = "~/Desktop/pdam_control_vs_mid_DEG.csv")

# Compare C and high
DEG_control_vs_high <- results(DEG.int, contrast = c("Treatment", "control", "high"))
DEG_control_vs_high
DEG_control_vs_high.sig.num <- sum(DEG_control_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_high.sig.num
# 462 DEGs
DEG_control_vs_high.sig <- subset(DEG_control_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_high.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_high.sig.list$contrast <- as_factor(c("Treatment_high_vs_control")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_control_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_control_vs_high.sig.list <- as.data.frame(counts(DEG_control_vs_high.sig.list))
DEG_control_vs_high.sig.list["Treatment_Compare"] <- "CvsHigh" # adding treatment comparison column
write.csv(DEG_control_vs_high.sig.list, file = "~/Desktop/pdam_control_vs_high_DEG.csv")

# Compare mid and high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_mid_vs_high.sig.num
# 8 DEGs 
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subsey list of significant genes from original count data 
DEG_mid_vs_high.sig.list$contrast <- as_factor(c("Treatment_high_vs_control")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_mid_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list))
DEG_mid_vs_high.sig.list["Treatment_Compare"] <- "MidsHigh" # adding treatment comparison column
write.csv(DEG_mid_vs_high.sig.list, file = "~/Desktop/pdam_mid_vs_high_DEG.csv")

##### Unique genes from intersections of DEG 
DEGs_CvsMid <- as.data.frame(DEG_control_vs_mid.sig.list)
#colnames(DEGs_CvsMid) <- "DEGs"
DEGs_CvsHigh <- as.data.frame(DEG_control_vs_high.sig.list)
#colnames(DEGs_CvsHigh) <- "DEGs"
DEGs_MidvsHigh <- as.data.frame(DEG_mid_vs_high.sig.list)
#colnames(DEGs_MidvsHigh) <- "DEGs"

DEGs.all <- rbind(DEGs_CvsMid, DEGs_CvsHigh, DEGs_MidvsHigh)
write.csv(DEGs.all, file = "~/Desktop/pdam_DEGs.all_treatment.csv")
DEGs.all <- select(DEGs.all, -Treatment_Compare)
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all <- unique(DEGs.all) # 549 unique DEGs
#unique.sig.num <- length(t(unique(DEGs.all)))

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


































































































## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(count_pdam, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pdam[gfilt,]  
dim(gkeep) # 20995 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(count_pdam[which(rownames(count_pdam) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(gcount_filt)
dim(gcount_filt) # 968 x 17 -- only 968 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_pdam, "~/Desktop/metadata_pdam_filtered.csv")
write.csv(gcount_filt, "~/Desktop/genecount_pdam_star_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(metadata_pdam) %in% colnames(gcount_filt))




## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
metadata_pdam$Treatment <- factor(metadata_pdam$Treatment, levels = c("control", "mid", "high", "unknown"))
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
gvst_pdam <- vst(gdds_pdam, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
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
gPCAdata <- plotPCA(gvst_pdam, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
pdam_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape=Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
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
ggsave("~/Desktop/pdam_heatmap_NCBI_star_PCA.pdf", pdam_heatmap_PCA, width = 8, height = 8, units = c("in"))













####### Erin chille code 





## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_pdam <- DESeq(gdds_pdam) #run differential expression test by group using the Wald model
res_DEG_pdam <- results(DEG_pdam)
res_DEG_pdam_Ordered <- res_DEG_pdam[order(res_DEG_pdam$pvalue),]
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing

# # Summary all results 
# #summary is just printing a table for you, you need to tell it what threshold you want
# help("summary",package="DESeq2")
# alpha <- 0.05 #set alpha to 0.05, this will control FDR
# summary(res_DEG_pdam) #default FDR is still 0.1
# summary(res_DEG_pdam, alpha) #no showing all genes with FRD < 0.05
# # To get significant genes, indepdent filtering in results() has alpha argument, usede to optimize a cutoff on mean normalized count
# res_DEG_pdam_05 <- results(DEG_pdam, alpha= alpha) #set FDR to 0.05 now
# res_DEG_pdam_05_Sig <- res_DEG_pdam[which(res_DEG_pdam$padj < alpha),]
# summary(res_DEG_pdam_05) #this is all the genes
# summary(res_DEG_pdam_05_Sig) #this is the significant ones!
# sum(res_DEG_pdam_05_Sig$padj < 0.05, na.rm=TRUE) #singificant genes 
# sig="significant"
# res_DEG_pdam_05_Sig$Significance <- sig
# resS4_05_nonSig <- res_DEG_pdam[which(res_DEG_pdam$padj > alpha),] #create list of nonsig
# nonsig <- "non-significant"
# 
# # Order results tables by adj pvalues
# head( res_DEG_pdam_05_Sig[ order( res_DEG_pdam_05_Sig$log2FoldChange ), ] ) #head for strongest downregulation
# tail( res_DEG_pdam_05_Sig[ order( res_DEG_pdam_05_Sig$log2FoldChange ), ] ) #tail for strongest up regulation





DEG_pdam$Treatment
resultsNames(DEG_pdam)
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"
# Why no mid v high? - bc its a contrast

# Explore significant p-values for treatments 
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
dim(pdam_DEGs.control_vs_mid) # 1114 by 8
write.csv(pdam_DEGs.control_vs_mid, "~/Desktop/pdam_DEGs.control_vs_mid.csv")
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
dim(pdam_DEGs.control_vs_high) # 1454 by 8
write.csv(pdam_DEGs.control_vs_high, "~/Desktop/pdam_DEGs.control_vs_high.csv")
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
dim(pdam_DEGs.mid_vs_high) # 33 by 8
write.csv(pdam_DEGs.mid_vs_high, "~/Desktop/pdam_DEGs.mid_vs_high.csv")
pdam_sig.num.mid_vs_high

DEGs_pdam_contrast_all_treatments <- bind_rows(pdam_DEGs.control_vs_mid, pdam_DEGs.control_vs_high, pdam_DEGs.mid_vs_high) # bind mid_vs_control results, high_vs_control results, and mid_vs_high results together by row  
dim(DEGs_pdam_contrast_all_treatments) # 2601 by 8
DEG_pdam.sg.num <- sum(DEGs_pdam_contrast_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_pdam.sg.num # 2601 adjusted p-values 
summary(DEGs_pdam_contrast_all_treatments)
write.csv(DEGs_pdam_contrast_all_treatments, file="~/Desktop/DEGs_pdam_contrast_all_treatments.csv")
# I thought the values in this file would be the numbers that went into the venn diagram but I guess not 



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
dim(DEGvst) # 1176 by 12 
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

# Save results
# Save CSV of differentially expressed genes, heatmap, and PCAplot
write.csv(counts(DEG_pdam_list), file="~/Desktop/pdam_DEGlist.csv")

pdam_DEGs_heatmap_PCA <- grid.arrange(pdam_DEGPCAplot, pdam_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/pdam_DEGs_heatmap_PCA.pdf", pdam_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))



## Venn diagrams of pairwise comparison DEGs in pdam
pdam_control_vs_mid_venn <- as.data.frame(pdam_DEGs.control_vs_mid$gene_id) # list of genes in control_vs_mid
colnames(pdam_control_vs_mid_venn) <- c("gene_id") # label col names 
pdam_control_vs_mid_venn$contrast <- as.character(c("control vs mid")) # set contrast as a character
head(pdam_control_vs_mid_venn)

pdam_control_vs_high_venn <- as.data.frame(pdam_DEGs.control_vs_high$gene_id) # list of genes in control_vs_high
colnames(pdam_control_vs_high_venn) <- c("gene_id") # label col names
pdam_control_vs_high_venn$contrast <- as.character(c("control vs high")) # set contrast as character
head(pdam_control_vs_high_venn)

pdam_mid_vs_high_venn <- as.data.frame(pdam_DEGs.mid_vs_high$gene_id) # list of genes in mid_vs_high
colnames(pdam_mid_vs_high_venn) <- c("gene_id") # label col names
pdam_mid_vs_high_venn$contrast <- as.character(c("mid vs high")) # set contrast as character
head(pdam_mid_vs_high_venn)

# Bind all venn dataframes together
pdam_DEG_all_venn <- bind_rows(pdam_control_vs_mid_venn, pdam_control_vs_high_venn, pdam_mid_vs_high_venn)
head(pdam_DEG_all_venn)

# Prepare a palette of 3 colors with R colorbrewer:
ann_colors <- c("cadetblue3", "palevioletred","darkgreen") 

# Build venn diagram 
venn.diagram(
  x = list(
    pdam_DEG_all_venn %>% filter(contrast=="control vs mid") %>% select(gene_id) %>% unlist() , 
    pdam_DEG_all_venn %>% filter(contrast=="control vs high") %>% select(gene_id) %>% unlist() , 
    pdam_DEG_all_venn %>% filter(contrast=="mid vs high") %>% select(gene_id) %>% unlist()
  ),
  category.names = c("control vs mid" , "control vs high" , "mid vs high"),
  filename = '~/Desktop/pdam_DEGs_venn.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = ann_colors,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
# not sure if this looks correct, but it is interesting 
















# NEED TO WORK ON THIS AS OF 8/10/20
### Gene ontology analysis of planula genes

# Prepare data
# Obtain names of all expressed genes (poverA = 0.85,5), and all differentially expressed genes (p<0.05)
gcount_filt <- gcount_filt %>% as.data.frame()
gcount_filt <- cbind(gene_id = rownames(gcount_filt), gcount_filt)
rownames(gcount_filt) <- NULL
dim(gcount_filt) # 21369 by 13
head(gcount_filt)

pdam_DEG_IDs <- as.data.frame(counts(DEG_pdam_list))
dim(pdam_DEG_IDs) # 1136 by 12
pdam_DEG_IDs <- cbind(gene_id = rownames(pdam_DEG_IDs), pdam_DEG_IDs)
rownames(pdam_DEG_IDs) <- NULL
colnames(pdam_DEG_IDs) <-colnames(gcount_filt)
head(pdam_DEG_IDs)
dim(pdam_DEG_IDs)
head(gcount_filt)
dim(gcount_filt)


# Import merged annotation file - ONLY pdam merged in this file 
map <- read.csv(file="~/Desktop/stringTie_subset_pdam_merged.gtf", header=FALSE, sep="\t", skip=2) #load sample info
map <- subset(map, V3=="transcript")
map <- map[,c(1,4,5,9)]
map <- separate(map, V9, into = c("gene_id", "transcript_id", "gene_name"), sep=";")
map$gene_id <- gsub("gene_id ","",map$gene_id) #remove extra characters
map$gene_id <- gsub(" ","",map$gene_id) #remove extra characters
map$transcript_id <- gsub("transcript_id ","",map$transcript_id) #remove extra characters
map$transcript_id <- gsub(" ","",map$transcript_id) #remove extra characters
map$gene_name <- gsub("gene_name ","",map$gene_name) #remove extra characters
# map$gene_name <- gsub("gene_name ","",map$gene_name) #remove extra characters
map$gene_name <- gsub(" ","",map$gene_name) #remove extra characters
colnames(map) <- c("scaffold", "start", "stop", "gene_id", "transcript_id", "gene_name")
write.csv(map, file="~/Desktop/map_annotated_.csv")
dim(map) # 51884 by 6
head(map) # supposed to be filtering by gene id or name comparing with gcount_filt, but not working well



dim(gcount_filt) # 21369 by 13

# Not sure if the below commented lines are correct...need to filter the map by gcounts
# gcount_filt_sep <- separate(gcount_filt, gene_id, into = c("gene_id", "gene_name"))
# head(gcount_filt_sep)

# map_pdam <- filter(map, gene_name %in% gcount_filt_sep$gene_id)
# erin: map_pln <- filter(map, gene_name %in% gcount_filt_pln$gene_id) #Should be 25,358
# dim(map_pdam)

# ######## Vingette code with deseq2 objects from above 
DEG_pdam # DESeq2 object 
# Get results of DESEq2 object 
res_DEG_pdam <- results(DEG_pdam)
# order by smallest pvalues
resOrdered_DEG_pdam <- DEG_pdam[order(DEG_pdam$pvalue),]










# 
# ######## Vingette code 
# ### Above, I loaded in data, made sure cols of countdata == rows of metadata
# 
# # Make count data matrix for DESeq2
# cts <- as.matrix(countdata_pdam)
# head(cts)
# # Make Treatment a factor in metadata
# metadata_pdam$Treatment <- factor(metadata_pdam$Treatment)
# 
# # contruct DESEq2 data set with metadata and gene counts
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = metadata_pdam, 
#                               design = ~ Treatment)
# dds
# 
# # Prefilter from dds object 
# keep <- rowSums(counts(dds) >= 10)
# dds <- dds[keep,]
# 
# # set reference level so that DESeq2 knows what to compare against 
# dds$Treatment <- relevel(dds$Treatment, ref = "control")
# 
# # Run DESeq analysis and get results
# dds <- DESeq(dds)
# res <- results(dds)
# res
# 
# # Specify coefficient or contrast to build results table for 
# res_mid_vs_control <- results(dds, name = "Treatment_mid_vs_control")
# res_mid_vs_high <- results(dds, contrast = c("Treatment", "mid", "high"))
# # difference between using name vs contrast in results? Maybe use contrast only when comparing two experimental treatment groups (as opposed comparing to control)
# 
# # Log fold change - essentially seeing how the gene changed in comparison between two treatments 
# # Shrinkage of effect size (LFC estimates - quantitative measure of magnitude of some phenomenon) hrlps with visualizing and ranking of genes 
# # To shrink, we give dds object and name of coefficients to shrink 
# resultsNames(dds) # high vs mid not here because it is a contrast
# resLFC_mid_vs_control <- lfcShrink(dds, coef = "Treatment_mid_vs_control", type = "apeglm") # using apeglm method from Zhu, A., Ibrahim, J.G., Love, M.I. (2018)
# resLFC_mid_vs_control # only chaning the LFC in object
# 
# ## Looking at p-values to evaluate significance 
# # Order results table by smallest pvalue 
# resOrdered <- res[order(res$pvalue),]
# resOrdered # lowest pvalue around 0.03
# sum(res$padj < 0.1, na.rm=TRUE)
# 
# ## Mine / Connelly cod
# # Convert countdata to matrix 
# countdata_pdam_matrix <- as.matrix(countdata_pdam)
# dim(countdata_pdam_matrix)
# 
# # Check to make sure rownames in metadata match colnames in matrix
# all(rownames(metadata_pdam) %in% colnames(countdata_pdam_matrix))
# 
# # Create deseq2 data object 
# dds_pdam_sediment <- DESeqDataSetFromMatrix(countData = countdata_pdam_matrix,
#                                             colData = metadata_pdam,
#                                             design =~Treatment)
# # will give this error: some variables in design formula are characters, converting to factors
# relevel(dds_pdam_sediment$Treatment, ref = "control")
# head(dds_pdam_sediment)
# 
# # Perform DESeq2 analysis 
# dds_pdam_sediment <- DESeq(dds_pdam_sediment)
# dds_pdam_sediment
# dim(dds_pdam_sediment)
# as_data_frame(colData(dds_pdam_sediment))
# 
# # Remove low counts
# keep <- rowSums(counts(dds_pdam_sediment), na.rm = TRUE) >= 10
# dds_pdam_sediment <- dds_pdam_sediment[keep,]
# dds_pdam_sediment
# dim(dds_pdam_sediment)
# 
# # Obtain DESeq2 results 
# res_dds_pdam_sediment <- results(dds_pdam_sediment)
# head(res_dds_pdam_sediment)
# 
# summary(res_dds_pdam_sediment)
# 
# plotDispEsts(dds_pdam_sediment)
# plotMA(res_dds_pdam_sediment, ylim = c(-10, 10))
# summary(res_dds_pdam_sediment)
# 
# ## Visualize gene count data
# 
# # Log transform count data
# # First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# # Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# #Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# # To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# # In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
# 
# 
# SF.dds_pdam_sediment <- estimateSizeFactors( dds_pdam_sediment ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
# print(sizeFactors(SF.dds_pdam_sediment))
# 
# # Size factors all less than 4, so we can use VST
# gvst_pdam_sediment <- vst(dds_pdam_sediment, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
# head(assay(gvst_pdam_sediment))
# 
# # Plot heatmap of sample-to-sample distances 
# gsampleDists_pdam_sediment <- dist(t(assay(gvst_pdam_sediment))) #calculate distance matix
# gsampleDistMatrix_pdam_sediment <- as.matrix(gsampleDists_pdam_sediment) #distance matrix
# rownames(gsampleDistMatrix_pdam_sediment) <- colnames(gvst_pdam_sediment) #assign row names
# colnames(gsampleDistMatrix_pdam_sediment) <- NULL #assign col names
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistMatrix_pdam_sediment, #plot matrix
#          clustering_distance_rows=gsampleDists_pdam_sediment, #cluster rows
#          clustering_distance_cols=gsampleDists_pdam_sediment, #cluster columns
#          col=colors) #set colors
# 
# # PCA plot of samples
# gPCAdata_pdam_sediment <- plotPCA(gvst_pdam_sediment, intgroup = c("condition"), returnData=TRUE)
# 
# rld <- rlog(dds_pdam_sediment)
# plotPCA(rld)
# colData(dds_pdam_sediment)
# 













