# Title: DESeq2 with Ofav samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/16/20

# Code for Francois sedimentation data. O. fav only samples analyzed here aligned against O. fav. STAR was read aligner with gff annotation file from NCBI

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
ofav_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/gene_count_ofav_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(ofav_counts) # 30180 x 15
for ( col in 1:ncol(ofav_counts)){
  colnames(ofav_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(ofav_counts)[col])
}
for ( col in 1:ncol(ofav_counts)){
  colnames(ofav_counts)[col] <-  gsub("X", "", colnames(ofav_counts)[col])
}
# In Hollie's code, she read a functional annotation file in here 
annot <- read.csv("~/Desktop/GFFs/GCF_002042975.1_ofav_dov_v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
# annot$gene <- annot$attr
annot <- annot[!grepl("##", annot$scaffold),]
annot$gene <-gsub(";.*", "", annot$attr)
annot$gene <-gsub("ID=", "", annot$gene)
# regmatches(annot$gene, gregexpr("(?<=gene=).*", annot$gene, perl = TRUE))
# annot$gene <- gsub(";.*", "", annot$gene) # removing everything after LOC term

# Load metadata
metadata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata.csv", header = TRUE)
dim(metadata) # 45 by 12
head(metadata)
# Selecting only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Acerv species only 
ofav_metadata <- subset(metadata, Species=="Obicella faveolata")
# Renaming treatments
ofav_metadata$Treatment <- gsub("Ctl", "control", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T1", "Treatment1", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T2", "Treatment2", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T3", "Treatment3", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T4", "Treatment4", ofav_metadata$Treatment)
# Removing unwanted text from SampleID
ofav_metadata$SampleID <- gsub(".txt.gz", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- gsub(";.*", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- gsub(".fastq.gz", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- sub("\\.", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- sub("X", "", ofav_metadata$SampleID)
# Making sampleID as rownames in metadata 
rownames(ofav_metadata) <- ofav_metadata$SampleID

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(ofav_counts, filt) # create filter for counts data 
keep <- ofav_counts[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
ofav_counts_filt <- as.matrix(ofav_counts[which(rownames(ofav_counts) %in% gn.keep),]) 
storage.mode(ofav_counts_filt) <- "integer" # stores count data as integer 
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(ofav_metadata) %in% colnames(ofav_counts_filt)) # must come out TRUE
# Set Treatment as a factor
ofav_metadata$Treatment <- factor(ofav_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
data <- DESeqDataSetFromMatrix(countData = ofav_counts_filt, colData = ofav_metadata, design = ~ Treatment)

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
DEG.int.res <- results(DEG.int) # save DE results ; why does it say 'Wald test p-value: Treatment Treatment4 vs control' for DEG.int.res? Is it only looking at treatment 4 and control? In DESeq object created above, it says that design is Treatment
resultsNames(DEG.int) # view DE results 

DEG_control_vs_T1 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment1"))
DEG_control_vs_T1
DEG_control_vs_T1.sig.num <- sum(DEG_control_vs_T1$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T1.sig.num
# 8 DEGs
DEG_control_vs_T1.sig <- subset(DEG_control_vs_T1, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T1.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T1.sig)),] # subsey list of significant genes from original count data 
DEG_control_vs_T1.sig.list$contrast <- as_factor(c("Treatment1_vs_control")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_control_vs_T1.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_T1.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T1.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning message:
#   In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#               Estimated rdf < 1.0; not estimating variance
write.csv(counts(DEG_control_vs_T1.sig.list), file = "~/Desktop/ofav_control_vs_T1_DEG.csv")

DEG_control_vs_T2 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment2"))
DEG_control_vs_T2
DEG_control_vs_T2.sig.num <- sum(DEG_control_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T2.sig.num
# 8 DEGs
DEG_control_vs_T2.sig <- subset(DEG_control_vs_T2, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T2.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T2.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_control_vs_T2.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_T2.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T2.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning message:
#   In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#               Estimated rdf < 1.0; not estimating variance
write.csv(counts(DEG_control_vs_T2.sig.list), file = "~/Desktop/ofav_control_vs_T2_DEG.csv")

DEG_control_vs_T3 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment3"))
DEG_control_vs_T3
DEG_control_vs_T3.sig.num <- sum(DEG_control_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T3.sig.num
# 12 DEGs
DEG_control_vs_T3.sig <- subset(DEG_control_vs_T3, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T3.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T3.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_control_vs_T3.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_T3.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T3.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(DEG_control_vs_T3.sig.list), file = "~/Desktop/ofav_control_vs_T3_DEG.csv")

DEG_control_vs_T4 <- results(DEG.int, contrast = c("Treatment", "control", "Treatment4"))
DEG_control_vs_T4
DEG_control_vs_T4.sig.num <- sum(DEG_control_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_T4.sig.num
# 10 DEGs
DEG_control_vs_T4.sig <- subset(DEG_control_vs_T4, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_T4.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_T4.sig)),] # subsey list of significant genes from original count data 
SFtest <- estimateSizeFactors(DEG_control_vs_T4.sig.list)
print(sizeFactors(SFtest))
DEG_control_vs_T4.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_T4.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# Warning message:
#   In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
#               Estimated rdf < 1.0; not estimating variance
write.csv(counts(DEG_control_vs_T4.sig.list), file = "~/Desktop/acerv_control_vs_T4_DEG.csv")

DEG_T1_vs_T2 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment2"))
DEG_T1_vs_T2
DEG_T1_vs_T2.sig.num <- sum(DEG_T1_vs_T2$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T2.sig.num
# 0 DEGs

DEG_T1_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment3"))
DEG_T1_vs_T3
DEG_T1_vs_T3.sig.num <- sum(DEG_T1_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T3.sig.num
# 0 DEGs

DEG_T1_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment1", "Treatment4"))
DEG_T1_vs_T4
DEG_T1_vs_T4.sig.num <- sum(DEG_T1_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T1_vs_T4.sig.num
# 0 DEGs

DEG_T2_vs_T3 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment3"))
DEG_T2_vs_T3
DEG_T2_vs_T3.sig.num <- sum(DEG_T2_vs_T3$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T3.sig.num
# 0 DEGs

DEG_T2_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment2", "Treatment4"))
DEG_T2_vs_T4
DEG_T2_vs_T4.sig.num <- sum(DEG_T2_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T2_vs_T4.sig.num
# 0 DEGs

DEG_T3_vs_T4 <- results(DEG.int, contrast = c("Treatment", "Treatment3", "Treatment4"))
DEG_T3_vs_T4
DEG_T3_vs_T4.sig.num <- sum(DEG_T3_vs_T4$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_T3_vs_T4.sig.num
# 0 DEGs

##### Unique genes from intersections of DEG in CvsT1, CvsT2, CvsT3, CvsT4, T1vsT2, T1vsT3, T1vsT4, T2vsT3, T2vsT4, T3vsT4
DEGs_CvsT1 <- as.data.frame(rownames(DEG_control_vs_T1.sig.list))
colnames(DEGs_CvsT1) <- "DEGs"
DEGs_CvsT2 <- as.data.frame(rownames(DEG_control_vs_T2.sig.list))
colnames(DEGs_CvsT2) <- "DEGs"
DEGs_CvsT3 <- as.data.frame(rownames(DEG_control_vs_T3.sig.list))
colnames(DEGs_CvsT3) <- "DEGs"
DEGs_CvsT4 <- as.data.frame(rownames(DEG_control_vs_T4.sig.list))
colnames(DEGs_CvsT4) <- "DEGs"

DEGs.all <- rbind(DEGs_CvsT1, DEGs_CvsT2, DEGs_CvsT3, DEGs_CvsT4)
DEGs.all <- unique(DEGs.all)
#DEGs.all$DEGs <- gsub("\\|.*", "", DEGs.all$DEGs)
#DEGs.all$DEGs <- gsub(".*-", "", DEGs.all$DEGs)


# unique.sig.num <- length(t(unique(DEGs.all)))

unique.sig.list <- data[which(rownames(data) %in% DEGs.all$DEGs), ] # subset list of sig transcripts from original count data
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# - note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(unique.sig.list), file = "~/Desktop/ofav_unique.sig.list.csv")

# PCA.plot <- plotPCA(unique.rsig, intgroup = "Treatment") # plot PCA of all samples for DEG only 
# PCA.plot
# PC.info <- PCA.plot$data

# PCA plot of diff-expressed genes 
ofav_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_ofav <- round(100*attr(ofav_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
ofav_DEGPCAplot <- ggplot(ofav_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_ofav[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_ofav[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="skyblue1", Treatment2="skyblue2", Treatment3="skyblue3", Treatment4="skyblue4")) +
  #scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
ofav_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- ofav_DEGPCAplot$data
ggsave("~/Desktop/ofav_DEGs_PCA.pdf", ofav_DEGPCAplot)

df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
colnames(df) <- "Treatment"
colnames(ofav_counts)
col.order <- c("17_ctl2_Of_ZTH_1",
               "18_T33_Of_VLL",
               "23_ctl1_Of_CT_1",
               "26_T12_Of_WCL",
               "30_T23_Of_RPG",
               "32_T22_Of_EVR",
               "36_T43_Of_JJN",
               "40_T13_Of_GWS",
               "43_ctl3_Of_JVP_1",
               "44_T41_Of_PVT_1",
               "48_T31_Of_JNO",
               "50_T21_Of_YZB",
                "51_T42_Of_UOF",
               "59_T11_Of_TQP",
               "60_T32_Of_WY")

# Removing excess and isolating gene name              
unique.DEG.annot <- as.data.frame(counts(unique.sig.list))
unique.DEG.annot$gene <- rownames(unique.DEG.annot)
unique.DEG.annot$gene <- gsub("\\|.*", "", unique.DEG.annot$gene)

unique.DEG.annot <- merge(unique.DEG.annot, annot, by = "gene")
# unique.DEG.annot <- unique.DEG.annot[!duplicated(unique.DEG.annot$gene),]
rownames(unique.DEG.annot) <- unique.DEG.annot$gene
write.csv(unique.DEG.annot, file = "~/Desktop/ofav_unique_DEG_annotated.csv")

unique.DEG.annot <- unique.DEG.annot[,2:16]
rownames(df) <- colnames(unique.DEG.annot)
# unique.DEG.annot <- unique.DEG.annot[,-16]
mat <- as.matrix(unique.DEG.annot)

mat <- mat[,col.order]
#dev.off()
#pdf(file = "~/Desktop/Unique_Heatmap.DEG_Annotated.pdf")
ofav_heatmap <- pheatmap(mat, 
         annotation_col = df,
         annotation_colors = ann_colors,
         scale = "row",
         show_rownames = T,
         fontsize_row = 4,
         cluster_cols = T,
         show_colnames = T)
#dev.off()
# plot has all treatment comparisons 
ggsave("~/Desktop/ofav_DEGs_heatmap.pdf", ofav_heatmap)


## Connecting IPS annotations to full annotations
# Read in ofav annot file
annot <- annot[!grepl("##", annot$scaffold),]
ofav_GO <- filter(annot, grepl("XP_", attr)) # only want rows with proteins
ofav_GO$prot <-gsub(";.*", "", ofav_GO$attr)
ofav_GO$prot <-gsub(".*-", "", ofav_GO$prot)

# read in interproscan file 
ofav_IPS <- read.csv("~/Desktop/ofav.interpro.gff3",header = FALSE, sep="\t", skip=4)
length(unique(ofav_IPS$V1)) # 1245302
colnames(ofav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
ofav_IPS_GO <- filter(ofav_IPS, grepl("GO:", attr)) # select only rows with GO terms

# merge annot and interproscan file by protein
ofav_merge <- merge(ofav_GO, ofav_IPS_GO, by = "prot")
ofav_merge <- na.omit(ofav_merge)

# subset bu Pfam predictor
pfam_ofav <- subset(ofav_merge, Predict == "Pfam")
# Isolate the gene id
## need to manually extract gene
pfam_ofav$gene <- regmatches(pfam_ofav$attr.x, gregexpr("(?<=gene=).*", pfam_ofav$attr.x, perl = TRUE)) #removing everything up to LOC
pfam_ofav$gene <- gsub(";.*", "", pfam_ofav$gene) # removing everything after LOC term

# Use gene id to merge with DEGs file with full annot file -- will give final annotation of DEGs with GO terms, etc
colnames(DEGs.all) <- "gene"
DEGs.all$gene <- gsub("\\|.*", "", DEGs.all$gene)
DEGs.all$gene <- gsub(".*-", "", DEGs.all$gene)

ofav_full_annot <- merge(DEGs.all, pfam_ofav, by = "gene", all.x = TRUE)
write.csv(ofav_full_annot, file = "~/Desktop/ofav_full_annot.csv")

































































































































############################################################################################################################################################


# Load gene count matrix
ofav_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_ofav_only_matrix.csv", header = TRUE, row.names = "gene_id")
dim(ofav_counts) # 30180 x 15
head(ofav_counts)
for ( col in 1:ncol(ofav_counts)){
  colnames(ofav_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(ofav_counts)[col])
}
for ( col in 1:ncol(ofav_counts)){
  colnames(ofav_counts)[col] <-  gsub("X", "", colnames(ofav_counts)[col])
}
# what would happen if I removed sample 23? It's a potential outlier based on its alignment to the Ofav genome 
ofav_counts <- ofav_counts[, -3] 

# Load metadata
metadata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata.csv", header = TRUE)
dim(metadata) # 45 by 11
head(metadata)
# Selecting only the columns I need for analyses 
metadata <- select(metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Select Acerv species only 
ofav_metadata <- subset(metadata, Species=="Obicella faveolata")
# Renaming treatments
ofav_metadata$Treatment <- gsub("Ctl", "control", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T1", "Treatment1", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T2", "Treatment2", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T3", "Treatment3", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T4", "Treatment4", ofav_metadata$Treatment)
# Removing unwanted text from SampleID
ofav_metadata$SampleID <- gsub(".txt.gz", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- gsub(";.*", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- gsub(".fastq.gz", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- sub("\\.", "", ofav_metadata$SampleID)
ofav_metadata$SampleID <- sub("X", "", ofav_metadata$SampleID)

# Removing sample 23 to see if outlier affects results 
ofav_metadata <- ofav_metadata[!grepl("23_ctl1_Of_CT_1", ofav_metadata$SampleID),]

# Making sampleID as rownames in metadata 
rownames(ofav_metadata) <- ofav_metadata$SampleID

# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(ofav_metadata) %in% colnames(ofav_counts)) # must come out TRUE





## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(ofav_counts, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- ofav_counts[gfilt,]  
dim(gkeep) # 18815 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
ofav_gcount_filt <- as.data.frame(ofav_counts[which(rownames(ofav_counts) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(ofav_gcount_filt)
dim(ofav_gcount_filt) # 18815 x 15 -- only 18815 genes kept after filtering 
# without sample23:
dim(ofav_gcount_filt) # 19704 x 14 -- only 19704 genes kept after filtering 


# Write acerv metadata and counts tables with corrected column and row names and filtered gene counts
write.csv(ofav_metadata, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_ofav_filtered.csv")
write.csv(ofav_gcount_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_ofav_star_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(ofav_metadata) %in% colnames(ofav_gcount_filt))





## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
ofav_metadata$Treatment <- factor(ofav_metadata$Treatment, levels = c("control", "Treatment1", "Treatment2", "Treatment3", "Treatment4"))
head(ofav_metadata)

ofav_gcount_matrix <- as.matrix(ofav_gcount_filt)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_ofav <- DESeqDataSetFromMatrix(countData = ofav_gcount_matrix,
                                    colData = ofav_metadata,
                                    design = ~Treatment)
gdds_ofav





## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_ofav <- estimateSizeFactors(gdds_ofav) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_ofav
print(sizeFactors(SF.gdds_ofav)) #view size factors

# size factors all less than 4, can use VST
# vst not working, using rlog 
gvst_ofav <- vst(gdds_ofav, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_ofav))
dim(gvst_ofav)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_ofav))) # calculate distance matrix, t returns transpose of assay(gvst_ofav)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_ofav) # assign row names 
colnames(gsampleDistsMatrix) <- NULL # assign col names 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
ofav_heatmap <- pheatmap(gsampleDistsMatrix, # plot matrix
                         clustering_distance_rows = gsampleDists, # cluster rows
                         clustering_distance_cols = gsampleDists, # cluster cols
                         col = colors) # set colors 
## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_ofav, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
ofav_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  #geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() + 
  ggtitle("O. fav (star)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
ofav_heatmap_PCA <- grid.arrange(ofav_PCAplot, ofav_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/ofav_star_heatmap_PCA.pdf", ofav_heatmap_PCA, width = 8, height = 8, units = c("in"))
# without sample 23
ggsave("~/Desktop/ofav_star_heatmap_PCA_outlier.pdf", ofav_heatmap_PCA, width = 8, height = 8, units = c("in"))






## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_ofav <- DESeq(gdds_ofav) #run differential expression test by group using the Wald model
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
res_DEG_ofav <- results(DEG_ofav)
res_DEG_ofav_Ordered <- res_DEG_ofav[order(res_DEG_ofav$pvalue),]
DEG_ofav$Treatment
resultsNames(DEG_ofav)
# [1] "Intercept"                       "Treatment_Treatment1_vs_control"
# [3] "Treatment_Treatment2_vs_control" "Treatment_Treatment3_vs_control"
# [5] "Treatment_Treatment4_vs_control"

# Explore significant p-values for treatments 
# Control vs treatment1
DEG_ofav_results_control_vs_T1 <- results(DEG_ofav, name = "Treatment_Treatment1_vs_control") # results only for control vs treatment1 treatments 
results_ordered_DEG_ofav_results_control_vs_T1 <- DEG_ofav_results_control_vs_T1[order(DEG_ofav_results_control_vs_T1$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_control_vs_T1) # view summary of results with adj p < 0.1
ofav_sig.num.control_vs_T1 <- sum(DEG_ofav_results_control_vs_T1$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 8 significantly differentially expressed genes between control and T1 that are less than 0.05
# When sample 23 removed, only 3 significantly differentially expressed genes between control and T1 that are less than 0.05
ofav_DEGs.control_vs_T1 <- subset(DEG_ofav_results_control_vs_T1, padj<0.05) # subset only <0.05 padj values
ofav_DEGs.control_vs_T1 <- as.data.frame(ofav_DEGs.control_vs_T1) # make df
ofav_DEGs.control_vs_T1$contrast <- as_factor(c("Treatment_Treatment1_vs_control")) # set contrast as a factor 
ofav_DEGs.control_vs_T1 <- cbind(gene_id = rownames(ofav_DEGs.control_vs_T1), ofav_DEGs.control_vs_T1) # make gene id a row and bind it to the rest of the df
rownames(ofav_DEGs.control_vs_T1) <- NULL # remove row names 
ofav_DEGs.control_vs_T1
dim(ofav_DEGs.control_vs_T1) # 8 by 8
write.csv(ofav_DEGs.control_vs_T1, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/ofav_DEGs.control_vs_T1.csv")
ofav_sig.num.control_vs_T1

# Control vs treatment2
DEG_ofav_results_control_vs_T2 <- results(DEG_ofav, name = "Treatment_Treatment2_vs_control") # results only for control vs treatment2 treatments 
results_ordered_DEG_ofav_results_control_vs_T2 <- DEG_ofav_results_control_vs_T2[order(DEG_ofav_results_control_vs_T2$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_control_vs_T2) # view summary of results with adj p < 0.1
ofav_sig.num.control_vs_T2 <- sum(DEG_ofav_results_control_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 8 significantly differentially expressed genes between control and T1 that are less than 0.05
# When sample 23 removed, only 2 significantly differentially expressed genes between control and T1 that are less than 0.05
ofav_DEGs.control_vs_T2 <- subset(DEG_ofav_results_control_vs_T2, padj<0.05) # subset only <0.05 padj values
ofav_DEGs.control_vs_T2 <- as.data.frame(ofav_DEGs.control_vs_T2) # make df
ofav_DEGs.control_vs_T2$contrast <- as_factor(c("Treatment_Treatment2_vs_control")) # set contrast as a factor 
ofav_DEGs.control_vs_T2 <- cbind(gene_id = rownames(ofav_DEGs.control_vs_T2), ofav_DEGs.control_vs_T2) # make gene id a row and bind it to the rest of the df
rownames(ofav_DEGs.control_vs_T2) <- NULL # remove row names 
ofav_DEGs.control_vs_T2
dim(ofav_DEGs.control_vs_T2) # 8 by 8
write.csv(ofav_DEGs.control_vs_T2, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/ofav_DEGs.control_vs_T2.csv")
ofav_sig.num.control_vs_T2

# Control vs treatment3
DEG_ofav_results_control_vs_T3 <- results(DEG_ofav, name = "Treatment_Treatment3_vs_control") # results only for control vs treatment3 treatments 
results_ordered_DEG_ofav_results_control_vs_T3 <- DEG_ofav_results_control_vs_T3[order(DEG_ofav_results_control_vs_T3$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_control_vs_T3) # view summary of results with adj p < 0.1
ofav_sig.num.control_vs_T3 <- sum(DEG_ofav_results_control_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 12 significantly differentially expressed genes between control and T1 that are less than 0.05
ofav_DEGs.control_vs_T3 <- subset(DEG_ofav_results_control_vs_T3, padj<0.05) # subset only <0.05 padj values
ofav_DEGs.control_vs_T3 <- as.data.frame(ofav_DEGs.control_vs_T3) # make df
ofav_DEGs.control_vs_T3$contrast <- as_factor(c("Treatment_Treatment3_vs_control")) # set contrast as a factor 
ofav_DEGs.control_vs_T3 <- cbind(gene_id = rownames(ofav_DEGs.control_vs_T3), ofav_DEGs.control_vs_T3) # make gene id a row and bind it to the rest of the df
rownames(ofav_DEGs.control_vs_T3) <- NULL # remove row names 
ofav_DEGs.control_vs_T3
dim(ofav_DEGs.control_vs_T3) # 12 by 8
write.csv(ofav_DEGs.control_vs_T3, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/ofav_DEGs.control_vs_T3.csv")
ofav_sig.num.control_vs_T3

# Control vs treatment4
DEG_ofav_results_control_vs_T4 <- results(DEG_ofav, name = "Treatment_Treatment4_vs_control") # results only for control vs treatment4 treatments 
results_ordered_DEG_ofav_results_control_vs_T4 <- DEG_ofav_results_control_vs_T4[order(DEG_ofav_results_control_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_control_vs_T4) # view summary of results with adj p < 0.1
ofav_sig.num.control_vs_T4 <- sum(DEG_ofav_results_control_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 10 significantly differentially expressed genes between control and T1 that are less than 0.05
ofav_DEGs.control_vs_T4 <- subset(DEG_ofav_results_control_vs_T4, padj<0.05) # subset only <0.05 padj values
ofav_DEGs.control_vs_T4 <- as.data.frame(ofav_DEGs.control_vs_T4) # make df
ofav_DEGs.control_vs_T4$contrast <- as_factor(c("Treatment_Treatment4_vs_control")) # set contrast as a factor 
ofav_DEGs.control_vs_T4 <- cbind(gene_id = rownames(ofav_DEGs.control_vs_T4), ofav_DEGs.control_vs_T4) # make gene id a row and bind it to the rest of the df
rownames(ofav_DEGs.control_vs_T4) <- NULL # remove row names 
ofav_DEGs.control_vs_T4
dim(ofav_DEGs.control_vs_T4) # 10 by 8
write.csv(ofav_DEGs.control_vs_T4, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/ofav_DEGs.control_vs_T4.csv")
ofav_sig.num.control_vs_T4

# treatment1 vs treatment2
DEG_ofav_results_T1_vs_T2 <- results(DEG_ofav, contrast = c("Treatment", "Treatment1", "Treatment2")) # results only for treatment1 vs treatment2 treatments 
results_ordered_DEG_ofav_results_T1_vs_T2 <- DEG_ofav_results_T1_vs_T2[order(DEG_ofav_results_T1_vs_T2$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T1_vs_T2) # view summary of results with adj p < 0.1
ofav_sig.num.T1_vs_T2 <- sum(DEG_ofav_results_T1_vs_T2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between T1 and T2 that are less than 0.05

# treatment1 vs treatment3
DEG_ofav_results_T1_vs_T3 <- results(DEG_ofav, contrast = c("Treatment", "Treatment1", "Treatment3")) # results only for treatment1 vs treatment3 treatments 
results_ordered_DEG_ofav_results_T1_vs_T3 <- DEG_ofav_results_T1_vs_T3[order(DEG_ofav_results_T1_vs_T3$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T1_vs_T3) # view summary of results with adj p < 0.1
ofav_sig.num.T1_vs_T3 <- sum(DEG_ofav_results_T1_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between control and T1 that are less than 0.05

# treatment1 vs treatment4
DEG_ofav_results_T1_vs_T4 <- results(DEG_ofav, contrast = c("Treatment", "Treatment1", "Treatment4")) # results only for treatment1 vs treatment4 treatments 
results_ordered_DEG_ofav_results_T1_vs_T4 <- DEG_ofav_results_T1_vs_T4[order(DEG_ofav_results_T1_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T1_vs_T4) # view summary of results with adj p < 0.1
ofav_sig.num.T1_vs_T4 <- sum(DEG_ofav_results_T1_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between control and T1 that are less than 0.05
# weird, seems like the T1 comparisons are the same results with T2, T3, T4

# treatment2 vs treatment3
DEG_ofav_results_T2_vs_T3 <- results(DEG_ofav, contrast = c("Treatment", "Treatment2", "Treatment3")) # results only for treatment2 vs treatment3 treatments 
results_ordered_DEG_ofav_results_T2_vs_T3 <- DEG_ofav_results_T2_vs_T3[order(DEG_ofav_results_T2_vs_T3$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T2_vs_T3) # view summary of results with adj p < 0.1
ofav_sig.num.T2_vs_T3 <- sum(DEG_ofav_results_T2_vs_T3$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between control and T1 that are less than 0.05

# treatment2 vs treatment4
DEG_ofav_results_T2_vs_T4 <- results(DEG_ofav, contrast = c("Treatment", "Treatment2", "Treatment4")) # results only for treatment2 vs treatment4 treatments 
results_ordered_DEG_ofav_results_T2_vs_T4 <- DEG_ofav_results_T2_vs_T4[order(DEG_ofav_results_T2_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T2_vs_T4) # view summary of results with adj p < 0.1
ofav_sig.num.T2_vs_T4 <- sum(DEG_ofav_results_T2_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between control and T1 that are less than 0.05
## ???? same result still...not sure what's going on 

# treatment3 vs treatment4
DEG_ofav_results_T3_vs_T4 <- results(DEG_ofav, contrast = c("Treatment", "Treatment3", "Treatment4")) # results only for treatment3 vs treatment4 treatments 
results_ordered_DEG_ofav_results_T3_vs_T4 <- DEG_ofav_results_T3_vs_T4[order(DEG_ofav_results_T3_vs_T4$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_T3_vs_T4) # view summary of results with adj p < 0.1
ofav_sig.num.T3_vs_T4 <- sum(DEG_ofav_results_T3_vs_T4$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 0 significantly differentially expressed genes between T3 and T4 that are less than 0.05
# okay weird. these results are all the same. i need to look back at these DEGs

# Combining all DEGs among all treatment comparison
DEGs_ofav_all_treatments <- bind_rows(ofav_DEGs.control_vs_T1, ofav_DEGs.control_vs_T2, ofav_DEGs.control_vs_T3, ofav_DEGs.control_vs_T4) # bind CvsT1, CvsT2, CvsT3, CvsT4 results together by row - comparisons where there were DEGs
dim(DEGs_ofav_all_treatments) # 38 by 8
DEG_ofav.sg.num <- sum(DEGs_ofav_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_ofav.sg.num # 38 adjusted p-values 
summary(DEGs_ofav_all_treatments)
write.csv(DEGs_ofav_all_treatments, file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/DEGs_ofav_all_treatments.csv")





## Visualize diff-expressed genes

# Subset and log-transform count data 
# Subset list of genes by those which padj>0.
dim(DEGs_ofav_all_treatments)
DEGs_ofav <- DEGs_ofav_all_treatments$gene_id # list all gene names 
DEGs_ofav <- unique(DEGs_ofav) # select only unique gene names 
DEG_ofav_list <- gdds_ofav[which(rownames(gdds_ofav) %in% DEGs_ofav)] # filter gdds_acerv DESeq2 object by unique gene names
dim(DEG_ofav_list) # 15 x 15
print(counts(DEG_ofav_list))

# As determined above, size factors all less than 4, so proceed with VST
#apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_ofav_list, blind = FALSE, nsub = nrow(counts(DEG_ofav_list)))
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
dim(DEGvst) # 15 x 15
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
# Testing if the first two command in the heatmap is necessary given that we are clustering the columns anyways.
ofav_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_ofav.sg.num) #sort by decreasing sig ?
mat_ofav <- assay(DEGvst)[ofav_topVarGenes, ] #make an expression object
mat_ofav <- mat_ofav - rowMeans(mat_ofav) #diff_ofav in expression compared to average across all samples
dim(mat_ofav)
ann_colors <- list(Treatment= c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange"))
df_DEG_ofav <- as.data.frame(colData(DEGvst)[c("Treatment")]) #make dataframe for column naming and associated treatment
ofav_DEGheatmap <- pheatmap(mat_ofav, scale= "row", legend=TRUE, annotation_legend=TRUE, annotation_col=df_DEG_ofav, annotation_colors = ann_colors,
                            clustering_distance_rows="euclidean", clustering_method = "average",
                            show_rownames =FALSE,
                            show_colnames =TRUE,
                            cluster_cols = TRUE)
ofav_DEGheatmap 

# PCA plot of diff-expressed genes 
ofav_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_ofav <- round(100*attr(ofav_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
ofav_DEGPCAplot <- ggplot(ofav_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_ofav[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_ofav[2],"% variance")) +
  scale_color_manual(values = c(control="black", Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
ofav_DEGPCAplot
# PCA plot is of differentially expressed genes only

# Save results
write.csv(counts(DEG_ofav_list), file="~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/ofav_DEG_list_unique.csv")
ofav_DEGs_heatmap_PCA <- grid.arrange(ofav_DEGPCAplot, ofav_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/ofav_DEGs_heatmap_PCA.pdf", ofav_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))

## Can I make DEG PCAs with individual contrasts? eg control vs T1, T1 vs T4, etc 




