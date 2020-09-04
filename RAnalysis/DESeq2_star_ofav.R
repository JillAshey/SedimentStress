# Title: DESeq2 with ofav samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/03/20

# Code for Francois sedimentation data in FL. Ofav samples aligned against Ofav genome. STAR was read aligner with gff annotation file from NCBI.

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

# Load gene count matrix for Ofav
countdata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_ofav_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 30180 by 33
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(countdata)[col])
}
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub("X", "", colnames(countdata)[col])
}

# Load metadata
ofav_metadata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_FL_metadata_raw.csv", header = TRUE)
dim(ofav_metadata) # 45 by 12
head(ofav_metadata)
# Removing columns I don't need for the analyses 
ofav_metadata <- select(ofav_metadata, c(Rep, Species, Treatment.in.mg.L.of.sediment, Location, File.Name.fastq))
# Renaming cols
colnames(ofav_metadata) <-c("Replicate","Species", "Treatment", "Location", "SampleID")
# Renaming treatments
ofav_metadata$Treatment <- gsub("Ctl", "control", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T1", "Treatment1", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T2", "Treatment2", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T3", "Treatment3", ofav_metadata$Treatment)
ofav_metadata$Treatment <- gsub("T4", "Treatment4", ofav_metadata$Treatment)
# Removing txt.gz from SampleID
ofav_metadata$SampleID <- gsub(".txt.gz", "", ofav_metadata$SampleID)
# Select ofav species only 
ofav_metadata <- subset(ofav_metadata, Species=="Obicella faveolata")
# For the time being, remove samples that I don't have data for 
ofav_metadata <- ofav_metadata %>% drop_na()

# Select only Ofav samples in count data
id <- ofav_metadata$SampleID
id # as reference
col <- colnames(countdata)
col # as reference 
ofav_countdata <- select(countdata, c("18_T33_Of_VLL", "26_T12_Of_WCL", "30_T23_Of_RPG", "32_T22_Of_EVR", "36_T43_Of_JJN", "40_T13_Of_GWS", "48_T31_Of_JNO", "50_T21_Of_YZB", "51_T42_Of_UOF","59_T11_Of_TQP", "60_T32_Of_WY"))
col_ofav <- colnames(ofav_countdata)
col_ofav <- as.data.frame(col_ofav)
ofav_metadata <- cbind(ofav_metadata, col_ofav)
ofav_metadata <- select(ofav_metadata, -SampleID)
names(ofav_metadata)[names(ofav_metadata) == "col_ofav"] <- "SampleID" 
rownames(ofav_metadata) <- ofav_metadata$SampleID
ofav_ID <- ofav_metadata$SampleID
# ofav_countdata <- select(ofav_countdata, all_of(ofav_ID))
all(rownames(ofav_metadata) %in% colnames(ofav_countdata)) # must come out TRUE


## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(ofav_countdata, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- ofav_countdata[gfilt,]  
dim(gkeep) # 18874 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(ofav_countdata[which(rownames(ofav_countdata) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(gcount_filt)
dim(gcount_filt) # 18874 x 11 -- only 18874 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(ofav_metadata, "Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/metadata_ofav_raw_filtered.csv")
write.csv(gcount_filt, "Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_ofav_matrix_filtered.csv")     
          
#Checking again that all row and column names match. Must return "TRUE"
all(rownames(ofav_metadata) %in% colnames(gcount_filt))


## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
ofav_metadata$Treatment <- factor(ofav_metadata$Treatment, levels = c("Treatment1", "Treatment2", "Treatment3", "Treatment4"))
head(ofav_metadata)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_ofav <- DESeqDataSetFromMatrix(countData = gcount_filt,
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
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(Treatment1="cadetblue3", Treatment2="palevioletred", Treatment3="darkgreen", Treatment4="black")) +
  coord_fixed() + 
  ggtitle("O. fav (star)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
# PCA plot is very strange. Groups clumping but all treatments are in the different clumps 

# save plots 
ofav_heatmap_PCA <- grid.arrange(ofav_PCAplot, ofav_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Plots/ofav_star_heatmap_PCA.pdf", ofav_heatmap_PCA, width = 8, height = 8, units = c("in"))



## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_ofav <- DESeq(gdds_ofav) #run differential expression test by group using the Wald model
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing
res_DEG_ofav <- results(DEG_ofav)
res_DEG_ofav_Ordered <- res_DEG_ofav[order(res_DEG_ofav$pvalue),]
DEG_ofav$Treatment
resultsNames(DEG_ofav)
# [1] "Intercept"                          "Treatment_Treatment2_vs_Treatment1" "Treatment_Treatment3_vs_Treatment1" "Treatment_Treatment4_vs_Treatment1"
# i think it is assuming that Treatment 1 is my control. but alas it is not. 

# I'm only going to do contrasts now?? Not sure what the difference is, will look into it 
# Treatment1 vs Treatment 2
DEG_ofav_results_1.vs.2 <- results(DEG_ofav, contrast = c("Treatment", "Treatment1", "Treatment")) # results only for 1 vs 2 treatments 
# tried to run as contrast, did not work. It gave me the error: error in cleanContrast(object, contrast, expanded = isExpanded, listValues = listValues,as Treatment1 is the reference level, was expecting Treatment_Treatment_vs_Treatment1 to be present in 'resultsNames(object)'
# So they are using treatment 1 as ref level. I'll proceed but these results are probably rough
# W/ reference levek
DEG_ofav_results_1.vs.2 <- results(DEG_ofav, name = "Treatment_Treatment2_vs_Treatment1") # results only for 1 vs 1 treatments 
results_ordered_DEG_ofav_results_1.vs.2 <- DEG_ofav_results_1.vs.2[order(DEG_ofav_results_1.vs.2$pvalue),] # order from smallest pvalue 
summary(DEG_ofav_results_1.vs.2) # view summary of results with adj p < 0.1
# No LCF up or down, only 38 outliers
ofav_sig.num.1vs.2<- sum(DEG_ofav_results_1.vs.2$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05? -- 0
# There are 0 diff expressed genes between treatment 1 and 2 that are less than 0.05. Going to take a break from this analysis until I get control samples









