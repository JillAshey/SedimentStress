# Title: DESeq2 with mcap samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. M. cap only samples analyzed here aligned against M. cap. HISAT2 was read aligner.

# Read in count data for all samples
countdata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/hisat2/gene_count_mcap_hisat2_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 9352 x 64
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata)[col]) # removing .fastq.trim.fq.bam.merge.gtf from end of colnames
} 

# Load metadata 
metadata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Data/sediment_metadata_allSamples_raw.csv", header = TRUE)
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
rownames(metadata) <- metadata$SampleID # make sampleID the row names 
metadata <- metadata[-65,] # remove random blank space at the end
# Select mcap species only 
metadata_mcap <- subset(metadata, Species=="Montipora capitata")

# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
mcap_ID <- metadata_mcap$SampleID
count_mcap <- select(countdata, all_of(mcap_ID))
all(rownames(metadata_mcap) %in% colnames(count_mcap)) # must come out TRUE



## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(count_mcap, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_mcap[gfilt,]  
dim(gkeep) # only 2809 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt_mcap <- as.data.frame(count_mcap[which(rownames(count_mcap) %in% gn.keep),])
head(gcount_filt_mcap)
dim(gcount_filt_mcap) # 2809 x 14 -- only 2809 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_mcap, "~/Desktop/metadata_mcapSamples_raw_filtered.csv")
write.csv(gcount_filt_mcap, "~/Desktop/gene_count_mcap_hisat2_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(metadata_mcap) %in% colnames(gcount_filt_mcap))



## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
metadata_mcap$Treatment <- factor(metadata_mcap$Treatment, levels = c("control", "mid", "high", "unknown"))
head(metadata_mcap)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_mcap <- DESeqDataSetFromMatrix(countData = gcount_filt_mcap,
                                    colData = metadata_mcap,
                                    design = ~1) # because we don't know the treatment metadata for the mcap samples yet, there is only one treatment 'unknown'. So we can just set design as 1
gdds_mcap



## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_mcap <- estimateSizeFactors(gdds_mcap) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_mcap
print(sizeFactors(SF.gdds_mcap)) #view size factors

# size factors all less than 4, can use VST
gvst_mcap <- vst(gdds_mcap, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_mcap))
dim(gvst_mcap)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_mcap))) # calculate distance matrix, t returns transpose of assay(gvst)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_mcap) # assign row names
colnames(gsampleDistsMatrix) <- NULL # assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
mcap_heatmap<- pheatmap(gsampleDistsMatrix, # plot matrix
         clustering_distance_rows = gsampleDists, # cluster rows
         clustering_distance_cols = gsampleDists, # cluster cols
         col = colors) # set colors

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_mcap, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
mcap_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape =Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() +
  ggtitle("M. capitata (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
mcap_heatmap_PCA <- grid.arrange(mcap_PCAplot, mcap_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/mcap_heatmap_hisat2_PCA.pdf", mcap_heatmap_PCA, width = 8, height = 8, units = c("in"))









