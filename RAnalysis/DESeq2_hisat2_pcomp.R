# Title: DESeq2 with pcomp samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. Pcomp only samples analyzed here aligned against Pcomp. HISAT2 was read aligner.

# Read in count data for all samples
countdata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/hisat2/gene_count_pcomp_hisat2_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 74726 x 64
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
# Select mcap species only 
metadata_pcomp <- subset(metadata, Species=="Porites compressa")

# Subset count data for only pcomp samples based on SampleID and make sure rows of metadata = cols of count data
pcomp_ID <- metadata_pcomp$SampleID
count_pcomp <- select(countdata, all_of(pcomp_ID))
all(rownames(metadata_pcomp) %in% colnames(count_pcomp)) 



## Pre-filter gene counts

# Set filter values for PoverA, P=85% percent of the samples have counts over A=5. 
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(count_pcomp, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pcomp[gfilt,]  
dim(gkeep) # only 17950 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt_pcomp <- as.data.frame(count_pcomp[which(rownames(count_pcomp) %in% gn.keep),])
head(gcount_filt_pcomp)
dim(gcount_filt_pcomp) # 17950 x 17 -- only 17950 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_pcomp, "~/Desktop/metadata_pcomp_raw_filtered.csv")
write.csv(gcount_filt_pcomp, "~/Desktop/gene_count_pcomp_hisat2_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
all(rownames(metadata_pcomp) %in% colnames(gcount_filt_pcomp))



## Construct DESeq2 dataset 

# Set Treatment as a factor and give levels 
metadata_pcomp$Treatment <- factor(metadata_pcomp$Treatment, levels = c("control", "mid", "high", "unknown"))
head(metadata_pcomp)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
gdds_pcomp <- DESeqDataSetFromMatrix(countData = gcount_filt_pcomp,
                                    colData = metadata_pcomp,
                                    design = ~1) # because we don't know the treatment metadata for the mcap samples yet, there is only one treatment 'unknown'. So we can just set design as 1
gdds_pcomp



## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.gdds_pcomp <- estimateSizeFactors(gdds_pcomp) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.gdds_pcomp
print(sizeFactors(SF.gdds_pcomp)) #view size factors

# size factors all less than 4, can use VST
gvst_pcomp <- vst(gdds_pcomp, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(gvst_pcomp))
dim(gvst_pcomp)

# Using vst object, Plot heat-map of sample-to-sample distances
gsampleDists <- dist(t(assay(gvst_pcomp))) # calculate distance matrix, t returns transpose of assay(gvst)
gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
rownames(gsampleDistsMatrix) <- colnames(gvst_pcomp) # assign row names
colnames(gsampleDistsMatrix) <- NULL # assign col names
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
pcomp_heatmap<- pheatmap(gsampleDistsMatrix, # plot matrix
                        clustering_distance_rows = gsampleDists, # cluster rows
                        clustering_distance_cols = gsampleDists, # cluster cols
                        col = colors) # set colors

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(gvst_pcomp, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pcomp_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape =Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() +
  ggtitle("P. compressa (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background

# save plots 
pcomp_heatmap_PCA <- grid.arrange(pcomp_PCAplot, pcomp_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/pcomp_heatmap_hisat2_PCA.pdf", pcomp_heatmap_PCA, width = 8, height = 8, units = c("in"))
