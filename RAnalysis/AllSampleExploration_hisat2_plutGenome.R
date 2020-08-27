# Title: PCA all sample exploration
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Examining PCA plot with all samples, regardless of species. At the moment, Francois isn't totally sure which samples are which, but hope to look at PCA plot to see how certain samples are clustering
# Using data generated from HISAT2 with Plutea genome

# Read in count data for all samples
countdata_plob_allsamples <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/gene_count_plobGenome_allsamples_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata_plob_allsamples) # 31126 x 64
head(countdata_plob_allsamples)
for ( col in 1:ncol(countdata_plob_allsamples)){
  colnames(countdata_plob_allsamples)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata_plob_allsamples)[col])
}

# Load metadata 
metadata <- read.csv("~/Desktop/sediment_metadata_allSamples_raw.csv", header = TRUE)
head(metadata)
names(metadata)[names(metadata) == "File.Name.fastq"] <- "SampleID"
names(metadata)[names(metadata) == "Treatment.in.mg.L.of.sediment"] <- "Treatment"
names(metadata)[names(metadata) == "Time.point.in.days"] <- "Days"
metadata$Treatment <- gsub("400", "high", metadata$Treatment)
metadata$Treatment <- gsub("40", "mid", metadata$Treatment)
metadata$Treatment <- gsub("0", "control", metadata$Treatment)
metadata$SampleID <- gsub(".fastq.gz", "", metadata$SampleID)
metadata$SampleID <- paste0("X", metadata$SampleID)
metadata$Days <- gsub("7", "seven" ,metadata$Days)
metadata$Days <- gsub("4", "four" ,metadata$Days)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[-65,]


# Check to make sure rownames in metadata file match column names in countdata file - must return TRUE
col_order <- rownames(metadata)
countdata_plob_allsamples <- countdata_plob_allsamples[, col_order]
all(rownames(metadata) %in% colnames(countdata_plob_allsamples)) 


# Not sure if I should filter because it leave me with so few genes
# Pre-filter gene counts
# Set filter values for PoverA, P=85% percent of the samples have counts over A=5.
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
gfilt <- genefilter(countdata_plob_allsamples, filt)
gfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- countdata_plob_allsamples[gfilt,]
dim(gkeep) # new df of genes that passed filtering

# List names of genes that passed filtering
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names
gcount_filt <- as.data.frame(countdata_plob_allsamples[which(rownames(countdata_plob_allsamples) %in% gn.keep),])
head(gcount_filt)
dim(gcount_filt) # 21 by 64...only that many left after filtering?

# Write treatment, gene and transcript count files with corrected column and row headers
# write.csv(metadata_pdam, "~/Desktop/metadata_pdamSamples_filtered.csv")
# write.csv(gcount_filt, "~/Desktop/genecount_pdamSamples_filtered.csv")

#Checking again that all row and column names match. Must return "TRUE"
gcount_filt <- gcount_filt[, col_order]
all(rownames(metadata) %in% colnames(gcount_filt))




## Construct DESeq2 dataset 

# Set group as a factor and give levels 
metadata$Treatment <- factor(metadata$Treatment, levels = c("control", "mid", "high", "unknown"))
head(metadata)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
dds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                    colData = metadata,
                                    design = ~Treatment)
dds




## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.dds <- estimateSizeFactors(dds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.dds
print(sizeFactors(SF.dds)) #view size factors

# some definitely not less than 4? use rlog
grlog <- vst(dds, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(grlog))
dim(grlog)

# Using vst object, Plot heat-map of sample-to-sample distances
# gsampleDists <- dist(t(assay(grlog))) # calculate distance matrix, t returns transpose of assay(grlog)
# gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
# rownames(gsampleDistsMatrix) <- colnames(grlog) # assign row names 
# colnames(gsampleDistsMatrix) <- NULL # assign col names 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistsMatrix, # plot matrix
#          clustering_distance_rows = gsampleDists, # cluster rows
#          clustering_distance_cols = gsampleDists, # cluster cols
#          col = colors) # set colors 

## Using rlog object, make PCA plot of samples 
# one PCA with treatment and days
gPCAdata <- plotPCA(grlog, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_all_treatment_day <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape = Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() + 
  ggtitle("All species on P. lutea genome by treatment and day (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
# appears to be 3-4 groupings. 

# one pca with treatment and species 
gPCAdata <- plotPCA(grlog, intgroup = c("Treatment", "Species"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_all_treatment_species <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape = Species)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() + 
  ggtitle("All species on P. lutea genome by treatment and species (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background




## Will now plot PCA by species to get clearer idea of spread

## Pdam

# Subset data in both metadata and counts dfs so it is just pdam samples
metadata_pdam <- subset(metadata, Species=="Pocillopora damicornis")
pdam_ID <- metadata_pdam$SampleID
count_pdam <- select(countdata_plob_allsamples, all_of(pdam_ID))
all(rownames(metadata_pdam) %in% colnames(count_pdam)) 


# Not sure if I should filter because it leave me with so few genes, as I mapped pdam to plut genome
# Pre-filter gene counts
# Set filter values for PoverA, P=85% percent of the samples have counts over A=5.
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
pdamfilt <- genefilter(count_pdam, filt)
pdamfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pdam[pdamfilt,]
dim(gkeep) # new df of genes that passed filtering

# List names of genes that passed filtering
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names
pdamcount_filt <- as.data.frame(count_pdam[which(rownames(count_pdam) %in% gn.keep),])
head(pdamcount_filt)
dim(pdamcount_filt) # 35 by 17
all(rownames(metadata_pdam) %in% colnames(pdamcount_filt))


## Construct DESeq2 dataset 

# Set group as a factor and give levels 
metadata_pdam$Treatment <- factor(metadata_pdam$Treatment, levels = c("control", "mid", "high", "unknown"))
head(metadata_pdam)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
pdam_dds <- DESeqDataSetFromMatrix(countData = pdamcount_filt,
                              colData = metadata_pdam,
                              design = ~Treatment)
pdam_dds
  


## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.pdam_dds <- estimateSizeFactors(pdam_dds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.pdam_dds
print(sizeFactors(SF.pdam_dds)) #view size factors

# some definitely not less than 4? using rlog
pdam_rlog <- rlog(pdam_dds, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(pdam_rlog))
dim(pdam_rlog)

# Using vst object, Plot heat-map of sample-to-sample distances
# gsampleDists <- dist(t(assay(pdam_rlog))) # calculate distance matrix, t returns transpose of assay(gvst)
# gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
# rownames(gsampleDistsMatrix) <- colnames(pdam_rlog) # assign row names 
# colnames(gsampleDistsMatrix) <- NULL # assign col names 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistsMatrix, # plot matrix
#          clustering_distance_rows = gsampleDists, # cluster rows
#          clustering_distance_cols = gsampleDists, # cluster cols
#          col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(pdam_rlog, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_pdam <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape = Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() + 
  ggtitle("P. damicornis on P. lutea genome by treatment and day (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background



## Plob

# Subset data in both metadata and counts dfs so it is just pdam samples
metadata_plob <- subset(metadata, Species=="Porites lobata")
plob_ID <- metadata_plob$SampleID
count_plob <- select(countdata_plob_allsamples, all_of(plob_ID))
all(rownames(metadata_plob) %in% colnames(count_plob)) 


# Not sure if I should filter because it leave me with so few genes, as I mapped pdam to plut genome
# Pre-filter gene counts
# Set filter values for PoverA, P=85% percent of the samples have counts over A=5.
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
plobfilt <- genefilter(count_plob, filt)
plobfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_plob[plobfilt,]
dim(gkeep) # new df of genes that passed filtering

# List names of genes that passed filtering
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names
plobcount_filt <- as.data.frame(count_plob[which(rownames(count_plob) %in% gn.keep),])
head(plobcount_filt)
dim(plobcount_filt) # 12504 by 17
all(rownames(metadata_plob) %in% colnames(plobcount_filt))


## Construct DESeq2 dataset 

# Set group as a factor and give levels 
metadata_plob$Treatment <- factor(metadata_plob$Treatment, levels = c("control", "mid", "high", "unknown"))
head(metadata_plob)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
plob_dds <- DESeqDataSetFromMatrix(countData = plobcount_filt,
                                   colData = metadata_plob,
                                   design = ~Treatment)
plob_dds



## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.lob_dds <- estimateSizeFactors(plob_dds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.lob_dds
print(sizeFactors(SF.lob_dds)) #view size factors

# All less than 4, can use vst
plob_vst <- vst(plob_dds, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(plob_vst))
dim(plob_vst)

# Using vst object, Plot heat-map of sample-to-sample distances
# gsampleDists <- dist(t(assay(plob_vst))) # calculate distance matrix, t returns transpose of assay(gvst)
# gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
# rownames(gsampleDistsMatrix) <- colnames(plob_vst) # assign row names 
# colnames(gsampleDistsMatrix) <- NULL # assign col names 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistsMatrix, # plot matrix
#          clustering_distance_rows = gsampleDists, # cluster rows
#          clustering_distance_cols = gsampleDists, # cluster cols
#          col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(plob_vst, intgroup = c("Treatment", "Days"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_plob <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape =Days)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() +
  ggtitle("P. lobata on P. lutea genome by treatment and day (hisat2)") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background



## Pcomp

# Subset data in both metadata and counts dfs so it is just pdam samples
metadata_pcomp <- subset(metadata, Species=="Porites compressa")
pcomp_ID <- metadata_pcomp$SampleID
count_pcomp <- select(countdata_plob_allsamples, all_of(pcomp_ID))
all(rownames(metadata_pcomp) %in% colnames(count_pcomp)) 


# Not sure if I should filter because it leave me with so few genes, as I mapped pdam to plut genome
# Pre-filter gene counts
# Set filter values for PoverA, P=85% percent of the samples have counts over A=5.
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
pcompfilt <- genefilter(count_pcomp, filt)
pcompfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pcomp[pcompfilt,]
dim(gkeep) # new df of genes that passed filtering

# List names of genes that passed filtering
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names
pcompcount_filt <- as.data.frame(count_pcomp[which(rownames(count_pcomp) %in% gn.keep),])
head(pcompcount_filt)
dim(pcompcount_filt) # 14597 by 17
all(rownames(metadata_pcomp) %in% colnames(pcompcount_filt))


## Construct DESeq2 dataset 

# Set group as a factor and give levels 
# metadata_pcomp$Treatment <- factor(metadata_pcomp$Treatment, levels = c("control", "mid", "high", "unknown"))
# head(metadata_pcomp)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
pcomp_dds <- DESeqDataSetFromMatrix(countData = pcompcount_filt,
                                   colData = metadata_pcomp,
                                   design = ~1) # can set as 1 because treatments are all unknown
pcomp_dds



## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.comp_dds <- estimateSizeFactors(pcomp_dds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.comp_dds
print(sizeFactors(SF.comp_dds)) #view size factors

# All less than 4, can use vst
pcomp_vst <- vst(pcomp_dds, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(pcomp_vst))
dim(pcomp_vst)

# Using vst object, Plot heat-map of sample-to-sample distances
# gsampleDists <- dist(t(assay(pcomp_vst))) # calculate distance matrix, t returns transpose of assay(gvst)
# gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
# rownames(gsampleDistsMatrix) <- colnames(pcomp_vst) # assign row names 
# colnames(gsampleDistsMatrix) <- NULL # assign col names 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistsMatrix, # plot matrix
#          clustering_distance_rows = gsampleDists, # cluster rows
#          clustering_distance_cols = gsampleDists, # cluster cols
#          col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(pcomp_vst, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_pcomp <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() + 
  ggtitle("P. compressa on P. lutea genome by treatment (hisat2)") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background



## Mcap

# Subset data in both metadata and counts dfs so it is just pdam samples
metadata_mcap <- subset(metadata, Species=="Montipora capitata")
mcap_ID <- metadata_mcap$SampleID
count_mcap <- select(countdata_plob_allsamples, all_of(mcap_ID))
all(rownames(metadata_mcap) %in% colnames(count_mcap)) 


# Not sure if I should filter because it leave me with so few genes, as I mapped pdam to plut genome
# Pre-filter gene counts
# Set filter values for PoverA, P=85% percent of the samples have counts over A=5.
filt <- filterfun(pOverA(0.85,5)) # creating filter function

# Create filter for counts data
mcapfilt <- genefilter(count_mcap, filt)
mcapfilt # gives T or F for which genes have < 5 counts

# Id genes to keep by count filter
gkeep <- count_pcomp[mcapfilt,]
dim(gkeep) # new df of genes that passed filtering

# List names of genes that passed filtering
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names
mcapcount_filt <- as.data.frame(count_mcap[which(rownames(count_mcap) %in% gn.keep),])
head(mcapcount_filt)
dim(mcapcount_filt) # 117 by 14
all(rownames(metadata_mcap) %in% colnames(mcapcount_filt))


## Construct DESeq2 dataset 

# Set group as a factor and give levels 
# metadata_pcomp$Treatment <- factor(metadata_pcomp$Treatment, levels = c("control", "mid", "high", "unknown"))
# head(metadata_pcomp)

# Create a DESeqDataSet design from gene count matrix and labels. 
# Here we set the design to test for any differences in gene expression across treatments
mcap_dds <- DESeqDataSetFromMatrix(countData = mcapcount_filt,
                                    colData = metadata_mcap,
                                    design = ~1) # can set as 1 because treatments are all unknown
mcap_dds



## Visualize gene count data
# We're looking to see if the samples of the same treatments cluster -- PURELY FOR VISUALIZATION

# Log-transform the count data
# First we are going to log-transform the data using a variance stabilizing transforamtion (vst). This is only for visualization purposes. 
# Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). 
# Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. 
# In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4. Otherwise, there could be artefacts in our results.
SF.mcap_dds <- estimateSizeFactors(mcap_dds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 to use vst
SF.mcap_dds
print(sizeFactors(SF.mcap_dds)) #view size factors

# All less than 4, can use vst
# Giving me an error trying to use vst, so just going to use rlog
mcap_rlog <- rlog(mcap_dds, blind = FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
head(assay(mcap_rlog))
dim(mcap_rlog)

# Using vst object, Plot heat-map of sample-to-sample distances
# gsampleDists <- dist(t(assay(mcap_rlog))) # calculate distance matrix, t returns transpose of assay(gvst)
# gsampleDistsMatrix <- as.matrix(gsampleDists) # create distance matrix
# rownames(gsampleDistsMatrix) <- colnames(mcap_rlog) # assign row names 
# colnames(gsampleDistsMatrix) <- NULL # assign col names 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #assign colors
# pheatmap(gsampleDistsMatrix, # plot matrix
#          clustering_distance_rows = gsampleDists, # cluster rows
#          clustering_distance_cols = gsampleDists, # cluster cols
#          col = colors) # set colors 

## Using vst object, make PCA plot of samples 
gPCAdata <- plotPCA(mcap_rlog, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
pca_mcap <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  geom_text(aes(label=name),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen", unknown="black")) +
  coord_fixed() + 
  ggtitle("M. capitata on P. lutea genome by treatment (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background



## Make pdf with all pca plots 
ggsave(file = "~/Desktop/PCAs_AllSamples_plutGenome.pdf", height = 4, width = 4.25)
pca_all_treatment_day
pca_all_treatment_species
pca_mcap
pca_pcomp
pca_pdam
pca_plob


