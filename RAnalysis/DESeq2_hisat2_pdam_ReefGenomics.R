# Title: DESeq2 with pdam samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. dam only samples analyzed here aligned against P. dam. HISAT2 was read aligner. ReefGenomics fasta and gff files were used.

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
countdata <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/hisat2/gene_count_pdam_rgGFF_hisat2_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 1675 x 64
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata)[col])
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
metadata$Treatment <- gsub("<NA>", "unknown", metadata$Treatment)
rownames(metadata) <- metadata$SampleID # make sampleID the row names 
metadata <- metadata[-65,] # remove random blank space at the end
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
dim(gkeep) # 872 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(count_pdam[which(rownames(count_pdam) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(gcount_filt)
dim(gcount_filt) # 872 x 17 -- only 872 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
write.csv(metadata_pdam, "~/Desktop/metadata_pdam_raw_filtered.csv")
write.csv(gcount_filt, "~/Desktop/genecount_pdam_rgGFF_hisat2_filtered.csv")

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
  ggtitle("P. dam (hisat2)") +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
# control and mid clumped on PC2, high and control spread over PC1

# save plots 
pdam_heatmap_PCA <- grid.arrange(pdam_PCAplot, pdam_heatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/pdam_heatmap_hisat2_PCA.pdf", pdam_heatmap_PCA, width = 8, height = 8, units = c("in"))



















## Run DGE analysis

# After looking at plots to make sure all genes are looking okay, run DE analysis 
# Use Wald model 
DEG_pdam <- DESeq(gdds_pdam) #run differential expression test by group using the Wald model
res_DEG_pdam <- results(DEG_pdam) # obtain results
res_DEG_pdam_Ordered <- res_DEG_pdam[order(res_DEG_pdam$pvalue),] # order results by lowest pvalue
head(res_DEG_pdam_Ordered)

# Obtain results names in order to compare by treatment 
DEG_pdam$Treatment
resultsNames(DEG_pdam)
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"
# Why no mid v high? - bc its a contrast


## Explore significant p-values for treatments 
# Control vs mid
DEG_pdam_results_control_vs_mid <- results(DEG_pdam, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_mid <- DEG_pdam_results_control_vs_mid[order(DEG_pdam_results_control_vs_mid$pvalue),] # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_mid) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_mid <- sum(DEG_pdam_results_control_vs_mid$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 147 significantly differentially expressed genes between control and mid that are less than 0.05
pdam_DEGs.control_vs_mid <- subset(DEG_pdam_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_mid <- as.data.frame(pdam_DEGs.control_vs_mid) # make df
pdam_DEGs.control_vs_mid$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
pdam_DEGs.control_vs_mid <- cbind(gene_id = rownames(pdam_DEGs.control_vs_mid), pdam_DEGs.control_vs_mid) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_mid) <- NULL # remove row names 
pdam_DEGs.control_vs_mid
dim(pdam_DEGs.control_vs_mid) # 147 by 8
write.csv(pdam_DEGs.control_vs_mid, "~/Desktop/pdam_DEGs.control_vs_mid.csv")
pdam_sig.num.control_vs_mid

# Control vs high
DEG_pdam_results_control_vs_high <- results(DEG_pdam, name = "Treatment_high_vs_control") # results only for control vs high treatments 
results_ordered_DEG_pdam_results_control_vs_high <- DEG_pdam_results_control_vs_high[order(DEG_pdam_results_control_vs_high$pvalue),] # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_high <- sum(DEG_pdam_results_control_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 361 significantly differentially expressed genes between control and mid that are less than 0.05
pdam_DEGs.control_vs_high <- subset(DEG_pdam_results_control_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_high <- as.data.frame(pdam_DEGs.control_vs_high) # make df
pdam_DEGs.control_vs_high$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
pdam_DEGs.control_vs_high <- cbind(gene_id = rownames(pdam_DEGs.control_vs_high), pdam_DEGs.control_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_high) <- NULL # remove row names 
pdam_DEGs.control_vs_high
dim(pdam_DEGs.control_vs_high) # 361 by 8
write.csv(pdam_DEGs.control_vs_high, "~/Desktop/pdam_DEGs.control_vs_high.csv")
pdam_sig.num.control_vs_high

# Mid vs high
DEG_pdam_results_mid_vs_high <- results(DEG_pdam, contrast = c("Treatment", "mid", "high")) # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_mid_vs_high <- order(DEG_pdam_results_mid_vs_high$pvalue) #Order p-values by smallest value first
summary(DEG_pdam_results_mid_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.mid_vs_high <- sum(DEG_pdam_results_mid_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
pdam_sig.num.mid_vs_high # 7 significantly differentially expressed genes between high and mid 
pdam_DEGs.mid_vs_high <- subset(DEG_pdam_results_mid_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.mid_vs_high <- as.data.frame(pdam_DEGs.mid_vs_high) # make df
pdam_DEGs.mid_vs_high$contrast <- as_factor(c("mid_vs_high")) # set contrast as a factor
pdam_DEGs.mid_vs_high <- cbind(gene_id = rownames(pdam_DEGs.mid_vs_high), pdam_DEGs.mid_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.mid_vs_high) <- NULL # remove row names 
pdam_DEGs.mid_vs_high
dim(pdam_DEGs.mid_vs_high) # 7 by 8
write.csv(pdam_DEGs.mid_vs_high, "~/Desktop/pdam_DEGs.mid_vs_high.csv")
pdam_sig.num.mid_vs_high

# Combine treatment comparisons together 
DEGs_pdam_contrast_all_treatments <- bind_rows(pdam_DEGs.control_vs_mid, pdam_DEGs.control_vs_high, pdam_DEGs.mid_vs_high) # bind mid_vs_control results, high_vs_control results, and mid_vs_high results together by row  
dim(DEGs_pdam_contrast_all_treatments) # 515 by 8
DEG_pdam.sg.num <- sum(DEGs_pdam_contrast_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_pdam.sg.num # 515 adjusted p-values 
summary(DEGs_pdam_contrast_all_treatments)
write.csv(DEGs_pdam_contrast_all_treatments, file="~/Desktop/DEGs_pdamSamples_contrast_all_treatments.csv")




## Visualize diff-expressed genes

# get rid of redundant genes 
dim(DEGs_pdam_contrast_all_treatments)
DEGs_pdam <- DEGs_pdam_contrast_all_treatments$gene_id # list all gene names 
DEGs_pdam <- unique(DEGs_pdam) # select only unique gene names 
DEG_pdam_list <- gdds_pdam[which(rownames(gdds_pdam) %in% DEGs_pdam)] # filter gdds_plob DESeq2 object by unique gene names
dim(DEG_pdam_list) # 380 x 12
DEG_pdamSamples_list_save <- print(counts(DEG_pdam_list))
write.csv(DEG_pdamSamples_list_save, file="~/Desktop/DEG_pdamSamples_list_alltreatments_save.csv")


# As determined above, size factors all less than 4, so proceed with VST
# apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_pdam_list, blind = FALSE, nsub = nrow(counts(DEG_pdam_list)))
dim(DEGvst) # 380 by 12 
print(assay(DEGvst)) # look at vst-transformed gene count data 

# Plot heat map with diff expressed genes
pdam_topVarGenes <- head(order(rowVars(assay(DEGvst)),decreasing=TRUE), DEG_pdam.sg.num) #sort counts by decreasing sig ?
mat_pdam <- assay(DEGvst)[pdam_topVarGenes, ] #make an expression object
mat_pdam <- mat_pdam - rowMeans(mat_pdam) #diff in expression compared to average across all samples
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
pdam_DEGPCAdata <- plotPCA(DEGvst, intgroup = c("Treatment", "Days"), returnData=TRUE)
percentVar_pca_pdam <- round(100*attr(pdam_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
pdam_DEGPCAdata$Days <- gsub("7", "seven", pdam_DEGPCAdata$Days)
pdam_DEGPCAdata$Days <- gsub("4", "four", pdam_DEGPCAdata$Days)

pdam_DEGPCAplot <- ggplot(pdam_DEGPCAdata, aes(PC1, PC2, color=Treatment, fill = Treatment, shape = Days)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_pdam[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_pdam[2],"% variance")) +
  scale_color_manual(values = c(control="cadetblue3", mid="palevioletred", high="darkgreen")) +
  scale_shape_manual(values=c(16, 17)) + 
  coord_fixed() +
  ggtitle("P. damicornis samples on P. damicornis genome by Treatment and Day (hisat2)")
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
pdam_DEGPCAplot
# PCA plot is of differentially expressed genes only 

# Save results
ggsave("~/Desktop/pdamSamples_DEGs_PCA_Treatment_Days.pdf", pdam_DEGPCAplot, width = 8, height = 8, units = c("in"))







## GO analysis 
# Polina code below 


# Control vs mid 
DEG_pdam_results_control_vs_mid <- results(DEG_pdam, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
pdam_DEGs.control_vs_mid <- subset(DEG_pdam_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
df_con_mid <- as.data.frame(pdam_DEGs.control_vs_mid) %>% arrange(padj)
df_con_mid$name_t <- rownames(df_con_mid)
df_con_mid <- df_con_mid %>% separate(name_t, "go_name", sep='-', extra="drop")
#df_con_mid <- df_con_mid %>% mutate(go_annots = list_of_ontologies[go_name])
mid_sig_genes <- filter(df_con_mid, padj < 0.05)$go_name

# Control vs high
DEG_pdam_results_control_vs_high <- results(DEG_pdam, name = "Treatment_high_vs_control") # results only for control vs high treatments 
pdam_DEGs.control_vs_high <- subset(DEG_pdam_results_control_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_high$names  <- rownames(pdam_DEGs.control_vs_high)
df_con_high <- as.data.frame(pdam_DEGs.control_vs_high) %>% arrange(padj)
df_con_high$name_t <- rownames(df_con_high)
df_con_high <- df_con_high %>% separate(name_t, "go_name", sep='-', extra="drop")
#df_con_high <- df_con_high %>% mutate(go_annots = list_of_ontologies[go_name])
high_sig_genes <- filter(df_con_high, padj < 0.05)$go_name

# mid vs high
DEG_pdam_results_mid_vs_high <- results(DEG_pdam, contrast = c("Treatment", "mid", "high")) # results only for control vs mid treatments 
pdam_DEGs.mid_vs_high <- subset(DEG_pdam_results_mid_vs_high, padj<0.05) # subset only <0.05 padj values
df_mid_high <- as.data.frame(pdam_DEGs.mid_vs_high) %>% arrange(padj)
df_mid_high$name_t <- rownames(df_mid_high)
df_mid_high <- df_mid_high %>% separate(name_t, "go_name", sep='-', extra="drop")
mid.high_sig_genes <- filter(df_mid_high, padj < 0.05)$go_name

#length(high_sig_genes)
#length(mid_sig_genes)
#length(unique(c(high_sig_genes,mid_sig_genes)))

hi_OE_gv <- as.integer(genes%in%c(df_con_high$names))
mid_OE_gv <- as.integer(genes%in%c(mid_sig_genes))
names(hi_OE_gv) <- genes
names(mid_OE_gv) <- genes

## Make vectors for goseq
genes <- unlist(lapply(rownames(DEG_pdam), function(x) strsplit(x, '-')[[1]][1]))
genes <- rownames(DEG_pdam)
#length(genes)
#length(unique(genes))

# Get all of the GO Annotation information from the gff3
# How to make these files:
# ** for (tab) instead type ctrl + v and then tab. This enters a 'real' tab character. I don't know if this is
#   different for macs.
# tmp.gff3 -> grep "maker(tab)gene" pdam_annotation.gff3 > tmp.gff3
# go_only.gff3 -> grep "GO:" tmp.gff3 > go_only.gff3
GO_annot <- read.table("~/Desktop/go_only.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
list_of_annots <- lapply(GO_annot$Attr, function(x) strsplit(x,split=";")[[1]])
list_of_ontologies <- lapply(list_of_annots, grep, value=TRUE, pattern="Ontology_term")
names(list_of_ontologies) <- unlist(lapply(list_of_annots, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))
list_of_ontologies <- lapply(list_of_ontologies, function(x) strsplit(strsplit(x, "=")[[1]][-1],",")[[1]])

# Get all of the gene length information from another gff3 (a less filtered one)
all_annot <- read.table("~/Desktop/tmp.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
gene_lengths <- mutate(all_annot, len = End - Start)$len
list_of_all_annots <- lapply(all_annot$Attr, function(x) strsplit(x,split=";")[[1]])
names(gene_lengths) <- unlist(lapply(list_of_all_annots, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))

# ThisIsFine.jpeg
# - over the list of all ontology-gene pairs
# - go into the vector of go annots
# - create single-row dataframes containing the name of the protein and the annot for EACH annot
# - step back out to the second do.call, rbind them into a DF of all annots for 1 gene
# - repeat this for all genes
# - step back out to the first do.call, rbind all of these dfs into a df of all annots for all genes
# Yeah.
# This takes a little bit to run. Just like a couple of seconds (compared to instantaneous.)
# I'm sure it's possible to make it a little faster, but we're doing a LOT of rbind calls.
annotation_df <- do.call(
  rbind,
  lapply(
    seq_along(list_of_ontologies), 
    function(x) do.call(
      rbind,
      lapply(list_of_ontologies[[x]],
             function(y) data.frame(protein=names(list_of_ontologies)[[x]], annotation=y)
      )
    )
  )
)

hi_pwf <- nullp(hi_OE_gv, names(gene_lengths), bias.data=gene_lengths)

hi_pwf <- nullp(hi_OE_gv)


















## My code, trying GOseq

# Control vs mid
DEG_pdam_results_control_vs_mid <- results(DEG_pdam, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_mid <- DEG_pdam_results_control_vs_mid[order(DEG_pdam_results_control_vs_mid$pvalue),] # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_mid) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_mid <- sum(DEG_pdam_results_control_vs_mid$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 147 significantly differentially expressed genes between control and mid that are less than 0.05
pdam_DEGs.control_vs_mid <- subset(DEG_pdam_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_mid <- as.data.frame(pdam_DEGs.control_vs_mid) # make df
pdam_DEGs.control_vs_mid$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
pdam_DEGs.control_vs_mid <- cbind(gene_id = rownames(pdam_DEGs.control_vs_mid), pdam_DEGs.control_vs_mid) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_mid) <- NULL # remove row names 
pdam_DEGs.control_vs_mid
dim(pdam_DEGs.control_vs_mid) # 147 by 8

# reading data 
# reading DE genes for control vs mid
de.genes <- pdam_DEGs.control_vs_mid$gene_id
length(de.genes)
is.vector(de.genes)

# reading all genes analyzed 
assayed.genes<- rownames(countdata_pdam)
length(assayed.genes)
is.vector(assayed.genes)

# mark DE genes in all genes 
gene.vector <- as.integer(assayed.genes%in%de.genes)
names(gene.vector) <- assayed.genes
head(gene.vector)
length(gene.vector)
is.vector(gene.vector)
gene.vec_df <- as.data.frame(gene.vector) # need to remove everything before | 

# gff annotation file (from NCBI)
og_all_annot <- read.table("~/Desktop/GCF_003704095.1_ASM370409v1_genomic.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
# subset by gene type only 
annotate.genes <- subset(og_all_annot, Type == "gene")
dim(annotate.genes)
annotate.genes <- separate(annotate.genes, Attr, into = c("ID", "Dbxref", "Name", "xloc", "gbkey", "gene", "gene_biotype"), sep = ";")
annotate.genes$ID <- gsub("ID=", "", annotate.genes$ID)
annotate.genes$Dbxref <- gsub("Dbxref=", "", annotate.genes$Dbxref)
annotate.genes$Name <- gsub("Name=", "", annotate.genes$Name)
annotate.genes$gbkey <- gsub("gene=", "", annotate.genes$gbkey)
annotate.genes$gbkey <- gsub("gbkey=", "", annotate.genes$gbkey)
annotate.genes$xloc <- gsub("gbkey=", "", annotate.genes$xloc)
annotate.genes$xloc <- gsub("end_range=", "", annotate.genes$xloc)
annotate.genes$gene_biotype <- gsub("gene_biotype=", "", annotate.genes$gene_biotype)

# calculating gene lengths
gene_lengths <- mutate(annotate.genes, len = End - Start)$len
gene_ID <- annotate.genes$ID
length(gene_ID)

# binding together gene lengths and gene id
gene_length_ID <- cbind(gene_ID, gene_lengths) 
gene_length_ID <- as.data.frame(gene_length_ID)
dim(gene_length_ID) # 21837 by 2

length(gene.vector) # 23079 by 2

# Problem is that gene length and vector are not same vector length. Why? 
# May be an issue with subsetting by 'gene' above.











# hi_pwf <- nullp(gene.vector, names(gene_length_ID), bias.data=gene_length_ID)
# gene_lengths <- mutate(og_all_annot, len = End - Start)$len
# names(gene_lengths) <- unlist(lapply(og_all_annot, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))
# 
# annotate.genes <- subset(og_all_annot, Type == "gene")
# gene_lengths <- mutate(annotate.genes, len = End - Start)$len

# this was for NCBI with known pdam samples (before other metadata was gathered)
## Obtaining and tidying data 
# countdata_pdam <- read.csv("Desktop/PutnamLab/Repositories/Tufts_URI_RNAseq/Tufts_URI_CSM_RNASeq/Seneca/STAR_pipeline/Output/DESeq2/gene_count_pdamSamples_matrix.csv", header = TRUE, row.names = "gene_id")
# dim(countdata_pdam) # 23079 x 12
# head(countdata_pdam)
# for ( col in 1:ncol(countdata_pdam)){
#   colnames(countdata_pdam)[col] <-  sub(".fastq.trim.fq.bam.merge.gtf", "", colnames(countdata_pdam)[col])
# }


  