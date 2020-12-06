# Title: DESeq2 with plob samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. lob only samples analyzed here aligned against P. lutea. STAR was read aligner 

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
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_gene_count_matrix.csv")
dim(countdata) # 31126 x 17
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub("X", "", colnames(countdata)[col])
}
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  gsub(".fastq.trim.fq.noGFF.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(countdata)[col])
}
rownames(countdata) <- countdata$gene_id

# Load functional annotation file 
annot <- read.csv("~/Desktop/GFFs/Plut.GFFannotation.fixed_transcript.gff",header = FALSE, sep="\t", skip=6)
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
annot$gene_id <-gsub(";.*", "", annot$attr)
annot$gene_id <-gsub("ID=", "", annot$gene_id)
# subset by gene ? not sure 

# Load plob metadata
metadata_plob <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/plob_metadata_raw_filtered.csv")
metadata_plob <- na.omit(metadata_plob)
metadata_plob$SampleID <- gsub("X", "", metadata_plob$SampleID)
rownames(metadata_plob) <- metadata_plob$SampleID

# Subset count data for only plob samples with metadata based on SampleID and make sure rows of metadata = cols of count data
plob_ID <- metadata_plob$SampleID
count_plob <- select(countdata, all_of(plob_ID))
all(rownames(metadata_plob) %in% colnames(count_plob)) # must come out TRUE

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_plob, filt) # create filter for counts data 
keep <- count_plob[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
plob_counts_filt <- as.matrix(count_plob[which(rownames(count_plob) %in% gn.keep),]) 
dim(plob_counts_filt) # 15369 x 12
storage.mode(plob_counts_filt) <- "integer" # stores count data as integer 
write.csv(plob_counts_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_plob) %in% colnames(plob_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_plob$Treatment <- factor(metadata_plob$Treatment, levels = c("control", "mid", "high"))
data <- DESeqDataSetFromMatrix(countData = plob_counts_filt, colData = metadata_plob, design = ~ Treatment)

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
DEG_control_vs_mid <- as.data.frame(DEG_control_vs_mid) # make full results into a df
DEG_control_vs_mid["Treatment_Compare"] <- "CvsMid" # add treatment comparison col
write.csv(DEG_control_vs_mid, file = "~/Desktop/plob_control_vs_mid_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_mid.sig.num <- sum(DEG_control_vs_mid$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_mid.sig.num
# 109 DEGs
DEG_control_vs_mid.sig <- subset(DEG_control_vs_mid, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_mid.sig["Treatment_Compare"] <- "CvsMid" # adding treatment comparison column
DEG_control_vs_mid.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_mid.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_mid.sig.list <- as.data.frame(counts(DEG_control_vs_mid.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_mid.sig.list_full <- cbind(DEG_control_vs_mid.sig, DEG_control_vs_mid.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_mid.sig.list_full, file = "~/Desktop/plob_control_vs_mid_DEG_full.csv") # write out csv
DEG_control_vs_mid.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_mid.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare C and high
DEG_control_vs_high <- results(DEG.int, contrast = c("Treatment", "control", "high"))
DEG_control_vs_high
DEG_control_vs_high <- as.data.frame(DEG_control_vs_high) # make full results into a df
DEG_control_vs_high["Treatment_Compare"] <- "CvsHigh" # add treatment comparison col
write.csv(DEG_control_vs_high, file = "~/Desktop/plob_control_vs_high_all_genes.csv") # maybe include gene counts too?
DEG_control_vs_high.sig.num <- sum(DEG_control_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_control_vs_high.sig.num
# 92 DEGs
DEG_control_vs_high.sig <- subset(DEG_control_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_control_vs_high.sig["Treatment_Compare"] <- "CvsHigh" # adding treatment comparison column
DEG_control_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_control_vs_high.sig)),] # subset list of significant genes from original count data 
DEG_control_vs_high.sig.list <- as.data.frame(counts(DEG_control_vs_high.sig.list)) # make list of sig gene counts into a df
DEG_control_vs_high.sig.list_full <- cbind(DEG_control_vs_high.sig, DEG_control_vs_high.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_control_vs_high.sig.list_full, file = "~/Desktop/plob_control_vs_high_DEG_full.csv") # write out csv
DEG_control_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_control_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 

# Compare mid and high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high <- as.data.frame(DEG_mid_vs_high) # make full results into a df
DEG_mid_vs_high["Treatment_Compare"] <- "MidvsHigh" # add treatment comparison col
write.csv(DEG_mid_vs_high, file = "~/Desktop/plob_mid_vs_high_all_genes.csv") # maybe include gene counts too?
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_mid_vs_high.sig.num
# 18 DEGs
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig["Treatment_Compare"] <- "MidvsHigh" # adding treatment comparison column
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subset list of significant genes from original count data 
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list)) # make list of sig gene counts into a df
DEG_mid_vs_high.sig.list_full <- cbind(DEG_mid_vs_high.sig, DEG_mid_vs_high.sig.list) # bind results with gene counts for DEGs
write.csv(DEG_mid_vs_high.sig.list_full, file = "~/Desktop/plob_mid_vs_high_DEG_full.csv") # write out csv
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 


##### Unique genes from intersections of DEG in CvsMid, CvsHigh, MidvsHigh
DEGs.all <- rbind(DEG_control_vs_mid.sig.list_full, 
                  DEG_control_vs_high.sig.list_full,
                  DEG_mid_vs_high.sig.list_full)
write.csv(DEGs.all, file = "~/Desktop/plob_DEGs.all_treatment.csv")
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all_plob <- DEGs.all$DEGs
DEGs.all_plob <- unique(DEGs.all_plob)
DEGs.all_plob <- as.data.frame(DEGs.all_plob) # 219 unique DEGs among treatment comparisons

unique.sig.list <- data[which(rownames(data) %in% DEGs.all$DEGs), ] # subset list of sig transcripts from original count data
SFtest <- estimateSizeFactors(unique.sig.list)
print(sizeFactors(SFtest))
unique.vst.sig <- varianceStabilizingTransformation(unique.sig.list, blind = FALSE)
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
write.csv(counts(unique.sig.list), file = "~/Desktop/plob_unique.sig.list.csv")


# PCA plot of diff-expressed genes 
plob_DEGPCAdata <- plotPCA(unique.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_plob <- round(100*attr(plob_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
plob_DEGPCAplot <- ggplot(plob_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=8) +
  #geom_text(aes(label=name), hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar_pca_plob[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_plob[2],"% variance")) +
  scale_color_manual(values = c(control="gray", mid = "darksalmon", high = "darkred")) +
  coord_fixed() +
  #ggtitle("P. lobata") +
  theme_bw() + #Set background color
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size=25),
        #title = element_text(size=30),
        legend.position = "none",
        panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
plob_DEGPCAplot
# PCA plot is of differentially expressed genes only
#PC.info <- plob_DEGPCAplot$data
ggsave("~/Desktop/plob_DEGs_PCA.png", plob_DEGPCAplot, width = 30, height = 20,, units = "cm")





## Heatmap of DEGs
# This heatmap is going to be wild. Grouping columns by treatment, putting gene id on left hand side and GO term on right hand side
# Also taking out legend 

# Need a file with counts, gene names, GO IDs, term, and ontology

# To get a file with gene names associated with GO terms, I have to look at pdam_GOterms_ByGene
# This file is pdam genes with GO terms associated with them. One GO term per line, so multiple pdam gene names sometimes if genes have multiple GO terms
# This file includes all gene names and GO IDs
plob_sig <- read.csv("~/Desktop/plob_GOterms_ByGene.csv", header = TRUE)
plob_sig <- select(plob_sig, -X)
colnames(plob_sig)[1] <-"gene"
head(plob_sig)

# Getting col order to order unique counts data
list(metadata_plob)
list <- metadata_plob[order(metadata_plob$Treatment),] # need to order them so it will group by treatment in plot
list(list$SampleID) # look at sample IDs and use that list to make col.order
col.order <- c("6_1",
               "9_1",
               "25_1",
               "26_1",
               "8_1",
               "21_1",
               "29_1",
               "34_1",
               "7_1",
               "22_1",
               "23_1",
               "27_1")
  
# Now I will order the counts data so the samples will group by treatment
unique.DEG.annot <- as.data.frame(counts(unique.sig.list)) # make df of sig genes with counts and sample IDs
list(colnames(unique.DEG.annot))
unique.DEG.annot2 <- unique.DEG.annot[, col.order]
unique.DEG.annot2$gene <- rownames(unique.DEG.annot2)

# Now we will take the unique.DEG.annot2 and merge it with acerv_sig
# The unique.DEG.annot2 file includes gene names for DEGs and counts data
test_merge <- merge(unique.DEG.annot2, plob_sig, by = "gene", all.x = TRUE)
# test_merge now holds gene names for DEGs, counts data, and GO.IDs

# Now we need info about term and ontology 
GO_all <- read.csv("~/Desktop/plob_GO_ALL.csv", header = TRUE) 
GO_all <- select(GO_all, -X)
colnames(GO_all)[1] <-"GO.ID"
GO_merge <- merge(test_merge, GO_all, by = "GO.ID", all.x = T)
# Great! GO_merge now contains GO IDs, gene names for DEGs, counts data, over and under represented pvalue, numCat, term, and ontology.
# All of the genes are in there (sometimes duplicated because multple GO/term/ontology info per gene) and some don't have any GO/term/ontology info. 
# *** I could also just plot GO_merge, but have duplicate genes for some that have multiple multple GO/term/ontology info per gene. Probably not the best idea tho

# Trying to aggregate based on GO terms. Hopefully this works because I also want the term and ontology to also aggregate but i think they may just go.
# Well maybe I could aggregate multiple times and then bind them? Lets see
agg_GO <- aggregate(GO_merge$GO.ID, list(GO_merge$gene), paste, collapse = ",") # aggregate GO terms 
colnames(agg_GO) <- c("gene", "GO.ID")
agg_term <- aggregate(GO_merge$term, list(GO_merge$gene), paste, collapse = ",") # aggregate term
colnames(agg_term) <- c("gene", "term")
agg_ont <- aggregate(GO_merge$ontology, list(GO_merge$gene), paste, collapse = ",") # aggregate ontology
colnames(agg_ont) <- c("gene", "ontology")
agg_over <- aggregate(GO_merge$over_represented_pvalue, list(GO_merge$gene), paste, collapse = ",")
colnames(agg_over) <- c("gene", "over_represented_pvalue")

# Now I'll merge them all together!
merge_all <- merge(agg_GO, agg_term, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, agg_ont, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, agg_over, by = "gene", all.x = TRUE)
merge_all <- merge(merge_all, unique.DEG.annot2, by = "gene", all.x = TRUE)
write.csv(merge_all, file = "~/Desktop/plob_GO_DEG.csv") # maybe include gene counts too?



## Hooray! Now I have a lovely file with counts, gene names, GO IDs, term, and ontology 
# Now I must put it in the heatmap...........

# First, lets make a matrix of gene counts 
rownames(merge_all) <- merge_all$gene
mat <- select(merge_all, -c("gene", "GO.ID", "term", "ontology", "over_represented_pvalue"))
mat <- as.matrix(mat) 

# Now lets make df of only treatment and sample ID
#df <- as.data.frame(colData(unique.vst.sig) [, c("Treatment")])
#colnames(df) <- "Treatment"
#df <- df[order(df$Treatment),]
#df <- as.data.frame(df)
#colnames(df) <- "Treatment"
df <- select(metadata_plob, c("Treatment"))
#df <- df[order(df$Treatment),]
# probably just easier to take treatment info straight from metadata file

# Now lets make a df of only gene names 
df_gene <- as.data.frame(merge_all$gene)
colnames(df_gene) <- "DEG"
rownames(df_gene) <- df_gene$DEG

# Some genes have multiple terms, so I am going to select the first term for every gene 
merge_all$term2 <- merge_all$term
merge_all$term <- gsub(",.*", "", merge_all$term)

# Some genes have NAs, so subbing blank for NA to see the actual terms
merge_all[is.na(merge_all$term)] <- " "

#Set colors for treatment
ann_colors <- list(Treatment = c(control="gray", mid = "darksalmon", high = "darkred"))


## Plot heatmap
plob_heatmap <- pheatmap(mat, 
                         annotation_col = df,
                         #annotation_row = df_gene,
                         annotation_colors = ann_colors,
                         annotation_legend = F,
                         cluster_rows = F,
                         show_rownames = T,
                         cluster_cols = F,
                         show_colnames = T,
                         scale = "row",
                         fontsize_row = 8,
                         labels_row = merge_all$term)
plob_heatmap
ggsave("~/Desktop/plob_heatmap.png", plob_heatmap, width = 30, height = 20,, units = "cm")






















## Going to remove the control from the list and see how DEG PCA looks without it--can look at spread of mid and high 
# Subset count data for only mcap samples based on SampleID and make sure rows of metadata = cols of count data
metadata_plob_treatment <- subset(metadata_plob, Treatment=="mid" | Treatment=="high")
plob_ID_treatment <- metadata_plob_treatment$SampleID
count_plob_treatment <- select(countdata, all_of(plob_ID_treatment))

# Filter reads by proportion of samples containing cutoff value
filt <- filterfun(pOverA(0.85, 5)) # set filter values for P over A; I used 0.85 and 5
tfil <- genefilter(count_plob_treatment, filt) # create filter for counts data 
keep <- count_plob_treatment[tfil,] # identify genes to keep based on filter
gn.keep <- rownames(keep)
plob_treatment_counts_filt <- as.matrix(count_plob_treatment[which(rownames(count_plob_treatment) %in% gn.keep),]) 
storage.mode(plob_treatment_counts_filt) <- "integer" # stores count data as integer 
#write.csv(plob_treatment_counts_filt, "~/Desktop/plob_treatment_counts_filt.csv")
# Checking to make sure rownames in metadata == colnames in counts data 
all(rownames(metadata_plob_treatment) %in% colnames(plob_treatment_counts_filt)) # must come out TRUE
# Set Treatment as a factor
metadata_plob_treatment$Treatment <- factor(metadata_plob_treatment$Treatment, levels = c("mid", "high"))
data <- DESeqDataSetFromMatrix(countData = plob_treatment_counts_filt, colData = metadata_plob_treatment, design = ~ Treatment)

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
#[1] "Intercept"             "Treatment_high_vs_mid"

# Compare mid and high
DEG_mid_vs_high <- results(DEG.int, contrast = c("Treatment", "mid", "high"))
DEG_mid_vs_high
DEG_mid_vs_high.sig.num <- sum(DEG_mid_vs_high$padj <0.05, na.rm = T) # identify # of significant pvalues with 10%FDR (padj<0.1) -  jk using 0.05
DEG_mid_vs_high.sig.num
# 14 DEGs 
# why is it different than when control is included?
DEG_mid_vs_high.sig <- subset(DEG_mid_vs_high, padj <0.05) # identify and subset significant pvalues
DEG_mid_vs_high.sig.list <- data[which(rownames(data) %in% rownames(DEG_mid_vs_high.sig)),] # subsey list of significant genes from original count data 
DEG_mid_vs_high.sig.list$contrast <- as_factor(c("Treatment_high_vs_mid")) # set contrast as a factor 
SFtest <- estimateSizeFactors(DEG_mid_vs_high.sig.list)
print(sizeFactors(SFtest))
DEG_mid_vs_high.vst.sig <- varianceStabilizingTransformation(DEG_mid_vs_high.sig.list, blind = FALSE) # apply a regularized log transformation to minimize effects of small counts and normalize wrt library 
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
DEG_mid_vs_high.sig.list <- as.data.frame(counts(DEG_mid_vs_high.sig.list))
DEG_mid_vs_high.sig.list["Treatment_Compare"] <- "MidsHigh" # adding treatment comparison column
#write.csv(DEG_mid_vs_high.sig.list, file = "~/Desktop/pdam_mid_vs_high_DEG.csv")

#write.csv(DEGs.all, file = "~/Desktop/pdam_DEGs.all_treatment.csv")
DEGs.all <- select(DEG_mid_vs_high.sig.list, -Treatment_Compare)
DEGs.all$DEGs <- rownames(DEGs.all)
DEGs.all <- unique(DEGs.all) # 4 unique DEGs
#unique.sig.num <- length(t(unique(DEGs.all)))

# PCA plot of diff-expressed genes 
plob_DEGPCAdata <- plotPCA(DEG_mid_vs_high.vst.sig, intgroup = c("Treatment"), returnData=TRUE)
percentVar_pca_plob <- round(100*attr(plob_DEGPCAdata, "percentVar")) #plot PCA of samples with all data
plob_DEGPCAplot <- ggplot(plob_DEGPCAdata, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_pca_plob[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_pca_plob[2],"% variance")) +
  scale_color_manual(values = c(mid = "pink", high = "darkgreen")) +
  coord_fixed() +
  ggtitle("Plob") + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background
plob_DEGPCAplot
# PCA plot is of differentially expressed genes only
PC.info <- plob_DEGPCAplot$data
ggsave("~/Desktop/plob_treatment_DEGs_PCA.pdf", plob_DEGPCAplot)




























