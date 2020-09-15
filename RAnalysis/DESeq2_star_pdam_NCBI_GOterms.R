# Title: DESeq2 with pdam samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. dam only samples analyzed here aligned against P. dam. STAR was read aligner with gff annotation file from NCBI.
# I edited NCBI file to include GO terms

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
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_pdam_GOterms_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 26336 x 64
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
# metadata$Treatment <- gsub("<NA>", "unknown", metadata$Treatment)
rownames(metadata) <- metadata$SampleID # make sampleID the row names 
metadata <- metadata[-65,] # remove random blank space at the end
metadata <- na.omit(metadata) # removing rows with NAs
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
dim(gkeep) # 17326 genes left after filtering

# List names of genes that passed filtering 
gn.keep <- rownames(gkeep)

# gene count data that was filtered in PoverA (P percent of samples that have counts over A) + gene names 
gcount_filt <- as.data.frame(count_pdam[which(rownames(count_pdam) %in% gn.keep),]) # only keep gene names that are in gn.keep
head(gcount_filt)
dim(gcount_filt) # 17326 x 12 -- only 17326 genes kept after filtering 

# Write treatment, gene and transcript count files with corrected column and row headers
# write.csv(metadata_pdam, "~/Desktop/metadata_pdam_filtered.csv")
write.csv(gcount_filt, "~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/genecount_pdam_GOterms_star_filtered.csv")

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
gPCAdata <- plotPCA(gvst_pdam, intgroup = c("Treatment"), returnData=TRUE) # create PCA loadings ?
percentVar <- round(100*attr(gPCAdata, "percentVar")) # calculating % variance for PCA axis titles ??
#plot PCA of samples with all data
pdam_PCAplot <- ggplot(gPCAdata, aes(PC1, PC2, color=Treatment, shape=Days)) + 
  geom_point(size=3) +
  # geom_text(aes(label=name),hjust=0, vjust=0) +
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
#pdam_heatmap_PCA <- grid.arrange(pdam_PCAplot, pdam_heatmap[[4]], nrow=2, clip="off")
#ggsave("~/Desktop/pdam_heatmap_NCBI_star_PCA.pdf", pdam_heatmap_PCA, width = 8, height = 8, units = c("in"))



## DGE analysis 

# Run DE analysis 
# Use Wald model 
DEG_pdam <- DESeq(gdds_pdam) #run differential expression test by group using the Wald model
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
res_DEG_pdam <- results(DEG_pdam)
res_DEG_pdam_Ordered <- res_DEG_pdam[order(res_DEG_pdam$pvalue),]

DEG_pdam$Treatment
resultsNames(DEG_pdam)
# [1] "Intercept"                 "Treatment_mid_vs_control"  "Treatment_high_vs_control"

# Explore significant p-values for treatments 
# Control vs mid
DEG_pdam_results_control_vs_mid <- results(DEG_pdam, name = "Treatment_mid_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_mid <- DEG_pdam_results_control_vs_mid[order(DEG_pdam_results_control_vs_mid$pvalue),] # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_mid) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_mid <- sum(DEG_pdam_results_control_vs_mid$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
# 168 significantly differentially expressed genes between control and mid that are less than 0.05
pdam_DEGs.control_vs_mid <- subset(DEG_pdam_results_control_vs_mid, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_mid <- as.data.frame(pdam_DEGs.control_vs_mid) # make df
pdam_DEGs.control_vs_mid$contrast <- as_factor(c("control_vs_mid")) # set contrast as a factor 
pdam_DEGs.control_vs_mid <- cbind(gene_id = rownames(pdam_DEGs.control_vs_mid), pdam_DEGs.control_vs_mid) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_mid) <- NULL # remove row names 
pdam_DEGs.control_vs_mid
dim(pdam_DEGs.control_vs_mid) # 168 by 8
write.csv(pdam_DEGs.control_vs_mid, "~/Desktop/pdam_DEGs.control_vs_mid.csv")
pdam_sig.num.control_vs_mid

# Control vs high 
DEG_pdam_results_control_vs_high <- results(DEG_pdam, name = "Treatment_high_vs_control") # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_control_vs_high <- order(DEG_pdam_results_control_vs_high$pvalue) # order from smallest pvalue 
summary(DEG_pdam_results_control_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.control_vs_high <- sum(DEG_pdam_results_control_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
pdam_sig.num.control_vs_high # 328 significantly differentially expressed genes between control and mid 
pdam_DEGs.control_vs_high <- subset(DEG_pdam_results_control_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.control_vs_high <- as.data.frame(pdam_DEGs.control_vs_high) # make df
pdam_DEGs.control_vs_high$contrast <- as_factor(c("control_vs_high")) # set contrast as a factor
pdam_DEGs.control_vs_high <- cbind(gene_id = rownames(pdam_DEGs.control_vs_high), pdam_DEGs.control_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.control_vs_high) <- NULL # remove row names 
pdam_DEGs.control_vs_high
dim(pdam_DEGs.control_vs_high) # 328 by 8
write.csv(pdam_DEGs.control_vs_high, "~/Desktop/pdam_DEGs.control_vs_high.csv")
pdam_sig.num.control_vs_high

# Mid vs high
DEG_pdam_results_mid_vs_high <- results(DEG_pdam, contrast = c("Treatment", "mid", "high")) # results only for control vs mid treatments 
results_ordered_DEG_pdam_results_mid_vs_high <- order(DEG_pdam_results_mid_vs_high$pvalue) #Order p-values by smallest value first
summary(DEG_pdam_results_mid_vs_high) # view summary of results with adj p < 0.1
pdam_sig.num.mid_vs_high <- sum(DEG_pdam_results_mid_vs_high$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
pdam_sig.num.mid_vs_high # 3 significantly differentially expressed genes between control and mid 
pdam_DEGs.mid_vs_high <- subset(DEG_pdam_results_mid_vs_high, padj<0.05) # subset only <0.05 padj values
pdam_DEGs.mid_vs_high <- as.data.frame(pdam_DEGs.mid_vs_high) # make df
pdam_DEGs.mid_vs_high$contrast <- as_factor(c("mid_vs_high")) # set contrast as a factor
pdam_DEGs.mid_vs_high <- cbind(gene_id = rownames(pdam_DEGs.mid_vs_high), pdam_DEGs.mid_vs_high) # make gene id a row and bind it to the rest of the df
rownames(pdam_DEGs.mid_vs_high) <- NULL # remove row names 
pdam_DEGs.mid_vs_high
dim(pdam_DEGs.mid_vs_high) # 3 by 8
write.csv(pdam_DEGs.mid_vs_high, "~/Desktop/pdam_DEGs.mid_vs_high.csv")
pdam_sig.num.mid_vs_high

DEGs_pdam_contrast_all_treatments <- bind_rows(pdam_DEGs.control_vs_mid, pdam_DEGs.control_vs_high, pdam_DEGs.mid_vs_high) # bind mid_vs_control results, high_vs_control results, and mid_vs_high results together by row  
dim(DEGs_pdam_contrast_all_treatments) # 499 by 8
DEG_pdam.sg.num <- sum(DEGs_pdam_contrast_all_treatments$padj < 0.05, na.rm=TRUE) #How many adjusted p-values were less than 0.05?
DEG_pdam.sg.num # 499 adjusted p-values 
summary(DEGs_pdam_contrast_all_treatments)
write.csv(DEGs_pdam_contrast_all_treatments, file="~/Desktop/DEGs_pdam_contrast_all_treatments.csv")


## Visualize diff-expressed genes

# Subset and log-transform count data 
# Subset list of genes by those which padj>0.
dim(DEGs_pdam_contrast_all_treatments)
DEGs_pdam <- DEGs_pdam_contrast_all_treatments$gene_id # list all gene names 
DEGs_pdam <- unique(DEGs_pdam) # select only unique gene names 
DEG_pdam_list <- gdds_pdam[which(rownames(gdds_pdam) %in% DEGs_pdam)] # filter gdds_pdam DESeq2 object by unique gene names
dim(DEG_pdam_list) # 368 x 12
print(counts(DEG_pdam_list))

# As determined above, size factors all less than 4, so proceed with VST
#apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size
DEGvst <- vst(DEG_pdam_list, blind = FALSE, nsub = nrow(counts(DEG_pdam_list)))
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
dim(DEGvst) # 368 by 12 
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
write.csv(counts(DEG_pdam_list), file="~/Desktop/pdam_DEG_list_unique.csv")

pdam_DEGs_heatmap_PCA <- grid.arrange(pdam_DEGPCAplot, pdam_DEGheatmap[[4]], nrow=2, clip="off")
ggsave("~/Desktop/pdam_DEGs_heatmap_PCA.pdf", pdam_DEGs_heatmap_PCA, width = 8, height = 8, units = c("in"))




## GOseq

# Make count tables into vectors so they can be read by goseq
gene.vector=as.vector(count_pdam)
#names(gene.vector)=assayed.genes
head(gene.vector)
dim(gene.vector)
DEG.vector <- as.vector(DEGs_pdam_contrast_all_treatments)
head(DEG.vector)
dim(DEG.vector)

# goseq has genomes on file, but no coral genomes so this info is not relevant to me
supportedOrganisms()
supportedGenomes()

# To run goseq, I first need 

# need list of GO terms, DEG in binary, all genes (names), 

# I need a list of all genes examined and DEGs for control vs mid 
control_vs_mid_names <- pdam_DEGs.control_vs_mid$gene_id # DEGs for control vs mid 

# gff annotations. need to get gene id from here 
all_annot <- read.table("~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
all_annot <- all_annot[-1,]
all_annot$Start <- as.numeric(all_annot$Start)
all_annot$End <- as.numeric(all_annot$End)
all_annot$gene_id <- sub(";.*", "", all_annot$Attr)
all_annot$gene_id <- gsub("ID=", "", all_annot$gene_id) #remove ID= 
all_annot <- mutate(all_annot, len = End - Start) # calculate gene length
gene_lengths <- subset(all_annot[,c("gene_id", "len")]) # make df with just gene id and length in it 

# now I need to edit DEGs so that gene_id names match - this means removing everything that comes after | in the rows that begin with gene or STRG
# But going to leave out STRG for now, as those are novel loci
rna <- pdam_DEGs.control_vs_mid %>%
  filter(!str_detect(gene_id, 'gene'))
rna <- rna %>%
  filter(!str_detect(gene_id, 'STRG')) # only rna, will join with gene once i remove the |
gene <- pdam_DEGs.control_vs_mid %>%
  filter(!str_detect(gene_id, 'rna'))
gene <- gene %>%
  filter(!str_detect(gene_id, 'STRG')) 
gene$gene_id <- gsub("\\|.*", "", gene$gene_id)
pdam_DEGs.control_vs_mid <- rbind(rna, gene) 



m <- merge(pdam_DEGs.control_vs_mid, gene_lengths, by = "gene_id", all=TRUE)
m$lfcSE <- gsub("0.*", "1", m$lfcSE)
m$lfcSE <- gsub("1.*", "1", m$lfcSE)
m[is.na(m)] <- 0
DEGs_mid_control <- select(m, c("gene_id", "lfcSE"))
colnames(DEGs_mid_control) <- c("gene_id", "DEgenes")

DEGs_mid_control$DEgenes <- as.numeric(DEGs_mid_control$DEgenes)
test<- nullp(DEGs_mid_control$DEgenes, DEGs_mid_control$gene_id, bias.data=gene_lengths$len)
# now ready to run go seq






finaltable$finaltable













  
GO_annot <- read.table("~/Desktop/GFFs/go_only.gff3", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))


hello <- GO_annot %>%
  filter(!str_detect(Attr, 'GO'))

list_of_annots <- lapply(GO_annot$Attr, function(x) strsplit(x,split=";")[[1]])
list_of_ontologies <- lapply(list_of_annots, grep, value=TRUE, pattern="Ontology_term")
names(list_of_ontologies) <- unlist(lapply(list_of_annots, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))
list_of_ontologies <- lapply(list_of_ontologies, function(x) strsplit(strsplit(x, "=")[[1]][-1],",")[[1]])





GO_annot$gene_id <- sub(";.*", "", GO_annot$Attr)
GO_annot$gene_id <- gsub("ID=", "", GO_annot$gene_id) #remove ID= 

hello <- GO_annot %>%
  select(!str_detect(Attr, 'GO'))

blahagain <- str_split(GO_annot$Attr, ";")
blahagain <- as.data.frame(blahagain)

list_of_all_annots <- lapply(GO_annot$Attr, function(x) strsplit(x,split=";")[[1]])

list_of_all_annots <- lapply(GO_annot$Attr, function(x) strsplit(x,split=";"))



# ID=exon-XM_027179716.1-4 ;
# Parent=rna-XM_027179716.1 ;
# Dbxref=GeneID:113664116,Genbank:XM_027179716.1 ;
# gbkey=mRNA ;
# gene=LOC113664116 ;
# product=cytochrome P450 4V2-like ;
# transcript_id=XM_027179716.1 ;
# pdam_00021920 ;
# Pocillopora damicornis ;
# 113664116 ;
# cytochrome P450 4V2-like ;
# GO:0005506,GO:0016705,GO:0020037,GO:0055114                                                                                                                                                                                                                                   

blah <- separate(GO_annot, Attr, into = c("id", "parent", "dbxref", "gbkey", "gene", "product", "transcript_id", "pdam", "species", "Number", "Kind", "GO"), sep = ";")





GO_annot$test <-sub("*.([^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;])", "\\1", GO_annot$Attr) 
GO_annot$test <-sub(".*([^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.[^;]*.)", "", GO_annot$Attr) 


#remove everything after the third . in the gene column


Plut.gff$transcript_id <-sub("^([^;]*.[^;]*.[^;]*).*", "\\1", Plut.gff$transcript_id) #remove everything after the third . in the gene column



LPS_goseq_res <- goseq(LPS_pwf, names(gene_lengths), gene2cat = annotation_df, method="Wallenius", use_genes_without_cat=TRUE)











# import gff with GO terms in it 
all_annot <- read.table("~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
all_annot <- all_annot[-1,]
all_annot$Start <- as.numeric(all_annot$Start)
all_annot$End <- as.numeric(all_annot$End)
gene_lengths <- mutate(all_annot, len = End - Start)$len
list_of_all_annots <- lapply(all_annot$Attr, function(x) strsplit(x,split=";")[[1]])
names(gene_lengths) <- unlist(lapply(list_of_all_annots, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))





head(count_pdam)
genes <- unlist(lapply(rownames(rawdat), function(x) strsplit(x, '-')[[1]][1])) # may need to take this part out



# create DESeq2 object 
# metadata_pdam$Treatment <- relevel(metadata_pdam$Treatment, ref="control")
count_pdam_matrix <- as.matrix(count_pdam)
dds <- DESeqDataSetFromMatrix(count_pdam_matrix, colData = metadata_pdam, design = ~Treatment)

## Retrieve out the original counts
rawdat <- counts(dds)
dim(rawdat)

test <- count(gdds_pdam)


## need to compare all gene names to the DE gene names from control vs mid 

















##### Polina Code - she did work with GO terms with Connelly data, so I am going to try to use her code to replicate that 

library("DESeq2")
library("tximport")
library("tidyr")
library("dplyr")
library("readr")
#library("ggpubr")
library("ggrepel")
library("genefilter")
library("goseq")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

# from above
# Metadata pdam
head(metadata_pdam)
dim(metadata_pdam)
# count data pdam - making sure i have gene_id as row names here
countdata <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_pdam_GOterms_matrix.csv", header = TRUE, row.names = "gene_id")
dim(countdata) # 26336 x 64
head(countdata)
for ( col in 1:ncol(countdata)){
  colnames(countdata)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(countdata)[col])
}
pdam_ID <- metadata_pdam$SampleID
count_pdam <- select(countdata, all_of(pdam_ID))
all(rownames(metadata_pdam) %in% colnames(count_pdam)) # must come out TRUE

# create DESeq2 object 
# metadata_pdam$Treatment <- relevel(metadata_pdam$Treatment, ref="control")
count_pdam_matrix <- as.matrix(count_pdam)
dds <- DESeqDataSetFromMatrix(count_pdam_matrix, colData = metadata_pdam, design = ~Treatment)

## Retrieve out the original counts
rawdat <- counts(dds)
dim(rawdat)

## Create filters
filt <- filterfun(pOverA(0.85,5))
gfilt <- genefilter(rawdat, filt)
dds <- dds[gfilt,]

# Run DESeq2 & look at results 
dds <- DESeq(dds)
dds
dim(dds)
# more filtering
keep <- rowSums(counts(dds), na.rm = TRUE) >= 10
dds <- dds[keep,]
dds

results <- results(dds)
results

resultsNames(dds)
summary(results)

# Manual size factor checking
SF.dds <- estimateSizeFactors(dds)
# use VST
gvst <- vst(dds, blind=FALSE)

# plot heatmap of sample-sample distances
sampleDists <- dist(t(assay(gvst)))
sampleDistsMatrix <- as.matrix(sampleDists)
rownames(sampleDistsMatrix) <- colnames(gvst)
colnames(sampleDistsMatrix) <- NULL
colors <- colorRampPalette (rev(brewer.pal(0, "Blues")) )(255)

# save to a pdf
pdf("connelly_pdam_heatmaps.pdf")
pheatmap(sampleDistsMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

### PCA plots
PCAdat <- plotPCA(gvst, intgroup = c("Treatment"), returnData = TRUE)
percentVar <- round(100*attr(PCAdat, "percentVar"))
pdf("connelly_pdam_PCA.pdf")


ggplot(PCAdat, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #  scale_color_manual(values = c(Control="cadetblue3", LPS="palevioletred", Antibiotics="darkgreen", Heat="purple", Antibiotics_Heat="red1", Antibiotics_Heat_LPS="orange3")) +
  coord_fixed() + 
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_rect(fill = "white")) #Set the plot background
dev.off()



### Actual differential expression analysis
dds <- DESeq(dds)
control_mid_res <- results(dds, contrast = c("Treatment", "control", "mid"))
control_high_res <- results(dds, contrast = c("Treatment", "control", "high"))
mid_high_res <- results(dds, contrast = c("Treatment", "mid", "high")) 
head(as.data.frame(control_mid_res) %>% arrange(padj) %>% filter(padj <= 0.05))

## Make vectors for goseq
genes <- unlist(lapply(rownames(rawdat), function(x) strsplit(x, '-')[[1]][1])) # may need to take this part out
length(genes)
length(unique(genes))


# Control vs mid DE genes
df_control_mid <- as.data.frame(control_mid_res) %>% arrange(padj)
df_control_mid$name_t <- rownames(df_control_mid)
df_control_mid <- df_control_mid %>% separate(name_t, "go_name", sep='-', extra="drop")
#df_con_high <- df_con_high %>% mutate(go_annots = list_of_ontologies[go_name])
control_mid_sig_genes <- filter(df_control_mid, padj < 0.05)$go_name

# Control vs high DE genes
df_control_high <- as.data.frame(control_high_res) %>% arrange(padj)
df_control_high$name_t <- rownames(df_control_high)
df_control_high <- df_control_high %>% separate(name_t, "go_name", sep='-', extra="drop")
#df_con_mid <- df_con_mid %>% mutate(go_annots = list_of_ontologies[go_name])
control_high_sig_genes <- filter(df_control_high, padj < 0.05)$go_name

# Mid vs high DE genes 
df_mid_high <- as.data.frame(mid_high_res) %>% arrange(padj)
df_mid_high$name_t <- rownames(df_mid_high)
df_mid_high <- df_mid_high %>% separate(name_t, "go_name", sep='-', extra="drop")
#df_con_mid <- df_con_mid %>% mutate(go_annots = list_of_ontologies[go_name])
mid_high_sig_genes <- filter(df_mid_high, padj < 0.05)$go_name

length(control_mid_sig_genes)
length(control_high_sig_genes)
length(unique(c(control_high_sig_genes,control_mid_sig_genes)))
length(mid_high_sig_genes)

control_mid_OE_gv <- as.integer(genes%in%c(control_mid_sig_genes))
control_high_OE_gv <- as.integer(genes%in%c(control_high_sig_genes))
mid_high_OE_gv <- as.integer(genes%in%c(mid_high_sig_genes))
names(control_mid_OE_gv) <- genes
names(control_high_OE_gv) <- genes
names(mid_high_OE_gv) <- genes

# Get all of the GO Annotation information from the gff3
# How to make these files:
# ** for (tab) instead type ctrl + v and then tab. This enters a 'real' tab character. I don't know if this is
#   different for macs.
# tmp.gff3 -> grep "maker(tab)gene" pdam_annotation.gff3 > tmp.gff3
# go_only.gff3 -> grep "GO:" tmp.gff3 > go_only.gff3
GO_annot <- read.table("~/Desktop/GFFs/go_only.gff3", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
list_of_annots <- lapply(GO_annot$Attr, function(x) strsplit(x,split=";")[[1]])
list_of_ontologies <- lapply(list_of_annots, grep, value=TRUE, pattern="Ontology_term")
names(list_of_ontologies) <- unlist(lapply(list_of_annots, function(x) strsplit(grep("ID=",x, value=TRUE),"=")[[1]][2]))
# list_of_ontologies <- lapply(list_of_ontologies, function(x) strsplit(strsplit(x, "=")[[1]][-1],",")[[1]])

# Get all of the gene length information from another gff3 (a less filtered one)
all_annot <- read.table("~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms.gff", sep="\t", quote="", stringsAsFactors=FALSE, col.names=c("ID","Source","Type","Start","End","Score","Strand","Phase","Attr"))
all_annot <- all_annot[-1,]
all_annot$Start <- as.numeric(all_annot$Start)
all_annot$End <- as.numeric(all_annot$End)
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

















