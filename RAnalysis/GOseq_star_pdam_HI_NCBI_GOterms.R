# Title: GO with pdam Hawaii samples
# Project: Sedimentation GO-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. dam only samples analyzed here aligned against P. dam. STAR was read aligner with gff annotation file from NCBI.
# I edited NCBI file to include GO terms

## Goal: determine GO terms for DEGs

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

## Read data in 
# all expressed pdam genes filtered (PoverA 0.85 5)
pdam_gene_all <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/genecount_pdam_GOterms_star_filtered.csv", header = TRUE)
# all differentially expressed unique pdam genes 
pdam_DEG <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/pdam_DEG_list_unique.csv", header = TRUE)
# differentially expressed genes - control vs mid
pdam_DEG_control_vs_mid <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/pdam_DEGs.control_vs_mid.csv", header = TRUE)
# differentially expressed genes - control vs high
pdam_DEG_control_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/pdam_DEGs.control_vs_high.csv", header = TRUE)
# differentially expressed genes - mid vs high
pdam_DEG_mid_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/pdam_DEGs.mid_vs_high.csv", header = TRUE)
# NCBI annotation file edited to include pdam_xx terms and GO terms
# pdam_annotations <- read.table("~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms_sepcol.csv", header = TRUE) need to find way to save final table raw with GO terms so it actually splits by col
finaltable <- read.csv("~/Desktop/GFFs/finaltable.csv", header = TRUE)

# control vs mid -- focus on this set first 
## To run GOseq, I first need to calculate probability weighting function - gives certain 'weights' to set of genes based on biased data (gene lengths) and gene status (DEG y or n)
## I calculate probability weighting function with nullp, which requires list of all genes examined with the DEGs labelled as 1 and other genes as 0
## I need a list of all genes examined and DEGs for control vs mid 
pdam_DEG_control_vs_mid_names <- pdam_DEG_control_vs_mid # DEGs names for control vs mid 

# I need list of all genes that I looked at. I can get that from the ID= in finaltable
finaltable <- mutate(finaltable, len = gene.stop - gene.start) # calculate gene length
finaltable$gene_id <- sub(";.*", "", finaltable$gene) # isolate ID=
finaltable$gene_id <- gsub("ID=", "", finaltable$gene_id) #remove ID= 
gene_lengths <- subset(finaltable[,c("gene_id", "len")]) # make df with just gene id and length in it 
dim(gene_lengths)

# Need to edit the gene_id in the pdam_DEG_control_vs_mid_names. Some of them have a '|' and the LOC following. The ones that have this pattern are gene- and STRG (novel loci detected by STAR). Given I do not have annotations for these novel 
# loci, I am going to exclude STRG from GO analysis. Now I need to just remove the |LOCxxx after gene-. 
rna <- pdam_DEG_control_vs_mid_names %>% # removing all gene_ids that begin with gene
  filter(!str_detect(gene_id, 'gene'))
rna <- rna %>% # removing all gene_ids that begin with STRG
  filter(!str_detect(gene_id, 'STRG')) # only rna, will join with gene once i remove the | character from gene. 
gene <- pdam_DEG_control_vs_mid_names %>% # removing all gene_ids that begin with rna
  filter(!str_detect(gene_id, 'rna')) 
gene <- gene %>% # removing all gene_ids that begin with STRG
  filter(!str_detect(gene_id, 'STRG')) 
gene$gene_id <- gsub("\\|.*", "", gene$gene_id) # removing '|' and everything after it 
pdam_DEG_control_vs_mid_names <- rbind(rna, gene) # binding gene and rna names together again
# Now I have list of DEGs with proper gene_ids that can match to those in gene_lengths 

# Compare gene_ids in pdam_DEG_control_vs_mid_names to those in gene_lengths 
m <- merge(pdam_DEG_control_vs_mid_names, gene_lengths, by = "gene_id", all=TRUE) # merged them by gene_id, but there are still the other cols from pdam_DEG_control_vs_mid_names that show me which genes are diff expressed. 
# there is probably an easier and cleaner way to do this 
unique(m$lfcSE) 
# [1] NA           "1"          "1.22926851" "1.1"       
# because the lfcSE came with the DEG gene ids, I can use this to create binary for DEGs in all genes examined
# Everything labelled NA is not a diff expressed gene and gets a 0; everything labelled anything other than NA is a diff expressed gene and gets a 1
m$lfcSE <- gsub("1.*", "1", m$lfcSE)
m[is.na(m)] <- 0
DEGs_mid_control <- select(m, c("gene_id", "lfcSE", "len")) # select gene id and column with binary for gene_ids
colnames(DEGs_mid_control) <- c("gene_id", "DEgenes", "len") # rename so I know its represented DEGs binary 
dim(DEGs_mid_control)
dim(gene_lengths)
head(DEGs_mid_control)

# Now I can calculate pwf. Hooray
DEGs_mid_control$DEgenes <- as.numeric(DEGs_mid_control$DEgenes) # change to numeric 
pdam_mid_vs_control_pwf <- nullp(DEGs_mid_control$DEgenes, DEGs_mid_control$gene_id, bias.data=gene_lengths$len)

pdam_mid_vs_control_pwf <- nullp(DEGs_mid_control, gene_lengths$gene_id, bias.data=gene_lengths$len)

plotPWF(pdam_mid_vs_control_pwf) # doesnt look great...but not sure how to normalize 
dim(pdam_mid_vs_control_pwf)

# I now believe I can move onto goseq. I need the pwf file I just created and id names for genes. 
# also need go terms with associated gene names
go_annots <- subset(finaltable[,c("GOterms", "gene_id")])
dim(go_annots)
head(go_annots)
go_annots <- go_annots[,c(2,1)] # switching so that gene ids are first col

dim(go_annots)
dim(pdam_mid_vs_control_pwf)
dim(DEGs_mid_control)

testing <- goseq(pdam_mid_vs_control_pwf, DEGs_mid_control$gene_id, gene2cat = go_annots, method="Wallenius", use_genes_without_cat=TRUE)
# use getgo function?
# Using manually entered categories.
# Error in `[.default`(summary(map), , 1) : incorrect number of dimensions
# In addition: Warning message:
#   In goseq(pdam_mid_vs_control_pwf, DEGs_mid_control$gene_id, gene2cat = go_annots,  :
#              Gene column could not be identified in gene2cat conclusively, using the one headed gene_id

# Maybe I need to subset only DEGs with GO terms associated with them...?
go_annots[go_annots == ""] <- NA
go_annots_filtered <- na.omit(go_annots)

# Merge go_annots_filtered and gene_lengths by gene_ids
m_annot <- merge(go_annots_filtered, gene_lengths, by = "gene_id", all=TRUE) # merged by gene_id
m_annot_filtered <- na.omit(m_annot)
# Merge m_annot_filtered by DEGs_mid_control by gene_ids
full_annot_filtered <- merge(m_annot_filtered, DEGs_mid_control, by = "gene_id", all=TRUE) # merged by gene_id










LPS_goseq_res <- goseq(LPS_pwf, names(gene_lengths), gene2cat = annotation_df, method="Wallenius", use_genes_without_cat=TRUE)


## from manual 
# Need to obtain two vectors, one containing all genes assayed and one containing all DE genes 
de.genes <- pdam_DEG_control_vs_mid_names[,2]
assayed.genes <- gene_lengths[,1] # make df with just gene id and length in it 

gene.vector <- as.integer(assayed.genes%in%de.genes)
names(gene.vector) = assayed.genes
head(gene.vector)

# First, need to quantify length bias with probability weighting function

















