# Title: DESeq2 with Mcap samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/27/20

# Code for Francois sedimentation data. Mcap only samples analyzed here aligned against Mcap STAR was read aligner with annotation file. 
# Used gffreads to turn gff into gtf so it could be properly read by star

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
mcap_counts <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/star/gene_count_mcap_gtf_matrix.csv", header = TRUE, row.names = "gene_id")
dim(mcap_counts) # 63227 x 11
for ( col in 1:ncol(mcap_counts)){
  colnames(mcap_counts)[col] <-  sub(".fastq.trim.fq.Aligned.sortedByCoord.out.bam.merge.gtf", "", colnames(mcap_counts)[col])
}
for ( col in 1:ncol(mcap_counts)){
  colnames(mcap_counts)[col] <-  gsub("X", "", colnames(mcap_counts)[col])
}
annot <- read.csv("~/Desktop/GFFs/mcap.annotation.gtf", header = FALSE, sep = "\t")
colnames(annot) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
unique(annot$id) # [1] "transcript" "exon"      
annot <- annot %>% separate(attr, c("transcript_id", "gene"), sep = ";")
annot$gene <-gsub("gene_id", "", annot$gene)
annot$gene <-gsub(" ", "", annot$gene)

# Load metadata
metadata <- read.csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/sediment_HI_metadata_raw.csv", header = TRUE)
dim(metadata) # 65 x 9
head(metadata)

