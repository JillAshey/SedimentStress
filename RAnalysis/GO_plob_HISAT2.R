# Title: GO with P. lobata samples
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Code for Francois sedimentation data. P. lobata only samples analyzed here. HISAT2 was read aligner.
# GO analysis

# Load libraries 
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
library("forcats")
library("gridExtra")

# Obtain all expressed genes (PoverA = 0.85, 5)
countdata_plob <- read.csv("~/Desktop/genecount_plob_filtered.csv", header=TRUE, sep = ",") 
dim(countdata_plob) # should be 12074
colnames(countdata_plob) <- c("gene_id",     "X6_1",  "X7_1",  "X8_1",  "X9_1",  "X21_1", "X22_1", "X23_1", "X25_1", "X26_1", "X27_1", "X29_1", "X34_1")
head(countdata_plob)

# Obtain all differentially expressed genes 
DEG_plob <- read.csv("~/Desktop/DEG_plob_list_alltreatments_save.csv", header = TRUE)
dim(DEG_plob) # should be 369
colnames(DEG_plob) <- c("gene_id",     "X6_1",  "X7_1",  "X8_1",  "X9_1",  "X21_1", "X22_1", "X23_1", "X25_1", "X26_1", "X27_1", "X29_1", "X34_1")
head(DEG_plob)

# Input merged annotated gtf file 
##### USE MERGED GTF FILE THAT ONLY CONTAINS THE KNOWN PLOB SAMPLES - GO BACK IN BLUEWAVES AND DO THIS 
map_plob <- read.csv(file="~/Desktop/plobSamples_merged.annotated.gtf", header=FALSE, sep="\t", skip=2) #load sample info
map_plob <- separate(map_plob, V9, into = c("transcript_id", "gene_id", "gene_name", "xloc", "ref_gene_id", "cmp_ref", "class_code", "tss_id"), sep = ";")
map_plob <- subset(map_plob, V3=="transcript")
map_plob <- map_plob[,c(1,4,5,9:11)]
map_plob$gene_id <- gsub("gene_id ","",map_plob$gene_id) #remove extra characters
map_plob$gene_id <- gsub(" ","",map_plob$gene_id) #remove extra characters
map_plob$transcript_id <- gsub("transcript_id ","",map_plob$transcript_id) #remove extra characters
map_plob$gene_name <- gsub("gene_name ","",map_plob$gene_name) #remove extra characters
map_plob$gene_name <- gsub("gene_name ","",map_plob$gene_name) #remove extra characters
map_plob$gene_name <- gsub(" ","",map_plob$gene_name) #remove extra characters
map_plob$gene_name <- gsub("xlocXLOC_[0-9][0-9][0-9][0-9][0-9][0-9]", "unknown", map_plob$gene_name)
colnames(map_plob) <- c("scaffold", "start", "stop", "transcript_id", "gene_id", "gene_name")

#Remove duplicate rows
map_plob <- map_plob %>% mutate_all(na_if," ")
map_plob <- map_plob[,c(5,6)]
map_unique <- unique(map_plob)
dim(map_unique)
tail(map_unique)

# Import the reference annotation file
ref_annotations <- read.table("~/Desktop/GFFs/Plut.GFFannotation.fixed_transcript.gff", header = FALSE)
ref_annotations <- ref_annotations %>%
  filter(V3 == 'gene')
ref_annotations <- ref_annotations[,c(1,4,5,9)]
colnames(ref_annotations) <- c("scaffold", "start", "stop", "gene_name")
ref_annotations$Length <- ref_annotations$stop  - ref_annotations$start
ref_annotations$gene_name <- gsub("ID=","",ref_annotations$gene_name) #remove extra characters
ref_annotations$gene_name<-gsub(";.*", "", ref_annotations$gene_name) #from gene name remove everything after the first ;
dim(ref_annotations)
head(ref_annotations)



## Build GOSeq dataframe 

# Build a dataframe that links the gene IDs of expressed  genes (poverA = 0.85,5), the gene names of those genes (from the gene map), and the gene lengths (from the annotation file)
# Find gene names associated with the gene IDs of expressed planuala genes (poverA = 0.85,5)
genes.map_unique <- merge(countdata_plob, map_unique, by.x="gene_id") # not all samples coming in here
dim(genes.map_unique) # 221 by 14
genes.map_unique <- genes.map_unique[,-c(2:10)]
genes.map_unique$gene_name <- replace_na(genes.map_unique$gene_name, "unknown")
tail(genes.map_unique)

# Find gene positions in ref corresponding to expressed planuale genes (poverA = 0.85,5)
# that only has the expressed planuale genes (poverA = 0.85,5)
genes.map_unique.ref <- merge(genes.map_unique, ref_annotations, by.x="gene_name") 
dim(genes.map_unique.ref) # 221 by 9
#removed all genes with unknown names
tail(genes.map_unique.ref)

# make gene IDs unique 
genes.map_unique.ref$gene_id <- make.unique(as.character(genes.map_unique.ref$gene_id), sep = ".")




#### Build GOSEQ vector 

# GOseq requires a vector of all genes and all differentially expressed genes. 
# Make gene vector
DEG <- filter(genes.map_unique.ref, gene_id%in%DEG_plob$gene_id) #make vector of differentially expressed genes
dim(DEG) # 15 by 9 ???
DEG_names <- as.vector(DEG$gene_id)

## Not very successful here, may be because it does not have GO terms in annotation file?

