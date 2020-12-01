# Title: P. damicornis GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Using the GFF file from reef genomics. fule is not properly running in stringTie after HISAT2, so I'm going to make some adjustments to see if that can fix it
# I'll be adding gene_id= to 'gene' column


#Load libraries
# library(tidyverse)
# 
# pdam.gff <- read.csv(file="~/Desktop/pdam_annotation.gff3", header=FALSE, sep="\t", skip=2) 
# 
# #rename columns
# colnames(pdam.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
# 
# # seeing what kinds of components are in annotation
# unique(pdam.gff$id)
# 
# # Making ref_gene_id, moving gene id into separate column
# pdam.gff$gene_id <- sub(";.*", "", pdam.gff$gene)
# pdam.gff$gene_id <- gsub("ID=", "transcript_id=", pdam.gff$gene_id) #remove ID= and replace with transcript_id=
# # or could do: pdam.gff$gene_id <- gsub("ID=", "gene_id=", pdam.gff$gene_id) #remove ID= and replace with gene_id=
# head(pdam.gff)
# 
# # paste gene id back into gene column
# pdam.gff$gene <- paste(pdam.gff$gene, pdam.gff$gene_id)
# head(pdam.gff$gene)
# 
# # remove white spaces
# pdam.gff$gene <- gsub(" ", "", pdam.gff$gene, fixed = TRUE) #remove white space
# head(pdam.gff)
# 
# # Remove last column 
# pdam.gff <- pdam.gff[,-10]
# 
# #save file
# write.table(pdam.gff, file="~/Desktop/pdam_annotation_AddTranscript_id_fixed.gtf", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
# #write.table(pdam.gff, file="~/Desktop/pdam_annotation_Addgene_id_fixed.gtf", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)


## May not need this info anymore because I realized my mistake in HISAT2. Because I am now using the reef genomics gff file, I have to rebuild hisat2 index
# with the reef genomics fasta file


























