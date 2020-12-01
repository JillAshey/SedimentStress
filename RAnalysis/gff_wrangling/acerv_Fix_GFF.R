# Title: A. cervicornis GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/01/20

# Need to do some acerv gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="~/Desktop/GFFs/Acerv_assembly_v1.0.gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 

# Checking what kinds of ids are in gff
unique(Acerv.gff$id)
# [1] "gene"        "mRNA"        "exon"        "CDS"         "start_codon" "stop_codon"  "tRNA"       

#If id == mRNA, exon, start_codon, stop_codon, CDS, tRNA, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id),  paste0(gene)))
head(Acerv.gff)

# Remove last col
Acerv.gff <- Acerv.gff[,-10]
head(Acerv.gff)  

#save file
write.table(Acerv.gff, file="~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
