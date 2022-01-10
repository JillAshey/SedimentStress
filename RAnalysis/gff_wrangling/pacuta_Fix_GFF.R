# Title: P. acuta GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 01/10/2022

# Need to do some pacuta gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Pacuta.gff <- read.csv(file="~/Desktop/GFFs/pacuta/braker_v1.gff3", header=FALSE, sep="\t") 

#rename columns
colnames(Pacuta.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Pacuta.gff)

# Make col for transcript ids
Pacuta.gff$transcript_id <- sub(";.*", "", Pacuta.gff$gene)
Pacuta.gff$transcript_id <- gsub("ID=", "", Pacuta.gff$transcript_id) #remove ID= 

#remove everything after the second . in the gene column
Pacuta.gff$transcript_id <-sub("^([^.]*.[^.]*).*", "\\1", Pacuta.gff$transcript_id) #remove everything after the third . in the gene column
head(Pacuta.gff) 

# Check what kinds of ids are in gff
unique(Pacuta.gff$id)
# "start_codon" "CDS"         "exon"        "intron"      "stop_codon" 
# no ids == gene...so the transcript id identifier will be added to all ids?

#If id == start_codon, CDS, exon, intron, stop_codon, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Pacuta.gff <- Pacuta.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, "transcript_id=", Pacuta.gff$transcript_id),  paste0(gene)))
head(Pacuta.gff)

# Remove transcript_id col
Pacuta.gff <- Pacuta.gff[,-10]
head(Pacuta.gff) 

#save file
write.table(Pacuta.gff, file="~/Desktop/GFFs/pacuta/Pacuta.gff.annotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)




