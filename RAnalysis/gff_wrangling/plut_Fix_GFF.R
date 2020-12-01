# Title: P. lutea GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Need to do some mcap gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column

#Load libraries
library(tidyverse)

#Load  gene gff
Plut.gff <- read.csv(file="/Users/hputnam/Desktop/plut2v1.1.genes.gff3", header=FALSE, sep="\t", skip=1, nrows = 10) 

#rename columns
colnames(Plut.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Plut.gff)

Plut.gff$transcript_id <- sub(";.*", "", Plut.gff$gene)
Plut.gff$transcript_id <- gsub("ID=", "", Plut.gff$transcript_id) #remove ID= 

Plut.gff$transcript_id <-sub("^([^.]*.[^.]*.[^.]*).*", "\\1", Plut.gff$transcript_id) #remove everything after the third . in the gene column
head(Plut.gff)  

#If id == mRNA, exon, five_prime_UTR, three_prime_UTR, CDS, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Plut.gff$transcript_id),  paste0(gene)))

head(Plut.gff)

Plut.gff <- Plut.gff[,-10]
head(Plut.gff)     

#save file
write.table(Plut.gff, file="~/Desktop/Plut.GFFannotation.fixed_transcript.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
