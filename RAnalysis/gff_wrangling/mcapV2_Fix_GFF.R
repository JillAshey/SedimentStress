# Title: Mcap GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 02/17/2022

# Need to do some mcap v2 gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Mcap.gff <- read.csv(file="~/Desktop/GFFs/mcap/V2/Montipora_capitata_HIv1.genesNoCopies.gff3_polished", header=FALSE, sep="\t") 

#rename columns
colnames(Mcap.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Mcap.gff)

# This gff is a bit annoying to work with because it doesn't have consistent organization in the 'gene' column. for example, it
# has gene = xxxxx, exon = xxx.t1.exon, cds = cds.xxx.t1. So I'm going to mutate it so that I can pull out the relevant info for each row
#make transcript id col that is the same as gene col
Mcap.gff$transcript_id <- Mcap.gff$gene
#for everything other than gene, subset everything after 'Parent='. if gene, subset everything before the first ';' 
Mcap.gff <- Mcap.gff %>% 
  mutate(transcript_id = ifelse(id != "gene", sub(".*Parent=", "", Mcap.gff$transcript_id), sub(";.*", "", Mcap.gff$transcript_id)))
#for everything other than gene, subset everything after ';'. if gene, remove 'ID='
Mcap.gff <- Mcap.gff %>% 
  mutate(transcript_id = ifelse(id != "gene", sub(";.*", "", Mcap.gff$transcript_id), sub("ID=", "", Mcap.gff$transcript_id)))
Mcap.gff$transcript_id <- gsub(".t1", "", Mcap.gff$transcript_id) #remove .t1



Mcap.gff <- Mcap.gff %>% 
  mutate(gene = ifelse(id != "gene", sub(".*Parent=", "", Mcap.gff$gene), sub(";.*", "", Mcap.gff$gene)))

                       






#remove everything after the xxx in transcript_id col
## not sure if to remove anything yet...

# Check what kinds of ids are in gff
unique(Mcap.gff$id)
# "gene" "mRNA" "exon" "CDS" 

#If id == gene, mRNA, exon, or CDS, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Mcap.gff <- Mcap.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Mcap.gff$transcript_id),  paste0(gene)))
head(Mcap.gff)

# Remove transcript_id col
Mcap.gff <- Mcap.gff[,-10]
head(Mcap.gff) 

#save file
write.table(Mcap.gff, file="~/Desktop/GFFs/mcap/V2/Mcap.gff.annotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)




