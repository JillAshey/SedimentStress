# Title: M. cavernosa GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/01/20

# Need to do some mcav gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Mcav.gff <- read.csv(file="~/Desktop/GFFs/Mcavernosa.maker.coding.gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Mcav.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Mcav.gff)

Mcav.gff <- read.csv(file="~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3", header=FALSE, sep="\t", skip=1) 

# Creating transcript id
Mcav.gff$transcript_id <- sub(";.*", "", Mcav.gff$gene)
Mcav.gff$transcript_id <- gsub("ID=", "", Mcav.gff$transcript_id) #remove ID= 

# Checking what kinds of ids are in gff
unique(Mcav.gff$id)
# [1] "gene"            "mRNA"            "exon"            "CDS"             ""                "five_prime_UTR"  "three_prime_UTR"

#If id == mRNA, exon, five_prime_UTR, three_prime_UTR, CDS, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Mcav.gff <- Mcav.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Mcav.gff$transcript_id),  paste0(gene)))
head(Mcav.gff)
# may want to edit row 9 to tidy up transcript_id= (ie take off the :exon, etc)

# Remove last col
Mcav.gff <- Mcav.gff[,-10]
head(Mcav.gff)  

#save file
write.table(Mcav.gff, file="~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)











Mcav.gff <- read.csv(file="~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3", header=FALSE, sep="\t", skip=1) 
colnames(Mcav.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Mcav.gff)
dim(Mcav.gff) # 370403 by 9

# Remove empty rows
Mcav.gff <- Mcav.gff %>% 
  filter(!str_detect(scaffold, "##"))

unique(Mcav.gff$id)
# [1] "mRNA"            "exon"            "CDS"             "gene"            "five_prime_UTR"  "three_prime_UTR"

# Find out many of each unique part is in the annotation file
mrna_Mcav.gff <- subset(Mcav.gff, id=="mRNA")
dim(mrna_Mcav.gff) # 25142 rows
exon_Mcav.gff <- subset(Mcav.gff, id=="exon")
dim(exon_Mcav.gff) # 151448 rows
CDS_Mcav.gff <- subset(Mcav.gff, id=="CDS")
dim(CDS_Mcav.gff) # 144270 rows
gene_Mcav.gff <- subset(Mcav.gff, id=="gene")
dim(gene_Mcav.gff) # 25141 rows
five_prime_UTR_Mcav.gff <- subset(Mcav.gff, id=="five_prime_UTR")
dim(five_prime_UTR_Mcav.gff) # 11575 rows
three_prime_UTR_Mcav.gff <- subset(Mcav.gff, id=="three_prime_UTR")
dim(three_prime_UTR_Mcav.gff) # 12827 rows

# calculating gene.dff = how many bp each section of genome is 
Mcav.gff <- Mcav.gff %>%
  mutate(Mcav.gff, gene.diff = gene.stop - gene.start)

# See which gene.diff are negative and which are positive
resultPos <- Mcav.gff %>% # df that is only positive gene diffs
  filter(gene.diff > 0)
dim(resultPos) #369640 by 10
resultNeg <- Mcav.gff %>% # df that is only negative gene diffs
  filter(gene.diff < 0)
dim(resultNeg) #650 by 10
# There shouldn't be any negative values in gene.diff because that makes no sense. How could a section be - bp?

# not a huge difference between the df with all genes and the df with genes that have + gene.diff
gene_resultPos <- subset(resultPos, id=="gene")
dim(gene_resultPos) # 24815 rows

# I'm going to save resultPos and use it for Mcav to see how it goes 
write.csv(resultPos, file="~/Desktop/GFFs/Mcav.gff.annotations.fixed_gene.diff.Pos.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE) 








