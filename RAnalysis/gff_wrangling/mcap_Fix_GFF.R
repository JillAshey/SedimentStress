# Title: M. capitata GFF adjustments
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/19/20

# Need to do some mcap gff adjustments so it can run properly in STAR. Here, I'll be adding gene_id= and transcript_id=, and replacing 'id' with exon only 

#Load libraries
library(tidyverse)

#Load  gene gff
Mcap.gff <- read.csv(file="~/Desktop/Mcap.GFFannotation.gff.1", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Mcap.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")

# seeing what kinds of components are in annotation
unique(Mcap.gff$id)


#If id ==CDS add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse(Gene.Predict %in% c("AUGUSTUS") &
                         id == "CDS" ,
                       gsub("transcript_id ", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))
head(Mcap.gff)

#If id ==intron add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse(Gene.Predict %in% c("AUGUSTUS") &
                         id == "intron" ,
                       gsub("transcript_id ", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))
head(Mcap.gff)

#If id ==CDS add ;gene_id= <gene line ID > stopping at first ; , else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse(Gene.Predict %in% c("AUGUSTUS") &
                         id == "CDS" ,
                       gsub("gene_id ", "gene_id=", gene, fixed = TRUE), gsub("", "", gene)))
head(Mcap.gff)

#If id ==intron add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse(Gene.Predict %in% c("AUGUSTUS") &
                         id == "intron" ,
                       gsub("transcript_id ", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))
head(Mcap.gff)

#If id ==intron, remove last character in gene string, else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse( id == "intron" ,
                        substr(Mcap.gff$gene, 1, nchar(Mcap.gff$gene)-1), gsub("", "", gene)))
head(Mcap.gff)

#If id ==CDS, remove last character in gene string, else replace with nothing
Mcap.gff <- Mcap.gff %>%
  mutate(gene = ifelse( id == "CDS" ,
                        substr(Mcap.gff$gene, 1, nchar(Mcap.gff$gene)-1), gsub("", "", gene)))
head(Mcap.gff)

# Remove space in between transcript_id and gene_id
Mcap.gff$gene <- gsub(" ", "", Mcap.gff$gene)

#save file
write.table(Mcap.gff, file="~/Desktop/Mcap.GFFannotation.fixed_no_spaces.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)

# The original Mcap annotation file only had gene, CDS and intron in id column. In order to run properly, STAR needs at least some exons in the id column. Fix from online said to change all ids to exon
# But when I changed all ids to exon only in shell, the program also changed the gene_id in the gene column to exon_id and STAR needs gene_id to run (this is my assumption at the moment)
# So I am changing all ids in id column to exon in R, which will leave the gene column as is
unique(Mcap.gff$id)
Mcap.gff$id <- gsub("CDS", "exon", Mcap.gff$id)
Mcap.gff$id <- gsub("gene", "exon", Mcap.gff$id)
Mcap.gff$id <- gsub("intron", "exon", Mcap.gff$id)
unique(Mcap.gff$id)

#save file
write.table(Mcap.gff, file="~/Desktop/Mcap.GFFannotation.fixed_all_exons_no_spaces.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)






# 
# # new column for transcript id
# Mcap.gff$transcript_id <- sub(";.*", "", Mcap.gff$gene)
# head(Mcap.gff)
# Mcap.gff$transcript_id <- gsub("transcript_id", "", Mcap.gff$transcript_id)

# # Trying to get = in between space between id and string
# # Mcap.gff$gene_temp_blah <- sub("transcript_id.", "", Mcap.gff$gene)
# # Mcap.gff$gene_temp <- sub("exon_id", "", Mcap.gff$gene_temp_blah)
# # Mcap.gff$gene_temp_blah <- sub("*.;", "", Mcap.gff$gene_temp_blah)
# Mcap.gff$transcript_id <- gsub("transcript_id", "", Mcap.gff$transcript_id) #remove transcript_id
# head(str(Mcap.gff$gene, "^id"))
# #If id ==exon add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
# Mcap.gff <- Mcap.gff %>% 
#   mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Mcap.gff$transcript_id),  paste0(gene)))
