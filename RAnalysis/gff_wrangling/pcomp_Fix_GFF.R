#Load libraries
library(tidyverse)

#Load  gene gff
Pcomp.gff <- read.csv(file="~/Desktop/Porites_compressa_genemodels_v1.0_201905.gff", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Pcomp.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")

Pcomp.gff$transcript_id <- sub(";.*", "", Pcomp.gff$gene)
Pcomp.gff$transcript_id <- gsub("ID=", "", Pcomp.gff$transcript_id) #remove ID= 

Pcomp.gff$transcript_id <-sub("^([^.]*.[^.]*.[^.]*).*", "\\1", Pcomp.gff$transcript_id) #remove everything after the third . in the gene column
head(Pcomp.gff)  

#If id == mRNA, exon, five_prime_UTR, three_prime_UTR, CDS, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Pcomp.gff <- Pcomp.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Pcomp.gff$transcript_id),  paste0(gene)))

head(Pcomp.gff)

Pcomp.gff <- Pcomp.gff[,-10]
head(Pcomp.gff)     

#save file
write.table(Pcomp.gff, file="~/Desktop/Pcomp.GFFannotation.fixed_transcript.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
