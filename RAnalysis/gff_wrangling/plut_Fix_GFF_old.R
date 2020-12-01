#Load libraries
library(tidyverse)

#Load  gene gff
Plut.gff <- read.csv(file="~/Desktop/plut2v1.1.genes.gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Plut.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Plut.gff)

Plut.gff$transcript_id <- sub(";.*", "", Plut.gff$gene)
Plut.gff$transcript_id <- gsub("ID=", "", Plut.gff$transcript_id) #remove ID= 
head(Plut.gff)

Plut.gff$transcript_id <-sub("^([^.]*.[^.]*.[^.]*).*", "\\1", Plut.gff$transcript_id) #remove everything after the third . in the gene column
head(Plut.gff)  

#If id == mRNA, exon, five_prime_UTR, three_prime_UTR, CDS, add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with original gene
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Plut.gff$transcript_id),  paste0(gene)))

head(Plut.gff)

Plut.gff <- Plut.gff[,-10]
head(Plut.gff)  









Plut.gff$transcript_id <- sub(";.*", "", Plut.gff$gene)
Plut.gff$transcript_id <- gsub("ID=", "", Plut.gff$transcript_id)
head(Plut.gff)

Plut.gff$transcript_id <- sub("(.[^]+).*", "\\1", Plut.gff$transcript_id)
head(Plut.gff)




#If id ==mRNA add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "mRNA" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))


#If id ==exon add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "exon" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))

#If id ==five_prime_UTR add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "five_prime_UTR" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))

#If id ==three_prime_UTR add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "three_prime_UTR" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))

#If id ==CDS add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "CDS" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))

#If id ==gene add ;transcript_id= <gene line ID without ID= stopping at first ; , else replace with nothing
Plut.gff <- Plut.gff %>% 
  mutate(gene = ifelse(Gene.Predict %in% c("maker") & 
                         id == "gene" ,  
                       gsub("ID=", "transcript_id=", gene, fixed = TRUE), gsub("", "", gene)))

#save file
write.table(Plut.gff, file="~/Desktop/Plut.GFFannotation.fixed.gff", sep="\t", 
            col.names = FALSE, row.names=FALSE, quote=FALSE)
