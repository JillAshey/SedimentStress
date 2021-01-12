# Title: Sediment stress - plob file compilation
# Author: Jill Ashey
# date: 11/30/20

library("tidyverse")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("gridExtra")
library("unpivotr")

# Compiling files by treatment comparison so all the info (GO terms, counts, pvalues, etc) is in one place 


## C vs Mid
# DEGs with info about pvalues, counts
DEG_control_vs_mid <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_control_vs_mid_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_mid)[1] <- "gene"

# GO terms with gene names 
GO_plob <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/plob_GOterms.csv", header = TRUE)
GO_plob <- select(GO_plob, -c(X, Predict))
GO_plob <- unique(GO_plob)
colnames(GO_plob)[1] <- "gene"
GO_plob$gene <- gsub(".m1", "", GO_plob$gene)
GO_plob$gene <- gsub("model", "TU", GO_plob$gene)

# Import reference annotation file 
ref <- read.csv("~/Desktop/GFFs/Plut.GFFannotation.fixed_transcript.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- subset(ref, id == "gene") 
ref$gene <- gsub(";.*", "", ref$attr)
ref$gene <- gsub("ID=", "", ref$gene)
ref <- select(ref, c(scaffold, gene.start, gene.stop, gene))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref) 

# Merge prot names from ref annotation file with GO_pdam
test <- merge(GO_plob, ref, by = "gene", all.x = TRUE)
test <- unique(test)

# Merge this ref/GO info with DEG file
merge <- merge(DEG_control_vs_mid, test, by = "gene", all.x = TRUE)
merge <- unique(merge)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_mid_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_mid_all <- select(DEG_control_vs_mid_all, -GO_term)
final_DEG_control_vs_mid_all <- unique(DEG_control_vs_mid_all)
write.csv(final_DEG_control_vs_mid_all, file = "~/Desktop/plob_control_vs_mid_DEG_GO_all.csv")





## C vs High
# DEGs with info about pvalues, counts
DEG_control_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_control_vs_high_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_high)[1] <- "gene"

# GO terms with gene names 
# Read in above

# read in reference annotation gff file
# Read in above

# Merge prot names from ref annotation file with GO_pdam
# Done above - merged into variable 'test'

# Merge this ref/GO info with DEG file
merge <- merge(DEG_control_vs_high, test, by = "gene", all.x = TRUE)
merge <- unique(merge)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_high_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_high_all <- select(DEG_control_vs_high_all, -GO_term)
final_DEG_control_vs_high_all <- unique(DEG_control_vs_high_all)
write.csv(final_DEG_control_vs_high_all, file = "~/Desktop/plob_control_vs_high_DEG_GO_all.csv")





## C vs High
# DEGs with info about pvalues, counts
DEG_mid_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_mid_vs_high_DEG_full.csv", header = TRUE)
colnames(DEG_mid_vs_high)[1] <- "gene"

# GO terms with gene names 
# Read in above

# read in reference annotation gff file
# Read in above

# Merge prot names from ref annotation file with GO_pdam
# Done above - merged into variable 'test'

# Merge this ref/GO info with DEG file
merge <- merge(DEG_mid_vs_high, test, by = "gene", all.x = TRUE)
merge <- unique(merge)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_mid_vs_high_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_mid_vs_high_all <- select(DEG_mid_vs_high_all, -GO_term)
final_DEG_mid_vs_high_all <- unique(DEG_mid_vs_high_all)
write.csv(final_DEG_mid_vs_high_all, file = "~/Desktop/plob_mid_vs_high_DEG_GO_all.csv")




