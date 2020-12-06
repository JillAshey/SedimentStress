# Title: Acerv functional annotation
# Author: Jill Ashey
# date: 12/05/20



library("tidyverse")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("gridExtra")
library("unpivotr")

# Compiling files by treatment comparison so all the info (GO terms, counts, pvalues, etc) is in one place 



## C vs T1
# DEGs with info about pvalues, counts
DEG_control_vs_T1 <- read.csv("~/Desktop/acerv_control_vs_T1_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T1)[1] <- "gene"

# GO terms with gene names 
GO_acerv <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv_GOterms.csv", header = TRUE)
GO_acerv <- select(GO_acerv, -c(X, Predict))
GO_acerv <- unique(GO_acerv)
colnames(GO_acerv)[1] <- "gene"
GO_acerv$gene <- gsub("model", "TU", GO_acerv$gene)

# Merge this information
merge <- merge(DEG_control_vs_T1, GO_acerv, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_T1_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_T1_all <- select(DEG_control_vs_T1_all, -GO_term)
DEG_control_vs_T1_all <- unique(DEG_control_vs_T1_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
ref <- read.table("~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", sep = "\t", header = FALSE)
ref <- subset(ref, V3 == "gene")
colnames(ref) <- c("scaffold", "gene.predict", "id", "start", "stop", "pos1", "pos2", "pos3", "attr", "gene")
ref <- select(ref, c(scaffold, start, stop, gene))
#ref$transcript_id <- gsub(".*;", "", ref$transcript_id)
#ref$transcript_id <- gsub("Name=", "", ref$transcript_id)
ref <- ref %>% mutate(ref, length = stop - start)

# Merge ref and DEG_control_vs_T1_all 
final_DEG_CvsT1 <- merge(DEG_control_vs_T1_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT1, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")




## C vs T2
# DEGs with info about pvalues, counts
DEG_control_vs_T2 <- read.csv("~/Desktop/acerv_control_vs_T2_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T2)[1] <- "gene"

# Read in GO terms with gene names 
# Use GO_acerv that was read in above 
GO_acerv

# Merge this information
merge <- merge(DEG_control_vs_T2, GO_acerv, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_T2_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_T2_all <- select(DEG_control_vs_T2_all, -GO_term)
DEG_control_vs_T2_all <- unique(DEG_control_vs_T2_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use GO_acerv that was read in above 
ref

# Merge ref and DEG_control_vs_T2_all 
final_DEG_CvsT2 <- merge(DEG_control_vs_T2_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT2, file = "~/Desktop/acerv_control_vs_T2_DEG_GO_all.csv")




## C vs T3
# DEGs with info about pvalues, counts
DEG_control_vs_T3 <- read.csv("~/Desktop/acerv_control_vs_T3_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T3)[1] <- "gene"

# Read in GO terms with gene names 
# Use GO_acerv that was read in above 
GO_acerv

# Merge this information
merge <- merge(DEG_control_vs_T3, GO_acerv, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_T3_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_T3_all <- select(DEG_control_vs_T3_all, -GO_term)
DEG_control_vs_T3_all <- unique(DEG_control_vs_T3_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use GO_acerv that was read in above 
ref

# Merge ref and DEG_control_vs_T3_all 
final_DEG_CvsT3 <- merge(DEG_control_vs_T3_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT3, file = "~/Desktop/acerv_control_vs_T3_DEG_GO_all.csv")




## C vs T4
# DEGs with info about pvalues, counts
DEG_control_vs_T4 <- read.csv("~/Desktop/acerv_control_vs_T4_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T4)[1] <- "gene"

# Read in GO terms with gene names 
# Use GO_acerv that was read in above 
GO_acerv

# Merge this information
merge <- merge(DEG_control_vs_T4, GO_acerv, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_T4_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_T4_all <- select(DEG_control_vs_T4_all, -GO_term)
DEG_control_vs_T4_all <- unique(DEG_control_vs_T4_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use GO_acerv that was read in above 
ref

# Merge ref and DEG_control_vs_T4_all 
final_DEG_CvsT4 <- merge(DEG_control_vs_T4_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT4, file = "~/Desktop/acerv_control_vs_T4_DEG_GO_all.csv")




## T1 vs T4
# DEGs with info about pvalues, counts
DEG_T1_vs_T4 <- read.csv("~/Desktop/acerv_T1_vs_T4_DEG_full.csv", header = TRUE)
colnames(DEG_T1_vs_T4)[1] <- "gene"

# Read in GO terms with gene names 
# Use GO_acerv that was read in above 
GO_acerv

# Merge this information
merge <- merge(DEG_T1_vs_T4, GO_acerv, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_T1_vs_T4_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_T1_vs_T4_all <- select(DEG_T1_vs_T4_all, -GO_term)
DEG_T1_vs_T4_all <- unique(DEG_T1_vs_T4_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use GO_acerv that was read in above 
ref

# Merge ref and DEG_control_vs_T4_all 
final_DEG_T1vsT4 <- merge(DEG_T1_vs_T4_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_T1vsT4, file = "~/Desktop/acerv_T1_vs_T4_DEG_GO_all.csv")

























# # Merge GO term info with GO.terms
# go_full <- merge(GO.terms, enriched_GO, by = "category", all.x=TRUE)
# colnames(go_full)[2] <- "gene"
# 
# # Merge all go term info with deg
# all_plz <- merge(go_full, deg, by = "gene", all.x = TRUE)
# all_plz <- select(all_plz, -c(X))
# 
# # Trying to get it down to only the 30 DEGs
# # Going to try to aggregate by term
# agg2 <- aggregate(all_plz$term, list(all_plz$gene), paste, collapse = ",")
# colnames(agg2) <- c("gene", "term")
# 
# # Merge the all_plz df and the agg2 df
# all_plz2 <- merge(agg2, all_plz, by = "gene", all.x = TRUE)
# all_plz2 <- select(all_plz2, -c(term.y, category))
# all_plz2 <- unique(all_plz2)
# 
# # split GO terms in deg file
# split_GO <- strsplit(as.character(deg$GO_terms), ",")
# GO.terms <- data.frame(v1 = rep.int(deg$gene, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
# colnames(GO.terms) <- c("gene_id", "category")




