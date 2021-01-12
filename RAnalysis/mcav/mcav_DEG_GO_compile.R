# Title: Sediment stress - mcav file compilation
# Author: Jill Ashey
# date: 11/18/20

library("tidyverse")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("gridExtra")
library("unpivotr")


# Compiling files by treatment comparison so all the info (GO terms, counts, pvalues, etc) is in one place 



## C vs T1
# DEGs with info about pvalues, counts
DEG_control_vs_T1 <- read.csv("~/Desktop/mcav_control_vs_T1_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T1)[1] <- "gene"

# GO terms with gene names 
GO_mcav <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/mcav_GOterms.csv", header = TRUE)
GO_mcav <- select(GO_mcav, -c(X, Predict))
GO_mcav <- unique(GO_mcav)
colnames(GO_mcav)[1] <- "gene"
GO_mcav$gene <- gsub("-RA", "", GO_mcav$gene)

# Merge this information
merge <- merge(DEG_control_vs_T1, GO_mcav, by = "gene", all.x = TRUE)
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
ref <- read.table("~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3", sep = "\t", header = FALSE)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- subset(ref, id == "gene") 
ref$gene <- gsub(";.*", "", ref$attr)
ref$gene <- gsub("ID=", "", ref$gene)
ref <- select(ref, c(scaffold, gene.start, gene.stop, gene))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref) # 25142 x 5

# Merge ref and DEG_control_vs_T1_all 
final_DEG_CvsT1 <- merge(DEG_control_vs_T1_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT1, file = "~/Desktop/mcav_control_vs_T1_DEG_GO_all.csv")





## C vs T2
# DEGs with info about pvalues, counts
DEG_control_vs_T2 <- read.csv("~/Desktop/mcav_control_vs_T2_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T2)[1] <- "gene"

# GO terms with gene names 
# Use GO_mcav that was read in above 
GO_mcav

# Merge this information
merge <- merge(DEG_control_vs_T2, GO_mcav, by = "gene", all.x = TRUE)
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
# Use ref that was read in above 
ref

# Merge ref and DEG_control_vs_T2_all 
final_DEG_CvsT2 <- merge(DEG_control_vs_T2_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT2, file = "~/Desktop/mcav_control_vs_T2_DEG_GO_all.csv")





## C vs T3
# DEGs with info about pvalues, counts
DEG_control_vs_T3 <- read.csv("~/Desktop/mcav_control_vs_T3_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T3)[1] <- "gene"

# GO terms with gene names 
# Use GO_mcav that was read in above 
GO_mcav

# Merge this information
merge <- merge(DEG_control_vs_T3, GO_mcav, by = "gene", all.x = TRUE)
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
# Use ref that was read in above 
ref

# Merge ref and DEG_control_vs_T2_all 
final_DEG_CvsT3 <- merge(DEG_control_vs_T3_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT3, file = "~/Desktop/mcav_control_vs_T3_DEG_GO_all.csv")





## C vs T4
# DEGs with info about pvalues, counts
DEG_control_vs_T4 <- read.csv("~/Desktop/mcav_control_vs_T4_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_T4)[1] <- "gene"

# GO terms with gene names 
# Use GO_mcav that was read in above 
GO_mcav

# Merge this information
merge <- merge(DEG_control_vs_T4, GO_mcav, by = "gene", all.x = TRUE)
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
# Use ref that was read in above 
ref

# Merge ref and DEG_control_vs_T3_all 
final_DEG_CvsT4 <- merge(DEG_control_vs_T4_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_CvsT4, file = "~/Desktop/mcav_control_vs_T4_DEG_GO_all.csv")




## T1 vs T3
# DEGs with info about pvalues, counts
DEG_T1_vs_T3 <- read.csv("~/Desktop/mcav_T1_vs_T3_DEG_full.csv", header = TRUE)
colnames(DEG_T1_vs_T3)[1] <- "gene"

# GO terms with gene names 
# Use GO_mcav that was read in above 
GO_mcav

# Merge this information
merge <- merge(DEG_T1_vs_T3, GO_mcav, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_T1_vs_T3_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_T1_vs_T3_all <- select(DEG_T1_vs_T3_all, -GO_term)
DEG_T1_vs_T3_all <- unique(DEG_T1_vs_T3_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use ref that was read in above 
ref

# Merge ref and DEG_control_vs_T2_all 
final_DEG_T1vsT3 <- merge(DEG_T1_vs_T3_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_T1vsT3, file = "~/Desktop/mcav_T1_vs_T3_DEG_GO_all.csv")





## T2 vs T3
# DEGs with info about pvalues, counts
DEG_T2_vs_T3 <- read.csv("~/Desktop/mcav_T2_vs_T3_DEG_full.csv", header = TRUE)
colnames(DEG_T2_vs_T3)[1] <- "gene"

# GO terms with gene names 
# Use GO_mcav that was read in above 
GO_mcav

# Merge this information
merge <- merge(DEG_T2_vs_T3, GO_mcav, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_T2_vs_T3_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_T2_vs_T3_all <- select(DEG_T2_vs_T3_all, -GO_term)
DEG_T2_vs_T3_all <- unique(DEG_T2_vs_T3_all)
#write.csv(DEG_control_vs_T1_all, file = "~/Desktop/acerv_control_vs_T1_DEG_GO_all.csv")

# Read in ref file 
# Use ref that was read in above 
ref

# Merge ref and DEG_control_vs_T2_all 
final_DEG_T2vsT3 <- merge(DEG_T2_vs_T3_all, ref, by = "gene", all.x = TRUE)
write.csv(final_DEG_T2vsT3, file = "~/Desktop/mcav_T2_vs_T3_DEG_GO_all.csv")


