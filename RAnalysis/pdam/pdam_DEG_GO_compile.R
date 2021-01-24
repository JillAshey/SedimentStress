# Title: Sediment stress - pdam file compilation
# Author: Jill Ashey
# date: 11/19/20

library("tidyverse")
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("gridExtra")
library("unpivotr")


# Compiling files by treatment comparison so all the info (GO terms, counts, pvalues, etc) is in one place 


## C vs Mid
# DEGs with info about pvalues, counts
DEG_control_vs_mid <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_control_vs_mid_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_mid)[1] <- "gene"

# GO terms with gene names 
GO_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/pdam_GOterms.csv", header = TRUE)
GO_pdam <- select(GO_pdam, -c(X, Predict))
GO_pdam <- unique(GO_pdam)
colnames(GO_pdam)[1] <- "prot"
# need to link the XP protein name with the LOC gene name

# read in reference annotation gff file
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
ref <- ref[grep("XP", ref$attr), ] # isolate XP protein name
ref$gene <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$gene <- gsub(";.*", "", ref$gene) # remove everything after ;
ref$prot <- gsub(";.*", "", ref$attr)
ref$prot <- gsub("ID=", "", ref$prot)
ref$prot <- gsub(".*-", "", ref$prot)

# Merge prot names from ref annotation file with GO_pdam
test <- merge(GO_pdam, ref, by = "prot", all.x = TRUE)
test <- unique(test)

# Merge this ref/GO info with DEG file
merge <- merge(DEG_control_vs_mid, test, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

# Because IPS used different databases to generate GO terms, there are multiple rows of the same gene name with different 
# GO terms (because they were generated from different db). 
agg <- aggregate(merge$GO_term, list(merge$gene), paste, collapse = ",")
colnames(agg) <- c("gene", "GO_terms")

# Merge the merge df and the agg df
DEG_control_vs_mid_all <- merge(merge, agg, by = "gene", all.x = TRUE)
DEG_control_vs_mid_all <- select(DEG_control_vs_mid_all, -GO_term)
DEG_control_vs_mid_all <- unique(DEG_control_vs_mid_all)
write.csv(DEG_control_vs_mid_all, file = "~/Desktop/pdam_control_vs_mid_DEG_GO_all.csv")

# Calculate length for genes
final_DEG_control_vs_mid_all <- DEG_control_vs_mid_all %>% mutate(DEG_control_vs_mid_all, length = gene.stop - gene.start)
write.csv(final_DEG_control_vs_mid_all, file = "~/Desktop/pdam_control_vs_mid_DEG_GO_all.csv")





## C vs High
# DEGs with info about pvalues, counts
DEG_control_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_control_vs_high_DEG_full.csv", header = TRUE)
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
DEG_control_vs_high_all <- unique(DEG_control_vs_high_all)
write.csv(DEG_control_vs_high_all, file = "~/Desktop/pdam_control_vs_high_DEG_GO_all.csv")

# Calculate length for genes
final_DEG_control_vs_high_all <- DEG_control_vs_high_all %>% mutate(DEG_control_vs_high_all, length = gene.stop - gene.start)
write.csv(final_DEG_control_vs_high_all, file = "~/Desktop/pdam_control_vs_high_DEG_GO_all.csv")





## C vs High
# DEGs with info about pvalues, counts
DEG_mid_vs_high <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_mid_vs_high_DEG_full.csv", header = TRUE)
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
DEG_mid_vs_high_all <- unique(DEG_mid_vs_high_all)
write.csv(DEG_mid_vs_high_all, file = "~/Desktop/pdam_mid_vs_high_DEG_GO_all.csv")

# Calculate length for genes
final_DEG_mid_vs_high_all <- DEG_mid_vs_high_all %>% mutate(DEG_mid_vs_high_all, length = gene.stop - gene.start)
write.csv(final_DEG_mid_vs_high_all, file = "~/Desktop/pdam_mid_vs_high_DEG_GO_all.csv")


