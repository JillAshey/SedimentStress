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
DEG_control_vs_mid <- read.csv("~/Desktop/pdam_control_vs_mid_DEG_full.csv", header = TRUE)
colnames(DEG_control_vs_mid)[1] <- "gene"

# GO terms with gene names 
GO_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/pdam_GOterms.csv", header = TRUE)
GO_pdam <- select(GO_pdam, -c(X, Predict))
GO_pdam <- unique(GO_pdam)
colnames(GO_pdam)[1] <- "gene"
#GO_mcav$gene <- gsub("-RA", "", GO_mcav$gene)

# Merge this information
merge <- merge(DEG_control_vs_mid, GO_pdam, by = "gene", all.x = TRUE)
merge <- unique(merge)
#merge <- select(merge, -GO_term)

