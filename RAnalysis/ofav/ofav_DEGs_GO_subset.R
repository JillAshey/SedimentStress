# Title: Ofav GO terms
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 09/27/20

# Code for Francois sedimentation data. Used interproscan to get GO terms for Ofav proteins. Now subsetting GO terms 
# and add those terms to annotation file

# Exploring how to merge ofav annotation gff and ofav interproscan gff

# Read in ofav annot file
ofav <- read.csv("~/Desktop/GFFs/GCF_002042975.1_ofav_dov_v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ofav) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ofav <- ofav[!grepl("##", ofav$scaffold),]
ofav_GO <- filter(ofav, grepl("XP_", attr)) # only want rows with proteins
ofav_GO$prot <-gsub(";.*", "", ofav_GO$attr)
ofav_GO$prot <-gsub(".*-", "", ofav_GO$prot)

# read in interproscan file 
ofav_IPS <- read.csv("~/Desktop/ofav.interpro.gff3",header = FALSE, sep="\t", skip=4)
length(unique(ofav_IPS$V1)) # 1245302
colnames(ofav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
ofav_IPS_GO <- filter(ofav_IPS, grepl("GO:", attr)) # select only rows with GO terms

# merge annot and interproscan file by protein
test <- merge(ofav_GO, ofav_IPS_GO, by = "prot")
test <- na.omit(test)

# subset bu Pfam predictor
pfam_test <- subset(test, Predict == "Pfam")
# Isolate the gene id
## need to manually extract gene
pfam_test$gene <- regmatches(pfam_test$attr.x, gregexpr("(?<=gene=).*", pfam_test$attr.x, perl = TRUE)) #removing everything up to LOC
pfam_test$gene <- gsub(";.*", "", pfam_test$gene) # removing everything after LOC term

# Use gene id to merge with DEGs file with full annot file -- will give final annotation of DEGs with GO terms, etc
colnames(DEGs.all) <- "gene"
full_annot <- merge(DEGs.all, pfam_test, by = "gene", all.x = TRUE)
write.csv(full_annot, file = "~/Desktop/full_annot_ofav.csv")

# Should I include the DEGsall file with the gene counts as well?

test_annot <- unique(full_annot$attr.y)


