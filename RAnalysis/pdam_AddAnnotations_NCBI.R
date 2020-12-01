# Title: P. damicornis NCBI
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/27/20

# Using NCBI annotation file and NCBI gene results to better annotate the NCBI gff. It appears to have correct number of genes and what not, but no GO terms, which stops me from doing GOSeq
# Connelly github includes pdam and GO terms, but does not say where they came from. I'm going to use them to see if I can figure out how the NCBI gff file is working

# Goal: Use data from NCBI and Connelly to check and assign GO terms to NCBI gff

# Links where I downloaded data
# https://www.ncbi.nlm.nih.gov/gene?LinkName=bioproject_gene&from_uid=454489 -- NCBI pdam gene data. Don't seem to have GO terms, but has pdam_xx terms
# https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/data/pdam_genome_genesGO.txt -- Connelly github list of pdam_xx terms with associated GO terms
# https://github.com/michaeltconnelly/EAPSI_Pocillopora_LPS/blob/master/data/pdam_gene_universe.txt -- Connelly github list of pdam_xx terms

# Load packages
library("tidyverse")
library("dplyr")
library("genefilter")
library("ggplot2")
library("spdep") 
library("clusterProfiler")
library("DataCombine")

## NCBI gff comparing to NCBI Pdam gene data
pdamgff3_NCBI <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff", header=FALSE, sep="\t", skip=6) # read in gff from NCBI
colnames(pdamgff3_NCBI) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene") # name cols
dim(pdamgff3_NCBI) # 516,693 x 9
# pdamgff3_NCBI <- na.omit(pdamgff3_NCBI) # remove NAs
# dim(pdamgff3_NCBI) # 507906 x 9

# Isolate LOC term from gene col into its own col 
pdamgff3_NCBI$Symbol <- pdamgff3_NCBI$gene # copying gene column into new column called Symbol
pdamgff3_NCBI$Symbol <- regmatches(pdamgff3_NCBI$Symbol, gregexpr("(?<=gene=).*", pdamgff3_NCBI$Symbol, perl = TRUE)) #removing everything in Symbol col up to LOC
pdamgff3_NCBI$Symbol <- gsub(";.*", "", pdamgff3_NCBI$Symbol) # removing everything after LOC term
length(unique(pdamgff3_NCBI$Symbol)) # 22800 unique LOC terms in gff file

# I have isolated the LOC terms in the NCBI gff file into Symbol colymn. Now I'm going to compare those LOC terms with the LOC terms in the NCBI pdam gene data
geneResult <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/gene_result_pdam_NCBI.csv", header=TRUE) # read in NCBI pdam gene data
geneResult <- select(geneResult, c(Org_name, GeneID, Symbol, Aliases, description)) # getting rid of some extra cols that im not interested in rn
dim(geneResult) # 22968 x 5
length(unique(geneResult$Symbol)) # 22702 unique LOC terms in geneResult file

# Now that I have isolated the LOC terms, I will compare LOC terms in pdamgff3_NCBI and geneResult to see if there are any matches.
LOC_compare <- pdamgff3_NCBI[pdamgff3_NCBI$Symbol %in% geneResult$Symbol,] # selecting LOC terms that match in both lists
# These are also viable options for comparing the lists by Symbol: 
# test <-subset(pdamgff3_NCBI, Symbol %in% geneResult$Symbol)
# library("data.frame")
# setDT(pdamgff3_NCBI)[Symbol %chin% geneResult$Symbol]
dim(LOC_compare) # 483,931 x 11
length(unique(LOC_compare$Symbol)) # 22702 unique LOC terms
# I have df with LOC terms that match in both pdamgff3_NCBI and geneResult. Not all pdamgff3_NCBI LOC terms were in geneResult, so they were dropped

# Now I have a list of unique LOC terms in geneResult and LOC_compare. There are pdam_xx terms associated with the unique LOC terms in the geneResult file. 
geneResult_pdam <- geneResult[str_detect(geneResult$Aliases, "pdam"), ]  # Extract rows with pdam terms with str_detect
dim(geneResult_pdam) # 15689 x 5 -- so 15689 LOC terms associated with a pdam term
length(unique(geneResult_pdam$Aliases)) # checking to make sure all pdam terms were unique 
length(unique(geneResult_pdam$Symbol)) # checking to make sure all LOC terms were unique 

# Compare LOC terms in pdamgff3_NCBI and geneResult_pdam by LOC to find similarities
pdam_LOC_compare <- pdamgff3_NCBI[pdamgff3_NCBI$Symbol %in% geneResult_pdam$Symbol,]
dim(pdam_LOC_compare) # 400523 x 11
length(unique(pdam_LOC_compare$Symbol)) # 15689 unique LOC terms
# I guess comparing pdamgff3_NCBI and geneResult above was redundant. Well, I now have df with pdam-term filtered LOC terms that match both pdamgff3_NCBI and geneResult_pdam


# Okay so I have df with the LOC terms in it that match both the NCBI gff file and the gene results file. Time to deal with pdam terms 

# I got both pdam_universal and pdam_gene_GOterms lists from Connelly github. now sure how he generated them, but he were using pdam terms 
pdam_universal <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/universal_pdam_Connelly.csv", header = FALSE) # list of all pdam terms that Connelly used in his analysis (?)
dim(pdam_universal) # 12702 x 1
colnames(pdam_universal) <- "Aliases"
pdam_gene_GOterms <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Data/gene_GOterms_pdam_Connelly.csv", header = FALSE) # List of pdam terms that Connelly used in analysis and their associated GO terms. Again, not sure how he got this data
dim(pdam_gene_GOterms) # 12702 x 15
colnames(pdam_gene_GOterms) <- c("Aliases", "GO1", "GO2", "GO3", "GO4", "GO5", "GO6", "GO7", "GO8", "GO9", "GO10", "GO11", "GO12", "GO13", "GO14")
list_compare <- pdam_universal[pdam_universal$Aliases %in% pdam_gene_GOterms$Aliases,] # check to see if they matched one another, they do

# Compare pdam_xx terms in geneResult with those in pdam_gene_GOterms to find similarities
pdam_term_compare <- geneResult_pdam[geneResult_pdam$Aliases %in% pdam_gene_GOterms$Aliases,]
dim(pdam_term_compare) # 10050 x 5
length(unique(pdam_term_compare$Symbol)) # 10050 unique LOC terms
length(unique(pdam_term_compare$Aliases)) # 10050 unique pdam terms
# Now I have a list of the unique LOC terms associated with pdam terms that also have GO terms!

# I have to add the LOC-associated pdam GOterms into pdamgff3_NCBI file...
# First, I am merging the NCBI gene results with pdamGO terms by Aliases
newtable <- merge(geneResult, pdam_gene_GOterms, by = "Aliases", all.x=TRUE)
# Then I can merge the pdamgff3_NCBI with newtable by Symbol
finaltable <- merge(pdamgff3_NCBI, newtable, by = "Symbol", all.x=TRUE) # yay!! final table! ugly, but info is there
write.table(finaltable, file="~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms_raw.gff", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)

# Cleaning up table so it fits gff format
finaltable <- subset(finaltable, select = -Symbol) # removing symbol column
finaltable <- finaltable %>% 
  filter(!str_detect(scaffold, '##')) # removing empty rows
finaltable <- finaltable %>% 
  unite(finaltable, Aliases:description, sep = ";", remove = TRUE, na.rm = TRUE) # unite columns with ; separator 
finaltable <- finaltable %>% mutate_all(na_if,"")
finaltable <- finaltable %>% 
  unite(finaltable, gene:finaltable, sep = ";", remove = TRUE, na.rm = TRUE) # unite cols with ; separator 
names(finaltable)[names(finaltable) == "finaltable"] <- "gene"
finaltable <- finaltable %>% 
  unite(finaltable, GO1:GO14, sep = "", remove = TRUE, na.rm = TRUE)
names(finaltable)[names(finaltable) == "finaltable"] <- "GOterms"
write.table(finaltable, file="~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms_sepcol.gff", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
write.table(finaltable, file="~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms_sepcol.csv")
finaltable <- finaltable %>% 
  unite(finaltable, gene:finaltable, sep = ";", remove = TRUE, na.rm = TRUE)

write.table(finaltable, file="~/Desktop/GFFs/pdam_NCBI_annotation_fixed_GOterms.gff", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)
write.table(finaltable, file="~/Desktop/GFFs/pdam_NCBI_annotation_fixed.txt", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)

dim(finaltable)


write.table(finaltable, file="~/Desktop/GFFs/test.txt", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)



# The gff file from Reef Genomics 
pdamgff3_rg <- read.csv("Desktop/GFFs/pdam_annotation.gff3", header=FALSE, sep="\t", skip=2)
colnames(pdamgff3_rg) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(pdamgff3_rg) # 835,524 x 9
pdamgff3_rg$Aliases <- sub(";.*", "", pdamgff3_rg$gene)
pdamgff3_rg$Aliases <- gsub("-.*", "", pdamgff3_rg$Aliases) # removing everything after '-'

compare_gffs <- finaltable[finaltable$Aliases %in% pdamgff3_rg$Aliases,]
unique(compare_gffs$Aliases) # so no actual matches

