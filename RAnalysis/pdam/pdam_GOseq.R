# Title: Sediment stress - pdam GOseq
# Author: Jill Ashey
# date: 9/30/20

### Gene Ontology Analysis of Differentially Expressed Genes in pdam - NCBI gff


# Load packages 
library("DESeq2")
library("tidyverse")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("ggplot2")
library("gplots")
library("limma")
library("spdep") 
library("adegenet") 
library("goseq")
library("forcats")
library("gridExtra")
library("GO.db")

## The gene ids in gcounts and DEG files are LOC terms. I need to use the ref file to assign the LOC terms to XM terms. Then I can use the length_pdam file!
# Import reference annotation gff file. This is from NCBI
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6) # ref file from NCBI
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$gene_id <- gsub(";.*", "", ref$gene_id) # remove everything after ;
ref <- ref[!grepl("character", ref$gene_id),] # remove rows that have 'character' in gene_id col
ref <- subset(ref, id == "mRNA") # select only genes -- maybe i do need to do this step??
ref$XM <- gsub(";.*", "", ref$attr)
ref$XM <- gsub("ID=", "", ref$XM)
ref$XM <- gsub(".*-", "", ref$XM)

# Obtain names of all expressed ofav genes (poverA = 0.85,5), and all differentially expressed planuala genes (p<0.05)
gcounts_filt_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam/pdam_counts_filt.csv", header = TRUE) # read data in
dim(gcounts_filt_pdam) # 13881 rows x 13
#for ( col in 1:ncol(gcounts_filt_pdam)){
# colnames(gcounts_filt_pdam)[col] <-  gsub("X", "", colnames(gcounts_filt_pdam)[col]) # remove X in front of col names
#}
colnames(gcounts_filt_pdam)[1] <-"gene_id" # make colnames a true column called gene_id
head(gcounts_filt_pdam)
length(unique(gcounts_filt_pdam$gene_id)) # 13881 total unique gene ids
# Merge ref and gcounts_filt_pdam to assign LOC terms in gcounts_filt_pdam to XM terms
gcounts_filt_pdam_XM <- merge(gcounts_filt_pdam, ref, by = "gene_id")
gcounts_filt_pdam_XM <- unique(gcounts_filt_pdam_XM)

dim(gcounts_filt_pdam_XM) # 17119

# Import file with DEGs for pdam
DEG_pdam <- read.csv("~/Desktop/pdam_unique.sig.list_20210326.csv", header = TRUE) # read in list of significant pdam genes
dim(DEG_pdam) # 549 x 13
#for ( col in 1:ncol(DEG_pdam)){
#   colnames(DEG_pdam)[col] <-  gsub("X", "", colnames(DEG_pdam)[col]) # remove X in front of col names
# }
colnames(DEG_pdam)[1] <-"gene_id" # make colnames a true column called gene_id
head(DEG_pdam)
length(unique(DEG_pdam$gene_id)) # 549 total unique gene ids 
# Merge ref and gcounts_filt_pdam to assign LOC terms in gcounts_filt_pdam to XM terms
DEG_pdam_XM <- merge(DEG_pdam, ref, by = "gene_id")
DEG_pdam_XM <- unique(DEG_pdam_XM)

# Read in length data (calculated directly from transcripts) and merge with filtered counts 
length_Pdam <- read.csv("~/Desktop/length.mRNA_Pdam.csv", header = F)
length_Pdam <- length_Pdam[grepl("XM", length_Pdam$V1),]
colnames(length_Pdam) <- c("gene_id", "length")

#length_merge <- merge(length_Plob, gcounts_filt_plob, by = "gene_id")











# Import reference annotation gff file. This is from NCBI
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6) # ref file from NCBI
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
ref$LOC <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$LOC <- gsub(";.*", "", ref$LOC) # remove everything after ;
ref <- ref[!grepl("character", ref$LOC),] # remove rows that have 'character' in gene_id col
ref <- subset(ref, id == "mRNA") # select only genes -- maybe i do need to do this step??
ref$gene_id <- gsub(";.*", "", ref$attr)
ref$gene_id <- gsub("ID=", "", ref$gene_id)
ref$gene_id <- gsub(".*-", "", ref$gene_id)

test <- merge(length_Pdam, ref, by = "gene_id")





#ref <- subset(ref, id == "gene") # select only genes -- skipping this for now because I am having some issues with using nullp and lengths of vectors 
#ref <- ref[grep("XP", ref$attr), ] # isolate XP protein name
#ref$prot <- gsub(";.*", "", ref$attr)
#ref$prot <- gsub("ID=", "", ref$prot)
#ref$prot <- gsub(".*-", "", ref$prot)
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$gene_id <- gsub(";.*", "", ref$gene_id) # remove everything after ;
ref <- ref[!grepl("character", ref$gene_id),] # remove rows that have 'character' in gene_id col
ref <- select(ref, c(scaffold, gene.start, gene.stop, id, gene_id)) # select only certain cols
ref <- ref %>% mutate(ref, length = gene.stop - gene.start) # calculate gene length
dim(ref) # 484440 x 6
length(unique(ref$gene_id)) # 22799 
ref <- unique(ref) 
dim(ref) # 388598 x 6
ref <- subset(ref, id == "gene") # select only genes -- maybe i do need to do this step??
dim(ref) # 21837 x 6
length(unique(ref$gene_id)) # 21562 
# in this ref file, there are 21562 unique gene ids 
# Strange that the map file (created by me with stringTie) has more unique gene ids that the ref file (from NCBI). It could be because star/stringTie created 'novel' genes based on 
# exons splicing in different ways.
# I only need gene id and length...going to keep just those and throw out scaffold and gene start/stop
ref <- select(ref, c(length, gene_id))



#### Build GOSEQ vector 
#GOseq requires a vector of all genes, all differentially expressed genes, and gene lengths
# Make DEG/gene vector 
DEG <- filter(ref, gene_id%in%DEG_pdam$gene_id) #make vector of differentially expressed genes
dim(DEG) # 593 rows 
DEG_names <- as.vector(DEG$gene_id)
gene_vector=as.integer(ref$gene_id%in%DEG_names)
names(gene_vector)=ref$gene_id

# Make ID vector 
IDvector <- ref$gene_id

# Make length vector
lengthVector <- ref$length

# Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, IDvector, bias.data=lengthVector) #weight vector by length of gene
# plot looks weird...does this mean that there is not much length bias in the dataset?








































### I dont think i actually need to use the merged file at all...all the info I need to merge with gcounts_filt is in ref (gene id and length)
#Import merged annotated gtf file. This was created with stringTie and merged all info about genes expressed in pdam for all samples
map <- read.csv("~/Desktop/GFFs/Pdam.merged.annotated.gtf", header=FALSE, sep="\t") # merged annotated pdam file created with stringTie
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name columns
map <- map[!grepl("#", map$scaffold),] # remove rows that have a # in scaffold col
map <- map[grep("LOC", map$attr), ] # select all rows that have LOC in attr col
map$gene_id <- regmatches(map$attr, gregexpr("(?<=gene_name).*", map$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and making new col called gene_id
map$gene_id <- gsub(";.*", "", map$gene_id) # remove everything after ;
map$gene_id <- gsub(" ", "", map$gene_id) # remove any blank spaces
map <- subset(map, id=="transcript") # select only transcripts
map <- select(map, c(gene_id)) # keep only gene id col
map_unique <- unique(map) # keep unique rows only
dim(map_unique) # 26139 x 4
length(unique(map_unique$gene_id)) # 22734 unique geneids
# # in this map file, there are 22734 unique gene ids. 

# Import reference annotation gff file. This is from NCBI
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6) # ref file from NCBI
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##sequence", ref$scaffold),] # remove rows that have a # in scaffold col
#ref <- subset(ref, id == "gene") # select only genes -- skipping this for now because I am having some issues with using nullp and lengths of vectors 
#ref <- ref[grep("XP", ref$attr), ] # isolate XP protein name
#ref$prot <- gsub(";.*", "", ref$attr)
#ref$prot <- gsub("ID=", "", ref$prot)
#ref$prot <- gsub(".*-", "", ref$prot)
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$gene_id <- gsub(";.*", "", ref$gene_id) # remove everything after ;
ref <- ref[!grepl("character", ref$gene_id),] # remove rows that have 'character' in gene_id col
ref <- dpylr::select(ref, c(scaffold, gene.start, gene.stop, id, gene_id)) # select only certain cols
ref <- ref %>% mutate(ref, length = gene.stop - gene.start) # calculate gene length
dim(ref) # 484440 x 6
length(unique(ref$gene_id)) # 22799 
ref <- unique(ref) 
dim(ref) # 388598 x 6
ref <- subset(ref, id == "gene") # select only genes -- maybe i do need to do this step??
dim(ref) # 21837 x 6
length(unique(ref$gene_id)) # 21562 
# in this ref file, there are 21562 unique gene ids 
# Strange that the map file (created by me with stringTie) has more unique gene ids that the ref file (from NCBI). It could be because star/stringTie created 'novel' genes based on 
# exons splicing in different ways.
# I only need gene id and length...going to keep just those and throw out scaffold and gene start/stop
ref <- select(ref, c(length, gene_id)) # select only certain cols

# Build a df that merges gene IDs of expressed genes (poverA) from gcounts_filt_pdam and gene length from ref
# pdam_filt_ref <- merge(ref, gcounts_filt_pdam, by = "gene_id", all.x = TRUE)
# pdam_filt_ref <- na.omit(pdam_filt_ref)
# dim(pdam_filt_ref) # 13568 x 14
# I thought it supposed to be = to the number of expressed genes (13381), but some genes had no length (and so had NAs)
# Not going to keep them around because goseq needs gene length to do proper calculations
# Let's see if its okay that pdam_filt_ref = 13568 unique gene ids

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5; gcounts_filt_pdam) and the gene IDs of those genes in the merged annotation file (from the gene map; map_unique)
pdam_filt.map_unique <- merge(gcounts_filt_pdam, map_unique, by = "gene_id", all.x = TRUE, all.y = TRUE) # merge gene counts and map by gene_id
pdam_filt.map_unique <- na.omit(pdam_filt.map_unique)
length(unique(pdam_filt.map_unique$gene_id))
dim(pdam_filt.map_unique)
# This seems to merely confirm that there are 13881 DEGs that were actually identified in stringTie

#Find gene positions in ref corresponding to expressed genes 
pdam_filt.map_unique.ref <- merge(pdam_filt.map_unique, ref, by = "gene_id", all.x= TRUE) # merge pdam_filt.map_unique and ref by gene_idpdam_filt.map_unique.ref <- select(pdam_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x)) # remove specified cols
#pdam_filt.map_unique.ref <- na.omit(pdam_filt.map_unique.ref) 
#not going to omit NAs right now because some genes do not have a length (and thus have NAs)
#even though the genes that not have lengths will most likely not be used, still want to keep them around for now
dim(pdam_filt.map_unique.ref) # 13881 x 14

#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(pdam_filt.map_unique.ref, gene_id%in%DEG_pdam$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 549 
DEG <- DEG$gene_id # I believe I only need the gene ids here 
DEG <- unique(DEG) # okay unique DEGs # is 549
DEG_names <- as.vector(DEG) # turn values in vector to be used in goseq

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(pdam_filt.map_unique.ref$gene_id%in%DEG_names) # assign 0 and 1s
names(gene_vector)=unique(pdam_filt.map_unique.ref$gene_id) # give gene ids as names in veector
length(unique(names(gene_vector))) # 13881
head(gene_vector)

#Make ID vector
ID_vector <- pdam_filt.map_unique.ref$gene_id # isolate gene ids
head(ID_vector)
ID_vector <- unique(ID_vector) # make sure unique gene ids only 

#Make length vector
length_vector <- pdam_filt.map_unique.ref$length # isolate gene lengths 
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector (w/ gene names) by length of gene
# what should plot look like?





### Prepare GO term dataframe 
# Import GO terms
annot_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/pdam/pdam_GO_20210125.csv", header=TRUE) # this is the file that contains the pdam genes that were assigned GO terms by InterProScan
annot_GO <- select(annot_GO, -X) # remove random X col
annot_GO$GO.ID <- gsub(",", ";", annot_GO$GO.ID)
splitted <- strsplit(as.character(annot_GO$GO.ID), ";") #split into multiple GO ids
GO.terms_OG <- data.frame(v1 = rep.int(annot_GO$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms_OG) <- c("gene_id", "GO.ID")
GO.terms_OG[GO.terms_OG == "NA"] <- NA

head(GO.terms_OG)
tail(GO.terms_OG)
GO.terms_OG$GO.ID<- as.character(GO.terms_OG$GO.ID)
GO.terms_OG$GO.ID <- replace_na(GO.terms_OG$GO.ID, "unknown")
GO.terms_OG$GO.ID <- as.factor(GO.terms_OG$GO.ID)
GO.terms_OG$gene_id <- as.factor(GO.terms_OG$gene_id)
GO.terms_OG <- unique(GO.terms_OG)

dim(GO.terms_OG) # 31327 x 2
head(GO.terms_OG, 10)
tail(GO.terms_OG, 10)



### Perform GOseq
# Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# should I include this: test.cats=c("GO:CC", "GO:BP", "GO:MF") ?
GO.wall<-goseq(DEG.pwf, ID_vector, gene2cat=GO.terms_OG, method="Wallenius", use_genes_without_cat=TRUE)
# Using manually entered categories.
# Calculating the p-values...
# 'select()' returned 1:1 mapping between keys and columns
write.csv(GO.wall, file = "~/Desktop/plob_GO_ALL_20210327.csv")


#Subset enriched GO terms by category and save as csv
#How many enriched GO terms do we have?
class(GO.wall)
head(GO.wall)
tail(GO.wall)
nrow(GO.wall)

#Find only enriched GO terms that are statistically significant at cutoff - 0.05
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
enriched.GO.05 <- enriched.GO.05[order(-enriched.GO.05$numDEInCat),]
enriched.GO.05$term <- as.factor(enriched.GO.05$term)
head(enriched.GO.05)
MF <- subset(enriched.GO.05, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),] # 16
CC <- subset(enriched.GO.05, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),] # 2
BP <- subset(enriched.GO.05, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),] # 15
write.csv(MF, file = "~/Desktop/pdam_MF_Sig_Enriched_GO.05_20210327.csv")
write.csv(CC, file = "~/Desktop/pdam_CC_Sig_Enriched_GO.05_20210327.csv")
write.csv(BP, file = "~/Desktop/pdam_BP_Sig_Enriched_GO.05_20210327.csv")
write.csv(enriched.GO.05, file = "~/Desktop/pdam_Sig_Enriched_GO.05_ALL_20210327.csv")

# Merge GO terms and enriched list 
colnames(GO.terms_OG) <- c("gene_id", "category")
enriched_GO.terms <- merge(enriched.GO.05, GO.terms_OG, by = "category", all.x = TRUE)
DEG_treatment <- read.csv("~/Desktop/pdam_DEGs.all_treatment_20210326.csv", header = TRUE)
DEG_treatment <- select(DEG_treatment, -X)
#colnames(DEG_treatment)[1] <-"gene_id"
ByTreatment_GO.terms <- merge(DEG_treatment, enriched_GO.terms, by = "gene_id", all.x = TRUE)
ByTreatment_GO.terms <- na.omit(ByTreatment_GO.terms)
# now I have a df with DEGs gene names, treatment comparisons, GO terms, and term info
MF <- subset(ByTreatment_GO.terms, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),]
CC <- subset(ByTreatment_GO.terms, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),]
BP <- subset(ByTreatment_GO.terms, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),]
write.csv(ByTreatment_GO.terms, file = "~/Desktop/plob_ByTreatment_GO.terms_20210326.csv")
write.csv(MF, file = "~/Desktop/plob_MF_Sig_Enriched_ByTreatment_GO.05_20210326.csv")
write.csv(CC, file = "~/Desktop/plob_CC_Sig_Enriched_ByTreatment_GO.05_20210326.csv")
write.csv(BP, file = "~/Desktop/plob_BP_Sig_Enriched_ByTreatment_GO.05_20210326.csv")

# Plot terms by number of differentially expressed functions
MFplot <- MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
MFplot

CCplot <- CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
CCplot

BPplot <- BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
BPplot

GOplot <- grid.arrange(MFplot, CCplot, BPplot, ncol=3, clip="off")
ggsave("~/Desktop/pdam_GOplot_05_20210328.pdf", GOplot, width = 21, height = 21, units = c("in"))

# Combining all and ordering by pvalue
GOplot2_pvalue <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  coord_flip() +
  ylim(0,0.05) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("p-value") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2_pvalue
ggsave("~/Desktop/pdam_GOplot2_pvalue_05_20210328.pdf", GOplot2_pvalue, width = 28, height = 28, units = c("in"))







### GO slim
GoSlims <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv", header=TRUE)
colnames(GoSlims) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims <- merge(Gene.GO.IDs, GoSlims, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims <- Gene.GO.IDs.slims[-1,]
















































# Build a df that shows individual GO terms assigned with gene (ie one GO term per row)
split_GO <- strsplit(as.character(annot_GO$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID") 
# write csv to save gene names with one GO name per row. All GO terms still there, but listed individually 
write.csv(GO.terms, file = "~/Desktop/pdam_GOterms_ByGeneOnly.csv")
# merge pdam_filt.map_unique.ref and GO terms so that there will be a df with individual go terms, gene lengths, gene names, etc
GO.terms <- merge(pdam_filt_ref, GO.terms, by.x = "gene_id")
dim(GO.terms) # 1072 x 19
# Write csv to save file with gene names, one GO term per row, gene length, and gene counts
write.csv(GO.terms, file = "~/Desktop/pdam_GOterms_ByGene.csv")
GO.terms <- select(GO.terms, c(gene_id, GO.ID)) # select only gene id and GO.ID
GO.terms$GO.ID<- as.character(GO.terms$GO.ID) # make GO.ID a character variable
dim(GO.terms) 
GO.terms[GO.terms == 0] <- "unknown" # make terms that = 0 are replaced with the word 'unknown'
GO.terms <- unique(GO.terms) 
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown") # replace word 'unknown' with NAs
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID) # make GO.ID a factor variable
GO.terms$gene_id <- as.factor(GO.terms$gene_id) # make gene id a factor variable
head(GO.terms, 10)
tail(GO.terms, 10)
dim(GO.terms) # 686 x 2



gcount <- select(gcounts_filt_pdam, "gene_id")
merge <- merge(gcount, annot_GO, by = "gene_id", all.x = T)
write.csv(merge, file = "~/Desktop/merge_pdam_check.csv")


















### Perform GOseq
# Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# should I include this: test.cats=c("GO:CC", "GO:BP", "GO:MF") ?
GO.wall<-goseq(DEG.pwf, ID_vector, gene2cat=GO.terms, method="Wallenius", use_genes_without_cat=TRUE)
# Using manually entered categories.
# Calculating the p-values...
# 'select()' returned 1:1 mapping between keys and columns
#Subset enriched GO terms by category and save as csv
write.csv(GO.wall, file = "~/Desktop/pdam_GO_ALL.csv")
class(GO.wall) #data.frame
head(GO.wall)
tail(GO.wall)
#How many enriched GO terms?
nrow(GO.wall) # 253

#Find only enriched GO terms that are statistically significant at cutoff - 0.05
enriched.GO.05.a<-GO.wall$category[GO.wall$over_represented_pvalue<.05]
enriched.GO.05<-data.frame(enriched.GO.05.a)
colnames(enriched.GO.05) <- c("category")
enriched.GO.05 <- merge(enriched.GO.05, GO.wall, by="category")
enriched.GO.05 <- enriched.GO.05[order(-enriched.GO.05$numDEInCat),]
enriched.GO.05$term <- as.factor(enriched.GO.05$term)
head(enriched.GO.05)
MF <- subset(enriched.GO.05, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),]
CC <- subset(enriched.GO.05, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),]
BP <- subset(enriched.GO.05, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),]
write.csv(MF, file = "~/Desktop/pdam_MF_Sig_Enriched_GO.05.csv")
write.csv(CC, file = "~/Desktop/pdam_CC_Sig_Enriched_GO.05.csv")
write.csv(BP, file = "~/Desktop/pdam_BP_Sig_Enriched_GO.05.csv")
write.csv(enriched.GO.05, file = "~/Desktop/pdam_Sig_Enriched_GO.05_ALL.csv")

# Merge GO terms and enriched list 
colnames(GO.terms) <- c("gene_id", "category")
enriched_GO.terms <- merge(enriched.GO.05, GO.terms, by = "category", all.x = TRUE)
DEG_treatment <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_DEGs.all_treatment.csv", header = TRUE)
colnames(DEG_treatment)[1] <-"gene_id"
ByTreatment_GO.terms <- merge(DEG_treatment, enriched_GO.terms, by = "gene_id", all.x = TRUE)
ByTreatment_GO.terms <- na.omit(ByTreatment_GO.terms)
# now I have a df with DEGs gene names, treatment comparisons, GO terms, and term info
MF <- subset(ByTreatment_GO.terms, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),]
CC <- subset(ByTreatment_GO.terms, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),]
BP <- subset(ByTreatment_GO.terms, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),]
write.csv(ByTreatment_GO.terms, file = "~/Desktop/pdam_ByTreatment_GO.terms.csv")
write.csv(MF, file = "~/Desktop/pdam_MF_Sig_Enriched_ByTreatment_GO.05.csv")
write.csv(CC, file = "~/Desktop/pdam_CC_Sig_Enriched_ByTreatment_GO.05.csv")
write.csv(BP, file = "~/Desktop/pdam_BP_Sig_Enriched_ByTreatment_GO.05.csv")

# Use MF, CC, BP from enriched.GO.05 (ie not by treatment)
# Plot terms by number of differentially expressed functions
MFplot <- MF %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Molecular Function") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
MFplot

CCplot <- CC %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Cellular Component") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
CCplot

BPplot <- BP %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") +
  xlab("") +
  ylab("") +
  ggtitle("Biological Process") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank())#Set the plot background
BPplot

# save MF, CC, and BP in grid format 
GOplot <- grid.arrange(MFplot, CCplot, BPplot, ncol=3, clip="off")
ggsave("~/Desktop/pdam_GOplot_05.pdf", GOplot, width = 21, height = 21, units = c("in"))

# Combining all into one plot
GOplot2 <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  #geom_text(aes(label = over_represented_pvalue), hjust = 0, vjust = 0, size = 1) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2
ggsave("~/Desktop/pdam_GOplot2_05.pdf", GOplot2, width = 28, height = 28, units = c("in"))

# Combining all and ordering by pvalue
GOplot2_pvalue <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("p-value") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2_pvalue
ggsave("~/Desktop/pdam_GOplot2_pvalue_05.pdf", GOplot2_pvalue, width = 28, height = 32, units = c("in"))


# Combining all and ordering by pvalue -- BP
GOplot2_pvalue <- BP %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("p-value") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2_pvalue
ggsave("~/Desktop/pdam_BP_GOplot2_pvalue_05.pdf", GOplot2_pvalue, width = 28, height = 32, units = c("in"))

# Combining all and ordering by pvalue -- CC
GOplot2_pvalue <- CC %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("p-value") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2_pvalue
ggsave("~/Desktop/pdam_CC_GOplot2_pvalue_05.pdf", GOplot2_pvalue, width = 28, height = 32, units = c("in"))

# Combining all and ordering by pvalue -- MF
GOplot2_pvalue <- MF %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("p-value") +
  theme_bw() + #Set background color 
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) #Set the plot background #set title attributes
GOplot2_pvalue
ggsave("~/Desktop/pdam_MF_GOplot2_pvalue_05.pdf", GOplot2_pvalue, width = 28, height = 32, units = c("in"))







for(go in enriched.GO.05.a[1:25]){
  print(GOTERM[[go]]) 
  }





