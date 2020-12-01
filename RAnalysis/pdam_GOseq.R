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


# Obtain names of all expressed ofav genes (poverA = 0.85,5), and all differentially expressed planuala genes (p<0.05)
gcounts_filt_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_counts_filt.csv", header = TRUE) # read data in
dim(gcounts_filt_pdam) # 13881 rows x 13
for ( col in 1:ncol(gcounts_filt_pdam)){
  colnames(gcounts_filt_pdam)[col] <-  gsub("X", "", colnames(gcounts_filt_pdam)[col]) # remove X in front of col names
}
colnames(gcounts_filt_pdam)[1] <-"gene_id" # make colnames a true column called gene_id
head(gcounts_filt_pdam)
length(unique(gcounts_filt_pdam$gene_id)) # 13881 total unique gene ids

DEG_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_unique.sig.list.csv", header = TRUE) # read in list of significant pdam genes
dim(DEG_pdam) # 549 x 13
for ( col in 1:ncol(DEG_pdam)){
  colnames(DEG_pdam)[col] <-  gsub("X", "", colnames(DEG_pdam)[col]) # remove X in front of col names
}
colnames(DEG_pdam)[1] <-"gene_id" # make colnames a true column called gene_id
head(DEG_pdam)
length(unique(DEG_pdam$gene_id)) # 549 total unique gene ids 

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Pdam.merged.annotated.gtf", header=FALSE, sep="\t") # read in merged annotated pdam file created with stringTie
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name columns
map <- map[!grepl("#", map$scaffold),] # remove rows that have a # in scaffold col
map <- map[grep("LOC", map$attr), ] # select all rows that have LOC in attr col
map$gene_id <- regmatches(map$attr, gregexpr("(?<=gene_name).*", map$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and making new col called gene_id
map$gene_id <- gsub(";.*", "", map$gene_id) # remove everything after ;
map$gene_id <- gsub(" ", "", map$gene_id) # remove any blank spaces
map <- subset(map, id=="transcript") # select only transcripts
map <- select(map, c(scaffold, gene.start, gene.stop, gene_id)) # keep only specified cols
map_unique <- unique(map) 
dim(map_unique) # 26139 x 4
length(unique(map_unique$gene_id)) # 22734 unique geneids


# read in reference annotation gff file
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
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
ref <- subset(ref, id == "gene") # select only genes -- maybe i do need to do this step??
dim(ref) # 21837 x 6
#ref <- unique(ref)

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
pdam_filt.map_unique <- merge(map_unique, gcounts_filt_pdam, by = "gene_id", all.x = TRUE) # merge gcounts_filt_pdam and map_unique by gene_id
dim(pdam_filt.map_unique) # should be same # of rows as gcounts_filt_pdam ?? but it is not... dim is 26139 x 16
pdam_filt.map_unique <- na.omit(pdam_filt.map_unique)
dim(pdam_filt.map_unique) # 16516 rows...still not same as gcounts_filt_pdam
length(unique(pdam_filt.map_unique$gene_id)) # so there are correct number of ids...maybe I have to remove scaffold?
pdam_filt.map_unique <- select(pdam_filt.map_unique, -scaffold)
pdam_filt.map_unique <- unique(pdam_filt.map_unique)
dim(pdam_filt.map_unique) # still the same...16516
pdam_filt.map_unique <- select(pdam_filt.map_unique, -c(gene.start, gene.stop)) # remove gene start and stop? idk trying to get number down to match gcounts_filt_pdam
pdam_filt.map_unique <- unique(pdam_filt.map_unique) # okay finally at 13881! But i will probably need gene lengths......
# well maybe I wont need gene lengths because i calculate gene lengths in ref above 

#Find gene positions in ref corresponding to expressed genes 
pdam_filt.map_unique.ref <- merge(pdam_filt.map_unique, ref, by = "gene_id", all.x= TRUE) # merge pdam_filt.map_unique and ref by gene_idpdam_filt.map_unique.ref <- select(pdam_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x)) # remove specified cols
pdam_filt.map_unique.ref <- select(pdam_filt.map_unique.ref, -c(scaffold, gene.start, gene.stop))
dim(pdam_filt.map_unique.ref) # 13881 x 18
#pdam_filt.map_unique.ref <- na.omit(pdam_filt.map_unique.ref) # remove NAs - not sure
  

#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(pdam_filt.map_unique.ref, gene_id%in%DEG_pdam$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 549 - yes
DEG <- DEG$gene_id # I believe I only need the gene ids here 
DEG <- unique(DEG) # okay unique DEGs # is 549
DEG_names <- as.vector(DEG)

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(pdam_filt.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=unique(pdam_filt.map_unique.ref$gene_id)
length(unique(names(gene_vector))) # 13881
head(gene_vector)

#Make ID vector
ID_vector <- pdam_filt.map_unique.ref$gene_id
head(ID_vector)
ID_vector <- unique(ID_vector)

#Make length vector
length_vector <- pdam_filt.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene


### Prepare GO term dataframe 
# Import GO terms
annot_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/pdam_GOterms_DEGs.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
colnames(annot_GO)[1] <-"gene_id"
annot_GO.ref <- merge(annot_GO, ref, by = "gene_id", all.x= TRUE)

annot_GO.ref <- unique(annot_GO.ref)
annot_GO.ref <- select(annot_GO.ref, c(gene_id, GO_term))

split_GO <- strsplit(as.character(annot_GO.ref$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO.ref$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID") 
# save df here to get gene names with one GO name per row. All GO terms still there, but listed individually 
write.csv(GO.terms, file = "~/Desktop/pdam_GOterms_ByGene.csv")

GO.terms <- merge(pdam_filt.map_unique.ref, GO.terms, by.x = "gene_id")
dim(GO.terms) # 1072 x 19

GO.terms <- select(GO.terms, c(gene_id, GO.ID))
GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
dim(GO.terms)
GO.terms[GO.terms == 0] <- "unknown"
GO.terms <- unique(GO.terms)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
head(GO.terms, 10)
tail(GO.terms, 10)


### Perform GOseq
# Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# should I include this: test.cats=c("GO:CC", "GO:BP", "GO:MF") ?
GO.wall<-goseq(DEG.pwf, ID_vector, gene2cat=GO.terms, method="Wallenius", use_genes_without_cat=TRUE)
# Using manually entered categories.
# Calculating the p-values...
# 'select()' returned 1:1 mapping between keys and columns

#Subset enriched GO terms by category and save as csv
#How many enriched GO terms do we have
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














