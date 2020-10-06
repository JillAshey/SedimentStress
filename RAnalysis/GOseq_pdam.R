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
gcounts_filt_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_counts_filt.csv", header = TRUE)
dim(gcounts_filt_pdam) # 12873 rows x 16
for ( col in 1:ncol(gcounts_filt_pdam)){
  colnames(gcounts_filt_pdam)[col] <-  gsub("X", "", colnames(gcounts_filt_pdam)[col])
}
colnames(gcounts_filt_pdam)[1] <-"gene_id"
head(gcounts_filt_pdam)

DEG_pdam <- read.csv("~/Desktop/pdam_unique.sig.list.csv", header = TRUE)
dim(DEG_pdam)
for ( col in 1:ncol(DEG_pdam)){
  colnames(DEG_pdam)[col] <-  gsub("X", "", colnames(DEG_pdam)[col])
}
colnames(DEG_pdam)[1] <-"gene_id"
head(DEG_pdam)

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Pdam.merged.annotated.gtf", header=FALSE, sep="\t")
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
map <- map[!grepl("#", map$scaffold),]
map <- map[grep("LOC", map$attr), ]
map$gene_id <- regmatches(map$attr, gregexpr("(?<=gene_name).*", map$attr, perl = TRUE)) #removing everything in Symbol col up to LOC
map$gene_id <- gsub(";.*", "", map$gene_id)
map$gene_id <- gsub(" ", "", map$gene_id)
map <- subset(map, id=="transcript") # select only transcripts
map <- select(map, c(scaffold, gene.start, gene.stop, gene_id))
map_unique <- unique(map) 
dim(map_unique) 

ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- ref[!grepl("##", ref$scaffold),]
ref <- ref[grep("XP", ref$attr), ]
ref$prot <- gsub(";.*", "", ref$attr)
ref$prot <- gsub("ID=", "", ref$prot)
ref$prot <- gsub(".*-", "", ref$prot)
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC
ref$gene_id <- gsub(";.*", "", ref$gene_id)
ref <- select(ref, c(scaffold, gene.start, gene.stop, prot, gene_id))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref)
head(ref)

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
pdam_filt.map_unique <- merge(gcounts_filt_pdam, map_unique, by = "gene_id", all.x = TRUE)
dim(pdam_filt.map_unique) # should be same # of rows as gcounts_filt_pdam ??

#Find gene positions in ref corresponding to expressed genes 
pdam_filt.map_unique.ref <- merge(pdam_filt.map_unique, ref, by = "gene_id", all.x= TRUE)
pdam_filt.map_unique.ref <- select(pdam_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x))


#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(pdam_filt.map_unique.ref, gene_id%in%DEG_pdam$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 62
DEG_names <- as.vector(DEG$gene_id)

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(pdam_filt.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=pdam_filt.map_unique.ref$gene_id
head(gene_vector)

#Make ID vector
ID_vector <- pdam_filt.map_unique.ref$gene_id
head(ID_vector)

#Make length vector
length_vector <- pdam_filt.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene



























# Obtain names of all expressed ofav genes (poverA = 0.85,5), and all differentially expressed planuala genes (p<0.05)
gcounts_filt_pdam <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam_counts_filt.csv", header = TRUE)
dim(gcounts_filt_pdam) # 13881 rows x 13
for ( col in 1:ncol(gcounts_filt_pdam)){
  colnames(gcounts_filt_pdam)[col] <-  gsub("X", "", colnames(gcounts_filt_pdam)[col])
}
colnames(gcounts_filt_pdam)[1] <-"gene_id"
dim(gcounts_filt_pdam)
head(gcounts_filt_pdam)

DEG_pdam <- read.csv("~/Desktop/pdam_unique.sig.list.csv", header = TRUE)
dim(DEG_pdam) # 549 x 13
for ( col in 1:ncol(DEG_pdam)){
  colnames(DEG_pdam)[col] <-  gsub("X", "", colnames(DEG_pdam)[col])
}
colnames(DEG_pdam)[1] <-"gene_id"

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Pdam.merged.annotated.gtf", header=FALSE, sep="\t")
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
map <- map[!grepl("#", map$scaffold),]
map <- map[grep("LOC", map$attr), ]
map$gene_id <- regmatches(map$attr, gregexpr("(?<=gene_name).*", map$attr, perl = TRUE)) #removing everything in Symbol col up to LOC
map$gene_id <- gsub(";.*", "", map$gene_id)
map$gene_id <- gsub(" ", "", map$gene_id)
map <- subset(map, id=="transcript") # select only transcripts
map <- select(map, c(scaffold, gene.start, gene.stop, gene_id))
map_unique <- unique(map) 
dim(map_unique) # 26139 rows

# Import reference annotation file 
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- ref[!grepl("##", ref$scaffold),]
ref <- ref[grep("XP", ref$attr), ]
ref$prot <- gsub(";.*", "", ref$attr)
ref$prot <- gsub("ID=", "", ref$prot)
ref$prot <- gsub(".*-", "", ref$prot)
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC
ref$gene_id <- gsub(";.*", "", ref$gene_id)
ref <- select(ref, c(scaffold, gene.start, gene.stop, prot, gene_id))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref)
head(ref)

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
pdam_filt.map_unique <- merge(gcounts_filt_pdam, map_unique, by = "gene_id")
dim(pdam_filt.map_unique) # 16516

#Find gene positions in ref corresponding to expressed genes 
pdam_filt.map_unique.ref <- merge(pdam_filt.map_unique, ref, by = "gene_id", all.x= TRUE)
pdam_filt.map_unique.ref <- select(pdam_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x))
pdam_filt.map_unique.ref <- na.omit(pdam_filt.map_unique.ref)
pdam_filt.map_unique.ref <- unique(pdam_filt.map_unique.ref)

#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(pdam_filt.map_unique.ref, gene_id%in%DEG_pdam$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 549 ? but it is 599. not sure why
#DEG <- unique(DEG)
DEG_names <- as.vector(DEG$gene_id)
DEG_names <- unique(DEG_names)

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(pdam_filt.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=pdam_filt.map_unique.ref$gene_id
head(gene_vector)

#Make ID vector
ID_vector <- pdam_filt.map_unique.ref$gene_id
head(ID_vector)
#ID_vector <- unique(ID_vector)

#Make length vector
length_vector <- pdam_filt.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene














### Prepare GO term dataframe 
# Import GO terms
annot_GO <- read.csv("~/Desktop/pdam_GOterms.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
annot_GO <- merge(annot_GO, ref, by = "prot", all.x = TRUE)
annot_GO <- select(annot_GO, c(prot, Predict, GO_term, gene_id))
annot_GO <- unique(annot_GO)

annot_GO.ref <- merge(annot_GO, ref, by = "gene_id", all.x= TRUE)
annot_GO.ref <- unique(annot_GO.ref)
annot_GO.ref <- select(annot_GO.ref, c(gene_id, Predict, GO_term))

split_GO <- strsplit(as.character(annot_GO.ref$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO.ref$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")

GO.terms <- merge(pdam_filt.map_unique.ref, GO.terms, by.x = "gene_id")
dim(GO.terms) 

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
        #axis.text.y = element_text(size = 50),
                                  #face = 'bold') +
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

ggsave("~/Desktop/pdam_MFplot_05.pdf", MFplot, width = 21, height = 21, units = c("in"))
ggsave("~/Desktop/pdam_CCplot_05.pdf", CCplot, width = 21, height = 21, units = c("in"))
ggsave("~/Desktop/pdam_BPplot_05.pdf", BPplot, width = 21, height = 21, units = c("in"))

GOplot <- grid.arrange(MFplot, CCplot, BPplot, ncol=3, clip="off")
ggsave("~/Desktop/pdam_GOplot_05.pdf", GOplot, width = 21, height = 21, units = c("in"))

GOplot2 <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
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
ggsave("~/Desktop/pdam_GOplot2_05.pdf", GOplot2, width = 12, height = 12, units = c("in"))



# Not many GO terms that are statistically significant over 0.05. Try finding enriched GO terms are a statistically significant cutoff of 0.1
enriched.GO.1.a<-GO.wall$category[GO.wall$over_represented_pvalue<.1]
enriched.GO.1<-data.frame(enriched.GO.1.a)
colnames(enriched.GO.1) <- c("category")
enriched.GO.1 <- merge(enriched.GO.1, GO.wall, by="category")
enriched.GO.1 <- enriched.GO.1[order(-enriched.GO.1$numDEInCat),]
enriched.GO.1$term <- as.factor(enriched.GO.1$term)
head(enriched.GO.1)
MF <- subset(enriched.GO.1, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),]
CC <- subset(enriched.GO.1, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),]
BP <- subset(enriched.GO.1, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),]
write.csv(MF, file = "~/Desktop/pdam_MF_Sig_Enriched_GO.1.csv")
write.csv(CC, file = "~/Desktop/pdam_CC_Sig_Enriched_GO.1.csv")
write.csv(BP, file = "~/Desktop/pdam_BP_Sig_Enriched_GO.1.csv")
write.csv(enriched.GO.05, file = "~/Desktop/pdam_Sig_Enriched_GO.1_ALL.csv")

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
ggsave("~/Desktop/pdam_GOplot_1.pdf", GOplot, width = 21, height = 21, units = c("in"))

GOplot2 <- enriched.GO.1 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
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
ggsave("~/Desktop/pdam_GOplot2_1.pdf", GOplot2, width = 28, height = 28, units = c("in"))



