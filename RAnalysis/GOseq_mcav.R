# Title: Sediment stress - mcav GOseq
# Author: Jill Ashey
# date: 10/1/20

### Gene Ontology Analysis of Differentially Expressed Genes in Mcav


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
gcounts_filt_mcav <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/mcav_counts_filtered.csv", header = TRUE)
dim(gcounts_filt_mcav) # 12873 rows x 16
for ( col in 1:ncol(gcounts_filt_mcav)){
  colnames(gcounts_filt_mcav)[col] <-  gsub("X", "", colnames(gcounts_filt_mcav)[col])
}
colnames(gcounts_filt_mcav)[1] <-"gene_id"
head(gcounts_filt_mcav)

# Load DEGs
DEG_mcav <- read.csv("~/Desktop/mcav_unique.sig.list.csv", header = TRUE)
dim(DEG_mcav)
for ( col in 1:ncol(DEG_mcav)){
  colnames(DEG_mcav)[col] <-  gsub("X", "", colnames(DEG_mcav)[col])
}
colnames(DEG_mcav)[1] <-"gene_id"

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Mcav.merged.annotated.gtf", header=FALSE, sep="\t")
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
map <- subset(map, id=="transcript") # select only transcripts
dim(map) 
map$gene_id <- gsub(";.*", "", map$attr)
map$gene_id <- gsub("transcript_id", "", map$gene_id)
map$gene_id <- gsub("-.*", "", map$gene_id)
# Remove any duplicates
map <- map %>% mutate_all(na_if, " ")
map <- select(map, c(scaffold, gene.start, gene.stop, gene_id))
map_unique <- unique(map) 
dim(map_unique) # 25605 rows

# Import reference annotation file 
ref <- read.table("~/Desktop/GFFs/Mcav.gff.annotations.fixed_transcript.gff3", sep = "\t", header = FALSE)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- subset(ref, id == "gene") 
ref$gene_id <- gsub(";.*", "", ref$attr)
ref$gene_id <- gsub("ID=", "", ref$gene_id)
ref <- select(ref, c(scaffold, gene.start, gene.stop, gene_id))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref)

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
mcav_filt.map_unique <- merge(gcounts_filt_mcav, map_unique, by = "gene_id", all.x = TRUE)
dim(mcav_filt.map_unique) # should be same # of rows as gcounts_filt_mcav - it is!

#Find gene positions in ref corresponding to expressed genes 
mcav_filt.map_unique.ref <- merge(mcav_filt.map_unique, ref, by = "gene_id", all.x= TRUE)
mcav_filt.map_unique.ref <- select(mcav_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x))


#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(mcav_filt.map_unique.ref, gene_id%in%DEG_mcav$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 62
DEG_names <- as.vector(DEG$gene_id)

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(mcav_filt.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=mcav_filt.map_unique.ref$gene_id
head(gene_vector)

#Make ID vector
ID_vector <- mcav_filt.map_unique.ref$gene_id
head(ID_vector)

#Make length vector
length_vector <- mcav_filt.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene


### Prepare GO term dataframe 
# Import GO terms
annot_GO <- read.csv("~/Desktop/mcav_GOterms.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
colnames(annot_GO)[1] <-"gene_id"
annot_GO$gene_id <- gsub("-.*", "", annot_GO$gene_id)
annot_GO.ref <- merge(annot_GO, ref, by = "gene_id", all.x= TRUE)

annot_GO.ref <- unique(annot_GO.ref)
annot_GO.ref <- select(annot_GO.ref, c(gene_id, Predict, GO_term))

split_GO <- strsplit(as.character(annot_GO.ref$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO.ref$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
# save df here to get gene names with one GO name per row. All GO terms still there, but listed individually 

GO.terms <- merge(mcav_filt.map_unique.ref, GO.terms, by.x = "gene_id")
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
write.csv(MF, file = "~/Desktop/mcav_MF_Sig_Enriched_GO.05.csv")
write.csv(CC, file = "~/Desktop/mcav_CC_Sig_Enriched_GO.05.csv")
write.csv(BP, file = "~/Desktop/mcav_BP_Sig_Enriched_GO.05.csv")
write.csv(enriched.GO.05, file = "~/Desktop/mcav_Sig_Enriched_GO.05_ALL.csv")

# Merge GO terms and enriched list 
colnames(GO.terms) <- c("gene_id", "category")
enriched_GO.terms <- merge(enriched.GO.05, GO.terms, by = "category", all.x = TRUE)
DEG_treatment <- read.csv("~/Desktop/mcav_DEGs.all_treatment.csv", header = TRUE)
colnames(DEG_treatment)[1] <-"gene_id"
ByTreatment_GO.terms <- merge(DEG_treatment, enriched_GO.terms, by = "gene_id", all.x = TRUE)
ByTreatment_GO.terms <- na.omit(ByTreatment_GO.terms)
# now I have a df with DEGs gene names, treatment comparisons, GO terms, and term info
write.csv(ByTreatment_GO.terms, file = "~/Desktop/mcav_ByTreatment_GO.terms")
MF <- subset(ByTreatment_GO.terms, ontology=="MF")
MF <- MF[order(-MF$numDEInCat),]
CC <- subset(ByTreatment_GO.terms, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),]
BP <- subset(ByTreatment_GO.terms, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),]

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
ggsave("~/Desktop/mcav_GOplot_05.pdf", GOplot, width = 21, height = 21, units = c("in"))

# Combining all 
GOplot2 <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  geom_text(aes(label = over_represented_pvalue), hjust = 0, vjust = 0, size = 1) +
  coord_flip() +
  ylim(0,15) +
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
ggsave("~/Desktop/mcav_GOplot2_05.pdf", GOplot2, width = 28, height = 28, units = c("in"))








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
write.csv(MF, file = "~/Desktop/mcav_MF_Sig_Enriched_GO.1.csv")
write.csv(CC, file = "~/Desktop/mcav_CC_Sig_Enriched_GO.1.csv")
write.csv(BP, file = "~/Desktop/mcav_BP_Sig_Enriched_GO.1.csv")
write.csv(enriched.GO.05, file = "~/Desktop/mcav_Sig_Enriched_GO.1_ALL.csv")

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
ggsave("~/Desktop/mcav_GOplot_1.pdf", GOplot, width = 21, height = 21, units = c("in"))

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
ggsave("~/Desktop/mcav_GOplot2_1.pdf", GOplot2, width = 28, height = 28, units = c("in"))






