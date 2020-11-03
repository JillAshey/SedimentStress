# Title: Sediment stress - plob GOseq
# Author: Jill Ashey
# date: 10/16/20

### Gene Ontology Analysis of Differentially Expressed Genes in Plob, using Plut genome


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
gcounts_filt_plob <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_counts_filt.csv", header = TRUE)
dim(gcounts_filt_plob) # 15369 rows x 13
for ( col in 1:ncol(gcounts_filt_plob)){
  colnames(gcounts_filt_plob)[col] <-  gsub("X", "", colnames(gcounts_filt_plob)[col])
}
colnames(gcounts_filt_plob)[1] <-"gene_id"
head(gcounts_filt_plob)

# Load DEGs
DEG_plob <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob_unique.sig.list.csv", header = TRUE)
dim(DEG_plob) # 153 x 13
for ( col in 1:ncol(DEG_plob)){
  colnames(DEG_plob)[col] <-  gsub("X", "", colnames(DEG_plob)[col])
}
colnames(DEG_plob)[1] <-"gene_id"

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Plob.merged.annotated.gtf", header=FALSE, sep="\t")
colnames(map) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
map <- subset(map, id=="transcript") # select only transcripts
dim(map) # 31126
map$gene_id <- regmatches(map$attr, gregexpr("(?<=gene_name).*", map$attr, perl = TRUE))
map$gene_id <- gsub(";.*", "", map$gene_id)
# Remove any duplicates
map <- map %>% mutate_all(na_if, " ")
map <- select(map, c(scaffold, gene.start, gene.stop, gene_id))
map_unique <- unique(map) 
dim(map_unique) # 31126 rows

# Import reference annotation file 
ref <- read.csv("~/Desktop/GFFs/Plut.GFFannotation.fixed_transcript.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref <- subset(ref, id == "gene") 
ref$gene_id <- gsub(";.*", "", ref$attr)
ref$gene_id <- gsub("ID=", "", ref$gene_id)
ref <- select(ref, c(scaffold, gene.start, gene.stop, gene_id))
ref <- ref %>% mutate(ref, length = gene.stop - gene.start)
dim(ref) # 31125 x 5 --lost a row? where did it go

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
plob_filt.map_unique <- merge(gcounts_filt_plob, map_unique, by = "gene_id", all.x = TRUE) # merge gene counts and map by gene_id
dim(plob_filt.map_unique) # should be same # of rows as gcounts_filt_plob - it is!

#Find gene positions (start, stop, length) in ref corresponding to expressed genes 
plob_filt.map_unique.ref <- merge(plob_filt.map_unique, ref, by = "gene_id", all.x= TRUE) 
plob_filt.map_unique.ref <- select(plob_filt.map_unique.ref, -c(scaffold.x, gene.start.x, gene.stop.x))


#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(plob_filt.map_unique.ref, gene_id%in%DEG_plob$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 153--it is!
DEG_names <- as.vector(DEG$gene_id)

#Make vector of all expressed genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(plob_filt.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=plob_filt.map_unique.ref$gene_id
tail(gene_vector)

#Make ID vector
ID_vector <- plob_filt.map_unique.ref$gene_id
head(ID_vector)

#Make length vector
length_vector <- plob_filt.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene
# plot looks really good compared to the other species that I've done


### Prepare GO term dataframe 
# Import GO terms
annot_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOseq/plob_GOterms.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
colnames(annot_GO)[1] <-"gene_id"
# to properly merge with ref, I need to remove the .m1 from gene_ids beginning with plut2 and remove model and replace with TU for gene_ids beginning with jamg
# hopefully it doesn't mess anything else up
annot_GO$gene_id <- gsub(".m1", "", annot_GO$gene_id)
annot_GO$gene_id <- gsub("model", "TU", annot_GO$gene_id)
annot_GO.ref <- merge(annot_GO, ref, by = "gene_id", all.x= TRUE)

annot_GO.ref <- unique(annot_GO.ref)
annot_GO.ref <- select(annot_GO.ref, c(gene_id, Predict, GO_term))

split_GO <- strsplit(as.character(annot_GO.ref$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO.ref$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
# save df here to get gene names with one GO name per row. All GO terms still there, but listed individually 

GO.terms <- merge(plob_filt.map_unique.ref, GO.terms, by.x = "gene_id")
dim(GO.terms) # 53388 x 18

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
#GO.terms$gene_id <- gsub(" ", "", GO.terms$gene_id) # getting rid of weird spacing in gene col


### Perform GOseq
# Find enriched GO terms, "selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes for RNA-seq data"
# should I include this: test.cats=c("GO:CC", "GO:BP", "GO:MF") ?
GO.wall<-goseq(DEG.pwf, ID_vector, gene2cat=GO.terms, method="Wallenius", use_genes_without_cat=TRUE)
# Using manually entered categories.
# Calculating the p-values...
# 'select()' returned 1:1 mapping between keys and columns

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
MF <- MF[order(-MF$numDEInCat),] # 10
CC <- subset(enriched.GO.05, ontology=="CC")
CC <- CC[order(-CC$numDEInCat),] # 3
BP <- subset(enriched.GO.05, ontology=="BP")
BP <- BP[order(-BP$numDEInCat),] # 7
write.csv(MF, file = "~/Desktop/plob_MF_Sig_Enriched_GO.05.csv")
write.csv(CC, file = "~/Desktop/plob_CC_Sig_Enriched_GO.05.csv")
write.csv(BP, file = "~/Desktop/plob_BP_Sig_Enriched_GO.05.csv")
write.csv(enriched.GO.05, file = "~/Desktop/plob_Sig_Enriched_GO.05_ALL.csv")

# Merge GO terms and enriched list 
colnames(GO.terms) <- c("gene_id", "category")
enriched_GO.terms <- merge(enriched.GO.05, GO.terms, by = "category", all.x = TRUE)
DEG_treatment <- read.csv("~/Desktop/plob_DEGs.all_treatment.csv", header = TRUE)
#DEG_treatment <- select(DEG_treatment, -X)
colnames(DEG_treatment)[1] <-"gene_id"
ByTreatment_GO.terms <- merge(DEG_treatment, enriched_GO.terms, by = "gene_id", all.x = TRUE)
ByTreatment_GO.terms <- na.omit(ByTreatment_GO.terms)
# now I have a df with DEGs gene names, treatment comparisons, GO terms, and term info
write.csv(ByTreatment_GO.terms, file = "~/Desktop/plob_ByTreatment_GO.terms")

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
ggsave("~/Desktop/plob_GOplot_05.pdf", GOplot, width = 21, height = 21, units = c("in"))

# Combining all 
GOplot2 <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_point(size=3, aes(colour = ontology)) +
  #geom_text(aes(label = over_represented_pvalue), hjust = 0, vjust = 0, size = 1) +
  coord_flip() +
  #ylim(0,25) +
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
# Warning messages:
#   1: Removed 1 rows containing missing values (geom_segment). 
# 2: Removed 1 rows containing missing values (geom_point). 
# 3: Removed 1 rows containing missing values (geom_text). 
# Got rid of my protein binding rows. thats okay for now because they werent all that interesting anyway
ggsave("~/Desktop/plot_GOplot2_05.pdf", GOplot2, width = 28, height = 28, units = c("in"))
