# Title: Sediment stress - acerv GOseq
# Author: Jill Ashey
# date: 9/24/20

### Gene Ontology Analysis of Differentially Expressed Genes in Acerv


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

# Obtain names of all expressed Acerv genes (poverA = 0.85,5), and all differentially expressed planuala genes (p<0.05)
gcounts_filt_acerv <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv_counts_filtered.csv", header = TRUE)
dim(gcounts_filt_acerv) # 9055 rows x 16
colnames(gcounts_filt_acerv) <- c("gene_id", "19_T33_Ac_WK", "24_T12_Ac_FM", "25_ctl1_Ac_GF_1", "27_ctl2_Ac_YG_1", 
                                  "31_T22_Ac_UV", "35_T43_Ac_MT", "37_T13_Ac_ML", "38_T23_Ac_IN", "41_ctl3_Ac_RN_1",
                                  "45_T41_Ac_SC_1", "47_T31_Ac_JB", "52_T11_Ac_II", "53_T21_Ac_NH", "54_T42_Ac_JQ", 
                                  "57_T32_Ac_NM")
head(gcounts_filt_acerv)
gcounts_filt_acerv$gene_id <- gsub("TU", "model", gcounts_filt_acerv$gene_id) # sub model for TU to match gene names in annot and interpro files
# gene_id <- as.data.frame(gcounts_filt_acerv$gene_id)

DEG_acerv <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv_unique.sig.list.csv", header = TRUE)
dim(DEG_acerv) # 50 x 16
colnames(DEG_acerv) <- c("gene_id", "19_T33_Ac_WK", "24_T12_Ac_FM", "25_ctl1_Ac_GF_1", "27_ctl2_Ac_YG_1", 
                                  "31_T22_Ac_UV", "35_T43_Ac_MT", "37_T13_Ac_ML", "38_T23_Ac_IN", "41_ctl3_Ac_RN_1",
                                  "45_T41_Ac_SC_1", "47_T31_Ac_JB", "52_T11_Ac_II", "53_T21_Ac_NH", "54_T42_Ac_JQ", 
                                  "57_T32_Ac_NM")
head(DEG_acerv)
# gene_id <- as.data.frame(DEG_acerv$gene_id)

#Import merged annotated gtf file
map <- read.csv("~/Desktop/GFFs/Acerv.merged.annotated.gtf", header=FALSE, sep="\t")
map <- subset(map, V3=="transcript") # select only transcripts
dim(map) # 48478 x 9
map <- separate(map, V9, into = c("transcript_id", "gene_id", "gene_name", "xloc", "cmp_ref", "class_code", "tss_id"), sep = ";")
# columns were not evenly separated, so some things are in different rows, but I don't think I'll be using those cols anyway
map <- select(map, c(V1, V4, V5, transcript_id, gene_id)) # might need to include gene name and/or xloc at some point
colnames(map) <- c("scaffold", "start", "stop", "transcript_id", "gene_id")
map$transcript_id <- gsub("transcript_id", "", map$transcript_id)
map$gene_id <- gsub("gene_id", "", map$gene_id)
map$gene_id <- gsub(" ","",map$gene_id) #remove blank spaces 
#map$gene_id <- gsub("XLOC_[0-9][0-9][0-9][0-9][0-9][0-9]", "unknown", map$gene_id)
#head(map)
#colnames(map) <- c("gene_id", "gene_id_test")
map$gene_id <- gsub("model", "TU", map$gene_id)
#test <- subset(map, gene_id %in% gcounts_filt_acerv$gene_id)

# Remove any duplicates
map <- map %>% mutate_all(na_if, " ")
map <- select(map, c(transcript_id, gene_id))
map_unique <- unique(map) # might have to subset the unique values by either transcript_id or gene_id
dim(map_unique) #48478 rows
# transcript id has model in the value, while gene id has TU in value name 

# Import reference annotation file 
ref <- read.table("~/Desktop/GFFs/Acerv.GFFannotations.fixed_transcript.gff3", sep = "\t", header = FALSE)
ref <- subset(ref, V3 == "gene")
colnames(ref) <- c("scaffold", "gene.predict", "id", "start", "stop", "pos1", "pos2", "pos3", "attr", "gene_id")
ref <- select(ref, c(scaffold, start, stop, gene_id))
#ref$transcript_id <- gsub(".*;", "", ref$transcript_id)
#ref$transcript_id <- gsub("Name=", "", ref$transcript_id)
ref <- ref %>% mutate(ref, length = stop - start)
dim(ref) # 15243 x 2

# Build a dataframe that links the gene IDs of expressed genes (poverA = 0.85,5), the gene ids of those genes (from the gene map), and the gene lengths (from the annotation file)
# Because the gene_ids in gcount_filt_acerv have TU instead of model in names, I'm going to replace model with TU in that file. I tried without but it only gave me 3 genes.
#gcounts_filt_acerv$gene_id <- gsub("TU", "model", gcounts_filt_acerv$gene_id)
#colnames(map_unique) <- c("gene_id", "gene_id_real") # need to rename transcript id to gene id to see if I can properly merge it with gcounts file
acerv_filt.map_unique <- merge(gcounts_filt_acerv, map_unique, by = "gene_id", all.x = TRUE)
dim(acerv_filt.map_unique) # should be same # of rows as gcounts_filt_acerv
acerv_filt.map_unique$gene_id <- gsub("model", "TU", acerv_filt.map_unique$gene_id)

#Find gene positions in ref corresponding to expressed genes 
acerv_genes.map_unique.ref <- merge(acerv_filt.map_unique, ref, by = "gene_id", all.x = TRUE)
#acerv_genes.map_unique.ref <- select(acerv_genes.map_unique.ref, -gene_id_test)

#### Build GOSEQ vector 
#GOseq requires a vector of all genes and all differentially expressed genes. 
#Make gene vector
DEG <- filter(acerv_genes.map_unique.ref, gene_id%in%DEG_acerv$gene_id) #make vector of differentially expressed genes
dim(DEG) #should be 50
DEG_names <- as.vector(DEG$gene_id)

#Make vector of all expressed planuale genes (poverA = 0.85,5) with
#non-differentially expressed genes as 0 and differentially expressed genes as 1
gene_vector=as.integer(acerv_genes.map_unique.ref$gene_id%in%DEG_names)
names(gene_vector)=acerv_genes.map_unique.ref$gene_id
head(gene_vector)

#Make ID vector
ID_vector <- acerv_genes.map_unique.ref$gene_id
head(ID_vector)

#Make length vector
length_vector <- acerv_genes.map_unique.ref$length
head(length_vector)

#Calculate Probability Weighting Function
DEG.pwf<-nullp(gene_vector, ID_vector, bias.data=length_vector) #weight vector by length of gene


### Prepare GO term dataframe 
# Import GO term
annot_GO <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv_GOterms.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
colnames(annot_GO)[1] <-"gene_id"
annot_GO$gene_id <- gsub("model", "TU", annot_GO$gene_id)

split_GO <- strsplit(as.character(annot_GO$GO_term), ",")
GO.terms <- data.frame(v1 = rep.int(annot_GO$gene_id, sapply(split_GO, length)), v2 = unlist(split_GO)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")

GO.terms <- merge(acerv_genes.map_unique.ref, GO.terms, by.x="gene_id")
dim(GO.terms)

GO.terms <- select(GO.terms, gene_id, GO.ID)
GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
dim(GO.terms)
GO.terms[GO.terms == 0] <- "unknown"
GO.terms <- unique(GO.terms)
write.csv(GO.terms, file = "~/Desktop/acerv_GOterms.unique.csv")
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
# Warning message:
#   In goseq(DEG.pwf, ID_vector, gene2cat = GO.terms, method = "Wallenius",  :
#              Missing length data for 56% of genes.  Accuarcy of GO test will be reduced.


#### Explore enriched GO terms

#Subset enriched GO terms by category and save as csv
#How many enriched GO terms do we have
class(GO.wall)
dim(GO.wall) # 1415 x 7
head(GO.wall)
tail(GO.wall)
#Find only enriched GO terms that are statistically significant at cutoff 
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
write.csv(MF, file = "~/Desktop/acerv_MF_Sig_Enriched_GO.05.csv")
write.csv(CC, file = "~/Desktop/acerv_CC_Sig_Enriched_GO.05.csv")
write.csv(BP, file = "~/Desktop/acerv_BP_Sig_Enriched_GO.05.csv")
write.csv(enriched.GO.05, file = "~/Desktop/acerv_Sig_Enriched_GO.05_ALL.csv")

# Merge GO terms with gene ids 
colnames(GO.terms) <- c("gene_id","category")
GO.terms_gene_ids <- merge(enriched.GO.05, GO.terms, by = "category", all.x = TRUE)
GO.terms_gene_ids <- unique(GO.terms_gene_ids)
GO.terms_gene_ids.acerv_genes.map_unique.ref <- merge(acerv_genes.map_unique.ref, GO.terms_gene_ids, by = "gene_id", all.x = TRUE)
GO.terms_gene_ids.acerv_genes.map_unique.ref <- na.omit(GO.terms_gene_ids.acerv_genes.map_unique.ref)
GO.terms_gene_ids.acerv_genes.map_unique.ref <- select(GO.terms_gene_ids.acerv_genes.map_unique.ref, -transcript_id)
GO.terms_gene_ids.acerv_genes.map_unique.ref <- na.omit(GO.terms_gene_ids.acerv_genes.map_unique.ref)
write.csv(GO.terms_gene_ids.acerv_genes.map_unique.ref, file = "~/Desktop/acerv_GOterms_gene_ids.csv")

# Merge GO terms and enriched list 
colnames(GO.terms) <- c("gene_id", "category")
enriched_GO.terms <- merge(enriched.GO.05, GO.terms, by = "category", all.x = TRUE)
DEG_treatment <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv_DEGs.all_treatment.csv", header = TRUE)
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
write.csv(ByTreatment_GO.terms, file = "~/Desktop/acerv_ByTreatment_GO.terms.csv")
write.csv(MF, file = "~/Desktop/acerv_MF_Sig_Enriched_ByTreatment_GO.05.csv")
write.csv(CC, file = "~/Desktop/acerv_CC_Sig_Enriched_ByTreatment_GO.05.csv")
write.csv(BP, file = "~/Desktop/acerv_BP_Sig_Enriched_ByTreatment_GO.05.csv")

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
ggsave("~/Desktop/acerv_GOplot_05.pdf", GOplot, width = 21, height = 21, units = c("in"))

enriched.GO.05$ontology <- factor(enriched.GO.05$ontology, levels = enriched.GO.05[order(enriched.GO.05$over_represented_pvalue), "ontology"])

# Combining into one plot
GOplot2 <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, numDEInCat)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=numDEInCat) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=numDEInCat), color="grey") +
  geom_text(aes(label = over_represented_pvalue), hjust = -1, vjust = 0, size = 2) +
  geom_point(size=3, aes(colour = ontology)) +
  coord_flip() +
  ylim(0,25) +
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
ggsave("~/Desktop/acerv_GOplot2_05.pdf", GOplot2, width = 28, height = 28, units = c("in"))

# Combining into one plot and order by pvalue 
GOplot2_pvalue <- enriched.GO.05 %>% drop_na(ontology) %>% mutate(term = fct_reorder(term, over_represented_pvalue)) %>%
  mutate(term = fct_reorder(term, ontology)) %>%
  ggplot( aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_text(aes(label = numDEInCat), hjust = -1, vjust = 0.5, size = 3) +
  geom_point(size=3, aes(colour = ontology)) +
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
ggsave("~/Desktop/acerv_GOplot2_pvalue_05.pdf", GOplot2_pvalue, width = 28, height = 28, units = c("in"))













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
write.csv(MF, file = "~/Desktop/acerv_MF_Sig_Enriched_GO.1.csv")
write.csv(CC, file = "~/Desktop/acerv_CC_Sig_Enriched_GO.1.csv")
write.csv(BP, file = "~/Desktop/acerv_BP_Sig_Enriched_GO.1.csv")
write.csv(enriched.GO.1, file = "~/Desktop/acerv_Sig_Enriched_GO.1_ALL.csv")

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
ggsave("~/Desktop/GOplot_1.pdf", GOplot, width = 21, height = 21, units = c("in"))

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
ggsave("~/Desktop/GOplot2_1.pdf", GOplot2, width = 28, height = 28, units = c("in"))
