# Title: P. damicornis annotation files
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 08/24/20

# Currently trying to align Francois data to pdam genome. Alignment goes fine, but there are two different annotation files for Pdam (one from NCBI, one from Reef Genomics)
# Even thought they are supposed to be similar, it seems that the GFF file from RG is better annotated (has GO terms)

# Goal: look for differences between the Pdam annotation files from NCBI and Reef Genomics. 

# Load packages
library("tidyverse")
library("dplyr")
library("genefilter")
library("ggplot2")
library("spdep") 
library("clusterProfiler")
library("DataCombine")

## Reef Genomics
# Load gff 
pdamgff3_rg <- read.csv("~/Desktop/GFFs/pdam_annotation.gff3", header=FALSE, sep="\t", skip=2)
colnames(pdamgff3_rg) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(pdamgff3_rg) # 835,524 x 9
file.size("~/Desktop/GFFs/pdam_annotation.gff3")

length(unique(pdamgff3_rg$scaffold))

# Identify the unique parts -- gene predict
unique(pdamgff3_rg$Gene.Predict)
# [1] "maker"           "snap_masked"     "augustus_masked" ""                "repeatmasker"    "blastn"         
# [7] "est2genome"      "tblastx"         "cdna2genome"     "blastx"          "protein2genome"  "."              
# [13] "repeatrunner" 
# Find out many of each unique part is in the annotation file
G.P_maker_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="maker")
dim(G.P_maker_pdamgff_rg) # 29372 rows
G.P_snap_masked_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="snap_masked")
dim(G.P_snap_masked_pdamgff_rg) # 17795 rows
G.P_augustus_masked_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="augustus_masked")
dim(G.P_augustus_masked_pdamgff_rg) # 11521 rows
G.P_blank_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="")
dim(G.P_blank_pdamgff_rg) # 500 rows
G.P_repeatmasker_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="repeatmasker")
dim(G.P_repeatmasker_pdamgff_rg) # 5850 rows
G.P_blastn_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="blastn")
dim(G.P_blastn_pdamgff_rg) # 29216 rows
G.P_est2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="est2genome")
dim(G.P_est2genome_pdamgff_rg) # 25691 rows
G.P_tblastx_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="tblastx")
dim(G.P_tblastx_pdamgff_rg) # 5049 rows
G.P_cdna2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="cdna2genome")
dim(G.P_cdna2genome_pdamgff_rg) # 1572 rows
G.P_blastx_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="blastx")
dim(G.P_blastx_pdamgff_rg) # 423427 rows
G.P_protein2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="protein2genome")
dim(G.P_protein2genome_pdamgff_rg) # 285177 rows
G.P_._pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict==".")
dim(G.P_._pdamgff_rg) # 336 rows
G.P_repeatrunner_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="repeatrunner")
dim(G.P_repeatrunner_pdamgff_rg) # 18 rows

# Identify the unique parts -- ids
unique(pdamgff3_rg$id)
# [1] "gene"                        "mRNA"                        "exon"                        "five_prime_UTR"             
# [5] "CDS"                         "three_prime_UTR"             "match"                       "match_part"                 
# [9] ""                            "expressed_sequence_match"    "translated_nucleotide_match" "protein_match"              
# [13] "contig"
# Find out many of each unique part is in the annotation file
gene_pdamgff_rg <- subset(pdamgff3_rg, id=="gene")
dim(gene_pdamgff_rg) # 1675 rows
mrna_pdamgff_rg <- subset(pdamgff3_rg, id=="mRNA")
dim(mrna_pdamgff_rg) # 1675 rows
exon_pdamgff_rg <- subset(pdamgff3_rg, id=="exon")
dim(exon_pdamgff_rg) # 12504 rows
five_prime_UTR_pdamgff_rg <- subset(pdamgff3_rg, id=="five_prime_UTR")
dim(five_prime_UTR_pdamgff_rg) # 721 rows
CDS_pdamgff_rg <- subset(pdamgff3_rg, id=="CDS")
dim(CDS_pdamgff_rg) # 12037 rows
three_prime_UTR_pdamgff_rg <- subset(pdamgff3_rg, id=="three_prime_UTR")
dim(three_prime_UTR_pdamgff_rg) # 760 rows
match_pdamgff_rg <- subset(pdamgff3_rg, id=="match")
dim(match_pdamgff_rg) # 6587 rows
match_part_pdamgff_rg <- subset(pdamgff3_rg, id=="match_part")
dim(match_part_pdamgff_rg) # 617,859 rows
blank_pdamgff_rg <- subset(pdamgff3_rg, id=="") # not sure what this stands for, just blank rows 
dim(blank_pdamgff_rg) # 500 rows
expressed_sequence_match_pdamgff_rg <- subset(pdamgff3_rg, id=="expressed_sequence_match")  
dim(expressed_sequence_match_pdamgff_rg) # 7165 rows
translated_nucleotide_match_pdamgff_rg <- subset(pdamgff3_rg, id=="translated_nucleotide_match")  
dim(translated_nucleotide_match_pdamgff_rg) # 181 rows
protein_match_pdamgff_rg <- subset(pdamgff3_rg, id=="protein_match")  
dim(protein_match_pdamgff_rg) # 173,524 rows
contig_pdamgff_rg <- subset(pdamgff3_rg, id=="contig")  
dim(contig_pdamgff_rg) # 336 rows

# Total count (just to make sure it matches with dim of pdamgff3_rg)
1675 + 1675 + 12504 + 721 + 12037 + 760 + 6587 + 617859 + 500 + 7165 + 181 + 173524 + 336 # 835524

# Find which genes have GO terms
GO_pdamgff_rg <- filter(pdamgff3_rg, grepl("GO:", gene)) # only contains genes and mRNAs
dim(GO_pdamgff_rg) # 1588 rows - so only 1588 genes/mRNAs out of the 1675 for each have GO terms. And the mRNAs seem to correspond to the genes so the number may be half that
gene_GO_pdamgff_rg <- subset(GO_pdamgff_rg, id=="gene")
dim(gene_GO_pdamgff_rg) # 794 rows - which is half of 1588. So only 794 genes have GO terms associated with them


## NCBI
# Load gff
pdamgff3_NCBI <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff", header=FALSE, sep="\t", skip=6)
colnames(pdamgff3_NCBI) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(pdamgff3_NCBI) # 516,693 x 9

length(unique(pdamgff3_NCBI$scaffold))

# Identify the unique parts -- gene predict 
unique(pdamgff3_NCBI$Gene.Predict)
# [1] ""            "RefSeq"      "Gnomon"      "tRNAscan-SE" "cmsearch"   
# Find out many of each unique part is in the annotation file
G.P_blank_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="")  
dim(G.P_blank_pdamgff_NCBI) # 8787 rows
G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="RefSeq")  
dim(G.P_RefSeq_pdamgff_NCBI) # 23494 rows
G.P_Gnomon_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="Gnomon")  
dim(G.P_Gnomon_pdamgff_NCBI) # 483246 rows
G.P_tRNAscanSE_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="tRNAscan-SE")  
dim(G.P_tRNAscanSE_pdamgff_NCBI) # 1013 rows
G.P_cmsearch_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="cmsearch")  
dim(G.P_cmsearch_pdamgff_NCBI) # 153 rows


# Identify the unique parts -- ids
unique(pdamgff3_NCBI$id)
# [1] ""           "region"     "gene"       "mRNA"       "exon"       "CDS"        "pseudogene" "lnc_RNA"    "tRNA"       "cDNA_match" "transcript"
# [12] "snRNA"      "snoRNA"     "guide_RNA"  "rRNA"       "intron"
# Find out many of each unique part is in the annotation file
blank_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="")  
dim(blank_pdamgff_NCBI) # 8787 rows
region_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="region")  
dim(region_pdamgff_NCBI) # 4393 rows
gene_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="gene")  
dim(gene_pdamgff_NCBI) # 21837 rows
mrna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="mRNA")  
dim(mrna_pdamgff_NCBI) # 25170 rows
exon_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="exon")  
dim(exon_pdamgff_NCBI) # 226,892 rows
CDS_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="CDS")  
dim(CDS_pdamgff_NCBI) # 206,863 rows
pseudogene_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="pseudogene")  
dim(pseudogene_pdamgff_NCBI) # 1237 rows
lnc_rna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="lnc_RNA")  
dim(lnc_rna_pdamgff_NCBI) # 1771 rows
trna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="tRNA")  
dim(trna_pdamgff_NCBI) # 329 rows
cDNA_match_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="cDNA_match")  
dim(cDNA_match_pdamgff_NCBI) # 19065 rows
transcript_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="transcript")  
dim(transcript_pdamgff_NCBI) # 295 rows
snRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="snRNA")  
dim(snRNA_pdamgff_NCBI) # 25 rows
snoRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="snoRNA")  
dim(snoRNA_pdamgff_NCBI) # 24 rows
guide_RNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="guide_RNA")  
dim(guide_RNA_pdamgff_NCBI) # 1 row
rRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="rRNA")  
dim(rRNA_pdamgff_NCBI) # 3 rows
intron_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="intron")  
dim(intron_pdamgff_NCBI) # 1 row

# Total count (just to make sure it matches with dim of pdamgff3_NCBI)
8787 + 4393 + 21837 + 25170 + 226892 + 206863 + 1237 + 1771 + 329 + 19065 + 295 + 25 + 24 + 1 + 3 + 1 # 516,693

# Find which genes have GO terms
GO_pdamgff_NCBI <- filter(pdamgff3_NCBI, grepl("GO:", gene)) # only contains genes and mRNAs
dim(GO_pdamgff_NCBI) # no GO terms 



## Summary 
# Removing blank rows to get accurate dimensions
pdamgff3_rg <- subset(pdamgff3_rg, id!="")
dim(pdamgff3_rg) # 835,024 rows
pdamgff3_NCBI <- subset(pdamgff3_NCBI, id!="")
dim(pdamgff3_NCBI) # 507,906 rows

# Top 3 ids for Reef Genomics gff:
  # 1) match_part = 617,859
  # 2) protein_match = 173,524
  # 3) exon = 12,504

# Top 3 ids for NCBI gff:
  # 1) exon = 226,892
  # 2) CDS = 206,863 
  # 3) mRNA = 25,170

# The two gffs shared these ids: gene, mRNA, exon, CDS
#               pdamgff3_rg             vs               pdamgff3_NCBI
#           gene = 1,675                                gene = 21,837
#           mRNA = 1,675                                mRNA = 25,170
#           exon = 12,504                               exon = 226,892

# only Reef Genomics gff had GO terms associated with genes 




### Looking at genome files for Pdam
pdamfasta_rg <- read.csv("~/Desktop/GFFs/pdam_scaffolds.fasta", header=FALSE, sep="\t")
dim(pdamfasta_rg) # 3,912,153 x 1
GO_pdamgff_rg <- filter(pdamgff3_rg, grepl("GO:", gene)) # only contains genes and mRNAs







## Messing around with Reef Genomics gff3 file to see if I can get same results as Polina for GO
# original file from RG, no transcript id added
pdam.gff <- read.csv(file="~/Desktop/GFFs/pdam_annotation.gff3", header=FALSE, sep="\t", skip=2)
colnames(pdam.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(pdam.gff) # 835524 x 9
unique(pdam.gff$id)
# x <- grep("maker", pdam.gff$Gene.Predict)
# length(x)
pdam.gff_maker <- subset(pdam.gff, Gene.Predict == "maker") 
dim(pdam.gff_maker) # 29372 x 9
pdam.gff_maker_genes_only <- subset(pdam.gff_maker, id=="gene")
dim(pdam.gff_maker_genes_only) # 1675 x 9
unique(pdam.gff_maker_genes_only$id)

# file from RG with transcript_id= added in last col
pdam.gff_fixed_transcript <- read.csv(file="~/Desktop/pdam_annotation_AddTranscript_id_fixed.gff3", header=FALSE, sep="\t")
colnames(pdam.gff_fixed_transcript) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(pdam.gff_fixed_transcript) # 835524 x 9
unique(pdam.gff_fixed_transcript$id)
pdam.gff_fixed_transcript_maker <- subset(pdam.gff_fixed_transcript, Gene.Predict == "maker") 
dim(pdam.gff_fixed_transcript_maker) # 29372 x 9
pdam.gff_fixed_transcript_genes_only <- subset(pdam.gff_fixed_transcript_maker, id=="gene")
dim(pdam.gff_fixed_transcript_genes_only) # 1675 x 9
unique(pdam.gff_fixed_transcript_genes_only$id)

write.table(pdam.gff_fixed_transcript_genes_only, file="~/Desktop/GFFs/pdam.gff_fixed_transcript_genes_only.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

# Getting GO terms 
pdam_original_GO <- pdam.gff_maker_genes_only %>% 
  filter(str_detect(gene, "GO:"))
dim(pdam_original_GO) # 794 x 9

pdam_fixed_GO <- pdam.gff_fixed_transcript_genes_only %>% 
  filter(str_detect(gene, "GO:"))
dim(pdam_fixed_GO) # 794 x 9


# This file was made in terminal. code was: grep "match(ctrl-v,tab)gene" pdam_annotation.gff3 > tmp_original.gff3. Not sure how grep actually worked here, but it changed all maker stuff to gene only 
# Trying to figure out what is up with it 
tmp_original.gff3 <- read.csv(file="~/Desktop/GFFs/tmp_original.gff3", header=FALSE, sep="\t")
colnames(tmp_original.gff3) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
dim(tmp_original.gff3) # 26077
unique(tmp_original.gff3$id)
pdam_original_GO_tmpfile <- tmp_original.gff3 %>% 
  filter(str_detect(gene, "GO:"))
dim(pdam_original_GO_tmpfile)
duplicated(pdam_original_GO_tmpfile$gene)
write.table(pdam.gff_fixed_transcript_genes_only, file="~/Desktop/GFFs/pdam.gff_fixed_transcript_genes_only.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)








### Looking at fasta files 
length(unique(pdamgff3_NCBI$scaffold)) # 8788
length(unique(pdamgff3_rg$scaffold)) # 338
# Both of these numbers dont really make sense to me. The # of unique scaffolds in NCBI is greater than the scaffolds in its genome and the # of unique scaffolds is really low in RG


## NCBI fasta scaffolds
scaffolds_NCBI <- read.csv(file="~/Desktop/scaffold_NCBI_name.txt", header=FALSE)
colnames(scaffolds_NCBI) <- c("scaffold", "RG_scaffold", "Kind")
head(scaffolds_NCBI)
dim(scaffolds_NCBI) # 4393 x 3
scaffolds_NCBI$scaffold <- gsub("Pocillopora damicornis isolate RSMAS unplaced genomic scaffold", "", scaffolds_NCBI$scaffold)
scaffolds_NCBI$length <- sub(":.*", "", scaffolds_NCBI$scaffold)
scaffolds_NCBI$scaffold <- gsub(".*>", "", scaffolds_NCBI$scaffold)
scaffolds_NCBI$RG_scaffold <- gsub("ASM370409v1 ", "", scaffolds_NCBI$RG_scaffold)
length(unique(scaffolds_NCBI$scaffold))
length(unique(scaffolds_NCBI$RG_scaffold))
# I have a list associating the NCBI scaffolds with the Reef Genomics scaffolds 

# Clean up NCBI gff
pdamgff3_NCBI <- pdamgff3_NCBI %>%
  filter(str_detect(scaffold, "##", negate = TRUE))
length(unique(pdamgff3_NCBI$scaffold)) # 4393 - aha! once I removed the blanks, the unique scaffold # in gff file now == unique scaffold # in fasta file  

merge_NCBI <- merge(scaffolds_NCBI, pdamgff3_NCBI, by = "scaffold", all=TRUE)



