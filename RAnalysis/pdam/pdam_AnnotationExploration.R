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
# pdamgff3_rg <- na.omit(pdamgff3_rg)
# dim(pdamgff3_rg) # 835024 x 9
# length(unique(pdamgff3_rg$scaffold))

# Identify the unique parts -- gene predict
unique(pdamgff3_rg$Gene.Predict)
# [1] "maker"           "snap_masked"     "augustus_masked" ""                "repeatmasker"    "blastn"         
# [7] "est2genome"      "tblastx"         "cdna2genome"     "blastx"          "protein2genome"  "."              
# [13] "repeatrunner" 
# Find out many of each unique part is in the annotation file
##Maker
G.P_maker_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="maker")
dim(G.P_maker_pdamgff_rg) # 29372 rows
unique(G.P_maker_pdamgff_rg$id)
# [1] "gene"            "mRNA"            "exon"            "five_prime_UTR"  "CDS"             "three_prime_UTR"
# Find out many of each unique part is in the maker subset of annotation file
gene_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="gene")
dim(gene_G.P_maker_pdamgff_rg) # 1675 rows
mRNA_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="mRNA")
dim(mRNA_G.P_maker_pdamgff_rg) # 1675 rows
exon_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="exon")
dim(exon_G.P_maker_pdamgff_rg) # 12504 rows
five_prime_UTR_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="five_prime_UTR")
dim(five_prime_UTR_G.P_maker_pdamgff_rg) # 721 rows
CDS_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="CDS")
dim(CDS_G.P_maker_pdamgff_rg) # 12037 rows
three_prime_UTR_G.P_maker_pdamgff_rg <- subset(G.P_maker_pdamgff_rg, id=="three_prime_UTR")
dim(three_prime_UTR_G.P_maker_pdamgff_rg) # 760 rows
##SNAP masked
G.P_snap_masked_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="snap_masked")
dim(G.P_snap_masked_pdamgff_rg) # 17795 rows
unique(G.P_snap_masked_pdamgff_rg$id)
# [1] "match"      "match_part"
match_G.P_snap_masked_pdamgff_rg <- subset(G.P_snap_masked_pdamgff_rg, id=="match")
dim(match_G.P_snap_masked_pdamgff_rg) # 2344 rows
match_part_G.P_snap_masked_pdamgff_rg <- subset(G.P_snap_masked_pdamgff_rg, id=="match_part")
dim(match_part_G.P_snap_masked_pdamgff_rg) # 15451 rows
##Augustus masked
G.P_augustus_masked_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="augustus_masked")
dim(G.P_augustus_masked_pdamgff_rg) # 11521 rows
unique(G.P_augustus_masked_pdamgff_rg$id)
# [1] "match"      "match_part"
match_G.P_augustus_masked_pdamgff_rg <- subset(G.P_augustus_masked_pdamgff_rg, id=="match")
dim(match_G.P_augustus_masked_pdamgff_rg) # 1318 rows
match_part_G.P_augustus_masked_pdamgff_rg <- subset(G.P_augustus_masked_pdamgff_rg, id=="match_part")
dim(match_part_G.P_augustus_masked_pdamgff_rg) # 10203 rows
##Blank
G.P_blank_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="")
dim(G.P_blank_pdamgff_rg) # 500 rows
unique(G.P_blank_pdamgff_rg$id)
##RepeatMasker
G.P_repeatmasker_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="repeatmasker")
dim(G.P_repeatmasker_pdamgff_rg) # 5850 rows
unique(G.P_repeatmasker_pdamgff_rg$id)
# [1] "match"      "match_part"
match_G.P_repeatmasker_pdamgff_rg <- subset(G.P_repeatmasker_pdamgff_rg, id=="match")
dim(match_G.P_repeatmasker_pdamgff_rg) # 2925 rows
match_part_G.P_repeatmasker_pdamgff_rg <- subset(G.P_repeatmasker_pdamgff_rg, id=="match_part")
dim(match_part_G.P_repeatmasker_pdamgff_rg) # 2925 rows
##BLASTn
G.P_blastn_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="blastn")
dim(G.P_blastn_pdamgff_rg) # 29216 rows
unique(G.P_blastn_pdamgff_rg$id)
# [1] "expressed_sequence_match" "match_part" 
expressed_sequence_match_G.P_blastn_pdamgff_rg <- subset(G.P_blastn_pdamgff_rg, id=="expressed_sequence_match")
dim(expressed_sequence_match_G.P_blastn_pdamgff_rg) # 3395 rows
match_part_G.P_blastn_pdamgff_rg <- subset(G.P_blastn_pdamgff_rg, id=="match_part")
dim(match_part_G.P_blastn_pdamgff_rg) # 25821 rows
##est2genome
G.P_est2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="est2genome")
dim(G.P_est2genome_pdamgff_rg) # 25691 rows
unique(G.P_est2genome_pdamgff_rg$id)
# [1] "expressed_sequence_match" "match_part"  
expressed_sequence_match_G.P_est2genome_pdamgff_rg <- subset(G.P_est2genome_pdamgff_rg, id=="expressed_sequence_match")
dim(expressed_sequence_match_G.P_est2genome_pdamgff_rg) # 3590 rows
match_part_G.P_est2genome_pdamgff_rg <- subset(G.P_est2genome_pdamgff_rg, id=="match_part")
dim(match_part_G.P_est2genome_pdamgff_rg) # 22101 rows
##tBLASTx
G.P_tblastx_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="tblastx")
dim(G.P_tblastx_pdamgff_rg) # 5049 rows
unique(G.P_tblastx_pdamgff_rg$id)
# [1] "translated_nucleotide_match" "match_part"
translated_nucleotide_match_G.P_tblastx_pdamgff_rg <- subset(G.P_tblastx_pdamgff_rg, id=="translated_nucleotide_match")
dim(translated_nucleotide_match_G.P_tblastx_pdamgff_rg) # 181 rows
match_part_G.P_tblastx_pdamgff_rg <- subset(G.P_tblastx_pdamgff_rg, id=="match_part")
dim(match_part_G.P_tblastx_pdamgff_rg) # 4868 rows
##cDNA2genome
G.P_cdna2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="cdna2genome")
dim(G.P_cdna2genome_pdamgff_rg) # 1572 rows
unique(G.P_cdna2genome_pdamgff_rg$id)
# [1] "expressed_sequence_match" "match_part"  
expressed_sequence_match_G.P_cdna2genome_pdamgff_rg <- subset(G.P_cdna2genome_pdamgff_rg, id=="expressed_sequence_match")
dim(expressed_sequence_match_G.P_cdna2genome_pdamgff_rg) # 180 rows
match_part_G.P_cdna2genome_pdamgff_rg <- subset(G.P_cdna2genome_pdamgff_rg, id=="match_part")
dim(match_part_G.P_cdna2genome_pdamgff_rg) # 1392 rows
##BLASTx
G.P_blastx_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="blastx")
dim(G.P_blastx_pdamgff_rg) # 423427 rows
unique(G.P_blastx_pdamgff_rg$id)
# [1] "protein_match" "match_part"
protein_match_G.P_cdna2genome_pdamgff_rg <- subset(G.P_blastx_pdamgff_rg, id=="protein_match")
dim(protein_match_G.P_cdna2genome_pdamgff_rg) # 81561 rows
match_part_G.P_blastx_pdamgff_rg <- subset(G.P_blastx_pdamgff_rg, id=="match_part")
dim(match_part_G.P_blastx_pdamgff_rg) # 341866 rows 
##protein2genome 
G.P_protein2genome_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="protein2genome")
dim(G.P_protein2genome_pdamgff_rg) # 285177 rows
unique(G.P_protein2genome_pdamgff_rg$id)
# [1] "protein_match" "match_part" 
protein_match_G.P_protein2genome_pdamgff_rg <- subset(G.P_protein2genome_pdamgff_rg, id=="protein_match")
dim(protein_match_G.P_protein2genome_pdamgff_rg) # 91954 rows
match_part_G.P_protein2genome_pdamgff_rg <- subset(G.P_protein2genome_pdamgff_rg, id=="match_part")
dim(match_part_G.P_protein2genome_pdamgff_rg) # 193223 rows 
##period
G.P_._pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict==".")
dim(G.P_._pdamgff_rg) # 336 rows
unique(G.P_._pdamgff_rg$id)
# [1] "contig"
##RepeatRunner
G.P_repeatrunner_pdamgff_rg <- subset(pdamgff3_rg, Gene.Predict=="repeatrunner")
dim(G.P_repeatrunner_pdamgff_rg) # 18 rows
unique(G.P_repeatrunner_pdamgff_rg$id)
# [1] "protein_match" "match_part" 
protein_match_G.P_repeatrunner_pdamgff_rg <- subset(G.P_repeatrunner_pdamgff_rg, id=="protein_match")
dim(protein_match_G.P_repeatrunner_pdamgff_rg) # 9 rows
match_part_G.P_repeatrunner_pdamgff_rg <- subset(G.P_repeatrunner_pdamgff_rg, id=="match_part")
dim(match_part_G.P_repeatrunner_pdamgff_rg) # 9 rows 


# Identify the unique parts -- ids
unique(pdamgff3_rg$id)
# [1] "gene"                        "mRNA"                        "exon"                        "five_prime_UTR"             
# [5] "CDS"                         "three_prime_UTR"             "match"                       "match_part"                 
# [9] ""                            "expressed_sequence_match"    "translated_nucleotide_match" "protein_match"              
# [13] "contig"
# Find out many of each unique part is in the annotation file
##gene
gene_pdamgff_rg <- subset(pdamgff3_rg, id=="gene")
dim(gene_pdamgff_rg) # 1675 rows
unique(gene_pdamgff_rg$Gene.Predict)
# [1] "maker"
##mRNA
mrna_pdamgff_rg <- subset(pdamgff3_rg, id=="mRNA")
dim(mrna_pdamgff_rg) # 1675 rows
unique(mrna_pdamgff_rg$Gene.Predict)
# [1] "maker"
##Exon
exon_pdamgff_rg <- subset(pdamgff3_rg, id=="exon")
dim(exon_pdamgff_rg) # 12504 rows
unique(exon_pdamgff_rg$Gene.Predict)
# [1] "maker"
##five prime UTR
five_prime_UTR_pdamgff_rg <- subset(pdamgff3_rg, id=="five_prime_UTR")
dim(five_prime_UTR_pdamgff_rg) # 721 rows
unique(five_prime_UTR_pdamgff_rg$Gene.Predict)
# [1] "maker"
##CDS
CDS_pdamgff_rg <- subset(pdamgff3_rg, id=="CDS")
dim(CDS_pdamgff_rg) # 12037 rows
unique(CDS_pdamgff_rg$Gene.Predict)
# [1] "maker"
##three prime UTR
three_prime_UTR_pdamgff_rg <- subset(pdamgff3_rg, id=="three_prime_UTR")
dim(three_prime_UTR_pdamgff_rg) # 760 rows
unique(three_prime_UTR_pdamgff_rg$Gene.Predict)
# [1] "maker"
##match
match_pdamgff_rg <- subset(pdamgff3_rg, id=="match")
dim(match_pdamgff_rg) # 6587 rows
unique(match_pdamgff_rg$Gene.Predict)
# [1] "snap_masked"     "augustus_masked" "repeatmasker" 
snap_masked_match_pdamgff_rg<- subset(match_pdamgff_rg, Gene.Predict=="snap_masked")
dim(snap_masked_match_pdamgff_rg) # 2344 rows
augustus_masked_match_pdamgff_rg<- subset(match_pdamgff_rg, Gene.Predict=="augustus_masked")
dim(augustus_masked_match_pdamgff_rg) # 1318 rows
repeatmasker_match_pdamgff_rg<- subset(match_pdamgff_rg, Gene.Predict=="repeatmasker")
dim(repeatmasker_match_pdamgff_rg) # 2925 rows
##Match part
match_part_pdamgff_rg <- subset(pdamgff3_rg, id=="match_part")
dim(match_part_pdamgff_rg) # 617,859 rows
unique(match_part_pdamgff_rg$Gene.Predict)
# [1] "snap_masked"     "augustus_masked" "repeatmasker"    "blastn"          "est2genome"      "tblastx"         "cdna2genome"     "blastx"          "protein2genome" 
# [10] "repeatrunner"
snap_masked_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="snap_masked")
dim(snap_masked_match_part_pdamgff_rg) # 15451 rows
augustus_masked_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="augustus_masked")
dim(augustus_masked_match_part_pdamgff_rg) # 10203 rows
repeatmasker_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="repeatmasker")
dim(repeatmasker_match_part_pdamgff_rg) # 2925 rows
blastn_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="blastn")
dim(blastn_match_part_pdamgff_rg) # 25821 rows
est2genome_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="est2genome")
dim(est2genome_match_part_pdamgff_rg) # 22101 rows
tblastx_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="tblastx")
dim(tblastx_match_part_pdamgff_rg) # 4868 rows
cdna2genome_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="cdna2genome")
dim(cdna2genome_match_part_pdamgff_rg) # 1392 rows
blastx_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="blastx")
dim(blastx_match_part_pdamgff_rg) # 341866 rows
protein2genome_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="protein2genome")
dim(protein2genome_match_part_pdamgff_rg) # 193223 rows
repeatrunner_match_part_pdamgff_rg<- subset(match_part_pdamgff_rg, Gene.Predict=="repeatrunner")
dim(repeatrunner_match_part_pdamgff_rg) # 9 rows
##blank
blank_pdamgff_rg <- subset(pdamgff3_rg, id=="") # not sure what this stands for, just blank rows 
dim(blank_pdamgff_rg) # 500 rows
unique(blank_pdamgff_rg$Gene.Predict)
##Expressed sequence match 
expressed_sequence_match_pdamgff_rg <- subset(pdamgff3_rg, id=="expressed_sequence_match")  
dim(expressed_sequence_match_pdamgff_rg) # 7165 rows
unique(expressed_sequence_match_pdamgff_rg$Gene.Predict)
# [1] "blastn"      "est2genome"  "cdna2genome"
blastn_expressed_sequence_match_pdamgff_rg <- subset(expressed_sequence_match_pdamgff_rg, Gene.Predict=="blastn")
dim(blastn_expressed_sequence_match_pdamgff_rg) # 3395 rows
est2genome_expressed_sequence_match_pdamgff_rg <- subset(expressed_sequence_match_pdamgff_rg, Gene.Predict=="est2genome")
dim(est2genome_expressed_sequence_match_pdamgff_rg) # 3590 rows
cdna2genome_expressed_sequence_match_pdamgff_rg <- subset(expressed_sequence_match_pdamgff_rg, Gene.Predict=="cdna2genome")
dim(cdna2genome_expressed_sequence_match_pdamgff_rg) # 180 rows
##Translated nucleotide match
translated_nucleotide_match_pdamgff_rg <- subset(pdamgff3_rg, id=="translated_nucleotide_match")  
dim(translated_nucleotide_match_pdamgff_rg) # 181 rows
unique(translated_nucleotide_match_pdamgff_rg$Gene.Predict)
# [1] "tblastx"
##Protein match 
protein_match_pdamgff_rg <- subset(pdamgff3_rg, id=="protein_match")  
dim(protein_match_pdamgff_rg) # 173,524 rows
unique(protein_match_pdamgff_rg$Gene.Predict)
# [1] "blastx"         "protein2genome" "repeatrunner"
blastx_protein_match_pdamgff_rg <- subset(protein_match_pdamgff_rg, Gene.Predict=="blastx")
dim(blastx_protein_match_pdamgff_rg) # 81561 rows
protein2genome_protein_match_pdamgff_rg <- subset(protein_match_pdamgff_rg, Gene.Predict=="protein2genome")
dim(protein2genome_protein_match_pdamgff_rg) # 91954 rows
repeatrunner_protein_match_pdamgff_rg <- subset(protein_match_pdamgff_rg, Gene.Predict=="repeatrunner")
dim(repeatrunner_protein_match_pdamgff_rg) # 9 rows
##Contig
contig_pdamgff_rg <- subset(pdamgff3_rg, id=="contig")  
dim(contig_pdamgff_rg) # 336 rows
unique(contig_pdamgff_rg$Gene.Predict)
# [1] "."

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
# pdamgff3_NCBI <- na.omit(pdamgff3_NCBI)
# dim(pdamgff3_NCBI) # 507906 x 9
# length(unique(pdamgff3_NCBI$scaffold))

# Identify the unique parts -- gene predict 
unique(pdamgff3_NCBI$Gene.Predict)
# [1] ""            "RefSeq"      "Gnomon"      "tRNAscan-SE" "cmsearch"   
# Find out many of each unique part is in the annotation file
##Blank
G.P_blank_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="")  
dim(G.P_blank_pdamgff_NCBI) # 8787 rows
unique(G.P_blank_pdamgff_NCBI$id)
#[1] ""
##RefSeq
G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="RefSeq")  
dim(G.P_RefSeq_pdamgff_NCBI) # 23494 rows
unique(G.P_RefSeq_pdamgff_NCBI$id)
# [1] "region"     "cDNA_match" "tRNA"       "exon"       "rRNA"       "gene"       "CDS"        "intron"    
region_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="region")  
dim(region_G.P_RefSeq_pdamgff_NCBI) # 4393 rows
cDNA_match_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="cDNA_match")  
dim(cDNA_match_G.P_RefSeq_pdamgff_NCBI) # 19065 rows
tRNA_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="tRNA")  
dim(tRNA_G.P_RefSeq_pdamgff_NCBI) # 329 rows
exon_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="exon")  
dim(exon_G.P_RefSeq_pdamgff_NCBI) # 226892 rows
rRNA_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="rRNA")  
dim(rRNA_G.P_RefSeq_pdamgff_NCBI) # 3 rows
gene_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="gene")  
dim(gene_G.P_RefSeq_pdamgff_NCBI) # 21837 rows
CDS_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="CDS")  
dim(CDS_G.P_RefSeq_pdamgff_NCBI) # 206863 rows
intron_G.P_RefSeq_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="intron")  
dim(intron_G.P_RefSeq_pdamgff_NCBI) # 1 rows
##Gnomon
G.P_Gnomon_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="Gnomon")  
dim(G.P_Gnomon_pdamgff_NCBI) # 483246 rows
unique(G.P_Gnomon_pdamgff_NCBI$id)
# [1] "gene"       "mRNA"       "exon"       "CDS"        "pseudogene" "lnc_RNA"    "transcript"
gene_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="gene")  
dim(gene_G.P_Gnomon_pdamgff_NCBI) # 21446 rows
mRNA_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="mRNA")  
dim(mRNA_G.P_Gnomon_pdamgff_NCBI) # 25170 rows
exon_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="exon")  
dim(exon_G.P_Gnomon_pdamgff_NCBI) # 226478 rows
CDS_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="CDS")  
dim(CDS_G.P_Gnomon_pdamgff_NCBI) # 206849 rows
pseudogene_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="pseudogene")  
dim(pseudogene_G.P_Gnomon_pdamgff_NCBI) # 1237 rows
lnc_RNA_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="lnc_RNA")  
dim(lnc_RNA_G.P_Gnomon_pdamgff_NCBI) # 1771 rows
transcript_G.P_Gnomon_pdamgff_NCBI <- subset(G.P_Gnomon_pdamgff_NCBI, id=="transcript")  
dim(transcript_G.P_Gnomon_pdamgff_NCBI) # 295 rows
##tRNAscan-SE
G.P_tRNAscanSE_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="tRNAscan-SE")  
dim(G.P_tRNAscanSE_pdamgff_NCBI) # 1013 rows
unique(G.P_tRNAscanSE_pdamgff_NCBI$id)
# [1] "gene" "tRNA" "exon"
gene_G.P_tRNAscanSE_pdamgff_NCBI <- subset(G.P_tRNAscanSE_pdamgff_NCBI, id=="gene")  
dim(gene_G.P_tRNAscanSE_pdamgff_NCBI) # 327 rows
tRNA_G.P_tRNAscanSE_pdamgff_NCBI <- subset(G.P_tRNAscanSE_pdamgff_NCBI, id=="tRNA")  
dim(tRNA_G.P_tRNAscanSE_pdamgff_NCBI) # 327 rows
exon_G.P_tRNAscanSE_pdamgff_NCBI <- subset(G.P_tRNAscanSE_pdamgff_NCBI, id=="exon")  
dim(exon_G.P_tRNAscanSE_pdamgff_NCBI) # 359 rows
##cmsearch
G.P_cmsearch_pdamgff_NCBI <- subset(pdamgff3_NCBI, Gene.Predict=="cmsearch")  
dim(G.P_cmsearch_pdamgff_NCBI) # 153 rows
unique(G.P_cmsearch_pdamgff_NCBI$id)
# [1] "gene"      "snRNA"     "exon"      "snoRNA"    "guide_RNA" "rRNA"     
gene_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="gene")  
dim(gene_G.P_cmsearch_pdamgff_NCBI) # 51 rows
snRNA_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="snRNA")  
dim(snRNA_G.P_cmsearch_pdamgff_NCBI) # 25 rows
exon_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="exon")  
dim(exon_G.P_cmsearch_pdamgff_NCBI) # 51 rows
snoRNA_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="snoRNA")  
dim(snoRNA_G.P_cmsearch_pdamgff_NCBI) # 24 rows
guide_RNA_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="guide_RNA")  
dim(guide_RNA_G.P_cmsearch_pdamgff_NCBI) # 1 row
rRNA_G.P_cmsearch_pdamgff_NCBI <- subset(G.P_cmsearch_pdamgff_NCBI, id=="rRNA")  
dim(rRNA_G.P_cmsearch_pdamgff_NCBI) # 1 row

# Identify the unique parts -- ids
unique(pdamgff3_NCBI$id)
# [1] ""           "region"     "gene"       "mRNA"       "exon"       "CDS"        "pseudogene" "lnc_RNA"    "tRNA"       "cDNA_match" "transcript"
# [12] "snRNA"      "snoRNA"     "guide_RNA"  "rRNA"       "intron"
# Find out many of each unique part is in the annotation file
##Blank
blank_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="")  
dim(blank_pdamgff_NCBI) # 8787 rows
unique(blank_pdamgff_NCBI$Gene.Predict)
##region
region_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="region")  
dim(region_pdamgff_NCBI) # 4393 rows
unique(region_pdamgff_NCBI$Gene.Predict)
#[1] "RefSeq"
##gene
gene_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="gene")  
dim(gene_pdamgff_NCBI) # 21837 rows
unique(gene_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"      "tRNAscan-SE" "cmsearch"    "RefSeq"     
Gnomon_gene_pdamgff_NCBI <- subset(gene_pdamgff_NCBI, Gene.Predict=="Gnomon")  
dim(Gnomon_gene_pdamgff_NCBI) # 21446 rows
tRNAscanSE_gene_pdamgff_NCBI <- subset(gene_pdamgff_NCBI, Gene.Predict=="tRNAscan-SE")  
dim(tRNAscanSE_gene_pdamgff_NCBI) # 327 rows
cmsearch_gene_pdamgff_NCBI <- subset(gene_pdamgff_NCBI, Gene.Predict=="cmsearch")  
dim(cmsearch_gene_pdamgff_NCBI) # 51 rows
RefSeq_gene_pdamgff_NCBI <- subset(gene_pdamgff_NCBI, Gene.Predict=="RefSeq")  
dim(RefSeq_gene_pdamgff_NCBI) # 13 rows
##mRNA
mrna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="mRNA")  
dim(mrna_pdamgff_NCBI) # 25170 rows
unique(mrna_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"
##Exon
exon_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="exon")  
dim(exon_pdamgff_NCBI) # 226,892 rows
unique(exon_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"      "tRNAscan-SE" "cmsearch"    "RefSeq"     
Gnomon_exon_pdamgff_NCBI <- subset(exon_pdamgff_NCBI, Gene.Predict=="Gnomon")  
dim(Gnomon_exon_pdamgff_NCBI) # 226478 rows
tRNAscanSE_exon_pdamgff_NCBI <- subset(exon_pdamgff_NCBI, Gene.Predict=="tRNAscan-SE")  
dim(tRNAscanSE_exon_pdamgff_NCBI) # 359 rows
cmsearch_exon_pdamgff_NCBI <- subset(exon_pdamgff_NCBI, Gene.Predict=="cmsearch")  
dim(cmsearch_exon_pdamgff_NCBI) # 51 rows
RefSeq_exon_pdamgff_NCBI <- subset(exon_pdamgff_NCBI, Gene.Predict=="RefSeq")  
dim(RefSeq_exon_pdamgff_NCBI) # 4 rows
##CDS
CDS_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="CDS")  
dim(CDS_pdamgff_NCBI) # 206,863 rows
unique(CDS_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon" "RefSeq"
Gnomon_CDS_pdamgff_NCBI <- subset(CDS_pdamgff_NCBI, Gene.Predict=="Gnomon")  
dim(Gnomon_CDS_pdamgff_NCBI) # 206849 rows
RefSeq_CDS_pdamgff_NCBI <- subset(CDS_pdamgff_NCBI, Gene.Predict=="RefSeq")  
dim(RefSeq_CDS_pdamgff_NCBI) # 14 rows
##Pseudogene
pseudogene_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="pseudogene")  
dim(pseudogene_pdamgff_NCBI) # 1237 rows
unique(pseudogene_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"
##lncRNA
lnc_rna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="lnc_RNA")  
dim(lnc_rna_pdamgff_NCBI) # 1771 rows
unique(lnc_rna_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"
##tRNA
trna_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="tRNA")  
dim(trna_pdamgff_NCBI) # 329 rows
unique(trna_pdamgff_NCBI$Gene.Predict)
# [1] "tRNAscan-SE" "RefSeq"    
tRNAscanSE_tRNA_pdamgff_NCBI <- subset(trna_pdamgff_NCBI, Gene.Predict=="tRNAscan-SE")  
dim(tRNAscanSE_tRNA_pdamgff_NCBI) # 327 rows
RefSeq_tRNA_pdamgff_NCBI <- subset(trna_pdamgff_NCBI, Gene.Predict=="RefSeq")  
dim(RefSeq_tRNA_pdamgff_NCBI) # 2 rows
##cDNA_match 
cDNA_match_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="cDNA_match")  
dim(cDNA_match_pdamgff_NCBI) # 19065 rows
unique(cDNA_match_pdamgff_NCBI$Gene.Predict)
# [1] "RefSeq"
##transcript
transcript_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="transcript")  
dim(transcript_pdamgff_NCBI) # 295 rows
unique(transcript_pdamgff_NCBI$Gene.Predict)
# [1] "Gnomon"
##snRNA
snRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="snRNA")  
dim(snRNA_pdamgff_NCBI) # 25 rows
unique(snRNA_pdamgff_NCBI$Gene.Predict)
# [1] "cmsearch"
##snoRNA
snoRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="snoRNA")  
dim(snoRNA_pdamgff_NCBI) # 24 rows
unique(snoRNA_pdamgff_NCBI$Gene.Predict)
# [1] "cmsearch"
##Guide RNA
guide_RNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="guide_RNA")  
dim(guide_RNA_pdamgff_NCBI) # 1 row
unique(guide_RNA_pdamgff_NCBI$Gene.Predict)
# [1] "cmsearch"
##rRNA
rRNA_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="rRNA")  
dim(rRNA_pdamgff_NCBI) # 3 rows
unique(rRNA_pdamgff_NCBI$Gene.Predict)
# [1] "cmsearch" "RefSeq"  
cmsearch_rRNA_pdamgff_NCBI <- subset(rRNA_pdamgff_NCBI, Gene.Predict=="cmsearch")  
dim(cmsearch_rRNA_pdamgff_NCBI) # 1 rows
RefSeq_rRNA_pdamgff_NCBI <- subset(rRNA_pdamgff_NCBI, Gene.Predict=="RefSeq")  
dim(RefSeq_rRNA_pdamgff_NCBI) # 2 rows
##intron
intron_pdamgff_NCBI <- subset(pdamgff3_NCBI, id=="intron")  
dim(intron_pdamgff_NCBI) # 1 row
unique(intron_pdamgff_NCBI$Gene.Predict)
# [1] "RefSeq"

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

write.table(merge_NCBI, file="~/Desktop/GFFs/scaffold_merge_NCBI.txt", sep="\t", col.names = TRUE, row.names=FALSE, quote=FALSE)



