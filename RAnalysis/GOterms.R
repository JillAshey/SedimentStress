# GO terms for species 



##### Florida

## GO terms for acerv
acerv_IPS <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/InterProScan/acerv.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(acerv_IPS$V1)) # 981372
colnames(acerv_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
acerv_IPS_GO <- filter(acerv_IPS, grepl("GO:", attr)) # select only rows with GO terms
acerv_IPS_GO$GO_term <- regmatches(acerv_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", acerv_IPS_GO$attr, perl = TRUE))
acerv_IPS_GO$GO_term <- gsub(";.*", "", acerv_IPS_GO$GO_term)
acerv_IPS_GO <- select(acerv_IPS_GO, c(prot, Predict, GO_term))
write.csv(acerv_IPS_GO, file = "~/Desktop/acerv_GOterms.csv")

## GO terms for mcav
mcav_IPS <- read.csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/InterProScan/mcav.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(mcav_IPS$V1)) # 716241
colnames(mcav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
mcav_IPS_GO <- filter(mcav_IPS, grepl("GO:", attr)) # select only rows with GO terms
mcav_IPS_GO$GO_term <- regmatches(mcav_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", mcav_IPS_GO$attr, perl = TRUE))
mcav_IPS_GO$GO_term <- gsub(";.*", "", mcav_IPS_GO$GO_term)
mcav_IPS_GO <- select(mcav_IPS_GO, c(prot, Predict, GO_term))
write.csv(mcav_IPS_GO, file = "~/Desktop/mcav_GOterms.csv")

## GO terms for ofav
ofav_IPS <- read.csv("~/Desktop/Interproscan/ofav.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(ofav_IPS$V1)) # 1245302
colnames(ofav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
ofav_IPS_GO <- filter(ofav_IPS, grepl("GO:", attr)) # select only rows with GO terms
ofav_IPS_GO$GO_term <- regmatches(ofav_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", ofav_IPS_GO$attr, perl = TRUE))
ofav_IPS_GO$GO_term <- gsub(";.*", "", ofav_IPS_GO$GO_term)
ofav_IPS_GO <- select(ofav_IPS_GO, c(prot, Predict, GO_term))
write.csv(ofav_IPS_GO, file = "~/Desktop/ofav_GOterms.csv")

# to get LOC term for ofav, need to merge annotion gff file and go term file by XP

# ref annotation file 
ref <- read.table("~/Desktop/GFFs/GCF_002042975.1_ofav_dov_v1_genomic.gff", sep = "\t", header = FALSE)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")
ref_XP <- ref[grep("XP", ref$attr), ]
ref_XP$prot <- sub(";.*", "", ref_XP$attr)
ref_XP$prot <- gsub(".*-", "", ref_XP$prot)
# ofav GO terms
annot_GO <- read.csv("~/Desktop/ofav_GOterms.csv", header=TRUE)
annot_GO <- select(annot_GO, -X)
# merge go terms and ref file
annot_GO.ref <- merge(annot_GO, ref_XP, by = "prot", all.x= TRUE)
# make gene id col
annot_GO.ref$gene_id <- regmatches(annot_GO.ref$attr, gregexpr("(?<=gene=).*", annot_GO.ref$attr, perl = TRUE))
annot_GO.ref$gene_id <- sub(";.*", "", annot_GO.ref$gene_id)
# select go term and gene id cols
annot_GO.ref <- select(annot_GO.ref, c(GO_term, gene_id))
# save go terms with associated gene ids 
write.csv(annot_GO.ref, "~/Desktop/ofav_GOterms_geneID.csv")

                                  
                                  
                              
##### Hawaii 

# GO terms for mcap
mcap_IPS <- read.csv("~/Desktop/InterProScan/mcap.interpro.gff3", header = FALSE, sep="\t", skip=4)
colnames(mcap_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
length(unique(mcap_IPS$prot)) # 1547527 unique protein ids
mcap_IPS <- mcap_IPS[!grepl("#", mcap_IPS$prot),] # remove rows that have a # in scaffold col
mcap_IPS_GO <- filter(mcap_IPS, grepl("GO:", attr)) # select only rows with GO terms
mcap_IPS_GO$GO_term <- regmatches(mcap_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", mcap_IPS_GO$attr, perl = TRUE))
mcap_IPS_GO$GO_term <- gsub(";.*", "", mcap_IPS_GO$GO_term)
mcap_IPS_GO <- select(mcap_IPS_GO, c(prot, Predict, GO_term))
write.csv(mcap_IPS_GO, file = "~/Desktop/mcap_GOterms.csv")


# GO terms for pcomp
pcomp_IPS <- read.csv("~/Desktop/InterProScan/pcomp.interpro.gff3", header = FALSE, sep="\t", skip=4)
colnames(pcomp_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
length(unique(pcomp_IPS$prot)) # 2086032 unique protein ids
pcomp_IPS <- pcomp_IPS[!grepl("#", pcomp_IPS$prot),] # remove rows that have a # in scaffold col
pcomp_IPS_GO <- filter(pcomp_IPS, grepl("GO:", attr)) # select only rows with GO terms
pcomp_IPS_GO$GO_term <- regmatches(pcomp_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", pcomp_IPS_GO$attr, perl = TRUE))
pcomp_IPS_GO$GO_term <- gsub(";.*", "", pcomp_IPS_GO$GO_term)
pcomp_IPS_GO <- select(pcomp_IPS_GO, c(prot, Predict, GO_term))
write.csv(pcomp_IPS_GO, file = "~/Desktop/pcomp_GOterms.csv")


# GO terms for pdam
pdam_IPS <- read.csv("~/Desktop/InterProScan/pdam.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(pdam_IPS$V1)) # 1100156
colnames(pdam_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
pdam_IPS_GO <- filter(pdam_IPS, grepl("GO:", attr)) # select only rows with GO terms
pdam_IPS_GO$GO_term <- regmatches(pdam_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", pdam_IPS_GO$attr, perl = TRUE))
pdam_IPS_GO$GO_term <- gsub(";.*", "", pdam_IPS_GO$GO_term)
pdam_IPS_GO <- select(pdam_IPS_GO, c(prot, Predict, GO_term))
write.csv(pdam_IPS_GO, file = "~/Desktop/pdam_GOterms.csv")

# GO terms for plob
plob_IPS <- read.csv("~/Desktop/InterProScan/plob.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(plob_IPS$V1)) # 1192653
colnames(plob_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
plob_IPS_GO <- filter(plob_IPS, grepl("GO:", attr)) # select only rows with GO terms
plob_IPS_GO$GO_term <- regmatches(plob_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", plob_IPS_GO$attr, perl = TRUE))
plob_IPS_GO$GO_term <- gsub(";.*", "", plob_IPS_GO$GO_term)
plob_IPS_GO <- select(plob_IPS_GO, c(prot, Predict, GO_term))
write.csv(plob_IPS_GO, file = "~/Desktop/plob_GOterms.csv")







## For Pdam and Ofav, I need to use gff annotation files to isolate XP and match with LOC term 
## check if mRNA and prot ids are 1:1 on bluewaves

# Pdam
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
ref <- ref[grep("XP", ref$attr), ] # isolate XP protein name
ref$prot <- gsub(";.*", "", ref$attr) # remove everything after ;
ref$prot <- gsub("ID=", "", ref$prot) # remove ID=
ref$prot <- gsub(".*-", "", ref$prot) # remove everything after -
ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
ref$gene_id <- gsub(";.*", "", ref$gene_id) # remove everything after ;

pdam_IPS_GO.ref <- merge(ref, pdam_IPS_GO, by = "prot", all.x = TRUE) # merge ref and pdam_IPS_GO by prot
pdam_IPS_GO.ref <- unique(pdam_IPS_GO.ref) # select only unique ones
pdam_IPS_GO.ref <- na.omit(pdam_IPS_GO.ref) # remove rows with NAs in them
pdam_IPS_GO.ref <- select(pdam_IPS_GO.ref, c(prot, gene_id, GO_term)) # select prot, gene_id, GO_term
pdam_IPS_GO.ref <- unique(pdam_IPS_GO.ref) # select only unique ones 

# deg_pdam from DEG pdam csv
DEG_pdam <- read.csv("~/Desktop/pdam_unique.sig.list.csv", header = TRUE) # read in pdam DEGs
colnames(DEG_pdam)[1] <-"gene_id" # name 1st col gene_id
DEG_test <- merge(pdam_IPS_GO.ref, DEG_pdam, by = "gene_id") # merge pdam_IPS_GO.ref and DEG_pdam by gene_id
write.csv(DEG_test, file = "~/Desktop/pdam_GOterms_DEGs.csv") # write out file with DEG counts, gene_id, prot, and Go terms

  
  
  