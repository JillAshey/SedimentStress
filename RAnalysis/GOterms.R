# GO terms for species 

# Breaking down annotations from Blast2GO, InterProScan, and Uniprot to get GO terms 



##### Florida

## GO terms for acerv

# InterProScan
acerv_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/acerv.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(acerv_IPS$V1)) # 981372
colnames(acerv_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
acerv_IPS_GO <- filter(acerv_IPS, grepl("GO:", attr)) # select only rows with GO terms
acerv_IPS_GO$GO_term <- regmatches(acerv_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", acerv_IPS_GO$attr, perl = TRUE))
acerv_IPS_GO$GO_term <- gsub(";.*", "", acerv_IPS_GO$GO_term)
acerv_IPS_GO <- select(acerv_IPS_GO, c(prot, GO_term))
colnames(acerv_IPS_GO) <- c("gene_id", "GO.ID")

# B2G + Uniprot
acerv_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/acerv_FuncAnn_UniP_B2G.csv", header = T)
length(unique(acerv_BU$gene_id)) # 29515
acerv_BU <- select(acerv_BU, c("gene_id", "GO.ID"))
acerv_BU$GO.ID <- gsub("F:", "", acerv_BU$GO.ID)
acerv_BU$GO.ID <- gsub("C:", "", acerv_BU$GO.ID)
acerv_BU$GO.ID <- gsub("P:", "", acerv_BU$GO.ID)
acerv_BU$GO.ID <- gsub(" ", "", acerv_BU$GO.ID)
acerv_BU_GO <- filter(acerv_BU, grepl("GO:", GO.ID)) # select only rows with GO terms

# Bind IPS and BU by row
acerv_GO <- rbind(acerv_BU_GO, acerv_IPS_GO)
acerv_GO <- unique(acerv_GO)
acerv_GO <- aggregate(acerv_GO$GO.ID, list(acerv_GO$gene_id), paste, collapse = ";")
colnames(acerv_GO) <- c("gene_id", "GO.ID")

# Save as a csv
write.csv(acerv_GO, file = "~/Desktop/acerv_GO_20210101.csv")



## GO terms for mcav

# InterProScan
mcav_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/mcav.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(mcav_IPS$V1)) # 716241
colnames(mcav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
mcav_IPS_GO <- filter(mcav_IPS, grepl("GO:", attr)) # select only rows with GO terms
mcav_IPS_GO$GO_term <- regmatches(mcav_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", mcav_IPS_GO$attr, perl = TRUE))
mcav_IPS_GO$GO_term <- gsub(";.*", "", mcav_IPS_GO$GO_term)
mcav_IPS_GO <- select(mcav_IPS_GO, c(prot, GO_term))
mcav_IPS_GO$GO_term <- gsub(",", ";", mcav_IPS_GO$GO_term)
colnames(mcav_IPS_GO) <- c("gene_id", "GO.ID")
dim(mcav_IPS_GO) # 55925 x 2

# B2G + Uniprot
mcav_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/mcav_FuncAnn_UniP_B2G.csv", header = T)
length(unique(mcav_BU$gene_id)) # 3191
mcav_BU <- select(mcav_BU, c("gene_id", "GO.ID"))
mcav_BU$GO.ID <- gsub("F:", "", mcav_BU$GO.ID)
mcav_BU$GO.ID <- gsub("C:", "", mcav_BU$GO.ID)
mcav_BU$GO.ID <- gsub("P:", "", mcav_BU$GO.ID)
mcav_BU$GO.ID <- gsub(" ", "", mcav_BU$GO.ID)
mcav_BU_GO <- filter(mcav_BU, grepl("GO:", GO.ID)) # select only rows with GO terms
dim(mcav_BU_GO) # 147 x 2

# Bind IPS and BU by row
mcav_GO <- rbind(mcav_BU_GO, mcav_IPS_GO)
mcav_GO <- unique(mcav_GO)
mcav_GO <- aggregate(mcav_GO$GO.ID, list(mcav_GO$gene_id), paste, collapse = ";")
colnames(mcav_GO) <- c("gene_id", "GO.ID")

# Save as a csv
write.csv(mcav_GO, file = "~/Desktop/mcav_GO_20210101.csv")



## GO terms for ofav

# InterProScan
ofav_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/ofav.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(ofav_IPS$V1)) # 1245302
colnames(ofav_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
ofav_IPS_GO <- filter(ofav_IPS, grepl("GO:", attr)) # select only rows with GO terms
ofav_IPS_GO$GO_term <- regmatches(ofav_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", ofav_IPS_GO$attr, perl = TRUE))
ofav_IPS_GO$GO_term <- gsub(";.*", "", ofav_IPS_GO$GO_term)
ofav_IPS_GO <- select(ofav_IPS_GO, c(prot, GO_term))
colnames(ofav_IPS_GO) <- c("gene_id", "GO.ID")

# B2G + Uniprot
ofav_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/ofav_FuncAnn_UniP_B2G.csv", header = T)
length(unique(ofav_BU$gene_id)) # 3191
ofav_BU <- select(ofav_BU, c("gene_id", "GO.ID"))
ofav_BU_GO <- filter(ofav_BU, grepl("GO:", GO.ID)) # select only rows with GO terms

# For gene_ids, I have ones that start with XP, XM, and XR. Need to look at reference annotation file to find a common gene id, which will be the LOC term in attr col
ref <- read.table("~/Desktop/GFFs/GCF_002042975.1_ofav_dov_v1_genomic.gff", sep = "\t", header = FALSE)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")

# Make separate df with only rows with XP 
ref_XP <- ref[grep("XP", ref$attr), ]
ref_XP$gene_id <- sub(";.*", "", ref_XP$attr)
ref_XP$gene_id <- gsub(".*-", "", ref_XP$gene_id)

# Merge ref_XP with ofav_BU_GO
ofav_XP <- merge(ref_XP, ofav_BU_GO, by = "gene_id")
# no matches 

# Merge ref_XP with ofav_IPS_GO
ofav_XP2 <- merge(ofav_IPS_GO, ref_XP, by = "gene_id")
ofav_XP2 <- unique(ofav_XP2)

# Isolate LOC terms 
ofav_XP2$LOC <- regmatches(ofav_XP2$attr, gregexpr("(?<=gene=).*", ofav_XP2$attr, perl = TRUE))
ofav_XP2$LOC <- sub(";.*", "", ofav_XP2$LOC)

# Subset LOC and GO terms
ofav_XP2 <- select(ofav_XP2, c("GO.ID", "LOC"))

# Make separate df with only rows with XM
ref_XM <- ref[grep("XM", ref$attr), ]
ref_XM$gene_id <- regmatches(ref_XM$attr, gregexpr("(?>XM).*", ref_XM$attr, perl = TRUE))
ref_XM$gene_id <- sub(";.*", "", ref_XM$gene_id)
ref_XM$gene_id <- sub("-.*", "", ref_XM$gene_id)

# Merge ref_XM with ofav_BU_GO
ofav_XM <- merge(ofav_BU_GO, ref_XM, by = "gene_id")
# 710 matches 

# Isolate LOC terms 
ofav_XM$LOC <- regmatches(ofav_XM$attr, gregexpr("(?<=gene=).*", ofav_XM$attr, perl = TRUE))
ofav_XM$LOC <- sub(";.*", "", ofav_XM$LOC)

# Subset LOC and GO terms
ofav_XM <- select(ofav_XM, c("GO.ID", "LOC"))

# Make separate df with only rows with XR
ref_XR <- ref[grep("XR", ref$attr), ]
ref_XR$gene_id <- regmatches(ref_XR$attr, gregexpr("(?>XR).*", ref_XR$attr, perl = TRUE))
ref_XR$gene_id <- sub(";.*", "", ref_XR$gene_id)
ref_XR$gene_id <- sub("-.*", "", ref_XR$gene_id)

# Merge ref_XR with ofav_BU_GO
ofav_XR <- merge(ofav_BU_GO, ref_XR, by = "gene_id")
# 298 matches 

# Isolate LOC terms 
ofav_XR$LOC <- regmatches(ofav_XR$attr, gregexpr("(?<=gene=).*", ofav_XR$attr, perl = TRUE))
ofav_XR$LOC <- sub(";.*", "", ofav_XR$LOC)

# Subset LOC and GO terms
ofav_XR <- select(ofav_XR, c("GO.ID", "LOC"))

# Bind ofav_XP2, ofav_XM, and ofav_XR together 
ofav_GO <- rbind(ofav_XP2, ofav_XM, ofav_XR)
ofav_GO$GO.ID <- gsub("F:", "", ofav_GO$GO.ID)
ofav_GO$GO.ID <- gsub("C:", "", ofav_GO$GO.ID)
ofav_GO$GO.ID <- gsub("P:", "", ofav_GO$GO.ID)
ofav_GO$GO.ID <- gsub(" ", "", ofav_GO$GO.ID)
ofav_GO <- unique(ofav_GO)
dim(ofav_GO) # 10246 x 2

# Rename and reorder columns 
colnames(ofav_GO)[2] <- "gene_id"
ofav_GO <- ofav_GO[c("gene_id", "GO.ID")]

# Save as a csv
write.csv(ofav_GO, file = "~/Desktop/ofav_GO_20210101.csv")









                                  
                              
##### Hawaii 

# GO terms for mcap



# GO terms for pcomp
# InterProScan
pcomp_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/pcomp.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(pcomp_IPS$V1)) # 2086032
colnames(pcomp_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
pcomp_IPS_GO <- filter(pcomp_IPS, grepl("GO:", attr)) # select only rows with GO terms
pcomp_IPS_GO$GO_term <- regmatches(pcomp_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", pcomp_IPS_GO$attr, perl = TRUE))
pcomp_IPS_GO$GO_term <- gsub(";.*", "", pcomp_IPS_GO$GO_term)
pcomp_IPS_GO <- select(pcomp_IPS_GO, c(prot, GO_term))
colnames(pcomp_IPS_GO) <- c("gene_id", "GO.ID")

# B2G + Uniprot
pcomp_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/pcomp_FuncAnn_UniP_B2G.csv", header = T)
length(unique(pcomp_BU$gene_id)) # 59049
pcomp_BU <- select(pcomp_BU, c("gene_id", "GO.ID"))
pcomp_BU$GO.ID <- gsub("F:", "", pcomp_BU$GO.ID)
pcomp_BU$GO.ID <- gsub("C:", "", pcomp_BU$GO.ID)
pcomp_BU$GO.ID <- gsub("P:", "", pcomp_BU$GO.ID)
pcomp_BU$GO.ID <- gsub(" ", "", pcomp_BU$GO.ID)
pcomp_BU_GO <- filter(pcomp_BU, grepl("GO:", GO.ID)) # select only rows with GO terms

# Bind IPS and BU by row
pcomp_GO <- rbind(pcomp_BU_GO, pcomp_IPS_GO)
pcomp_GO <- unique(pcomp_GO)
pcomp_GO <- aggregate(pcomp_GO$GO.ID, list(pcomp_GO$gene_id), paste, collapse = ";")
colnames(pcomp_GO) <- c("gene_id", "GO.ID")

# Save as a csv
write.csv(pcomp_GO, file = "~/Desktop/pcomp_GO_20210101.csv")



# GO terms for pcomp
# InterProScan
pdam_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/pdam.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(pdam_IPS$V1)) # 1100156
colnames(pdam_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
pdam_IPS_GO <- filter(pdam_IPS, grepl("GO:", attr)) # select only rows with GO terms
pdam_IPS_GO$GO_term <- regmatches(pdam_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", pdam_IPS_GO$attr, perl = TRUE))
pdam_IPS_GO$GO_term <- gsub(";.*", "", pdam_IPS_GO$GO_term)
pdam_IPS_GO <- select(pdam_IPS_GO, c(prot, GO_term))
colnames(pdam_IPS_GO) <- c("gene_id", "GO.ID")

# B2G + Uniprot
pdam_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/pdam_FuncAnn_UniP_B2G.csv", header = T)
length(unique(pdam_BU$gene_id)) # 3406
pdam_BU <- select(pdam_BU, c("gene_id", "GO.ID"))
pdam_BU$GO.ID <- gsub("F:", "", pdam_BU$GO.ID)
pdam_BU$GO.ID <- gsub("C:", "", pdam_BU$GO.ID)
pdam_BU$GO.ID <- gsub("P:", "", pdam_BU$GO.ID)
pdam_BU$GO.ID <- gsub(" ", "", pdam_BU$GO.ID)
pdam_BU_GO <- filter(pdam_BU, grepl("GO:", GO.ID)) # select only rows with GO terms

# For gene_ids, I have ones that start with XP, XM, and XR. Need to look at reference annotation file to find a common gene id, which will be the LOC term in attr col
ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr")

# Make separate df with only rows with XP 
ref_XP <- ref[grep("XP", ref$attr), ]
ref_XP$gene_id <- sub(";.*", "", ref_XP$attr)
ref_XP$gene_id <- gsub(".*-", "", ref_XP$gene_id)

# Merge ref_XP with pdam_BU_GO
pdam_XP <- merge(ref_XP, pdam_BU_GO, by = "gene_id")
# no matches 

# Merge ref_XP with ofav_IPS_GO
pdam_XP2 <- merge(pdam_IPS_GO, ref_XP, by = "gene_id")
pdam_XP2 <- unique(pdam_XP2)

# Isolate LOC terms 
pdam_XP2$LOC <- regmatches(pdam_XP2$attr, gregexpr("(?<=gene=).*", pdam_XP2$attr, perl = TRUE))
pdam_XP2$LOC <- sub(";.*", "", pdam_XP2$LOC)

# Subset LOC and GO terms
pdam_XP2 <- select(pdam_XP2, c("GO.ID", "LOC"))

# Make separate df with only rows with XM
ref_XM <- ref[grep("XM", ref$attr), ]
ref_XM$gene_id <- regmatches(ref_XM$attr, gregexpr("(?>XM).*", ref_XM$attr, perl = TRUE))
ref_XM$gene_id <- sub(";.*", "", ref_XM$gene_id)
ref_XM$gene_id <- sub("-.*", "", ref_XM$gene_id)

# Merge ref_XM with ofav_BU_GO
pdam_XM <- merge(pdam_BU_GO, ref_XM, by = "gene_id")
# 710 matches 

# Isolate LOC terms 
pdam_XM$LOC <- regmatches(pdam_XM$attr, gregexpr("(?<=gene=).*", pdam_XM$attr, perl = TRUE))
pdam_XM$LOC <- sub(";.*", "", pdam_XM$LOC)

# Subset LOC and GO terms
pdam_XM <- select(pdam_XM, c("GO.ID", "LOC"))

# Make separate df with only rows with XR
ref_XR <- ref[grep("XR", ref$attr), ]
ref_XR$gene_id <- regmatches(ref_XR$attr, gregexpr("(?>XR).*", ref_XR$attr, perl = TRUE))
ref_XR$gene_id <- sub(";.*", "", ref_XR$gene_id)
ref_XR$gene_id <- sub("-.*", "", ref_XR$gene_id)

# Merge ref_XR with ofav_BU_GO
pdam_XR <- merge(pdam_BU_GO, ref_XR, by = "gene_id")
# 1088 matches 

# Isolate LOC terms 
pdam_XR$LOC <- regmatches(pdam_XR$attr, gregexpr("(?<=gene=).*", pdam_XR$attr, perl = TRUE))
pdam_XR$LOC <- sub(";.*", "", pdam_XR$LOC)

# Subset LOC and GO terms
pdam_XR <- select(pdam_XR, c("GO.ID", "LOC"))

# Bind pdam_XP2, pdam_XM, and pdam_XR together 
pdam_GO <- rbind(pdam_XP2, pdam_XM, pdam_XR)
pdam_GO <- aggregate(pdam_GO$GO.ID, list(pdam_GO$gene_id), paste, collapse = ";")
pdam_GO <- unique(pdam_GO)

# Rename and reorder columns 
colnames(pdam_GO) <- c("gene_id", "GO.ID")
pdam_GO <- pdam_GO[c("gene_id", "GO.ID")]

# Save as a csv
write.csv(pdam_GO, file = "~/Desktop/pdam_GO_20210101.csv")



## GO terms for plob
# InterProScan
plob_IPS <- read.csv("~/Desktop/PutnamLab/FunctionalAnnotation/InterProScan/plob.interpro.gff3", header = FALSE, sep="\t", skip=4)
length(unique(plob_IPS$V1)) # 1192653
colnames(plob_IPS) <- c("prot", "Predict", "id", "start","stop", "pos1", "pos2","pos3", "attr")
plob_IPS_GO <- filter(plob_IPS, grepl("GO:", attr)) # select only rows with GO terms
plob_IPS_GO$GO_term <- regmatches(plob_IPS_GO$attr, gregexpr("(?<=Ontology_term=).*", plob_IPS_GO$attr, perl = TRUE))
plob_IPS_GO$GO_term <- gsub(";.*", "", plob_IPS_GO$GO_term)
plob_IPS_GO <- select(plob_IPS_GO, c(prot, GO_term))
plob_IPS_GO$GO_term <- gsub(",", ";", plob_IPS_GO$GO_term)
colnames(plob_IPS_GO) <- c("gene_id", "GO.ID")
dim(plob_IPS_GO) # 112162 x 2

# B2G + Uniprot
plob_BU <- read.csv("~/Desktop/PutnamLab/Repositories/FunctionalAnnotation/FunctionalAnnotation/final_annotations/plob_FuncAnn_UniP_B2G.csv", header = T)
length(unique(plob_BU$gene_id)) # 3191
plob_BU <- select(plob_BU, c("gene_id", "GO.ID"))
plob_BU$GO.ID <- gsub("F:", "", plob_BU$GO.ID)
plob_BU$GO.ID <- gsub("C:", "", plob_BU$GO.ID)
plob_BU$GO.ID <- gsub("P:", "", plob_BU$GO.ID)
plob_BU$GO.ID <- gsub(" ", "", plob_BU$GO.ID)
plob_BU_GO <- filter(plob_BU, grepl("GO:", GO.ID)) # select only rows with GO terms
dim(plob_BU_GO) # 147 x 2

# Bind IPS and BU by row
plob_GO <- rbind(plob_BU_GO, plob_IPS_GO)
plob_GO <- unique(plob_GO)
plob_GO <- aggregate(plob_GO$GO.ID, list(plob_GO$gene_id), paste, collapse = ";")
colnames(plob_GO) <- c("gene_id", "GO.ID")

# Save as a csv
write.csv(plob_GO, file = "~/Desktop/plob_GO_20210101.csv")










## For Pdam and Ofav, I need to use gff annotation files to isolate XP and match with LOC term 
## check if mRNA and prot ids are 1:1 on bluewaves

# Pdam
# ref <- read.csv("~/Desktop/GFFs/GCF_003704095.1_ASM370409v1_genomic.gff",header = FALSE, sep="\t", skip=6)
# colnames(ref) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "attr") # name cols
# ref <- ref[!grepl("##", ref$scaffold),] # remove rows that have a # in scaffold col
# ref <- ref[grep("XP", ref$attr), ] # isolate XP protein name
# ref$prot <- gsub(";.*", "", ref$attr) # remove everything after ;
# ref$prot <- gsub("ID=", "", ref$prot) # remove ID=
# ref$prot <- gsub(".*-", "", ref$prot) # remove everything after -
# ref$gene_id <- regmatches(ref$attr, gregexpr("(?<=gene=).*", ref$attr, perl = TRUE)) #removing everything in Symbol col up to LOC and creating new col called gene_id
# ref$gene_id <- gsub(";.*", "", ref$gene_id) # remove everything after ;
# 
# pdam_IPS_GO.ref <- merge(ref, pdam_IPS_GO, by = "prot", all.x = TRUE) # merge ref and pdam_IPS_GO by prot
# pdam_IPS_GO.ref <- unique(pdam_IPS_GO.ref) # select only unique ones
# pdam_IPS_GO.ref <- na.omit(pdam_IPS_GO.ref) # remove rows with NAs in them
# pdam_IPS_GO.ref <- select(pdam_IPS_GO.ref, c(prot, gene_id, GO_term)) # select prot, gene_id, GO_term
# pdam_IPS_GO.ref <- unique(pdam_IPS_GO.ref) # select only unique ones 
# 
# # deg_pdam from DEG pdam csv
# DEG_pdam <- read.csv("~/Desktop/pdam_unique.sig.list.csv", header = TRUE) # read in pdam DEGs
# colnames(DEG_pdam)[1] <-"gene_id" # name 1st col gene_id
# DEG_test <- merge(pdam_IPS_GO.ref, DEG_pdam, by = "gene_id") # merge pdam_IPS_GO.ref and DEG_pdam by gene_id
# write.csv(DEG_test, file = "~/Desktop/pdam_GOterms_DEGs.csv") # write out file with DEG counts, gene_id, prot, and Go terms

  
  
  