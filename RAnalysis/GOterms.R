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

                                  
                                  
                              
##### Hawaii 

# GO terms for mcap


# GO terms for pcomp


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










