
## Florida 

# Acerv
# acerv_list <- read.csv("~/Desktop//acerv_ByTreatment_GO.terms_20210124.csv", header = T)
# length(unique(acerv_list$gene_id)) # 9
# 
# acerv_GeneNames <- as.data.frame(unique(acerv_list$gene_id))
# colnames(acerv_GeneNames) <- "gene_id"
# acerv_GeneNames$gene_id <- gsub("TU", "model", acerv_GeneNames$gene_id)
# row.names(acerv_GeneNames) <- NULL
# 
# write.csv(acerv_GeneNames, "~/Desktop/acerv_GeneNames.csv")

# SUB ACERV
acerv_sub_list <- read.csv("~/Desktop/acerv_sub_ByTreatment_GO.terms_FullAnnot_20210222.csv", header = T)
acerv_sub_list <- na.omit(acerv_sub_list)
length(unique(acerv_sub_list$gene_id)) # 13

acerv_sub_GeneNames <- as.data.frame(unique(acerv_sub_list$gene_id))
colnames(acerv_sub_GeneNames) <- "gene_id"
#acerv_sub_GeneNames$gene_id <- gsub("TU", "model", acerv_sub_GeneNames$gene_id)
row.names(acerv_sub_GeneNames) <- NULL

write.csv(acerv_sub_GeneNames, "~/Desktop/acerv_sub_GeneNames.csv")

# Mcav
mcav_list <- read.csv("~/Desktop/mcav_ByTreatment_GO.terms_20210208.csv", header = T)
length(unique(mcav_list$gene_id)) # 22

mcav_GeneNames <- as.data.frame(unique(mcav_list$gene_id))
colnames(mcav_GeneNames) <- "gene_id"
row.names(mcav_GeneNames) <- NULL
write.csv(mcav_GeneNames, "~/Desktop/mcav_GeneNames.csv")



## Hawaii

# Pdam
pdam_list <- read.csv("~/Desktop/pdam_ByTreatment_GO.terms_20210208.csv", header = T)
length(unique(pdam_list$gene_id)) # 28
pdam_GeneNames <- as.data.frame(unique(pdam_list$gene_id))
colnames(pdam_GeneNames) <- "gene_id"
row.names(pdam_GeneNames) <- NULL
write.csv(pdam_GeneNames, "~/Desktop/pdam_GeneNames.csv")




# Plob
plob_list <- read.csv("~/Desktop/plob_ByTreatment_GO.terms_20210208.csv", header = T)
length(unique(plob_list$gene_id)) # 2916

jamg <- plob_list[c(1:102),]
jamg <- jamg$gene_id
jamg <- as.data.frame(jamg)
colnames(jamg) <- "gene_id"
jamg$gene_id <- gsub("TU", "model", jamg$gene_id)

plut <- plob_list[c(103:3038),]
plut <- plut$gene_id
plut <- as.data.frame(plut)
plut$plut <- paste0(plut$plut, ".m1")
colnames(plut) <- "gene_id"

plob_GeneNames <- rbind(jamg, plut)

write.csv(plob_GeneNames, "~/Desktop/plob_GeneNames.csv")
##### not working with extract_protein.py code 






