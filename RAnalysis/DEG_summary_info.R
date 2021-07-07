

## DEG summary info 

acerv <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/acerv/acerv_sub_DEGs.all_treatment_20210219.csv")
length(unique(acerv$gene_id)) # 215 unique DEGs
acerv_down <- acerv %>% filter(log2FoldChange < 0) 
length(count(acerv_down$gene_id)) # 154 unique DEGs downregulated
acerv_up <- acerv %>% filter(log2FoldChange > 0)
length(count(acerv_up$gene_id)) # 67 unique DEGs upregulated

mcav <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/mcav/mcav_DEGs.all_treatment_20210208.csv")
length(unique(mcav$gene_id)) # 62 unique DEGs
mcav_down <- mcav %>% filter(log2FoldChange < 0) 
length(unique(mcav_down$gene_id)) # 44 unique DEGs downregulated
mcav_up <- mcav %>% filter(log2FoldChange > 0)
length(unique(mcav_up$gene_id)) # 19 unique DEGs upregulated

pdam <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/pdam/pdam_DEGs.all_treatment_20210326.csv")
length(unique(pdam$gene_id)) # 549 unique DEGs
pdam_down <- pdam %>% filter(log2FoldChange < 0) 
length(unique(pdam_down$gene_id)) # 381 unique DEGs downregulated
pdam_up <- pdam %>% filter(log2FoldChange > 0)
length(unique(pdam_up$gene_id)) # 171 unique DEGs upregulated

plob <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/DESeq2/plob/plob_DEGs.all_treatment_20210326.csv")
length(unique(plob$gene_id)) # 153 unique DEGs
plob_down <- plob %>% filter(log2FoldChange < 0) 
length(unique(plob_down$gene_id)) # 128 unique DEGs downregulated
plob_up <- plob %>% filter(log2FoldChange > 0)
length(unique(plob_up$gene_id)) # 32 unique DEGs upregulated



## GO summary info 


acerv_go <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/acerv/acerv_sub_ByTreatment_GO.terms_20210327.csv")
length(unique(acerv_go$gene_id)) # 13 unique DEGs w/ GO terms
length(unique(acerv_go$category)) # 15 unique GO terms

mcav_go <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/mcav/mcav_ByTreatment_GO.terms_20210208.csv")
length(unique(mcav_go$gene_id)) # 22 unique DEGs w/ GO terms
length(unique(mcav_go$category)) # 33 unique GO terms

pdam_go <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/pdam/pdam_ByTreatment_GO.terms_20210508.csv")
length(unique(pdam_go$gene_id)) # 56 unique DEGs w/ GO terms
length(unique(pdam_go$category)) # 50 unique GO terms

plob_go <- read_csv("Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/plob/plob_ByTreatment_GO.terms_20210326.csv")
length(unique(plob_go$gene_id)) # 58 unique DEGs w/ GO terms
length(unique(plob_go$category)) # 20 unique GO terms



