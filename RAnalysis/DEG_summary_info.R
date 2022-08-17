

## DEG summary info 

acerv <- read_csv("Output/DESeq2/acerv/acerv_sub_DEGs.all_treatment_20210219.csv")
length(unique(acerv$gene_id)) # 215 unique DEGs
acerv_down <- acerv %>% filter(log2FoldChange < 0) 
length(count(acerv_down$gene_id)) # 154 unique DEGs downregulated
acerv_up <- acerv %>% filter(log2FoldChange > 0)
length(count(acerv_up$gene_id)) # 67 unique DEGs upregulated

mcav <- read_csv("Output/DESeq2/mcav/mcav_DEGs.all_treatment_20210208.csv")
length(unique(mcav$gene_id)) # 62 unique DEGs
mcav_down <- mcav %>% filter(log2FoldChange < 0) 
length(unique(mcav_down$gene_id)) # 44 unique DEGs downregulated
mcav_up <- mcav %>% filter(log2FoldChange > 0)
length(unique(mcav_up$gene_id)) # 19 unique DEGs upregulated

ofav <- read_csv("Output/DESeq2/ofav/ofav_DEGs_AllTreatments_20211104.csv")
length(unique(ofav$gene_id)) # 16 unique DEGs
ofav_down <- ofav %>% filter(log2FoldChange < 0) 
length(unique(ofav_down$gene_id)) # 3 unique DEGs downregulated
ofav_up <- ofav %>% filter(log2FoldChange > 0)
length(unique(ofav_up$gene_id)) # 13 unique DEGs upregulated


# pacuta

plob <- read_csv("Output/DESeq2/plob/plob_DEGs.all_treatment_20210326.csv")
length(unique(plob$gene_id)) # 153 unique DEGs
plob_down <- plob %>% filter(log2FoldChange < 0) 
length(unique(plob_down$gene_id)) # 128 unique DEGs downregulated
plob_up <- plob %>% filter(log2FoldChange > 0)
length(unique(plob_up$gene_id)) # 32 unique DEGs upregulated



## GO summary info 
acerv_go <- read_csv("Output/GOSeq/acerv/acerv_sub_ByTreatment_GO.terms_20211104.csv")
length(unique(acerv_go$gene_id)) # 83 unique DEGs w/ GO terms
length(unique(acerv_go$category)) # 366 unique GO terms

mcav_go <- read_csv("Output/GOSeq/mcav/mcav_ByTreatment_GO.terms_20211109.csv")
length(unique(mcav_go$gene_id)) # 35 unique DEGs w/ GO terms
length(unique(mcav_go$category)) # 220 unique GO terms

ofav_go <- read_csv("Output/GOSeq/ofav/ofav_ByTreatment_GO.terms_20211117.csv")
length(unique(ofav_go$gene_id)) # 7 unique DEGs w/ GO terms
length(unique(ofav_go$category)) # 20 unique GO terms

# pact

plob_go <- read_csv("Output/GOSeq/plob/plob_ByTreatment_GO.terms_20211106.csv")
length(unique(plob_go$gene_id)) # 83 unique DEGs w/ GO terms
length(unique(plob_go$category)) # 279 unique GO terms



