# Title: Evaluating shared GO terms
# Project: Sedimentation RNA-Seq
# Author: J. Ashey
# Date: 5/9/21


# should I be looking at the number of GO terms between species or the number of genes that have those GO terms between species?

# Load packages 
library("tidyverse")
library("readr")

## Read in GO sig enriched terms for each species

# Acerv
GO_acerv_sub <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Acervicornis_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
GO_acerv_sub <- select(GO_acerv_sub, c("category", "term", "ontology"))
length(unique(GO_acerv_sub$category)) # 15 GO terms

# Mcav
GO_mcav <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Mcavernosa_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
GO_mcav <- select(GO_mcav, c("category", "term", "ontology"))
length(unique(GO_mcav$category)) # 33 GO terms

# Pdam
GO_pdam <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Pdamicornis_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
GO_pdam <- select(GO_pdam, c("category", "term", "ontology"))
length(unique(GO_pdam$category)) # 50 GO terms

# Plob
GO_plob <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Plobata_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
GO_plob <- select(GO_plob, c("category", "term", "ontology"))
length(unique(GO_plob$category)) # 20 GO terms










# Merge and see what happens? Put larger df first 

# Acerv and Mcav
test_acerv_mcav <- inner_join(GO_mcav, GO_acerv_sub)

test <- full_join(GO_pdam, GO_acerv_sub, by = "term")






test_acerv_mcav <- merge(GO_mcav, GO_acerv_sub, by = "term", all.x = T)
unique(test_acerv_mcav$term.y) # 1 GO term in common: heme binding



test <- full_join(GO_mcav, GO_acerv_sub, by = "term")




# Acerv and Pdam
test_acerv_pdam <- merge(GO_pdam, GO_acerv_sub, by = "category", all = TRUE)



unique(test_acerv_pdam$term.y) # 1 GO term in common: heme binding
# Acerv and Plob
test_acerv_plob <- merge(GO_plob, GO_acerv_sub, by = "category", all.x = TRUE)
unique(test_acerv_plob$term.y) # 0 GO terms in common
# Mcav and Pdam
test_mcav_pdam <- merge(GO_pdam, GO_mcav, by = "category", all.x = TRUE)
unique(test_mcav_pdam$term.y) # 5 GO terms in common: heme binding, cell communication, membrane, metal ion binding, oxidation-reduction process
# Mcav and Plob
test_mcav_plob <- merge(GO_mcav, GO_plob, by = "category", all.x = TRUE)
unique(test_mcav_plob$term.y) # 0 GO terms in common
# Pdam and Plob
test_pdam_plob <- merge(GO_pdam, GO_plob, by = "category", all.x = TRUE)
unique(test_pdam_plob$term.y) # 9 GO terms in common: structural constituent of cytoskeleton, extracellular matrix structural constituent, calcium ion binding,
                                                    # protein binding, extracellular region, microtubule, microtubule-based process, chitin binding, cysteine-type peptidase activity
                     








