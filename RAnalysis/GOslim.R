## GO slim testing for sediment stress corals
# Sediment Stress
# 20210328
# J.Ashey


# Read in Packages
library(readr)
library(arsenal)
library(tidyverse)


#######Hawaii

## Pdam
# Read in pdam info
pdam.go <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Pdamicornis_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
pdam.go <- select(pdam.go, c("category", "term", "ontology", "over_represented_pvalue", "Treatment_Compare"))
colnames(pdam.go)[1] <-"GO.IDs" # rename column with GO terms to GO.IDs
length(unique(pdam.go$term)) # 50 GO terms


# Read in GOslim info
go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims.pdam <- merge(pdam.go, go.slim, by="GO.IDs", all = TRUE) # merge pdam info and GOslim
Gene.GO.IDs.slims.pdam <- na.omit(Gene.GO.IDs.slims.pdam)

# Plot pdam goslim
GO.pdam <- ggplot(Gene.GO.IDs.slims.pdam, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.pdam
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/pdam/pdam_GOslim_20210508.pdf", GO.pdam, width = 28, height = 28, units = c("in"))

## Plob
plob.go <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Plobata_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
plob.go <- select(plob.go, c("category", "term", "ontology", "over_represented_pvalue", "Treatment_Compare"))
colnames(plob.go)[1] <-"GO.IDs"
#plob.go <- plob.go[!plob.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims.plob <- merge(plob.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims.plob <- na.omit(Gene.GO.IDs.slims.plob)

GO.plob <- ggplot(Gene.GO.IDs.slims.plob, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.plob
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/plob/plob_GOslim_20210328.pdf", GO.plob, width = 28, height = 28, units = c("in"))




#######Florida

## Acerv
acerv.go <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Acervicornis_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
acerv.go <- select(acerv.go, c("category", "term", "ontology", "over_represented_pvalue", "Treatment_Compare"))
colnames(acerv.go)[1] <-"GO.IDs"
#plob.go <- plob.go[!plob.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims.acerv <- merge(acerv.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims.acerv <- na.omit(Gene.GO.IDs.slims.acerv)

GO.acerv <- ggplot(Gene.GO.IDs.slims.acerv, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.acerv
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/acerv/acerv_GOslim_20210328.pdf", GO.acerv, width = 28, height = 28, units = c("in"))


## Mcav
mcav.go <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/SuppTables/Mcavernosa_DEG_SuppTable.csv", locale = locale(encoding = "Latin1"))
mcav.go <- select(mcav.go, c("category", "term", "ontology", "over_represented_pvalue", "Treatment_Compare"))
colnames(mcav.go)[1] <-"GO.IDs"
#plob.go <- plob.go[!plob.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims.mcav <- merge(mcav.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims.mcav <- na.omit(Gene.GO.IDs.slims.mcav)

GO.mcav <- ggplot(Gene.GO.IDs.slims.mcav, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.mcav
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/mcav/mcav_GOslim_20210328.pdf", GO.mcav, width = 28, height = 28, units = c("in"))





### Compare GO slim terms and GO terms between species 
# Mcav & Acerv
MA_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.mcav$GO.Slim.Term, Gene.GO.IDs.slims.acerv$GO.Slim.Term) # "other metabolic processes",  "other molecular function"
MA_term <- intersect(Gene.GO.IDs.slims.mcav$term, Gene.GO.IDs.slims.acerv$term) # heme binding 

# Mcav & Pdam
MPd_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.mcav$GO.Slim.Term, Gene.GO.IDs.slims.pdam$GO.Slim.Term) # "kinase activity", "other biological processes", "other cellular component", "other membranes"
                          # "other metabolic processes",  "other molecular function",   "signal transduction", "transport", "transporter activity"
MPd_term <- intersect(Gene.GO.IDs.slims.mcav$term, Gene.GO.IDs.slims.pdam$term) # cell communication

# Mcav & Plob
MPl_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.mcav$GO.Slim.Term, Gene.GO.IDs.slims.plob$GO.Slim.Term) # "other biological processes", "other metabolic processes", "other molecular function", "signal transduction activity", "transporter activity"
MPl_term <- intersect(Gene.GO.IDs.slims.mcav$term, Gene.GO.IDs.slims.plob$term) # NA

# Plob & Pdam
PdPl_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.plob$GO.Slim.Term, Gene.GO.IDs.slims.pdam$GO.Slim.Term) # "cytoskeletal activity", "other biological processes", "other metabolic processes", "other molecular function", "transporter activity"
PdPl_term <- intersect(Gene.GO.IDs.slims.plob$term, Gene.GO.IDs.slims.pdam$term) # NA

# Plob & Acerv
APl_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.plob$GO.Slim.Term, Gene.GO.IDs.slims.acerv$GO.Slim.Term) # "other metabolic processes", "other molecular function"
APl_GO.term <- intersect(Gene.GO.IDs.slims.plob$term, Gene.GO.IDs.slims.acerv$term) # NA

# Pdam & Acerv
APd_GO.Slim.Term <- intersect(Gene.GO.IDs.slims.pdam$GO.Slim.Term, Gene.GO.IDs.slims.acerv$GO.Slim.Term) # "cell organization and biogenesis", "DNA metabolism", "other metabolic processes", "other molecular function", "protein metabolism", "RNA metabolism"
APd_GO.term <- intersect(Gene.GO.IDs.slims.pdam$term, Gene.GO.IDs.slims.acerv$term) # NA


### Plot all species on one GO slim plot 
# Add species name and site to each df
Gene.GO.IDs.slims.pdam$Species <- "P.damicornis"
Gene.GO.IDs.slims.pdam$Site <- "Florida"
Gene.GO.IDs.slims.plob$Species <- "P.lobata"
Gene.GO.IDs.slims.plob$Site <- "Florida"
Gene.GO.IDs.slims.acerv$Species <- "A.cervicornis"
Gene.GO.IDs.slims.acerv$Site <- "Hawaii"
Gene.GO.IDs.slims.mcav$Species <- "M.cavernosa"
Gene.GO.IDs.slims.mcav$Site <- "Hawaii"

# Bind all species df together for GO slim plot 
all.go.slim <- rbind(Gene.GO.IDs.slims.acerv, Gene.GO.IDs.slims.mcav, Gene.GO.IDs.slims.pdam, Gene.GO.IDs.slims.plob)

GO.all <- ggplot(all.go.slim, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.all
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/all_GOslim_20210509.pdf", GO.all, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/all_GOslim_20210509.jpeg", GO.all, width = 28, height = 28, units = c("in"))

# Subset and plot by MF
go.slim.MF <- all.go.slim %>% filter(ontology=="MF")

GO.subset.MF <- ggplot(go.slim.MF, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.subset.MF
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/MF_GOslim_20210705.pdf", GO.subset.MF, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/MF_GOslim_20210705.jpeg", GO.subset.MF, width = 28, height = 28, units = c("in"))

# Subset and plot by BP
go.slim.BP <- all.go.slim %>% filter(ontology=="BP")

GO.subset.BP <- ggplot(go.slim.BP, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.subset.BP
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/BP_GOslim_20210705.pdf", GO.subset.BP, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/BP_GOslim_20210705.jpeg", GO.subset.BP, width = 28, height = 28, units = c("in"))

# Subset and plot only 'other molecular function' & 'other metabolic processes'
go.slim.subset <- all.go.slim %>% filter(GO.Slim.Term=="other molecular function" | GO.Slim.Term=="other metabolic processes")

GO.subset <- ggplot(go.slim.subset, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.subset
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/MolecMetab.Subset_GOslim_20210509.pdf", GO.subset, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/MolecMetab.Subset_GOslim_20210509.jpeg", GO.subset, width = 28, height = 28, units = c("in"))

# Subset and plot by site 
## FL
go.FL.subset <- all.go.slim %>% filter(Site=="Florida")

GO.FL_subset <- ggplot(go.FL.subset, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.FL_subset
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/FL_GOslim_20210509.pdf", GO.FL_subset, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/FL_GOslim_20210509.jpeg", GO.FL_subset, width = 28, height = 28, units = c("in"))

## HI
go.HI.subset <- all.go.slim %>% filter(Site=="Hawaii")

GO.HI_subset <- ggplot(go.HI.subset, aes(x = Species, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO.HI_subset
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/HI_GOslim_20210509.pdf", GO.HI_subset, width = 28, height = 28, units = c("in"))
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/HI_GOslim_20210509.jpeg", GO.HI_subset, width = 28, height = 28, units = c("in"))














