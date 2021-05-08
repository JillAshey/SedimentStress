## GO slim testing for sediment stress corals
# Sediment Stress
# 20210328
# J.Ashey


# Read in Packages
library(readr)

## Pdam
# Read in pdam info
pdam.go <- read.csv("~/Desktop/pdam_ByTreatment_GO.terms_20210508.csv", header = T)
pdam.go <- pdam.go[, -1] # remove first column
colnames(pdam.go)[2] <-"GO.IDs" # rename column with GO terms to GO.IDs

# Read in GOslim info
go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims <- merge(pdam.go, go.slim, by="GO.IDs", all = TRUE) # merge pdam info and GOslim
Gene.GO.IDs.slims <- na.omit(Gene.GO.IDs.slims)


GO <- ggplot(Gene.GO.IDs.slims, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/pdam_GOslim_20210508.pdf", GO, width = 28, height = 28, units = c("in"))

## Plob
plob.go <- read.csv("~/Desktop/plob_ByTreatment_GO.terms_20210326.csv", header = T)
colnames(plob.go)[22] <-"GO.IDs"
plob.go <- plob.go[!plob.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims <- merge(plob.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims <- Gene.GO.IDs.slims[,-2]
Gene.GO.IDs.slims <- na.omit(Gene.GO.IDs.slims)


GO <- ggplot(Gene.GO.IDs.slims, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/plob_GOslim_20210328.pdf", GO, width = 28, height = 28, units = c("in"))



## Acerv
acerv.go <- read.csv("~/Desktop/acerv_sub_ByTreatment_GO.terms_20210327.csv", header = T)
colnames(acerv.go)[23] <-"GO.IDs"
acerv.go <- acerv.go[!acerv.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims <- merge(acerv.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims <- Gene.GO.IDs.slims[,-2]
Gene.GO.IDs.slims <- na.omit(Gene.GO.IDs.slims)


GO <- ggplot(Gene.GO.IDs.slims, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/acerv_GOslim_20210328.pdf", GO, width = 28, height = 28, units = c("in"))



## Mcav
mcav.go <- read.csv("~/Desktop/mcav_ByTreatment_GO.terms_20210208.csv", header = T)
colnames(mcav.go)[30] <-"GO.IDs"
mcav.go <- mcav.go[!mcav.go$ontology=='CC',]

go.slim <- read_csv("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/GOSeq/GO-GOslim.csv")
colnames(go.slim) <- c("GO.IDs", "GO.Term", "GO.Slim.Term", "Cat") #rename columns
Gene.GO.IDs.slims <- merge(mcav.go, go.slim, by="GO.IDs", all = TRUE)
Gene.GO.IDs.slims <- Gene.GO.IDs.slims[,-2]
Gene.GO.IDs.slims <- na.omit(Gene.GO.IDs.slims)


GO <- ggplot(Gene.GO.IDs.slims, aes(x = ontology, y = term)) + 
  geom_tile(aes(fill =over_represented_pvalue)) + 
  facet_grid(GO.Slim.Term ~ ., scales = "free_y", labeller = label_wrap_gen(width = 5, multi_line = TRUE))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     strip.text.y = element_text(angle=0, size = 25),
                     strip.text.x = element_text(size = 25),
                     axis.text = element_text(size = 15))
GO
ggsave("~/Desktop/PutnamLab/Repositories/SedimentStress/SedimentStress/Output/Figs/mcav_GOslim_20210328.pdf", GO, width = 28, height = 28, units = c("in"))

