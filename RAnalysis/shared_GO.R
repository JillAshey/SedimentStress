## Sediment stress - determining shared GO terms
## J. Ashey
## 20220125

# Load libraries
library(VennDiagram)
library(tidyverse)

## This code is showing the shared GO terms between species 

## read in BP GO data 
all.BP  <- read.csv("Output/GOSeq/all.go.BP_20220121.csv")
head(all.BP)
all.BP <- all.BP[,-c(1)] # remove weird first col

## subset by species 
acerv <- all.BP %>% filter(Species=="A.cervicornis")
length(unique(acerv$GO.IDs)) # 278 - there are some duplicates because certain terms matched to multiple GO slim terms 
mcav <- all.BP %>% filter(Species=="M.cavernosa")
length(unique(mcav$GO.IDs)) # 158 - there are some duplicates because certain terms matched to multiple GO slim terms 
ofav <- all.BP %>% filter(Species=="O.faveolata")
length(unique(ofav$GO.IDs)) # 27 - there are some duplicates because certain terms matched to multiple GO slim terms 
pacuta <- all.BP %>% filter(Species=="P.acuta") 
length(unique(pacuta$GO.IDs)) # 380 - there are some duplicates because certain terms matched to multiple GO slim terms 
plob <- all.BP %>% filter(Species=="P.lobata") 
length(unique(plob$GO.IDs)) # 198 - there are some duplicates because certain terms matched to multiple GO slim terms 


# Merge by species to find shared GO terms 

## acerv + mcav = 12 GO terms shared
AM <- merge(acerv, mcav, by = "GO.IDs")

## acerv + ofav = 2 GO terms shared
AO <- merge(acerv, ofav, by = "GO.IDs")

## acerv + pacuta = 5 GO terms shared
APa <- merge(acerv, pacuta, by = "GO.IDs")

## acerv + plob = 2 GO terms shared
APl <- merge(acerv, plob, by = "GO.IDs")

## mcav + ofav = 0 GO terms shared
MO <- merge(mcav, ofav, by = "GO.IDs")

## mcav + pacuta = 3 GO terms shared
MPa <- merge(mcav, pacuta, by = "GO.IDs")

## mcav + plob = 3 GO terms shared
MPl <- merge(mcav, plob, by = "GO.IDs")

## ofav + pacuta = 3 GO terms shared
OPa <- merge(ofav, pacuta, by = "GO.IDs")

## ofav + plob = 3 GO terms shared
OPl <- merge(ofav, plob, by = "GO.IDs")

## pacuta + plob = 3 GO terms shared
PaPl <- merge(pacuta, plob, by = "GO.IDs")




# Read in BP GO slim data 
all.BP.slim  <- read.csv("Output/GOSeq/all.go.BP.slim_20220121.csv")
all.BP.slim <- na.omit(all.BP.slim) # removing GO slim terms that did not have any GO ID matches in my dataset 

## subset by species 
acerv.slim <- all.BP.slim %>% filter(Species=="A.cervicornis")
length(unique(acerv.slim$GO.IDs)) # 184 - there are some duplicates because certain terms matched to multiple GO slim terms 
mcav.slim <- all.BP.slim %>% filter(Species=="M.cavernosa")
length(unique(mcav.slim$GO.IDs)) # 94 - there are some duplicates because certain terms matched to multiple GO slim terms 
ofav.slim <- all.BP.slim %>% filter(Species=="O.faveolata")
length(unique(ofav.slim$GO.IDs)) # 21 - there are some duplicates because certain terms matched to multiple GO slim terms 
pacuta.slim <- all.BP.slim %>% filter(Species=="P.acuta") 
length(unique(pacuta.slim$GO.IDs)) # 255 - there are some duplicates because certain terms matched to multiple GO slim terms 
plob.slim <- all.BP.slim %>% filter(Species=="P.lobata") 
length(unique(plob.slim$GO.IDs)) # 133 - there are some duplicates because certain terms matched to multiple GO slim terms 


# Merge by species to find shared GO terms / GO slim terms 

## acerv + mcav = 8 GO terms w/ GO slim info shared
AM.slim <- merge(acerv.slim, mcav.slim, by = "GO.IDs")

## acerv + ofav = 2 GO terms w/ GO slim info shared
AO.slim <- merge(acerv.slim, ofav.slim, by = "GO.IDs")

## acerv + pacuta = 4 GO terms w/ GO slim info shared
APa.slim <- merge(acerv.slim, pacuta.slim, by = "GO.IDs")

## acerv + plob = 5 GO terms w/ GO slim info shared
APl.slim <- merge(acerv.slim, plob.slim, by = "GO.IDs")

## mcav + ofav = 2 GO terms w/ GO slim info shared
MO.slim <- merge(mcav.slim, ofav.slim, by = "GO.IDs")

## mcav + pacuta = 2 GO terms w/ GO slim info shared
MPa.slim <- merge(mcav.slim, pacuta.slim, by = "GO.IDs")

## mcav + plob = 1 GO terms w/ GO slim info shared
MPl.slim <- merge(mcav.slim, plob.slim, by = "GO.IDs")

## ofav + pacuta = 1 GO terms w/ GO slim info shared
OPa.slim <- merge(ofav.slim, pacuta.slim, by = "GO.IDs")

## ofav + plob = 0 GO terms w/ GO slim info shared
OPl.slim <- merge(ofav.slim, plob.slim, by = "GO.IDs")

## pacuta + plob = 11 GO terms w/ GO slim info shared
PaPl.slim <- merge(pacuta.slim, plob.slim, by = "term")

