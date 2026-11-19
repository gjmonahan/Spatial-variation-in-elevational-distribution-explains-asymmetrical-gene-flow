setwd("C:/Users/gabri/Documents/Thom_Lab_Project/DB_PGLS")
library(ape)
library(phytools)
library(readxl)
library(tibble)
library(nlme)
library(dplyr)
##All elevation
unpruned_tree <- read.tree("1_22sp_tree.tre")
all_species <- unpruned_tree[["tip.label"]]
species_to_keep <- c("Stephanoxis_lalandi", "Poospiza_thoracica", "Piculus_aurulentus", "Drymophila_rubricollis", "Dysithamnus_xanthopterus",
                     "Synallaxis_cinerascens", "Syndactyla_rufosuperciliata", "Heliobletus_contaminatus", "Hemitriccus_obsoletus", "Phylloscartes_difficilis", "Saltator_maxillosus", "Poospiza_lateralis",
                     "Stephanophorus_diadematus", "Tangara_desmaresti")
pruned_tree <- keep.tip(unpruned_tree, species_to_keep)
print(pruned_tree, printlen=14)
write.tree(pruned_tree, file="1_14sp_tree.tre")
#import data sheet
data <- read_xlsx("PGLS_11_4.xlsx") 
data <- column_to_rownames(data, var = colnames(data)[1])

#Importing Variables
migration_ratio <-setNames(data[,"migration_ratio"],
                     rownames(data))

all(names(migration_ratio) %in% pruned_tree$tip.label) 

2_deg_elevation_ratio <-setNames(data[,"2_deg_elevation_ratio"],
                       rownames(data))
all(names(2_deg_elevation_ratio) %in% pruned_tree$tip.label) 

all_elevation_ratio <-setNames(data[,"all_elevation_ratio"],
                       rownames(data))
all(names(all_elevation_ratio) %in% pruned_tree$tip.label) 

####PGLS
spp<-rownames(data)
corBM<-corBrownian(phy=pruned_tree,form=~spp)
corBM

pgls_2_degree <-gls(log(2_deg_elevation_ratio)~log(migration_ratio),
                      data=data,correlation=corBM)
summary(pgls_2_degree)

pgls_all <-gls(log(all_elevation_ratio)~log(migration_ratio),
                      data=data,correlation=corBM)
summary(pgls_all)