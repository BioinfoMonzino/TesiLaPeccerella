rm(list = ls())
cat("\014")
graphics.off()

library(Seurat)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(limma)
library(DaMiRseq)
library(ggrepel)

setwd("C:/Users/mchiesa/Desktop/2.1.LIRA.vs.CTRL")

source("./0.Script/custom_seurat_functions.R")

# SALVATO QUI
load("analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.RData")
n_subset_cell <- 30000
###############################################################################
gc()

cd34_integ <- RenameIdents(cd34_integ,
                           `0` = "00_Monocyte", 
                           `1` = "01_Dendritic_Cell", 
                           `2` = "02_Bcell_Progen", 
                           `3` = "03_Neutrophil",
                           `4` = "04_Eosinophil",
                           `5` = "05_Neutrophil",
                           `6` = "06_Granul_Monoc_Progen", 
                           `7` = "07_Megakaryocytes_Platelet",
                           `8` = "08_Macrophages", 
                           `9` = "09_Basophil",
                           `10` = "10_Hematop_Stem_Cell", 
                           `11` = "11_Basophil", 
                           `12` = "12_Granul_Monoc_Progen",  
                           `13` = "13_Erythroblast", 
                           `14` = "14_Mast_cell", 
                           `15` = "15_Bcell_Progen", 
                           `16` = "16_Granul_Monoc_Progen",
                           `17` = "17_unknown",
                           `18` = "18_Granul_Monoc_Progen"
)


cd34_integ$cell_type <- Idents(cd34_integ)
cd34_integ$celltype.group <- paste(Idents(cd34_integ), 
                                   cd34_integ$drug, 
                                   sep = "_")
gc()




## subsampling
set.seed(123)
#cd34_integ.subsampled <- cd34_integ[, sample(colnames(cd34_integ), size = n_subset_cell, replace=F)]
cd34_integ.subsampled <- cd34_integ


metadata_seu<- cd34_integ.subsampled@meta.data
rm(cd34_integ)
gc()
##############
# take raw data and normalise it
ctrl_gr <- which(metadata_seu[,6] %in% "ctrl" & metadata_seu[,10] %in% c("00_Monocyte",
                                                                         "01_Dendritic_Cell",
                                                                         "03_Neutrophil",
                                                                         "04_Eosinophil",
                                                                         "05_Neutrophil",
                                                                         "08_Macrophages", 
                                                                         "09_Basophil",
                                                                         "11_Basophil", 
                                                                         "14_Mast_cell"
                                                                         )
                )

lira_gr <- which(metadata_seu[,6] %in% "lira" & metadata_seu[,10] %in% c("00_Monocyte",
                                                                         "01_Dendritic_Cell",
                                                                         "03_Neutrophil",
                                                                         "04_Eosinophil",
                                                                         "05_Neutrophil",
                                                                         "08_Macrophages", 
                                                                         "09_Basophil",
                                                                         "11_Basophil", 
                                                                         "14_Mast_cell"
                                                                         )
                 )




# LogCPM
# count_raw_CM_FB_DMD <- cd34_integ@assays$RNA@counts[, CM_FB_DMD] + 1 # +1 x zeri logaritmo
# count_raw_CM_FB_ISO <- cd34_integ@assays$RNA@counts[, CM_FB_ISO] + 1 # +1 x zeri logaritmo
# count_norm_CM_FB_DMD <- apply(count_raw_CM_FB_DMD, 2, function(x) log((x/sum(x))*10000,2))
# count_norm_CM_FB_ISO <- apply(count_raw_CM_FB_ISO, 2, function(x) log((x/sum(x))*10000,2))

# # CPM
# count_raw_CM_FB_DMD <- cd34_integ@assays$RNA@counts[, CM_FB_DMD] 
# count_raw_CM_FB_ISO <- cd34_integ@assays$RNA@counts[, CM_FB_ISO] 
# 
# cpm_norm_CM_FB_DMD <- apply(count_raw_CM_FB_DMD, 2, function(x) (x/sum(x))*10000)
# cpm_norm_CM_FB_ISO <- apply(count_raw_CM_FB_ISO, 2, function(x) (x/sum(x))*10000)

# Seurat Norm
seu_norm_ctrl <- cd34_integ.subsampled@assays[["RNA"]]@layers[["data"]][, ctrl_gr] 
seu_norm_lira <- cd34_integ.subsampled@assays[["RNA"]]@layers[["data"]][, lira_gr] 

metadata_seu_ctrl_gr <- metadata_seu[ctrl_gr,]
metadata_seu_lira_gr <- metadata_seu[lira_gr,]


DefaultAssay(cd34_integ.subsampled) <- "RNA"
rownames(seu_norm_ctrl) <- rownames(cd34_integ.subsampled)
colnames(seu_norm_ctrl) <- rownames(metadata_seu_ctrl_gr)

rownames(seu_norm_lira) <- rownames(cd34_integ.subsampled)
colnames(seu_norm_lira) <- rownames(metadata_seu_lira_gr)

rm(cd34_integ.subsampled)
gc()

####################################### Write Output
############ CTRL

########### Assicurati che i nomi delle cellule siano uguali nel dataset e nel metadata! 
#(spesso hanno caratteri diversi)

write.table(round(seu_norm_ctrl,2), # bisogna poi aggiungere 'Gene' come nome prima colonna
            "seu.norm.ctrl.sub.All.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

write.table(metadata_seu_ctrl_gr, # bisogna poi aggiungere 'Cell' come nome prima colonna
            "ctrl.Metadata.sub.All.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

############ Lira
write.table(round(seu_norm_lira,2), # bisogna poi aggiungere 'Gene' come nome prima colonna
            "seu.norm.lira.sub.All.txt", 
            sep="\t",
            quote=F,
            col.names = NA)

write.table(metadata_seu_lira_gr, 
            "lira.Metadata.sub.All.txt", 
            sep="\t",
            quote=F,
            col.names = NA) # aggiungere 'Cell' come nome prima colonna
