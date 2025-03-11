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
library(plotly)
library(scDblFinder)

setwd("C:/Users/mchiesa/Desktop/liraglut_4_samples")
source("./0.Script/custom_seurat_functions.R")

# ctrl_glp_neg <- "./scRNA-seq_0153/cellranger_1286/filtered_feature_bc_matrix.h5"
# ctrl_glp_pos <- "./scRNA-seq_0153/cellranger_1287/filtered_feature_bc_matrix.h5"
# lira_glp_neg <- "./scRNA-seq_0153/cellranger_1288/filtered_feature_bc_matrix.h5"
# lira_glp_pos <- "./scRNA-seq_0153/cellranger_1289/filtered_feature_bc_matrix.h5"

ctrl_glp_neg <- "./scRNA-seq_0157/cellranger_1286/filtered_feature_bc_matrix.h5"
ctrl_glp_pos <- "./scRNA-seq_0157/cellranger_1287/filtered_feature_bc_matrix.h5"
lira_glp_neg <- "./scRNA-seq_0157/cellranger_1288/filtered_feature_bc_matrix.h5"
lira_glp_pos <- "./scRNA-seq_0157/cellranger_1289/filtered_feature_bc_matrix.h5"

####################

npcs_redux <- 30
set.seed(1)

####### read h5
matrix_ctrl_glp_neg <- Read10X_h5(ctrl_glp_neg,use.names = T)
gc()
matrix_ctrl_glp_pos <- Read10X_h5(ctrl_glp_pos,use.names = T)
gc()
matrix_lira_glp_neg <- Read10X_h5(lira_glp_neg,use.names = T)
gc()
matrix_lira_glp_pos <- Read10X_h5(lira_glp_pos,use.names = T)
gc()

### crate Seu obj
srat_ctrl_glp_neg <- CreateSeuratObject(matrix_ctrl_glp_neg,
                                        project = "ctrl_glp_neg",
                                        min.cells = 30,
                                        min.features = 500) # questo impatta su # cellule
rm(matrix_ctrl_glp_neg)
gc()

srat_ctrl_glp_pos <- CreateSeuratObject(matrix_ctrl_glp_pos,
                                        project = "ctrl_glp_pos",
                                        min.cells = 30,
                                        min.features = 500) # questo impatta su # cellule
rm(matrix_ctrl_glp_pos)
gc()

srat_lira_glp_neg <- CreateSeuratObject(matrix_lira_glp_neg,
                                        project = "lira_glp_neg",
                                        min.cells = 30,
                                        min.features = 500) # questo impatta su # cellule
rm(matrix_lira_glp_neg)
gc()

srat_lira_glp_pos <- CreateSeuratObject(matrix_lira_glp_pos,
                                        project = "lira_glp_pos",
                                        min.cells = 30,
                                        min.features = 500) # questo impatta su # cellule
rm(matrix_lira_glp_pos)
gc()

###########     QC    #######
# calculate MT % for each cell
srat_ctrl_glp_neg[["percent.mt"]]  <- PercentageFeatureSet(srat_ctrl_glp_neg, pattern = "^MT")
srat_ctrl_glp_pos[["percent.mt"]]  <- PercentageFeatureSet(srat_ctrl_glp_pos, pattern = "^MT")
srat_lira_glp_neg[["percent.mt"]]  <- PercentageFeatureSet(srat_lira_glp_neg, pattern = "^MT")
srat_lira_glp_pos[["percent.mt"]]  <- PercentageFeatureSet(srat_lira_glp_pos, pattern = "^MT")

VlnPlot(srat_ctrl_glp_neg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_ctrl_glp_pos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_lira_glp_neg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_lira_glp_pos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## filtering cells
srat_ctrl_glp_neg <- subset(srat_ctrl_glp_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
srat_ctrl_glp_pos <- subset(srat_ctrl_glp_pos, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
srat_lira_glp_neg <- subset(srat_lira_glp_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
srat_lira_glp_pos <- subset(srat_lira_glp_pos, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

VlnPlot(srat_ctrl_glp_neg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_ctrl_glp_pos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_lira_glp_neg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(srat_lira_glp_pos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

gc()

#### Normalize separate slots (samples)
cd34_list <- list()

cd34_list[["srat_ctrl_glp_neg"]] <- srat_ctrl_glp_neg
rm(srat_ctrl_glp_neg)
gc()

cd34_list[["srat_ctrl_glp_pos"]] <- srat_ctrl_glp_pos
rm(srat_ctrl_glp_pos)
gc()

cd34_list[["srat_lira_glp_neg"]] <- srat_lira_glp_neg
rm(srat_lira_glp_neg)
gc()

cd34_list[["srat_lira_glp_pos"]] <- srat_lira_glp_pos
rm(srat_lira_glp_pos)
gc()

####
cd34_list <- lapply(X = cd34_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = F)
  x <- ScaleData(x, verbose = F, vars.to.regress = "percent.mt")
  x <- RunPCA(x, npcs = npcs_redux, verbose = F)
# x <- RunUMAP(x, reduction = "pca", dims = 1:npcs_redux, verbose = F)
  x <- RunUMAP(x, reduction = "pca", dims = 1:npcs_redux, verbose = F,n.components = 3)
# x <- RunTSNE(x, tsne.method = "Rtsne", dims = 1:npcs_redux, verbose = F)
})

gc()
################################ Data Integration
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = cd34_list)
cd34_anchors <- FindIntegrationAnchors(object.list = cd34_list,
                                         dims = 1:npcs_redux,
                                         anchor.features = features)
rm(cd34_list)
gc()

cd34_integ  <- IntegrateData(anchorset = cd34_anchors, dims = 1:npcs_redux)
rm(cd34_anchors)
gc()







################################################################################################
################################################################################################
#save.image("C:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.157.parz.RData")
load("C:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.RData")

################################################################################################




DefaultAssay(cd34_integ) <- "RNA"
cd34_integ <- JoinLayers(cd34_integ) # dalla versione 4 le conte sono divise per campione, con questa funzione si ricrea la matrice delle conte
gc()

## change Assay and integrate
DefaultAssay(cd34_integ) <- "integrated"
cd34_integ <- ScaleData(cd34_integ, verbose = F)
cd34_integ <- RunPCA(cd34_integ, npcs = npcs_redux, verbose = F)
gc()
#cd34_integ <- RunUMAP(cd34_integ, reduction = "pca", dims = 1:npcs_redux, verbose = F)
cd34_integ <- RunUMAP(cd34_integ, reduction = "pca", dims = 1:npcs_redux, verbose = F,n.components=3)
#cd34_integ <- RunTSNE(cd34_integ, tsne.method = "Rtsne", dims = 1:npcs_redux, verbose = F)
gc()


#################################### Add metadata informations
metadata_seu<- cd34_integ@meta.data
metadata_seu$receptor <- ifelse(metadata_seu$orig.ident %in% c("ctrl_glp_neg","lira_glp_neg"),"glp_neg","glp_pos")
metadata_seu$drug <- ifelse(metadata_seu$orig.ident %in% c("ctrl_glp_neg","ctrl_glp_pos"),"ctrl","lira")

cd34_integ$receptor <- metadata_seu$receptor
cd34_integ$drug <- metadata_seu$drug
#cd34_integ <- AddMetaData(object = cd34_integ, metadata = metadata_seu[,c("receptor","drug"),drop=F])


################################### Find clusters
ElbowPlot(cd34_integ,ndims = npcs_redux)
cd34_integ <- FindNeighbors(cd34_integ, dims = 1:npcs_redux, k.param = 30, verbose = F)
cd34_integ <- FindClusters(cd34_integ, verbose = F, resolution = 0.6) # prova resolution = xxx: Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
gc()

DimPlot(cd34_integ, reduction = "umap", split.by = "orig.ident",ncol = 2) + NoLegend()



#################################### identify doublets
# sce <- scDblFinder(cd34_integ@assays[["RNA"]]@layers[["data"]], clusters=Idents(cd34_integ))
# cd34_integ$scDblFinder.score <- sce$scDblFinder.score
# cd34_integ$scDblFinder.class <- sce$scDblFinder.class
# DimPlot(cd34_integ,
#         group.by="scDblFinder.class",
#         split.by="orig.ident",
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.box = T,
#         ncol = 2) + NoLegend()




####################################### 3D Plots
# cd34_integ@meta.data$seurat_clusters
# 
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# mycolor <- gg_color_hue(19) #,veri,colori,dei,cluster
# 
# plot_ly(x=cd34_integ@reductions[["umap"]]@cell.embeddings[,1],
#         y=cd34_integ@reductions[["umap"]]@cell.embeddings[,2], 
#         z=cd34_integ@reductions[["umap"]]@cell.embeddings[,3],
#         type="scatter3d",
#         mode="markers",
#         
#         color=cd34_integ@meta.data$seurat_clusters,
#         colors = mycolor,
#         
#         #colors=cd34_integ@assays$RNA["LYZ"],
#         
#         # color=cd34_integ@meta.data$orig.ident,
#         # colors = mycolor[c(1,9)],
#         size = 0.4
#         
#         )
# 
# 






DimPlot(cd34_integ,reduction = "umap",label = T,pt.size = 1,label.box = T) + NoLegend()
DimPlot(cd34_integ,label = T,pt.size = 0.5,split.by = "orig.ident") + NoLegend()
# DimPlot(cd34_integ,split.by = "new_groups",label = T) + NoLegend()
# RidgePlot(cd34_integ,features = features[15:20])


## ################### Cell ferquency by cluster ################
count_table <- as.data.frame(table(cd34_integ@meta.data$seurat_clusters,
                                   cd34_integ@meta.data$orig.ident))


# Stacked barplot with multiple groups

count_table_ctrl_glp_neg <- count_table[which(count_table$Var2 %in% "ctrl_glp_neg"),]
count_table_ctrl_glp_pos <- count_table[which(count_table$Var2 %in% "ctrl_glp_pos"),]
count_table_lira_glp_neg <- count_table[which(count_table$Var2 %in% "lira_glp_neg"),]
count_table_lira_glp_pos <- count_table[which(count_table$Var2 %in% "lira_glp_pos"),]



count_table_ctrl_glp_neg$Freq_Perc <- round(100*count_table_ctrl_glp_neg$Freq/sum(count_table_ctrl_glp_neg$Freq),1)
count_table_ctrl_glp_pos$Freq_Perc <- round(100*count_table_ctrl_glp_pos$Freq/sum(count_table_ctrl_glp_pos$Freq),1)
count_table_lira_glp_neg$Freq_Perc <- round(100*count_table_lira_glp_neg$Freq/sum(count_table_lira_glp_neg$Freq),1)
count_table_lira_glp_pos$Freq_Perc <- round(100*count_table_lira_glp_pos$Freq/sum(count_table_lira_glp_pos$Freq),1)




count_table <- rbind(count_table_ctrl_glp_neg,count_table_ctrl_glp_pos, count_table_lira_glp_neg, count_table_lira_glp_pos)

# ggplot(data=count_table, aes(x=Var1, y=Freq_Perc, fill=Var2)) +
#   geom_bar(stat="identity", position=position_dodge())


ggplot(data=count_table, aes(x=Var1, y=Freq_Perc, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~Var2)
  
  
  
# count_table
# 2 gruppi lira vs ctrl(k=30)
# cl 0
prop.test(x = c(8979, 6192), n = c(48092, 42069))
# cl 1
prop.test(x = c(7108, 4180), n = c(48092, 42069))
# cl 2
prop.test(x = c(6199, 3311), n = c(48092, 42069))
# cl 3
prop.test(x = c(3786, 2995), n = c(48092, 42069))
# cl 4
prop.test(x = c(734, 2592), n = c(48092, 42069))
# cl 5
prop.test(x = c(1486, 3292), n = c(48092, 42069))
# cl 6
prop.test(x = c(3112, 1665), n = c(48092, 42069))
# cl 7
prop.test(x = c(1524, 2964), n = c(48092, 42069))
# cl 8
prop.test(x = c(2452, 1988), n = c(48092, 42069))
# cl 9
prop.test(x = c(2373, 1750), n = c(48092, 42069))
# cl 10
prop.test(x = c(1956, 1974), n = c(48092, 42069))
# cl 11
prop.test(x = c(1982, 1624), n = c(48092, 42069))
# cl 12
prop.test(x = c(1671, 1535), n = c(48092, 42069))
# cl 13
prop.test(x = c(1332, 1622), n = c(48092, 42069))
# cl 14
prop.test(x = c(1079, 1265), n = c(48092, 42069))
# cl 15
prop.test(x = c(1084, 892), n = c(48092, 42069))
# cl 16
prop.test(x = c(883, 967), n = c(48092, 42069))
# cl 17
prop.test(x = c(108, 705), n = c(48092, 42069))
# cl 18
prop.test(x = c(244, 556), n = c(48092, 42069))


# save results
save.image("C:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.run157.RData")
