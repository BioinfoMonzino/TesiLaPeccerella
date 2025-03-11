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




setwd("C:/Users/mchiesa/Desktop/liraglut_4_samples/")
source("./0.Script/custom_seurat_functions.R")


# SALVATO QUI
load("c:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.RData")
gc()

#cd34_integ@meta.data$seurat_clusters





gc()

# k=30
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

DimPlot(cd34_integ,
        reduction = "umap",
        label = T,
        pt.size = 1.5,
        label.box = T,
        #group.by="active.ident",
        #split.by = "orig.ident"
) + NoLegend()
gc()

DimPlot(cd34_integ,
        reduction = "umap",
        #label = T,
        pt.size = 0.5,
        label.box = T,
        #group.by="active.ident",
        split.by = "orig.ident"
) 
gc()
gc()



metadata_seu<- cd34_integ@meta.data

metadata_seu$cell_macrogroups <- ifelse(metadata_seu$seurat_clusters %in% c("0","1","8"),"Monocytes",
                                        ifelse(metadata_seu$seurat_clusters %in% c("2","15"),"Bcell",
                                               ifelse(metadata_seu$seurat_clusters %in% c("3","4","5","9","11","14"),"Granulocytes",
                                                      ifelse(metadata_seu$seurat_clusters %in% c("12","16","18","6"),"GM_Prog",
                                                             ifelse(metadata_seu$seurat_clusters %in% c("7","13"),"Megakar_Erythr_Prog",
                                                                    ifelse(metadata_seu$seurat_clusters %in% c("10"),"Hematop_Stem_Cells",
                                                                           ifelse(metadata_seu$seurat_clusters %in% c("17"),"unknown",NA)
                                                                    )
                                                                    )
                                                             )
                                               )
                                        )
)

metadata_seu$cell_macrogroups2 <- ifelse(metadata_seu$seurat_clusters %in% c("2","15"),"Lymphoid",
                                        ifelse(metadata_seu$seurat_clusters %in% c("10"),"Hematop_Stem_Cells",
                                               ifelse(metadata_seu$seurat_clusters %in% c("17"),"unknown","Myeloid")
                                        )
)
 
metadata_seu$cell_macrogroups3 <- ifelse(metadata_seu$seurat_clusters %in% c("0","1","8"),"Mature",
                                        ifelse(metadata_seu$seurat_clusters %in% c("2","15"),"Progenitor",
                                               ifelse(metadata_seu$seurat_clusters %in% c("3","4","5","9","11","14"),"Mature",
                                                      ifelse(metadata_seu$seurat_clusters %in% c("6","12","16","18"),"Progenitor",
                                                             ifelse(metadata_seu$seurat_clusters %in% c("7","13"),"Progenitor",
                                                                    ifelse(metadata_seu$seurat_clusters %in% c("10"),"Hematop_Stem_Cells",
                                                                           ifelse(metadata_seu$seurat_clusters %in% c("17"),"unknown",NA)
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
)



cd34_integ <- AddMetaData(object = cd34_integ, metadata = metadata_seu[,10:12,drop=F])


DimPlot(cd34_integ,reduction = "umap",
        label = F,
        pt.size = 0.5,
        label.box = T,
        group.by="cell_macrogroups",
        #split.by = "orig.ident"
        ) 

DimPlot(cd34_integ,
        reduction = "umap",
        #label = T,
        pt.size = 1.5,
        label.box = T,
        #group.by="active.ident",
        #split.by = "orig.ident"
) 
gc()





gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycolor <- gg_color_hue(20) #,veri,colori,dei,cluster

plot_ly(x=cd34_integ@reductions[["umap"]]@cell.embeddings[,1],
        y=cd34_integ@reductions[["umap"]]@cell.embeddings[,2], 
        z=cd34_integ@reductions[["umap"]]@cell.embeddings[,3],
        type="scatter3d",
        mode="markers",
        
        color=cd34_integ@active.ident,
        colors = mycolor,
        
        #colors=cd34_integ@assays$RNA["LYZ"],
        
        # color=cd34_integ@meta.data$orig.ident,
        # colors = mycolor[c(1,9)],
        size = I(5)
        
)





##############################################?
#CD133 E CD34 :HSC
FeaturePlot(cd34_integ,features = "rna_PROM1",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_CD34",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#Cd45
FeaturePlot(cd34_integ,features = "rna_PTPRC",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


FeaturePlot(cd34_integ,features = "rna_CD38",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#Monocytes
FeaturePlot(cd34_integ,features = "rna_CD4",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#cd61: PLATELETS
FeaturePlot(cd34_integ,features = "rna_ITGB3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# CD123 GRANULOCITI E MONOCITI
FeaturePlot(cd34_integ,features = "rna_IL3RA",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# CXCR4 Cellule T
FeaturePlot(cd34_integ,features = "rna_CXCR4",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# rna_TNFRSF18 Cellule T
FeaturePlot(cd34_integ,features = "rna_TNFRSF18",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# rna_PRTN3 (GPA) GLOBULI ROSSI
FeaturePlot(cd34_integ,features = "rna_PRTN3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# rna_MAST CELL (CD203c)
FeaturePlot(cd34_integ,features = "rna_ENPP3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


# MEGAKARIOCYTE (PAPER PELLIN 2019)
FeaturePlot(cd34_integ,features = "rna_PF4",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_VWF",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#ERITROBLAST (PELLIN 2019)
FeaturePlot(cd34_integ,features = "rna_HBB",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#DENDRIITIC CELLS (PELLIN 2019)
FeaturePlot(cd34_integ,features = "rna_MPEG1",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


#GRANULOCITI(PELLIN 2019)/NEUTROFILI
FeaturePlot(cd34_integ,features = "rna_ELANE",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_MPO",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_LYZ",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_PRTN3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#LYNFOID (PELLIN 2019)
FeaturePlot(cd34_integ,features = "rna_DDIT4",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


#MAST CELL (HERAULT 2021)
FeaturePlot(cd34_integ,features = "rna_CPA3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_FCER1A",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#NETUROFILI (HERAULT 2021)
FeaturePlot(cd34_integ,features = "rna_MPO",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# LYMPHOID CELL 2 (HERAULT 2021)
FeaturePlot(cd34_integ,features = "rna_MZB1",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#LYMPHOID CELL 1 - T CELL (HERAULT 2021)
FeaturePlot(cd34_integ,features = "rna_TRBC2",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# MEGAKARIOCITES (HERAULT 2021))
FeaturePlot(cd34_integ,features = "rna_GP1BB",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


# DIV DIVISION (HERAULT 2021))
FeaturePlot(cd34_integ,features = "rna_HMGB2",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_UBE2C",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_BIRC5",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_CKS2",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_TOP2A",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# REP REPAIR (HERAULT 2021)) - HMMR2
FeaturePlot(cd34_integ,features = "rna_MFSD6",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_DUT",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_UNG",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


#EOSINOFILI (Kucinski 2024)
FeaturePlot(cd34_integ,features = "rna_PRG3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_EPX",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#ERITHROID PROGENITOR (KUCINSKI 2024)
FeaturePlot(cd34_integ,features = "rna_KLF1",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


#HSC PROGENITOR (KUCINSKI 2024)
FeaturePlot(cd34_integ,features = "rna_MECOM",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_CD34",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

#LYMPHOID PROGENITOR (KUCINSKI 2024)
FeaturePlot(cd34_integ,features = "rna_DNTT",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# t CELL PROG (KUCINSKI 2024)
FeaturePlot(cd34_integ,features = "rna_CD3E",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

# B CELL PROG (KUCINSKI 2024)
FeaturePlot(cd34_integ,features = "rna_VPREB3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


FeaturePlot(cd34_integ,features = "rna_IFITM3",pt.size = 0.5,cols = c("lightyellow","green", "blue"))

FeaturePlot(cd34_integ,features = "rna_CD3E",pt.size = 0.5,cols = c("lightyellow","green", "blue"))
FeaturePlot(cd34_integ,features = "rna_CD3E",pt.size = 0.5,cols = c("lightyellow","green", "blue"))


###############################


VlnPlot(cd34_integ,
        assay = "RNA",
        features = c("IFITM3","CD3E"),
        log = F,
        same.y.lims = T,
)



seu.obj.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> sign_marker
#top_n(n = 100, wt = p_val_adj) -> sign_marker
#write.table(sign_marker,"./TopMarker.list.per.Cluster.txt",sep="\t",quote=F)
gc()

DoHeatmap(cd34_integ, features = sign_marker$gene) + NoLegend()
gc()

#FeaturePlot(cd34_integ,pt.size = 0.1,features = c("Acta1","Myh8","Myh3","Myh1"),cols = c("lightgrey", "red"))

######## una volta identificati i tipi cellulari ? possibile rinominare i cluster coi nomi delle linee cellulare


# sign_marker
marker.to.plot <- sign_marker[which(sign_marker$p_val_adj < 10e-320),7,drop=T]
# DotPlot(cd34_integ, 
#         features = unique(marker.to.plot),
#         cols = c("yellow","purple4")) +
#   RotatedAxis() + 
#   NoLegend() +
#   theme(axis.text.x =element_text(size=8),
#   )

DotPlot(cd34_integ, 
        features = unique(marker.to.plot),
        cols = c("yellow","purple4"),assay = "RNA") +
  RotatedAxis() + 
  NoLegend() +
  theme(axis.text.x =element_text(size=4),
  )



## Carmeliet et al: selected best genes
#marker.to.plot <- c("FBLN5","HEY1","MECOM","CXCL12","RBP7","KDR","ENDOU","VCAM1","FMO1","FMO2","MGP","CFH","BGN","VWF","ISG15","IFIT3","IFI203","IFIT1","IFIT3B","COL4A2","APLN","SPARCL1","APLNR","TRP53I11","CCL21A","PRSS23","LYVE1","FXYD6","CP")
# marker.to.plot <- c("FBLN5", "STMN2", "8430408G22RIK", "GLUL", "FOS", "ID1", "HEY1", "SOX17", "GADD45G", "BTG2", "IER2", "KLF4", "ALPL", "JUNB", "DUSP1", "CYR61", "MECOM", "HSPA1A", "KLF2", "JUN", "CXCL12", "RBP7", "MGLL", "LY6C1", "AQP7", "BTNL9", "PDGFD", "SEMA7A", "CLEC2D", "UNC5B", "TPST1", "PLAT", "MAGIX", "PLSCR2", "FKBP3", "ESM1", "MCF2L", "SCGB3A1", "SLC26A10", "CABLES2", "KDR", "SSH2", "ENDOU", "FAM57B", "KIFC3", "NTF3", "TMEM182", "FAM212B", "4430402I18RIK", "CORO2B", "FAM117B", "LNX2", "3110062M04RIK", "INSIG1", "EBF3", "SLCO2B1", "MAP3K14", "AU020206", "DGKE", "2310040G24RIK", "VCAM1", "IER3", "PLTP", "PI16", "FMO2", "FMO1", "CALCRL", "EMP1", "AU021092", "BSG", "NFKBIA", "ELN", "BHLHE40", "ICAM1", "KLK8", "SOCS3", "SLFN2", "2200002D01RIK", "RND1", "EMP2", "RGCC", "LPL", "ISG15", "CAR4", "IFIT3", "RTP4", "IFI203", "RSAD2", "TCF15", "IFIT1", "LY6A", "MNDAL", "FABP4", "IIGP1", "AQP1", "IFIT3B", "AW112010", "TIMP4", "GPIHBP1", "GBP7", "TMSB10", "IFT122", "SPARC", "CLDN5", "COL4A2", "APLN", "VIM", "IGFBP7", "VWA1", "MEST", "SPARCL1", "ADM", "APLNR", "MEOX1", "TRP53I11", "PRNP", "NRP2", "MYCN", "TNFAIP8L1", "CCND1", "MGP", "CFH", "APOE", "CPE", "CYTL1", "BGN", "PLVAP", "DCN", "CTSH", "RBP1", "NPR3", "VWF", "H19", "TM4SF1", "IGFBP4", "FABP5", "TMEM108", "ID2", "CGNL1", "CLU", "CCL21A", "MMRN1", "FGL2", "PRSS23", "THY1", "IGFBP5", "FTH1", "LYVE1", "PRELP", "LCN2", "FXYD6", "NTS", "S100A6", "IFI27L2A", "CD63", "CD9", "PARD6G", "CP", "TIMP2", "LRG1")

marker.to.plot <- c("S100A9", "S100A8", "HLA-DRA", "VCAN", "FTH1", "CXCL8", "HLA-DRB1", "IFI30", "LYZ", "MRC1", "CD74", "FCN1", "CST3", "S100A10", "F13A1", "CD14", "TYROBP", "ANXA2", "NAMPT", "DPYD", "S100A4", "FUCA1", "HLA-DPA1", "S100A12", "AIF1", "HLA-DPB1", "HLA-DQB1", "S100B", "HLA-DQA1", "SLC8A1", "CCSER1", "JAML", "PKIB", "HLA-DMA", "HLA-DRB5", "CLEC10A", "ADAM28", "HLA-DMB", "SAMHD1", "GRN", "LMNA", "GSN", "HIST1H4C", "MKI67", "HMGB2", "H2AFZ", "TUBB", "TOP2A", "STMN1", "NUSAP1", "TUBA1B", "CENPF", "DIAPH3", "ATAD2", "AC020656.1", "PCLAF", "DUT", "SMC4", "RRM2", "TMPO", "TYMS", "DEK", "LINC01572", "PRTN3", "MPO", "ELANE", "AZU1", "SRGN", "CTSG", "PRSS57", "CST7", "SERPINB1", "NUCB2", "FNDC3B", "CALR", "KCNQ5", "P4HB", "PLAC8", "GYPC", "IGLL1", "HSP90B1", "PDE4D", "MS4A3", "SNHG29", "CDK6", "HSPA5", "PLCB1", "FAM107B", "PRG2", "CLC", "EPX", "RNASE2", "ACSM3", "PRG3", "IKZF2", "LTC4S", "SLC24A3", "DACH1", "CPA3", "STXBP5", "TNIK", "SLC39A11", "CD63", "THSD7A", "ANXA1", "PLIN2", "PKD2", "ALS2", "BPI", "RETN", "S100P", "MNDA", "CD24", "TSPO", "SLPI", "RPS26", "SERPINB10", "CSTA", "CFD", "TNFSF13B", "EPHA1-AS1", "WDR49", "SIPA1L1", "LYST", "SNHG5", "LRMDA", "NKG7", "MGST1", "PRKDC", "ADK", "ZEB2", "HBD", "GP1BB", "LTBP1", "PF4", "RAB27B", "DNM3", "PLEK", "ITGA2B", "PRKAR2B", "ABCC4", "MED12L", "THBS1", "INPP4B", "ANGPT1", "RAP1B", "ARHGAP6", "SOS1", "LAT", "PDLIM1", "TAGLN2", "TPM1", "LIMS1", "MEIS1", "CSTB", "CTSB", "SOD2", "GPNMB", "LGALS3", "FTL", "IGSF6", "CTSZ", "PLA2G7", "FCER1G", "MS4A7", "CD68", "NPC2", "CTSD", "CSF2RB", "PLPP1", "HDC", "CD52", "AKAP12", "ALOX5AP", "GCSAML", "RFLNB", "IFITM2", "ATP10D", "TCN1", "CLU", "RUNX1", "LGALS1", "NEGR1", "MSI2", "AFF3", "AUTS2", "SSBP2", "FAM30A", "NPM1", "C1QTNF4", "SOX4", "RNF220", "TCF4", "MIR181A1HG", "PEBP1", "COL24A1", "PDE7A", "RPL7A", "RPL14", "RPS4X", "IGFBP7", "SPINK2", "CD96", "BAALC", "HMGB1", "PDIA6", "HBB", "HBG2", "XACT", "BLVRB", "PRDX2", "HSP90AB1", "HSPD1", "APOC1", "HSPE1", "PTMA", "EIF5A", "TFRC", "MPC2", "RYR3", "NFIA", "UROD", "PIP5K1B", "ST8SIA6", "PVT1", "REXO2", "HPGD", "HPGDS", "LMO4", "GATA2", "RPS6KA5", "PRDX6", "RHEX", "ARHGAP18", "CTTNBP2", "PALM2-AKAP2", "KIT", "EXT1", "GRAP2", "SMYD3", "DENND1B", "IL18R1", "PTTG1", "CCNB1", "CDC20", "HMMR", "BIRC5", "CDKN3", "CCNB2", "DLGAP5", "HMGN2", "HSP90AA1", "NUCKS1", "RAN", "TUBB4B", "HSPA8", "NCL", "ARL6IP1", "AC079793.1", "ATP8B4", "MALAT1", "FHIT", "AP003086.1", "AP001011.1", "CPB2-AS1", "MBNL1-AS1", "ANKRD28", "CHST11", "DIAPH2", "MIR924HG", "FILIP1L", "RAB11A")


# DotPlot(cd34_integ, 
#         features = unique(marker.to.plot),
#         cols = c("yellow","purple4")) +
#   RotatedAxis() + 
#   NoLegend() +
#   theme(axis.text.x =element_text(size=8),
#   )

DotPlot(cd34_integ, 
        features = unique(marker.to.plot),
        cols = c("yellow","purple4"),assay = "RNA") +
  RotatedAxis() + 
  #NoLegend() +
  theme(axis.text.x =element_text(size=4),
  )

