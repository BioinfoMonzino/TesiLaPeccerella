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


load("C:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.RData")

gc()
#DimPlot(cd34_integ,reduction = "umap",label = T,pt.size = 0.5,label.box = T) + NoLegend()
# 
# seu.obj.markers <- FindAllMarkers(cd34_integ,
#                                   only.pos = TRUE,
#                                   min.pct = 0.2,
#                                   logfc.threshold = 0.2,
#                                   verbose = T,
# 
#                                   #min.diff.pct = 0.1,
#                                   test.use = "wilcox",
#                                   # test.use = "roc", #Roc 1 o 0 perfect classification
#                                   #test.use = "LR",
#                                   #test.use = "DESeq2"
# )
# 
# seu.obj.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 200, wt = avg_log2FC) -> sign_marker
#top_n(n = 100, wt = p_val_adj) -> sign_marker

# #DoHeatmap(cd34_integ, features = sign_marker$gene) + NoLegend()
# write.table(sign_marker,"./Data_HT/Results/0.Results.k.30.3D/TopMarker.list.per.Cluster.txt",sep="\t",quote=F)


gc()

##########################################################################?
###### differential analysis by condition
cd34_integ$celltype.group <- paste(Idents(cd34_integ), 
                                      cd34_integ$orig.ident, 
                                      sep = "_")
cd34_integ$celltype <- Idents(cd34_integ)
Idents(cd34_integ) <- "celltype.group"


df_data_counts <- cd34_integ@assays[["RNA"]]@layers[["counts"]] # genes x cells
DefaultAssay(cd34_integ) <- "RNA"
rownames(df_data_counts) <- rownames(cd34_integ)
colnames(df_data_counts) <- colnames(cd34_integ)
gc()

meta_cl <- cd34_integ@meta.data[,c("seurat_clusters","celltype.group","orig.ident")]
#df_cl <- as.data.frame(t(df_data_counts)) # cells x genes

rm(cd34_integ)
gc()

## select only ctrls cells to identify markers (otherwice PC crashes)
idx_ctrl_sample <- which(meta_cl$orig %in% c("ctrl_glp_neg","ctrl_glp_pos"))

meta_cl <- meta_cl[idx_ctrl_sample,]
df_data_counts <- df_data_counts[,idx_ctrl_sample]


##########################
# puoi ciclare sull'id del cluster!

for (i in 0:18){
  gc()
  cluster_id <- i
  
  meta_cl$class <- c()
  meta_cl$class <- ifelse(meta_cl$seurat_clusters == cluster_id ,"1_clust","0_other")
  
  

  gc()
  
  ###################
  SE<-DaMiR.makeSE(df_data_counts, meta_cl) #gene x cells
  #rm(data_cl_rnd)
  gc()
  data_norm <- DaMiR.normalization(SE, type = "logcpm" , minCounts=1,fSample = 0.05, hyper = "no")
  rm(SE)
  gc()
  
  # DE Analysis
  design <- model.matrix(~0+class, meta_cl)
  contrasts <- makeContrasts(clust.vs.other=class1_clust-class0_other,
                             levels=design)
  
  fit <- lmFit(assay(data_norm), design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  res_table <- topTable(fit2,coef = 1,adjust="BH",number=10000000000000000000)
  hist(res_table[,4], main="histogram of pValues", breaks=20, col="orange", xlab="pval")
  #head(res_table,20)
  rm(data_norm)
  gc()
  # FeaturePlot(cd34_integ, features = c("Ankrd1", "Actn2"), split.by = "orig.ident", max.cutoff = 3, 
  #             cols = c("yellow", "purple4"))
  # plots <- VlnPlot(cd34_integ, features = c("Ankrd1", "Actn2"), split.by = "orig.ident", group.by = "celltype", 
  #                  pt.size = 0, combine = FALSE,split.plot = TRUE)
  # CombinePlots(plots = plots, ncol = 1)
  
  
  res_table_plot <- res_table
  
  p_th <- 0.05 #10^-2
  adjp_th <- 10^-5
  Beta_th <- 0
  
  
  sign.lfc.p.adj.vect <- sign.lfc.p.vect <-sign.vect.adj.p <- sign.vect.p <- rep("NOT.SIGN", length(rownames(res_table_plot)))
  sign.vect.p[which(res_table_plot$P.Value<=p_th)] <- "SIGN"
  #sign.vect.adj.p[which(res_table_plot$adj.P.Val<=adjp_th )] <- "SIGN"
  sign.lfc.p.vect[which(res_table_plot$P.Value<=p_th & res_table_plot$logFC >= Beta_th)] <- "SIGN_UP"
  sign.lfc.p.vect[which(res_table_plot$P.Value<=p_th & res_table_plot$logFC <= -Beta_th)] <- "SIGN_DOWN"
  sign.lfc.p.vect[which(res_table_plot$adj.P.Val<=adjp_th & res_table_plot$logFC >= Beta_th)] <- "SIGN_UPUP"
  sign.lfc.p.vect[which(res_table_plot$adj.P.Val<=adjp_th & res_table_plot$logFC <= -Beta_th)] <- "SIGN_DOWNDOWN"
  res_table_plot$sign.p <- as.factor(sign.vect.p)
  res_table_plot$sign.adj.p <- as.factor(sign.vect.adj.p)
  res_table_plot$sign.lfc.p <- as.factor(sign.lfc.p.vect)
  res_table_plot$sign.lfc.adj.p <- as.factor(sign.lfc.p.adj.vect)
  res_table_plot$name <- rownames(res_table_plot)
  
  if(length(levels(res_table_plot$sign.lfc.p)) ==3){
    gg1 <- ggplot(aes(x=logFC, y=-log10(P.Value),color=sign.lfc.p.vect),data=res_table_plot, label=name) +
      geom_point(aes(alpha=as.numeric(sign.p)),size=5)+
      scale_color_manual(values=c("gray50","lightblue","salmon")) +
      #scale_color_manual(values=c("gray50","salmon","red")) +
      geom_hline(yintercept = -log10(p_th),color='gray50',linetype = "dashed") + 
      geom_vline(xintercept = Beta_th,color='gray50',linetype = "dashed") +
      geom_vline(xintercept = -Beta_th,color='gray50',linetype = "dashed") +
      #geom_hline(yintercept = -log10(res_table_plot$P.Value[max(which(res_table_plot$adj.P.Val<adjp_th))]),color='gray50',linetype = "dashed") + 
      theme(legend.position = "none") +
      #xlim(-max(abs(min(res_table_plot$logFC)),max(res_table_plot$logFC),0)-0.5,max(abs(min(res_table_plot$logFC)),max(res_table_plot$logFC))+0.5) +
      xlim(-3,3) +
      ylim(0,70) +
      xlab("log2 (FC)")
    
    
  }else{
    # plot Volcano plot 
    
    gg1 <- ggplot(aes(x=logFC, y=-log10(P.Value),color=sign.lfc.p.vect),data=res_table_plot, label=name) +
      geom_point(aes(alpha=as.numeric(sign.p)),size=5)+
      scale_color_manual(values=c("gray50","lightblue","blue","salmon","red")) +
      #scale_color_manual(values=c("gray50","salmon","red")) +
      geom_hline(yintercept = -log10(p_th),color='gray50',linetype = "dashed") + 
      geom_vline(xintercept = Beta_th,color='gray50',linetype = "dashed") +
      geom_vline(xintercept = -Beta_th,color='gray50',linetype = "dashed") +
      geom_hline(yintercept = -log10(res_table_plot$P.Value[max(which(res_table_plot$adj.P.Val<adjp_th))]),color='gray50',linetype = "dashed") + 
      theme(legend.position = "none") +
      #xlim(-max(abs(min(res_table_plot$logFC)),max(res_table_plot$logFC),0)-0.5,max(abs(min(res_table_plot$logFC)),max(res_table_plot$logFC))+0.5) +
      xlim(-4,4) +
      ylim(0,300) +
      xlab("log2 (FC)") #+
    # geom_text_repel(
    #   aes(x=logFC, y=-log10(P.Value)),
    #   label = ifelse(abs(res_table_plot$logFC) > 2.5, res_table_plot$name, ""),
    #   box.padding = unit(0.45, "lines"),
    #   hjust = 1,max.overlaps = 2000,
    #   color="black"
    # )
  }
  
  
  pdf(file = paste0("./Results/1.DE.cluster/clust.",cluster_id,".vs.other.pdf"),width = 10,height = 10)
  print(gg1)
  print(hist(res_table[,4], main="histogram of pValues", breaks=20, col="orange", xlab="p-value"))
  dev.off()
  
  write.table(res_table, paste0("./Results/1.DE.cluster/clust.",cluster_id,".vs.other.txt"),sep="\t",quote=F,col.names = NA)
  
}



