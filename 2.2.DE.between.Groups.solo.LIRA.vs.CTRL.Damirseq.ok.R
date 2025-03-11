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


## ################### Cell ferquency by cluster ################
# Stacked barplot with multiple groups
count_table <- as.data.frame(table(cd34_integ@meta.data$seurat_clusters,
                                   cd34_integ@meta.data$drug))

count_table_ctrl <- count_table[which(count_table$Var2 %in% "ctrl"),]
count_table_lira <- count_table[which(count_table$Var2 %in% "lira"),]



count_table_ctrl$Freq_Perc <- round(100*count_table_ctrl$Freq/sum(count_table_ctrl$Freq),1)
count_table_lira$Freq_Perc <- round(100*count_table_lira$Freq/sum(count_table_lira$Freq),1)




count_table <- rbind(count_table_ctrl, count_table_lira)

ggplot(data=count_table, aes(x=Var1, y=Freq_Perc, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge())


ggplot(data=count_table, aes(x=Var1, y=Freq_Perc, fill=Var2)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_wrap(~Var2)



###############################################################################
gc()

cd34_integ <- RenameIdents(cd34_integ,
                              `0` = "HSC_00",
                              `1` = "HSC_01",
                              `2` = "HSC_02",
                              `3` = "HSC_03",
                              `4` = "HSC_04",
                              `5` = "HSC_05",
                              `6` = "HSC_06",
                              `7` = "HSC_07",
                              `8` = "HSC_08",
                              `9` = "HSC_09",
                              `10` = "HSC_10",
                              `11` = "HSC_11",
                              `12` = "HSC_12",
                              `13` = "HSC_13",
                              `14` = "HSC_14",
                              `15` = "HSC_15",
                              `16` = "HSC_16",
                              `17` = "HSC_17",
                              `18` = "HSC_18"
                              #
                              
)


#DimPlot(cd34_integ,label = T,pt.size = 1,label.box = T) + NoLegend()

#cd34_integz <- RunUMAP(cd34_integ, reduction = "pca", dims = 1:npcs_redux, verbose = F,n.components = 3)
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# mycolor <- gg_color_hue(18) #,veri,colori,dei,cluster
# plot3d( 
#   x=cd34_integ@reductions[["tsne"]]@cell.embeddings[,1],
#   y=cd34_integ@reductions[["tsne"]]@cell.embeddings[,2], 
#   z=cd34_integ@reductions[["tsne"]]@cell.embeddings[,3], 
#   col = metadata_seu<- as.numeric(cd34_integ@meta.data[,7]), 
#   type = 's', 
#   radius = .1,
#   xlab="umap1", ylab="umap2", zlab="umap3")
# 


gc()

##########################################################################?
###### differential analysis by condition
cd34_integ$celltype.group <- paste(Idents(cd34_integ), 
                                      cd34_integ$drug, 
                                      sep = "_")
cd34_integ$celltype <- Idents(cd34_integ)
Idents(cd34_integ) <- "celltype.group"


df_data_counts <- cd34_integ@assays[["RNA"]]@layers[["counts"]]
#df_dati_scaled <- cd34_integ@assays$integrated@scale.data
df_metadata <- cd34_integ@meta.data


DefaultAssay(cd34_integ) <- "RNA"
rownames(df_data_counts) <- rownames(cd34_integ)
colnames(df_data_counts) <- colnames(cd34_integ)
rm(cd34_integ)
gc()



for (j in 0:18){
  cat(paste0(j,"\n"))
  gc()
  check_cluster <- j
  
  index_cl <- which(df_metadata$seurat_clusters == check_cluster)
  df_cl <- as.data.frame(t(df_data_counts[,index_cl]))
  meta_cl <- df_metadata[index_cl,c("celltype.group","drug")]
  
  data_meta <- cbind(meta_cl,df_cl)
  
  # # chose all poplulation (un-paired)
  meta_cl_rnd <- data_meta[,c(1,2)]
  data_cl_rnd <- data_meta[,-c(1,2)]
  colnames(meta_cl_rnd) <- c("group","class")
  
  
  
  
  ###################
  SE<-DaMiR.makeSE(t(data_cl_rnd), meta_cl_rnd)
  data_norm <- DaMiR.normalization(SE, type = "logcpm" , minCounts=1,fSample = 0.05, hyper = "no")
  
  # DE Analysis
  design <- model.matrix(~0+class, meta_cl_rnd)
  contrasts <- makeContrasts(lira.vs.ctrl=classlira-classctrl,
                             levels=design)
  
  fit <- lmFit(assay(data_norm), design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  res_table <- topTable(fit2,coef = 1,adjust="BH",number=10000000000000000000)
  
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
      #xlim(-3,3) +
      #ylim(0,70) +
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
     # xlim(-3,3) +
     # ylim(0,100) +
      xlab("log2 (FC)") #+
    # geom_text_repel(
    #   aes(x=logFC, y=-log10(P.Value)),
    #   label = ifelse(abs(res_table_plot$logFC) > 2.5, res_table_plot$name, ""),
    #   box.padding = unit(0.45, "lines"),
    #   hjust = 1,max.overlaps = 2000,
    #   color="black"
    # )
  }
  
  pdf(file = paste0("./Results/2.DE.group/2.1.LIRA.vs.CTRL/clust.",j,".lira.vs.ctrl.pdf"),width = 10,height = 10)
  print(gg1)
  print(hist(res_table[,4], main="histogram of pValues", breaks=20, col="orange", xlab="p-value"))
  dev.off()
  write.table(res_table, paste0("./Results/2.DE.group/2.1.LIRA.vs.CTRL//clust.",j,".lira.vs.ctrl.txt"),sep="\t",quote=F,col.names = NA)
}


