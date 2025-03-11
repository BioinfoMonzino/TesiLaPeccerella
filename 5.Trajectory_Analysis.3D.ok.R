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
library(slingshot)
library(umap)
#library(Rtsne)
library(slingshot)
library(Hmisc)
library(rgl)

setwd("C:/Users/mchiesa/Desktop/2.1.LIRA.vs.CTRL/")
source("./0.Script/custom_seurat_functions.R")


#load("C:/Users/mchiesa/Desktop/liraglut_4_samples/analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.RData")
load("C:/Users/mchiesa/Desktop/2.1.LIRA.vs.CTRL/analisi_Cri.4.gruppi.k.30.3D.reso.06.Final.RData")


dimred <- cbind(cd34_integ@reductions[["umap"]]@cell.embeddings[,1],
                cd34_integ@reductions[["umap"]]@cell.embeddings[,2],
                cd34_integ@reductions[["umap"]]@cell.embeddings[,3])
metadata_seu <- cd34_integ@meta.data
clustering <- as.factor(metadata_seu$seurat_clusters)

rm(cd34_integ)
gc()



set.seed(42)
pto_3D <- slingshot(data = dimred,
                    clusterLabels = clustering,
                    dist.method	= "mnn",
                    #dist.method	= "slingshot",
                    #end.clus = c("8","2","3","9","14","7"), 
                    start.clus = "10",
                    spar =1.1)

################# Plot inferred trajectories on 3D-UMAP:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycolor <- gg_color_hue(19) #,veri,colori,dei,cluster

color_vec <-c()
clust_num <- as.numeric(clustering)
for (i in 1:length(clust_num)){
  color_vec[i] <- mycolor[clust_num[i]+1]
}


plot3d(dimred,
       aspect = 'iso',
       size=1,
       col = color_vec
)

plot3d.SlingshotDataSet(SlingshotDataSet(pto_3D), add = TRUE, type='c', lwd=3)


# 
# plot3d(dimred,
#        aspect = 'iso',
#        size=1,
#        col = color_vec
# )
# 
# plot3d.SlingshotDataSet(SlingshotDataSet(pto_3D), add = TRUE, type='l', lwd=3)
# 
# 
# 
# plot3d(dimred,
#        aspect = 'iso',
#        size=1,
#        col = color_vec
# )
# 
# plot3d.SlingshotDataSet(SlingshotDataSet(pto_3D), add = TRUE, type='b', lwd=3)

# #############
# rifaccio con la distanza slingshot per la stima dello pseudotime: è + sensato
pto_3D_sl <- slingshot(data = dimred,
                       clusterLabels = clustering,
                       #dist.method	= "mnn",
                       dist.method	= "slingshot",
                       #end.clus = c("8","3","11","14","7","17"), 
                       start.clus = "10",
                       spar =1.1)



pseudotime_lineages <- as.data.frame(pto_3D_sl@assays@data@listData[["pseudotime"]])
pseudot<- as.data.frame(rowMeans(as.matrix(pseudotime_lineages),na.rm = T))
colnames(pseudot) <- "pseudot"

plot3d(dimred,
       aspect = 'iso',
       size=1,
       col = alpha(viridisLite::viridis(100)[cut(pseudot$pseudot, breaks=100)])
)
plot3d.SlingshotDataSet(SlingshotDataSet(pto_3D), add = TRUE, type='c', lwd=3)

xx <-SlingshotDataSet(pto_3D)

####################################################################################
# controlla i lineages e e poi media quelli consecutivi
i=1
#for (i in 1:length(colnames(pseudotime_lineages[,lineage_n]))){

lineage_n <- i

colori_lineage <- viridisLite::magma(100)[cut(pseudotime_lineages[,lineage_n], breaks=100)]
colori_lineage[which(is.na(colori_lineage))] <- "#008000"

plot3d(dimred,
       aspect = 'iso',
       size=1,
       col = colori_lineage
)
#}

idx_sampl_ok_lira <-  which(metadata_seu$drug %in% "lira")
idx_sampl_ok_ctrl <-  which(metadata_seu$drug %in% "ctrl")

pst_ctrls <- pseudotime_lineages[idx_sampl_ok_ctrl,]
pst_liras <- pseudotime_lineages[idx_sampl_ok_lira,]

pst_ctrls$drug <- "ctrl"
pst_liras$drug <- "lira"

df_pts <- rbind(pst_ctrls,pst_liras)


###########################
# lineage_macrofagi
df_pts_macrof <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(13,9,8)]),na.rm = T))
df_pts_macrof$drug <- df_pts$drug
df_pts_macrof[which(is.infinite(df_pts_macrof[,1])),1] <-NA
df_pts_macrof <- df_pts_macrof[complete.cases(df_pts_macrof),,drop=F]
colnames(df_pts_macrof)[1] <- "Lineage"

# lineage_neutro
df_pts_neutro <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(11,4,2,5)]),na.rm = T))
df_pts_neutro$drug <- df_pts$drug
df_pts_neutro[which(is.infinite(df_pts_neutro[,1])),1] <-NA
df_pts_neutro <- df_pts_neutro[complete.cases(df_pts_neutro),,drop=F]
colnames(df_pts_neutro)[1] <- "Lineage"

# lineage_meutro_baso
df_pts_neutro_baso <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(11,4,2,5,6,3,1)]),na.rm = T))
df_pts_neutro_baso$drug <- df_pts$drug
df_pts_neutro_baso[which(is.infinite(df_pts_neutro_baso[,1])),1] <-NA
df_pts_neutro_baso <- df_pts_neutro_baso[complete.cases(df_pts_neutro_baso),,drop=F]
colnames(df_pts_neutro_baso)[1] <- "Lineage"

# lineage_plt
df_pts_plt <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(12)]),na.rm = T))
df_pts_plt$drug <- df_pts$drug
df_pts_plt[which(is.infinite(df_pts_plt[,1])),1] <-NA
df_pts_plt <- df_pts_plt[complete.cases(df_pts_plt),,drop=F]
colnames(df_pts_plt)[1] <- "Lineage"


# 
# ###########################
# # lineage_macrofagi
# df_pts_macrof <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(13,11,9,8,10)]),na.rm = T))
# df_pts_macrof$drug <- df_pts$drug
# df_pts_macrof[which(is.infinite(df_pts_macrof[,1])),1] <-NA
# df_pts_macrof <- df_pts_macrof[complete.cases(df_pts_macrof),,drop=F]
# colnames(df_pts_macrof)[1] <- "Lineage"
# 
# # lineage_neutro
# df_pts_neutro <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(13,11,4,2,5)]),na.rm = T))
# df_pts_neutro$drug <- df_pts$drug
# df_pts_neutro[which(is.infinite(df_pts_neutro[,1])),1] <-NA
# df_pts_neutro <- df_pts_neutro[complete.cases(df_pts_neutro),,drop=F]
# colnames(df_pts_neutro)[1] <- "Lineage"
# 
# # lineage_meutro_baso
# df_pts_neutro_baso <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(13,11,4,2,5,6,3,1)]),na.rm = T))
# df_pts_neutro_baso$drug <- df_pts$drug
# df_pts_neutro_baso[which(is.infinite(df_pts_neutro_baso[,1])),1] <-NA
# df_pts_neutro_baso <- df_pts_neutro_baso[complete.cases(df_pts_neutro_baso),,drop=F]
# colnames(df_pts_neutro_baso)[1] <- "Lineage"
# 
# # lineage_plt
# df_pts_plt <- as.data.frame(rowMaxs( as.matrix(df_pts[,c(13,11,12)]),na.rm = T))
# df_pts_plt$drug <- df_pts$drug
# df_pts_plt[which(is.infinite(df_pts_plt[,1])),1] <-NA
# df_pts_plt <- df_pts_plt[complete.cases(df_pts_plt),,drop=F]
# colnames(df_pts_plt)[1] <- "Lineage"
# 





df_pts_plot <- df_pts_macrof
df_pts_plot <- df_pts_neutro
df_pts_plot <- df_pts_neutro_baso
df_pts_plot <- df_pts_plt



ggplot(data = df_pts_plot, aes(y = Lineage, x = drug)) + 
  geom_point(aes(y = Lineage, x = drug,color=drug),
             size=4,
             position = position_jitter(width =0.3),
             alpha=0.01) +
  geom_violin(aes( fill=drug), alpha=0.7,draw_quantiles = c(0.25,0.5,0.75),trim = F) +
  
  xlab("drug") +
  ylab("Pseudotime") +
  theme_bw()+
  theme(
    legend.position = "none"
  )
  # )+
  # ylim(0,30)
gc()
# 
# library(ggridges)
# # library(ggplot2)
# library(viridis)
# library(hrbrthemes)
# 
# ggplot(df_pts_plot, aes(x = Lineage, y = drug, fill = ..x..)) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   scale_fill_viridis(name = "Temp. [F]", option = "C") +
#   labs(title = 'Temperatures in Lincoln NE in 2016') +
#   theme_ipsum() +
#   theme(
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     strip.text.x = element_text(size = 8)
#   )
# 
# ggplot(df_pts_plot, aes(x = Lineage, y = drug, fill = drug)) +
#   geom_density_ridges(aes(
#     point_color = drug,
#     point_fill = drug),
#     #alpha = .2, 
#     jittered_points = TRUE,
#     alpha=0.5,point_alpha=0.05) +
#   scale_point_color_hue(l = 60) +
#   #annotation_logticks(scaled = TRUE,base = exp(1)) +
#   theme_ridges() + 
#   theme(
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     strip.text.x = element_text(size = 8)
#   ) 
# #stat_density_ridges(quantile_lines = TRUE, quantiles=c(0.025,0.5,0.975)) +



#### DE pseudotime
fit.2 <- lm(Lineage ~ drug, data = df_pts_plot)
summ_2 <- summary(fit.2)
summ_2

t.test(log(Lineage+1,2) ~ drug, data = df_pts_plot)
