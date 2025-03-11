rm(list = ls())
graphics.off()
cat("\014")
setwd("C:/Users/mchiesa/Desktop/2.1.LIRA.vs.CTRL/CCC")

library(ggplot2)
library(tibble)
library(tidyverse)


dati <- read.delim(file = "xBubblePlot.txt",na.strings = "")

pathw_unici <-  unique(dati$pathways)
dati_res <- as.data.frame(matrix(nrow = length(pathw_unici),
                                 ncol = 2*(ncol(dati)-1)))

rownames(dati_res) <- pathw_unici

colnames(dati_res)[c(TRUE,FALSE)] <- paste0(colnames(dati[2:ncol(dati)]),"_UP") #colonne dispari
colnames(dati_res)[c(FALSE,TRUE)] <- paste0(colnames(dati[2:ncol(dati)]),"_DOWN") #colonne dispari

dati_res[is.na(dati_res)] <- 0

kol <- 1
i=4 # Down #i=16 UP

for (j in 1:(ncol(dati)-1)){
  
  df_sel <- dati[,c(1,j+1)]
  
  for (i in 1:length(pathw_unici)){
    
    n_path_down <- length(which(df_sel$pathways %in% pathw_unici[i] & df_sel[,2] %in% "DOWN"))
    n_path_up <- length(which(df_sel$pathways %in% pathw_unici[i] & df_sel[,2] %in% "UP"))
    
    dati_res[i,kol] <- n_path_up
    dati_res[i,kol+1] <- n_path_down
  }
  kol <- kol +2
  
  
}

#Plot
df.res_plot<- dati_res %>%
  rownames_to_column(var = "id") %>%
  gather(key, counts, -id,na.rm = T)


df.res_plot %>%
  ggplot(aes(key, id)) +
  geom_point(aes(size = counts, colour=key, fill =key),alpha=0.5) +
  scale_size(range = c(-0.5,8)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        axis.text.x = element_text(face="bold", color="black",
                                   size=4,angle = 90),
        axis.text.y = element_text(face="bold", color="black",
                                   size=10))+
  scale_color_manual(values=rep(c("blue","red"),81))


