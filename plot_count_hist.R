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



dati <- read.delim("xCount_table.txt",stringsAsFactors = T)
#dati <- read.delim("xCount_table.inverted_x_test.txt",stringsAsFactors = T)

ggplot(data=dati, aes(x=cluster, y=perc, fill=cell_group)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~cell_group)
