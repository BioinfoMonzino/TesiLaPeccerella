rm(list = ls())
graphics.off()
cat("\014")

setwd("C:/Users/mchiesa/Desktop/2.1.LIRA.vs.CTRL/")

library(ggplot2)


dati_lolli <- read.delim("xLollipop.txt",stringsAsFactors = TRUE)

# ggplot(dati_lolli, aes(NES, goid, color = regulation)) +
#   geom_segment(aes(x = 0, y = goid, xend = NES, yend = goid),color="black") +
#   scale_size(range = c(1,5)) +
#   geom_point(aes(size = abs(NES),alpha = -log(p.value+0.0001,10))) +
#   scale_color_manual(values=c("blue","red")) +
#   facet_wrap(~cell, scales = "free_y") +
#   theme_bw() 
  
  
ggplot(dati_lolli, aes(NES, nome2, color = regulation)) +
  geom_segment(aes(x = 0, y = nome2, xend = NES, yend = nome2),color="black") +
  scale_size(range = c(1,5)) +
  geom_point(aes(size = abs(NES),alpha = -log(p.value+0.0001,10))) +
  scale_color_manual(values=c("blue","red")) +
  facet_wrap(~cell, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")+
  xlim(-4,4)
  

