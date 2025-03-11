library(readxl)
library(tidyverse)
library(ggplot2)
column_names <- c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9")
all_data <- list()
for (i in 1:9) {
  data <- read_excel("C:/Users/mchiesa/Desktop/filePlotCarla/Summary.GSEA.GO.xlsx", sheet = i)
  crescente_data<- arrange(data, NES)
  decrescente_data <- arrange(data, desc(NES))
  
  down_regolato <- crescente_data[1:10,]
  up_regolato <- decrescente_data[1:10,]
  merge_regulation <- rbind(up_regolato,down_regolato)
  merge_regulation$ID <- column_names[i] 
  all_data[[i]] <- merge_regulation
}


final_data <- bind_rows(all_data)

  
name_GO <- unique(final_data$GOID)
counter <- numeric(length(name_GO))
names(counter) <- name_GO 


for (i in seq_along(name_GO)) {
  count <- 0
  for (j in seq_along(final_data$GOID)) {
    if (name_GO[i] == final_data$GOID[j]) {
      count <- count + 1
    }
  }
  counter[i] <- count  
}


counter_data <- data.frame(GOID = name_GO, counter = counter)

counter_data <- arrange(counter_data,desc(counter))

ggplot(counter_data, aes(x= reorder(GOID, counter), y=counter)) + 
  geom_bar(aes(fill=as.factor(counter)), width=0.5,alpha=0.6,stat = "identity") +
  coord_flip() +
   expand_limits(x = 0, y = 0)+
  labs(title = "GO_ID Histogram", x = "GO_ID", y = "Count") +
  theme_minimal() +
theme(
  axis.text.y = element_text(size = 6),
  plot.title = element_text(hjust = 0.5, face = "bold") 
)

