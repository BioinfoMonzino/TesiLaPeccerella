library("xlsx")
##### LOAD DATA
data_ctrl <- read.delim("/Users/carla/TESI/CCC_All_cells/CTRL/pvalues.txt")
data_lira <- read.delim("/Users/carla/TESI/CCC_All_cells/LIRA/pvalues.txt")


compared_pvalue <- function(data_ctrl,data_lira) {
  df <- merge(data_ctrl,data_lira, by = c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","directionality","classification"), all = TRUE,suffixes = c(".ctrl", ".lira"))

  column_ctrl <- grep("\\.ctrl$", names(df), value = TRUE)
  column_lira <- grep("\\.lira$", names(df), value = TRUE)
  
  for(i in seq_along(column_ctrl)){
    col_ctrl <- column_ctrl[i]
    col_lira <- column_lira[i]
    col_result <- sub("\\.ctrl$", ".result", col_ctrl)
    col_regulation <- sub("\\.ctrl$", ".regulation", col_ctrl)
    
    df[[col_result]] <- ifelse(!is.na(df[[col_ctrl]]) & !is.na(df[[col_lira]]) &
                                 (df[[col_ctrl]] > 0.05 & df[[col_lira]] > 0.05), 0,
                               ifelse(!is.na(df[[col_ctrl]]) & !is.na(df[[col_lira]]) &
                                        (df[[col_ctrl]] < 0.05 & df[[col_lira]] <= 0.05), 0,
                                      ifelse(!is.na(df[[col_ctrl]]) & !is.na(df[[col_lira]]) &
                                               (df[[col_ctrl]] <= 0.05 & df[[col_lira]] > 0.05), 1,
                                             ifelse(!is.na(df[[col_ctrl]]) & !is.na(df[[col_lira]]) &
                                                      (df[[col_ctrl]] > 0.05 & df[[col_lira]] <= 0.05), 1, 0))))
     
    df[[col_regulation]] <- ifelse(df[[col_result]] == 1 & df[[col_ctrl]] <= 0.05, "DOWN",
                                  ifelse(df[[col_result]] == 1 & df[[col_lira]] <= 0.05, "UP", 0))
  } 
  
  
  result_columns <- grep("\\.result$", names(df), value = TRUE)
  regulation_columns <- grep("\\.regulation$", names(df), value = TRUE)
  other_columns <- setdiff(names(df), c(result_columns, regulation_columns))
  
  
  df <- df[, c(other_columns, result_columns, regulation_columns)]
  return(df)
}

df <- compared_pvalue(data_ctrl,data_lira)
write.csv(df, "/Users/carla/TESI/CCC_All_cells/df1.csv", row.names = FALSE)
write.xlsx(df, "/Users/carla/TESI/CCC_All_cells/df.xlsx")

