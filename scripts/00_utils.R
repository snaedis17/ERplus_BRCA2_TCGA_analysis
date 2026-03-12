library(dplyr)

norm12 <- function(x) substr(x, 1, 12)
norm16 <- function(x) substr(x, 1, 16)

clean_results <- function(res) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df %>% filter(!is.na(padj)) %>% arrange(padj)
  return(df)
}