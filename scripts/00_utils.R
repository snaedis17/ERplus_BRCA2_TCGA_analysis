library(dplyr)

# Normalize TCGA IDs
norm12 <- function(x) {
  substr(toupper(trimws(as.character(x))), 1, 12)
}

norm16 <- function(x) {
  substr(toupper(trimws(as.character(x))), 1, 16)
}

# Clean DESeq2 results
clean_results <- function(res) {

  df <- as.data.frame(res)

  df$gene <- rownames(df)

  df <- df %>%
    filter(!is.na(padj)) %>%
    arrange(padj)

  return(df)
}
