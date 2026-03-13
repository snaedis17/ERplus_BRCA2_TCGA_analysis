library(DESeq2)
library(dplyr)

source("scripts/00_utils.R")

# load processed data
counts <- readRDS("data_processed/counts_unique.rds")
sample_annot <- readRDS("data_processed/sample_annot.rds")

# build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_annot,
  design = ~ BRCA2_Status
)

# filter low expression
dds <- dds[rowSums(counts(dds)) > 20, ]

# run DESeq2
dds <- DESeq(dds)

res <- results(dds)

# clean results
res_clean <- clean_results(res)

# save
write.csv(res_clean,
          "results/deseq2_brca2_deg.csv",
          row.names = FALSE)
