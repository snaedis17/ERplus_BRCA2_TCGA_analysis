library(DESeq2)
library(dplyr)

source("scripts/00_utils.R")

counts <- readRDS("data_processed/counts_patient.rds")
annot <- readRDS("data_processed/sample_annotation.rds")

# match annotation
annot <- annot %>%
  distinct(patient_id, .keep_all = TRUE)

annot <- annot[match(colnames(counts), annot$patient_id), ]

# subset ER+
keep <- annot$ER_status == "ER+"

counts <- counts[, keep]
annot <- annot[keep, ]

annot$group <- annot$BRCA2_status

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = annot,
  design = ~ group
)

dds <- DESeq(dds)

res <- results(dds)

res_clean <- clean_results(res)

write.csv(
  res_clean,
  "results/deseq/BRCA2mut_vs_WT_ERpos.csv",
  row.names = FALSE
)




table(table(aliquot2pat$patient12))
