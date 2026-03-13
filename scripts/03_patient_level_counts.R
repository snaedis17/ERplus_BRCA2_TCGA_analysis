library(dplyr)

source("scripts/00_utils.R")

se <- readRDS("data_raw/tcga_brca_se.rds")
annot <- readRDS("data_processed/sample_annotation.rds")

counts <- assay(se)

sample_ids <- norm16(colnames(counts))
patient_ids <- norm12(sample_ids)

# collapse to patient level
counts_patient <- rowsum(
  t(counts),
  group = patient_ids
)

counts_patient <- t(counts_patient)

saveRDS(counts_patient, "data_processed/counts_patient.rds")