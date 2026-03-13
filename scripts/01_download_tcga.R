# ============================================================
# Download TCGA BRCA RNA-seq data
#
# Steps
# 1. Query TCGA BRCA RNA-seq STAR counts
# 2. Download data from GDC
# 3. Prepare SummarizedExperiment object
#
# Output
# data_raw/tcga_brca_se.rds
# ============================================================

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

se <- GDCprepare(query)

saveRDS(se, "data_raw/tcga_brca_se.rds")
