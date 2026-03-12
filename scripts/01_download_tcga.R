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