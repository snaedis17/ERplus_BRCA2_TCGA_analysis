library(dplyr)
library(SummarizedExperiment)

source("scripts/00_utils.R")

# Load TCGA object
se <- readRDS("data_raw/tcga_brca_se.rds")

coldata <- as.data.frame(colData(se))

# patient id
coldata$patient_id <- norm12(coldata$barcode)

# ER status
coldata$ER_status <- ifelse(
  coldata$paper_ER_status == "Positive",
  "ER+",
  "ER-"
)

# age
coldata$age <- coldata$age_at_index

# load BRCA2 mutation list
brca2 <- read.csv("data_raw/external/Maxwell_BRCA2uppl.csv")

brca2$patient_id <- norm12(brca2$Tumor_Sample_Barcode)

# assign BRCA2 status
coldata$BRCA2_status <- ifelse(
  coldata$patient_id %in% brca2$patient_id,
  "BRCA2mut",
  "BRCA2wt"
)

saveRDS(coldata, "data_processed/sample_annotation.rds")