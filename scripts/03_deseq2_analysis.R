# ============================================================
# Differential expression analysis (DESeq2)
# ER+ BRCA2 mutated vs wild-type tumors
#
# Steps
# 1. Load processed count matrix
# 2. Restrict to ER+ tumors
# 3. Build DESeq2 dataset
# 4. Filter low expression genes
# 5. Run DESeq2
# 6. Map ENSEMBL IDs to gene symbols
# 7. Save differential expression results
#
# Output
# results/deseq2_brca2_deg.csv
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
})

source("scripts/00_utils.R")

# ------------------------------------------------------------
# Load processed data
# ------------------------------------------------------------

counts <- readRDS("data_processed/counts_unique.rds")
sample_annot <- readRDS("data_processed/sample_annot.rds")

# Ensure matching order
stopifnot(all(colnames(counts) == sample_annot$aliquot))

# ------------------------------------------------------------
# Restrict to ER+ tumors
# ------------------------------------------------------------

sample_annot_er <- sample_annot %>%
  filter(ER == "ER+")

counts_er <- counts[, sample_annot_er$aliquot]

cat("ER+ samples:", ncol(counts_er), "\n")

# ------------------------------------------------------------
# Build DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_er,
  colData = sample_annot_er,
  design = ~ BRCA2_Status
)

# ------------------------------------------------------------
# Low expression filtering
# ------------------------------------------------------------

dds <- dds[rowSums(counts(dds)) > 20, ]

cat("Genes after filtering:", nrow(dds), "\n")

# ------------------------------------------------------------
# Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)

res <- results(dds)

# ------------------------------------------------------------
# Convert to dataframe
# ------------------------------------------------------------

res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)

# ------------------------------------------------------------
# Map ENSEMBL → gene symbols
# ------------------------------------------------------------

res_df$SYMBOL <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = res_df$ENSEMBL,
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

# ------------------------------------------------------------
# Clean results
# ------------------------------------------------------------

res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

dir.create("results", showWarnings = FALSE)

write_csv(
  res_df,
  "results/deseq2_brca2_mut_vs_wt_ERpos.csv"
)

saveRDS(
  res_df,
  "results/deseq2_brca2_mut_vs_wt_ERpos.rds"
)

cat("Results saved to results/ directory\n")
