# ============================================================
# Differential expression analysis (DESeq2)
# ER+ BRCA2 mutated vs wild-type tumors
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
# Low count filtering
# ------------------------------------------------------------

dds <- dds[rowSums(counts(dds)) > 10, ]

cat("Genes after count filtering:", nrow(dds), "\n")

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

# remove ENSG version numbers (important for annotation)
res_df$ENSEMBL <- sub("\\..*$", "", res_df$ENSEMBL)

# ------------------------------------------------------------
# Map gene symbols
# ------------------------------------------------------------

res_df$SYMBOL <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = res_df$ENSEMBL,
  column = "SYMBOL",
  keytype = "GENEID",
  multiVals = "first"
)

# ------------------------------------------------------------
# Map gene biotype
# ------------------------------------------------------------

res_df$biotype <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = res_df$ENSEMBL,
  column = "GENEBIOTYPE",
  keytype = "GENEID",
  multiVals = "first"
)

# ------------------------------------------------------------
# Clean results
# ------------------------------------------------------------

res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

cat("Genes with valid statistics:", nrow(res_df), "\n")

# ------------------------------------------------------------
# Protein-coding subset
# ------------------------------------------------------------

res_pc <- res_df %>%
  filter(biotype == "protein_coding")

cat("Protein-coding genes:", nrow(res_pc), "\n")

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

dir.create("results", showWarnings = FALSE)

write_csv(
  res_df,
  "results/deseq2_all_genes_ERpos_BRCA2.csv"
)

write_csv(
  res_pc,
  "results/deseq2_protein_coding_ERpos_BRCA2.csv"
)

saveRDS(
  res_df,
  "results/deseq2_all_genes_ERpos_BRCA2.rds"
)

saveRDS(
  res_pc,
  "results/deseq2_protein_coding_ERpos_BRCA2.rds"
)

cat("DESeq2 analysis complete\n")
