# ============================================================
# WT Young vs Old
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  library(AnnotationDbi)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})
filter  <- dplyr::filter
select  <- dplyr::select
arrange <- dplyr::arrange
# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/young"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
counts_one    <- readRDS("data/counts_ERpos_clean_final.rds")
annot_patient <- readRDS("data/annot_ERpos_clean_final.rds")

# ------------------------------------------------------------
# Keep only BRCA2 WT
# ------------------------------------------------------------
annot_patient <- annot_patient %>%
  filter(BRCA2_Status == "WT")

# ------------------------------------------------------------
# Define age groups
# ------------------------------------------------------------
annot_patient <- annot_patient %>%
  filter(!is.na(age_num)) %>%
  mutate(
    age_group = ifelse(age_num <= 40, "Young", "Old"),
    age_group = factor(age_group, levels = c("Old", "Young"))
  )

cat("Age groups:\n")
print(table(annot_patient$age_group))
#Old Young 
#663    41 
# ------------------------------------------------------------
# Sync counts
# ------------------------------------------------------------
counts_one <- counts_one[, annot_patient$patient12, drop = FALSE]
stopifnot(identical(colnames(counts_one), annot_patient$patient12))

# ------------------------------------------------------------
# Build coldata
# ------------------------------------------------------------
coldata <- data.frame(grp = annot_patient$age_group)
rownames(coldata) <- annot_patient$patient12

# ------------------------------------------------------------
# Filtering (edgeR best practice)
# ------------------------------------------------------------
dge <- DGEList(counts = counts_one)

keep <- filterByExpr(dge, group = coldata$grp)
counts_one <- counts_one[keep, , drop = FALSE]

cat("Genes after filtering:", nrow(counts_one), "\n")
#Genes after filtering: 17038 
# ------------------------------------------------------------
# DESeq2
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_one,
  colData   = coldata,
  design    = ~ grp
)

dds <- DESeq(dds)

# ------------------------------------------------------------
# Results
# ------------------------------------------------------------
res_raw <- results(dds, contrast = c("grp", "Young", "Old"))

res_shrunk <- lfcShrink(
  dds,
  coef = "grp_Young_vs_Old",
  type = "apeglm"
)

# ------------------------------------------------------------
# Annotate genes
# ------------------------------------------------------------
annotate_res <- function(res) {
  df <- as.data.frame(res)
  
  df$ENSEMBL <- gsub("\\..*", "", rownames(df))
  
  df$SYMBOL <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = df$ENSEMBL,
    keytype = "GENEID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  df %>% relocate(SYMBOL, ENSEMBL)
}

res_raw_df    <- annotate_res(res_raw)
res_shrunk_df <- annotate_res(res_shrunk)

# ------------------------------------------------------------
# COMBINE RAW + SHRUNK
# ------------------------------------------------------------
full_results <- data.frame(
  ENSEMBL = gsub("\\..*", "", rownames(res_raw)),
  SYMBOL  = res_raw_df$SYMBOL,
  
  baseMean = res_raw$baseMean,
  
  log2FC_raw   = res_raw$log2FoldChange,
  lfcSE        = res_raw$lfcSE,
  stat         = res_raw$stat,
  pvalue       = res_raw$pvalue,
  padj         = res_raw$padj,
  
  log2FC_shrunk = res_shrunk$log2FoldChange,
  
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(SYMBOL))

# ------------------------------------------------------------
# Significant genes (STRICT)
# ------------------------------------------------------------
sig_results <- full_results %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FC_shrunk) >= 1
  ) %>%
  arrange(padj)

cat("Significant genes:", nrow(sig_results), "\n")
#273
# ------------------------------------------------------------
# VOLCANO PLOT
# ------------------------------------------------------------
volcano_df <- full_results %>%
  filter(!is.na(padj), !is.na(log2FC_shrunk)) %>%
  mutate(
    negLogP = -log10(padj),
    category = case_when(
      padj < 0.05 & log2FC_shrunk > 1  ~ "Up",
      padj < 0.05 & log2FC_shrunk < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

top_genes <- volcano_df %>%
  filter(category != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 15)

p_volcano <- ggplot(volcano_df, aes(x = log2FC_shrunk, y = negLogP, color = category)) +
  geom_point(alpha = 0.7) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 20
  ) +
  
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NS" = "grey70"
  )) +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  theme_bw(base_size = 14) +
  
  labs(
    title = "WT Young vs Old",
    x = "log2 fold change (shrunk)",
    y = "-log10 adjusted p-value"
  )

print(p_volcano)

ggsave(
  file.path(results_dir, "volcano_young_vs_old.png"),
  p_volcano,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# PCA (VST)
# ------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

saveRDS(vst_mat, "data/vst_young_vs_old.rds")
saveRDS(coldata, "data/coldata_young_vs_old.rds")

pca <- prcomp(t(vst_mat), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = coldata$grp
)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA – WT Young vs Old",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)")
  )

print(p_pca)

ggsave(
  file.path(results_dir, "PCA_young_vs_old.png"),
  p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# SAVE ALL RESULTS
# ------------------------------------------------------------
saveRDS(res_raw_df, file.path(results_dir, "Young_vs_Old_raw.rds"))
saveRDS(res_shrunk_df, file.path(results_dir, "Young_vs_Old_shrunk.rds"))
saveRDS(full_results, file.path(results_dir, "Young_vs_Old_full.rds"))
saveRDS(dds, file.path(results_dir, "dds_young_vs_old.rds"))

write.csv(res_raw_df,
          file.path(results_dir, "Young_vs_Old_raw.csv"),
          row.names = FALSE)

write.csv(res_shrunk_df,
          file.path(results_dir, "Young_vs_Old_shrunk.csv"),
          row.names = FALSE)

write.csv(full_results,
          file.path(results_dir, "Young_vs_Old_full.csv"),
          row.names = FALSE)

write.csv(sig_results,
          file.path(results_dir, "Young_vs_Old_significant_strict.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# QC: raw vs shrunk comparison
# ------------------------------------------------------------
plot(res_raw$log2FoldChange, res_shrunk$log2FoldChange,
     xlab = "Raw LFC",
     ylab = "Shrunk LFC",
     main = "Shrinkage effect")
abline(0,1,col="red")
