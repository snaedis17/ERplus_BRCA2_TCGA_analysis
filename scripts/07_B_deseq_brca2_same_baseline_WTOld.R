# ============================================================
# deseq_brca2_same_baseline_WTOld.R
# PURPOSE:
# Use one common baseline (WT_Old) for all ER+ BRCA2/age comparisons
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
  library(ggplot2)
  library(ggrepel)
})

filter  <- dplyr::filter
select  <- dplyr::select
arrange <- dplyr::arrange
mutate  <- dplyr::mutate

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/same_baseline_WTOld"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
counts_one    <- readRDS("data/counts_ERpos_clean_final.rds")
annot_patient <- readRDS("data/annot_ERpos_clean_final.rds")

stopifnot(identical(colnames(counts_one), annot_patient$patient12))

# ------------------------------------------------------------
# Metadata
# IMPORTANT: WT_Old is the baseline
# ------------------------------------------------------------
annot_patient <- annot_patient %>%
  mutate(
    age = case_when(
      !is.na(age_num) ~ age_num,
      !is.na(age_at_index) ~ age_at_index,
      TRUE ~ NA_real_
    ),
    age = if_else(age > 90, 90, age),
    age_group = if_else(age <= 40, "Young", "Old"),
    age_group = factor(age_group, levels = c("Old", "Young")),
    BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut")),
    grp4 = paste(BRCA2_Status, age_group, sep = "_"),
    grp4 = factor(grp4, levels = c("WT_Old", "WT_Young", "Mut_Old", "Mut_Young"))
  ) %>%
  filter(!is.na(age), !is.na(grp4))

counts_one <- counts_one[, annot_patient$patient12, drop = FALSE]
stopifnot(identical(colnames(counts_one), annot_patient$patient12))

cat("\n================ GROUP SUMMARY ================\n")
print(table(annot_patient$grp4))
#WT_Old  WT_Young   Mut_Old Mut_Young 
#690        45        12         5 

cat("==============================================\n")

# ------------------------------------------------------------
# Filter low expression
# ------------------------------------------------------------
keep <- rowSums(counts_one >= 10) >= 5
counts_one <- counts_one[keep, , drop = FALSE]

cat("Genes kept:", nrow(counts_one), "\n")
#Genes kept: 17913 

# ------------------------------------------------------------
# DESeq2 with common baseline WT_Old
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_one,
  colData   = annot_patient,
  design    = ~ grp4
)

dds <- DESeq(dds)

saveRDS(dds, file.path(results_dir, "dds_grp4_WTOld_baseline.rds"))

cat("\nModel coefficients:\n")
print(resultsNames(dds))

# ------------------------------------------------------------
# Extract results
# All relative to WT_Old when possible
# ------------------------------------------------------------
res_list_raw <- list(
  WTYoung_vs_WTOld = results(dds, contrast = c("grp4", "WT_Young", "WT_Old")),
  MutOld_vs_WTOld  = results(dds, contrast = c("grp4", "Mut_Old",  "WT_Old")),
  MutYoung_vs_WTOld = results(dds, contrast = c("grp4", "Mut_Young", "WT_Old")),
  MutYoung_vs_WTYoung = results(dds, contrast = c("grp4", "Mut_Young", "WT_Young")),
  MutYoung_vs_MutOld  = results(dds, contrast = c("grp4", "Mut_Young", "Mut_Old"))
)

# ------------------------------------------------------------
# Shrinkage
# apeglm for direct coefficients
# ashr for arbitrary contrasts
# ------------------------------------------------------------
res_list_shrunk <- list(
  WTYoung_vs_WTOld = lfcShrink(
    dds,
    coef = "grp4_WT_Young_vs_WT_Old",
    type = "apeglm"
  ),
  
  MutOld_vs_WTOld = lfcShrink(
    dds,
    coef = "grp4_Mut_Old_vs_WT_Old",
    type = "apeglm"
  ),
  
  MutYoung_vs_WTOld = lfcShrink(
    dds,
    coef = "grp4_Mut_Young_vs_WT_Old",
    type = "apeglm"
  ),
  
  MutYoung_vs_WTYoung = lfcShrink(
    dds,
    contrast = c("grp4", "Mut_Young", "WT_Young"),
    type = "ashr"
  ),
  
  MutYoung_vs_MutOld = lfcShrink(
    dds,
    contrast = c("grp4", "Mut_Young", "Mut_Old"),
    type = "ashr"
  )
)

# ------------------------------------------------------------
# Annotation
# ------------------------------------------------------------
add_symbols <- function(res_obj) {
  df <- as.data.frame(res_obj)
  
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

res_list_raw    <- lapply(res_list_raw, add_symbols)
res_list_shrunk <- lapply(res_list_shrunk, add_symbols)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
for (nm in names(res_list_raw)) {
  saveRDS(
    res_list_raw[[nm]],
    file.path(results_dir, paste0("res_", nm, "_RAW.rds"))
  )
  
  saveRDS(
    res_list_shrunk[[nm]],
    file.path(results_dir, paste0("res_", nm, "_SHRUNK.rds"))
  )
  
  write.csv(
    res_list_raw[[nm]],
    file.path(results_dir, paste0("res_", nm, "_RAW.csv")),
    row.names = FALSE
  )
  
  write.csv(
    res_list_shrunk[[nm]],
    file.path(results_dir, paste0("res_", nm, "_SHRUNK.csv")),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# DEG count helper
# ------------------------------------------------------------
count_sig <- function(df) {
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(sig = padj < 0.05 & abs(log2FoldChange) >= 1) %>%
    summarise(
      total = sum(sig),
      up = sum(sig & log2FoldChange > 0),
      down = sum(sig & log2FoldChange < 0)
    )
}

cat("\n===== DEG COUNTS (SHRUNK) =====\n")
for (nm in names(res_list_shrunk)) {
  cat("\n", nm, "\n", sep = "")
  print(count_sig(res_list_shrunk[[nm]]))
}

#WTYoung_vs_WTOld
#total up down
#   272 54  218

#MutOld_vs_WTOld
#total up down
#    34 24   10

#MutYoung_vs_WTOld
#total up down
#    46 35   11

#MutYoung_vs_WTYoung
#total up down
#    25 19    6

#MutYoung_vs_MutOld
#total up down
#    23 16    7
# ------------------------------------------------------------
# Volcano plot helper
# ------------------------------------------------------------
plot_volcano <- function(res_df, title, filename) {
  df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neglog10 = -log10(padj),
      direction = case_when(
        padj < 0.05 & log2FoldChange > 1  ~ "Up",
        padj < 0.05 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  top_genes <- df %>%
    filter(direction != "NS", !is.na(SYMBOL), SYMBOL != "") %>%
    arrange(padj) %>%
    slice_head(n = 20)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10)) +
    geom_point(aes(color = direction), alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey80")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      data = top_genes,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = 20
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = title,
      x = "log2 fold change (shrunk)",
      y = "-log10 adjusted p-value"
    )
  
  print(p)
  
  ggsave(
    file.path(results_dir, filename),
    p,
    width = 6,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}


# ------------------------------------------------------------
# Plot key baseline comparisons
# ------------------------------------------------------------
plot_volcano(
  res_list_shrunk$WTYoung_vs_WTOld,
  "WT Young vs WT Old",
  "volcano_WTYoung_vs_WTOld.png"
)

plot_volcano(
  res_list_shrunk$MutOld_vs_WTOld,
  "Mut Old vs WT Old",
  "volcano_MutOld_vs_WTOld.png"
)

plot_volcano(
  res_list_shrunk$MutYoung_vs_WTOld,
  "Mut Young vs WT Old",
  "volcano_MutYoung_vs_WTOld.png"
)

plot_volcano(
  res_list_shrunk$MutYoung_vs_WTYoung,
  "Mut Young vs WT Young",
  "volcano_MutYoung_vs_WTYoung.png"
)

plot_volcano(
  res_list_shrunk$MutYoung_vs_MutOld,
  "Mut Young vs Mut Old",
  "volcano_MutYoung_vs_MutOld.png"
)

cat("\nDONE\n")
