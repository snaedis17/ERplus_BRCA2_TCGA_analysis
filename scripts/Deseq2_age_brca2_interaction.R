# ============================================================
Deseq2_age_brca2_interaction.R 
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

filter   <- dplyr::filter
select   <- dplyr::select
arrange  <- dplyr::arrange

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/age_interaction"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Helper: annotation
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

# ------------------------------------------------------------
# Volcano (clean)
# ------------------------------------------------------------
plot_volcano_clean <- function(df, title, file_name, use_shrunk = TRUE) {
  
  lfc_col <- ifelse(use_shrunk, "log2FC_shrunk", "log2FC_raw")
  
  plot_df <- df %>%
    filter(!is.na(padj), !is.na(.data[[lfc_col]])) %>%
    mutate(
      negLogP = -log10(padj),
      category = case_when(
        padj < 0.05 & .data[[lfc_col]] > 1  ~ "Up",
        padj < 0.05 & .data[[lfc_col]] < -1 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  top_genes <- plot_df %>%
    filter(category != "NS") %>%
    arrange(padj) %>%
    slice_head(n = 15)
  
  p <- ggplot(plot_df, aes(x = .data[[lfc_col]], y = negLogP, color = category)) +
    geom_point(alpha = 0.7, size = 2) +
    
    geom_text_repel(
      data = top_genes,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = 20
    ) +
    
    scale_color_manual(values = c(
      "Up" = "#d73027",
      "Down" = "#4575b4",
      "NS" = "grey80"
    )) +
    
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    
    theme_bw(base_size = 14) +
    
    labs(
      title = title,
      x = paste0("log2FC (", ifelse(use_shrunk, "shrunk", "raw"), ")"),
      y = "-log10 adjusted p-value"
    )
  
  print(p)
  
  ggsave(file.path(results_dir, file_name), p,
         width = 6, height = 5, dpi = 300, bg = "white")
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
counts_one <- readRDS("data/counts_ERpos_clean_final.rds")
annot_patient <- readRDS("data/annot_ERpos_clean_final.rds")

# ------------------------------------------------------------
# Preprocess
# ------------------------------------------------------------
annot_patient <- annot_patient %>%
  filter(!is.na(age_num)) %>%
  mutate(
    age_group = if_else(age_num <= 40, "Young", "Old"),
    age_group = factor(age_group, levels = c("Old", "Young")),
    BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut"))
  )

counts_one <- counts_one[, annot_patient$patient12]

# ------------------------------------------------------------
# Filtering
# ------------------------------------------------------------
dge <- DGEList(counts = counts_one)

keep <- filterByExpr(
  y = dge,
  group = interaction(annot_patient$age_group, annot_patient$BRCA2_Status)
)

counts_one <- counts_one[keep, ]
stopifnot(all(colnames(counts_one) == annot_patient$patient12))

cat("Genes kept:", nrow(counts_one), "\n")
#Genes kept: 17863 
# ------------------------------------------------------------
# DESeq2
# ------------------------------------------------------------

rownames(annot_patient) <- annot_patient$patient12

dds <- DESeqDataSetFromMatrix(
  countData = counts_one,
  colData   = annot_patient,
  design    = ~ age_group * BRCA2_Status
)

dds <- DESeq(dds)

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
res_age_wt      <- results(dds, name = "age_group_Young_vs_Old")
res_brca2_old   <- results(dds, name = "BRCA2_Status_Mut_vs_WT")
res_interaction <- results(dds, name = "age_groupYoung.BRCA2_StatusMut")

res_brca2_young <- results(
  dds,
  contrast = list(c("BRCA2_Status_Mut_vs_WT", "age_groupYoung.BRCA2_StatusMut"))
)

# ------------------------------------------------------------
# SHRINK (ONLY SIMPLE EFFECTS)
# ------------------------------------------------------------
res_age_wt_shrunk   <- lfcShrink(dds, coef = "age_group_Young_vs_Old", type = "apeglm")
res_brca2_old_shrunk <- lfcShrink(dds, coef = "BRCA2_Status_Mut_vs_WT", type = "apeglm")

# DO NOT TRUST SHRINK FOR INTERACTION
# res_interaction_shrunk <- lfcShrink(...)  ← not used

# ------------------------------------------------------------
# ANNOTATE
# ------------------------------------------------------------
res_age_wt_df      <- add_symbols(res_age_wt)
res_brca2_old_df   <- add_symbols(res_brca2_old)
res_interaction_df <- add_symbols(res_interaction)
res_brca2_young_df <- add_symbols(res_brca2_young)

res_age_wt_shrunk_df    <- add_symbols(res_age_wt_shrunk)
res_brca2_old_shrunk_df <- add_symbols(res_brca2_old_shrunk)

# ------------------------------------------------------------
# BUILD FULL TABLES (Adrian style)
# ------------------------------------------------------------
build_full <- function(raw, shrunk = NULL, df) {
  
  out <- data.frame(
    ENSEMBL = df$ENSEMBL,
    SYMBOL  = df$SYMBOL,
    baseMean = raw$baseMean,
    
    log2FC_raw = raw$log2FoldChange,
    lfcSE = raw$lfcSE,         
    pvalue = raw$pvalue,      
    padj = raw$padj
  )
  
  if (!is.null(shrunk)) {
    out$log2FC_shrunk <- shrunk$log2FoldChange
  }
  
  out
}

full_age        <- build_full(res_age_wt, res_age_wt_shrunk, res_age_wt_df)
full_brca2_old  <- build_full(res_brca2_old, res_brca2_old_shrunk, res_brca2_old_df)
full_interaction <- build_full(res_interaction, NULL, res_interaction_df)
full_brca2_young <- build_full(res_brca2_young, NULL, res_brca2_young_df)

# ------------------------------------------------------------
# VOLCANO
# ------------------------------------------------------------

plot_volcano_clean(full_interaction,
                   "Interaction: Age × BRCA2",
                   "volcano_interaction.png",
                   use_shrunk = FALSE)

plot_volcano_clean(full_brca2_old,
                   "BRCA2 effect in OLD",
                   "volcano_brca2_old.png",
                   use_shrunk = TRUE)

plot_volcano_clean(full_brca2_young,
                   "BRCA2 effect in YOUNG",
                   "volcano_brca2_young.png",
                   use_shrunk = FALSE)

plot_volcano_clean(full_age,
                   "Age effect in WT",
                   "volcano_age_wt.png",
                   use_shrunk = TRUE)

# ------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------
saveRDS(full_age, file.path(results_dir, "age_full.rds"))
saveRDS(full_brca2_old, file.path(results_dir, "brca2_old_full.rds"))
saveRDS(full_interaction, file.path(results_dir, "interaction_full.rds"))
saveRDS(full_brca2_young, file.path(results_dir, "brca2_young_full.rds"))

write.csv(full_age, file.path(results_dir, "age_full.csv"), row.names = FALSE)
write.csv(full_brca2_old, file.path(results_dir, "brca2_old_full.csv"), row.names = FALSE)
write.csv(full_interaction, file.path(results_dir, "interaction_full.csv"), row.names = FALSE)
write.csv(full_brca2_young, file.path(results_dir, "brca2_young_full.csv"), row.names = FALSE)

# ------------------------------------------------------------
# VST
# ------------------------------------------------------------
vst_mat <- assay(vst(dds, blind = TRUE))
saveRDS(vst_mat, file.path(results_dir, "vst_matrix.rds"))


#Shrinkage was applied to the main effects to improve stability of the estimated fold changes. Interaction terms were analyzed using the original values, as these reflect more complex relationships that are not well suited for shrinkage.

cat("\nDONE\n")


