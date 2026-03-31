# ============================================================
# 09_gsea_compare_key_contrasts.R
# ============================================================

# Clean workspace for reproducibility
rm(list = ls())

# Load required libraries (silencing startup messages)
suppressPackageStartupMessages({
  library(dplyr)     # data manipulation
  library(fgsea)     # fast GSEA implementation
  library(msigdbr)   # MSigDB gene sets
  library(ggplot2)   # plotting
})

# ------------------------------------------------------------
# Define output directories
# ------------------------------------------------------------
results_dir <- "results/within_age_brca2"
out_dir <- file.path(results_dir, "gsea_key_contrasts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load DESeq2 results
# ------------------------------------------------------------
# WT: age effect (Young vs Old in WT tumors)
# Mut: BRCA2 effect in young patients (Mut vs WT in young)
res_wt  <- readRDS(file.path(results_dir, "res_WTYoung_vs_WTOld_RAW.rds"))
res_mut <- readRDS(file.path(results_dir, "res_MutYoung_vs_WTYoung_RAW.rds"))

# Convert to data frames for easier handling
res_wt  <- as.data.frame(res_wt)
res_mut <- as.data.frame(res_mut)

# ------------------------------------------------------------
# Create ranked gene lists for GSEA
# ------------------------------------------------------------
# GSEA requires genes ranked by a statistic (here: Wald stat)
make_rank <- function(df) {
  df %>%
    # Remove genes with missing values or very low expression
    dplyr::filter(!is.na(stat), !is.na(SYMBOL), baseMean > 10) %>%
    
    # Ensure each gene appears only once
    dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
    
    # Create named vector: gene -> statistic
    { stats <- .$stat; names(stats) <- .$SYMBOL; stats } %>%
    
    # Sort from most upregulated to most downregulated
    sort(decreasing = TRUE)
}

# Generate ranked lists
r_wt  <- make_rank(res_wt)
r_mut <- make_rank(res_mut)

# ------------------------------------------------------------
# Load Hallmark gene sets (MSigDB)
# ------------------------------------------------------------
# Hallmark pathways summarize major biological processes
msig_h <- msigdbr(species = "Homo sapiens", category = "H")

# Convert into list format required by fgsea:
# pathway -> vector of genes
pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)

# ------------------------------------------------------------
# Run GSEA
# ------------------------------------------------------------
# Test whether genes in each pathway are enriched at the top/bottom of ranked lists
gsea_wt  <- fgsea(pathways_h, r_wt,  minSize = 15, maxSize = 500)
gsea_mut <- fgsea(pathways_h, r_mut, minSize = 15, maxSize = 500)

# Sort pathways by adjusted p-value (most significant first)
gsea_wt  <- as.data.frame(gsea_wt)  %>% arrange(padj)
gsea_mut <- as.data.frame(gsea_mut) %>% arrange(padj)

# ------------------------------------------------------------
# Save GSEA results
# ------------------------------------------------------------
saveRDS(gsea_wt,  file.path(out_dir, "gsea_wt.rds"))
saveRDS(gsea_mut, file.path(out_dir, "gsea_mut.rds"))

# ------------------------------------------------------------
# Prepare data for comparison
# ------------------------------------------------------------
# Keep only key columns and rename for clarity
prep <- function(df, label) {
  df %>%
    dplyr::select(pathway, NES, padj) %>%
    dplyr::rename(
      !!paste0("NES_", label)  := NES,
      !!paste0("padj_", label) := padj
    )
}

wt_clean  <- prep(gsea_wt, "WT")
mut_clean <- prep(gsea_mut, "Mut")

# Keep only pathways present in both analyses
cmp <- inner_join(wt_clean, mut_clean, by = "pathway")

# ------------------------------------------------------------
# Correlation analysis
# ------------------------------------------------------------
# Test whether pathway activity is concordant between WT and Mut
cor_res <- cor.test(cmp$NES_WT, cmp$NES_Mut, method = "spearman")
print(cor_res)

# ------------------------------------------------------------
# Plot pathway concordance
# ------------------------------------------------------------
# Each point = one pathway
# x-axis = aging effect (WT Young vs Old)
# y-axis = BRCA2 effect (Mut vs WT in young)
p <- ggplot(cmp, aes(NES_WT, NES_Mut)) +
  
  # Scatter points (each pathway)
  geom_point(alpha = 0.7) +
  
  # Linear trend (overall concordance)
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  
  # Reference lines (no enrichment)
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_bw() +
  
  labs(
    title = "Pathway concordance",
    subtitle = paste0("Spearman rho = ", round(cor_res$estimate, 3)),
    x = "WT Young vs Old (aging effect)",
    y = "Mut Young vs WT Young (BRCA2 effect)"
  )

# Save and display plot
ggsave(file.path(out_dir, "gsea_concordance.png"), p, width = 7, height = 5)
print(p)

# ------------------------------------------------------------
# Interpretation note
# ------------------------------------------------------------
# Pathways falling in opposite quadrants (e.g., +x / -y) indicate
# opposing biological effects between aging and BRCA2 mutation.

cat("\nScript 09 complete\n")
