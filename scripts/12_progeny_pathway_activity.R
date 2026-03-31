# ============================================================
# 12_progeny_pathway_activity.R
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(decoupleR)
  library(ggplot2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/within_age_brca2"
out_dir <- file.path(results_dir, "progeny_activity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load DESeq results (ALL contrasts)
# ------------------------------------------------------------
res_wt <- readRDS(file.path(results_dir, "res_WTYoung_vs_WTOld_RAW.rds"))
res_mut_young <- readRDS(file.path(results_dir, "res_MutYoung_vs_WTYoung_RAW.rds"))
res_mut_old   <- readRDS(file.path(results_dir, "res_MutOld_vs_WTOld_RAW.rds"))
res_mut_age   <- readRDS(file.path(results_dir, "res_MutYoung_vs_MutOld_RAW.rds"))

res_wt        <- as.data.frame(res_wt)
res_mut_young <- as.data.frame(res_mut_young)
res_mut_old   <- as.data.frame(res_mut_old)
res_mut_age   <- as.data.frame(res_mut_age)

# ------------------------------------------------------------
# Prepare matrix (STAT-based ranking)
# ------------------------------------------------------------
# ------------------------------------------------------------
# Prepare matrix (FIXED FILTER)
# ------------------------------------------------------------
prep_mat <- function(df, label) {
  df %>%
    dplyr::filter(!is.na(stat), !is.na(SYMBOL), baseMean > 10) %>%
    
    # velja strongest signal per gene
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(abs(stat), n = 1) %>%
    dplyr::ungroup() %>%
    
    {
      mat <- matrix(.$stat, ncol = 1)
      rownames(mat) <- .$SYMBOL
      colnames(mat) <- label
      mat
    }
}

mat_wt        <- prep_mat(res_wt,        "WT_age")
mat_mut_young <- prep_mat(res_mut_young, "Mut_young")
mat_mut_old   <- prep_mat(res_mut_old,   "Mut_old")
mat_mut_age   <- prep_mat(res_mut_age,   "Mut_age")

# ------------------------------------------------------------
# Load PROGENy model
# ------------------------------------------------------------
progeny_net <- decoupleR::get_progeny(organism = "human", top = 500)

# ------------------------------------------------------------
# Run PROGENy (ULM)
# ------------------------------------------------------------
run_prog <- function(mat) {
  res <- run_ulm(
    mat = mat,
    net = progeny_net,
    .source = "source",
    .target = "target",
    .mor = "weight"
  )
  
  res %>%
    dplyr::rename(
      pathway = source,
      NES = score,
      pval = p_value
    ) %>%
    dplyr::mutate(
      padj = p.adjust(pval, method = "fdr"),
      negLogFDR = -log10(padj)
    ) %>%
    dplyr::arrange(desc(abs(NES)))   
}
prog_wt        <- run_prog(mat_wt)
prog_mut_young <- run_prog(mat_mut_young)
prog_mut_old   <- run_prog(mat_mut_old)
prog_mut_age   <- run_prog(mat_mut_age)

# ------------------------------------------------------------
# Add contrast labels
# ------------------------------------------------------------
prog_wt$contrast        <- "WT Young vs Old"
prog_mut_young$contrast <- "Mut Young vs WT Young"
prog_mut_old$contrast   <- "Mut Old vs WT Old"
prog_mut_age$contrast   <- "Mut Young vs Mut Old"

# ------------------------------------------------------------
# Combine all
# ------------------------------------------------------------
prog_all <- dplyr::bind_rows(
  prog_wt,
  prog_mut_young,
  prog_mut_old,
  prog_mut_age
)

# ------------------------------------------------------------
# Select top pathways per contrast
# ------------------------------------------------------------
top_prog <- prog_all %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::group_by(contrast) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# Plot (FACET bubble plot)
# ------------------------------------------------------------
p <- ggplot(top_prog, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = negLogFDR, color = NES)) +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-5, 5),
    oob = scales::squish
  ) +
  
  scale_size(range = c(3, 10)) +
  
  facet_wrap(~contrast, scales = "free_y") +
  
  theme_bw() +
  labs(
    title = "PROGENy pathway activity across contrasts",
    x = "Pathway activity (NES)",
    y = "",
    size = "-log10(FDR)",
    color = "NES"
  )

ggsave(
  file.path(out_dir, "progeny_all_comparisons.png"),
  p,
  width = 10,
  height = 7,
  dpi = 300
)

print(p)

# ------------------------------------------------------------
# Save full results
# ------------------------------------------------------------
write.csv(
  prog_all,
  file.path(out_dir, "progeny_all_results.csv"),
  row.names = FALSE
)

cat("\nPROGENy analysis complete\n")
