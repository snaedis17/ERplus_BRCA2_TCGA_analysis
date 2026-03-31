# ============================================================
# 11_tf_activity_key_contrasts.R 
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(decoupleR)
  library(dorothea)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/within_age_brca2"
out_dir <- file.path(results_dir, "tf_activity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load DESeq results (FIXED NAMES)
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
# Prepare matrix (CRITICAL FILTER ADDED)
# ------------------------------------------------------------
prep_mat <- function(df, label) {
  df %>%
    dplyr::filter(!is.na(stat), !is.na(SYMBOL), baseMean > 10) %>%
    
    # velja eitt gildi per gene (sterkasta signal)
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

mat_wt        <- prep_mat(res_wt, "WT_age")
mat_mut_young <- prep_mat(res_mut_young, "Mut_young")
mat_mut_old   <- prep_mat(res_mut_old, "Mut_old")
mat_mut_age   <- prep_mat(res_mut_age, "Mut_age")

# ------------------------------------------------------------
# TF network
# ------------------------------------------------------------
data(dorothea_hs, package = "dorothea")

tf_net <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A","B","C")) %>%
  dplyr::group_by(tf) %>%
  dplyr::filter(n() >= 10) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# Run TF activity
# ------------------------------------------------------------
run_tf <- function(mat) {
  res <- run_ulm(
    mat = mat,
    net = tf_net,
    .source = "tf",
    .target = "target",
    .mor = "mor"
  )
  
  res %>%
    dplyr::rename(
      tf = source,
      NES = score,
      pval = p_value
    ) %>%
    dplyr::mutate(
      padj = p.adjust(pval, method = "fdr")
    ) %>%
    dplyr::arrange(desc(abs(NES)))   # þetta bætir clarity
}

tf_wt        <- run_tf(mat_wt)
tf_mut_young <- run_tf(mat_mut_young)
tf_mut_old   <- run_tf(mat_mut_old)
tf_mut_age   <- run_tf(mat_mut_age)

# ------------------------------------------------------------
# Concordance plot
# ------------------------------------------------------------
tf_cmp <- dplyr::inner_join(
  tf_wt %>% dplyr::select(tf, NES_WT = NES),
  tf_mut_young %>% dplyr::select(tf, NES_MUT = NES),
  by = "tf"
)

cor_res <- cor.test(tf_cmp$NES_WT, tf_cmp$NES_MUT, method = "spearman")

p <- ggplot(tf_cmp, aes(NES_WT, NES_MUT)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "TF activity concordance",
    subtitle = paste0("Spearman rho = ", round(cor_res$estimate, 3)),
    x = "WT Young vs Old",
    y = "Mut Young vs WT Young"
  )

ggsave(file.path(out_dir, "tf_concordance.png"), p, width = 7, height = 5)
print(p)

# ------------------------------------------------------------
# Bubble plot function
# ------------------------------------------------------------
plot_tf_bubble <- function(tf_df, title, filename) {
  
  tf_df <- tf_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(negLogFDR = -log10(padj))
  
  top <- tf_df %>%
    dplyr::filter(padj < 0.05) %>%                 
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 15)
  
  p <- ggplot(top, aes(x = NES, y = reorder(tf, NES))) +
    geom_point(aes(size = negLogFDR, color = NES)) +
    
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      limits = c(-5, 5),
      oob = scales::squish
    ) +
    
    scale_size(range = c(3, 10)) +
    
    theme_bw() +
    labs(
      title = title,
      x = "TF activity (NES)",
      y = "",
      size = "-log10(FDR)"
    )
  
  ggsave(file.path(out_dir, filename),
         p, width = 7, height = 5, dpi = 300)
  
  print(p)
}

# ------------------------------------------------------------
# Individual plots
# ------------------------------------------------------------
plot_tf_bubble(tf_wt,        "TF activity – WT Young vs Old", "tf_wt_age.png")
plot_tf_bubble(tf_mut_young, "TF activity – Mut Young vs WT Young", "tf_mut_young.png")
plot_tf_bubble(tf_mut_old,   "TF activity – Mut Old vs WT Old", "tf_mut_old.png")
plot_tf_bubble(tf_mut_age,   "TF activity – Mut Young vs Mut Old", "tf_mut_age.png")

cat("\nTF analysis complete\n")

# ------------------------------------------------------------
# Combine results (NO DOUBLE FDR)
# ------------------------------------------------------------
tf_wt$contrast        <- "WT Young vs Old"
tf_mut_young$contrast <- "Mut Young vs WT Young"
tf_mut_old$contrast   <- "Mut Old vs WT Old"
tf_mut_age$contrast   <- "Mut Young vs Mut Old"

tf_all <- dplyr::bind_rows(
  tf_wt,
  tf_mut_young,
  tf_mut_old,
  tf_mut_age
)

tf_all <- tf_all %>%
  dplyr::mutate(
    negLogFDR = -log10(padj)
  ) %>%
  dplyr::filter(!is.na(padj))

tf_top <- tf_all %>%
  dplyr::group_by(contrast) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# Combined plot
# ------------------------------------------------------------
p <- ggplot(tf_top, aes(x = NES, y = reorder(tf, NES))) +
  geom_point(aes(size = negLogFDR, color = NES)) +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-5, 5),
    oob = scales::squish
  ) +
  
  scale_size(range = c(3, 10)) +
  
  facet_wrap(~ contrast, scales = "free_y") +
  
  theme_bw() +
  labs(
    title = "Transcription factor activity across contrasts",
    x = "TF activity (NES)",
    y = "",
    size = "-log10(FDR)"
  )

ggsave(
  file.path(out_dir, "tf_all_contrasts_combined.png"),
  p,
  width = 12,
  height = 8,
  dpi = 300
)

print(p)
