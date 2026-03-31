# ============================================================
# 10_gsea_mutold_and_shared_pathways.R
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# ------------------------------------------------------------
# FIX namespace conflicts
# ------------------------------------------------------------
filter  <- dplyr::filter
select  <- dplyr::select
mutate  <- dplyr::mutate
arrange <- dplyr::arrange
rename  <- dplyr::rename

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
results_dir <- "results/within_age_brca2"
out_dir <- file.path(results_dir, "gsea_shared")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load DESeq results
# ------------------------------------------------------------
res_mut_old   <- readRDS(file.path(results_dir, "res_MutOld_vs_WTOld_RAW.rds"))
res_wt_age    <- readRDS(file.path(results_dir, "res_WTYoung_vs_WTOld_RAW.rds"))
res_mut_young <- readRDS(file.path(results_dir, "res_MutYoung_vs_WTYoung_RAW.rds"))
res_age_mut   <- readRDS(file.path(results_dir, "res_MutYoung_vs_MutOld_RAW.rds"))

res_mut_old   <- as.data.frame(res_mut_old)
res_wt_age    <- as.data.frame(res_wt_age)
res_mut_young <- as.data.frame(res_mut_young)
res_age_mut   <- as.data.frame(res_age_mut)

# ------------------------------------------------------------
# Ranked lists
# ------------------------------------------------------------
make_rank <- function(df) {
  df %>%
    filter(!is.na(stat), !is.na(SYMBOL), baseMean > 10) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    {
      stats <- .$stat
      names(stats) <- .$SYMBOL
      stats
    } %>%
    sort(decreasing = TRUE)
}

r_mut_old   <- make_rank(res_mut_old)
r_wt        <- make_rank(res_wt_age)
r_mut_young <- make_rank(res_mut_young)
r_mut_age   <- make_rank(res_age_mut)

# ------------------------------------------------------------
# Pathways (Hallmark)
# ------------------------------------------------------------
msig_h <- msigdbr(species = "Homo sapiens", category = "H")
pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)

# ------------------------------------------------------------
# Run GSEA
# ------------------------------------------------------------
gsea_mut_old   <- fgsea(pathways_h, r_mut_old,   minSize = 15, maxSize = 500)
gsea_wt        <- fgsea(pathways_h, r_wt,        minSize = 15, maxSize = 500)
gsea_mut_young <- fgsea(pathways_h, r_mut_young, minSize = 15, maxSize = 500)
gsea_mut_age   <- fgsea(pathways_h, r_mut_age,   minSize = 15, maxSize = 500)

gsea_mut_old   <- as.data.frame(gsea_mut_old)   %>% arrange(padj)
gsea_wt        <- as.data.frame(gsea_wt)        %>% arrange(padj)
gsea_mut_young <- as.data.frame(gsea_mut_young) %>% arrange(padj)
gsea_mut_age   <- as.data.frame(gsea_mut_age)   %>% arrange(padj)

# ------------------------------------------------------------
# Significant pathways (shared WT vs Mut Old)
# ------------------------------------------------------------
sig_wt  <- gsea_wt      %>% filter(!is.na(padj) & padj < 0.05)
sig_mut <- gsea_mut_old %>% filter(!is.na(padj) & padj < 0.05)

shared <- intersect(sig_wt$pathway, sig_mut$pathway)

cat("Shared pathways:", length(shared), "\n")

#12 shared pathways

# ------------------------------------------------------------
# Merge + direction
# ------------------------------------------------------------
shared_df <- sig_wt %>%
  select(pathway, NES_WT = NES, padj_WT = padj) %>%
  inner_join(
    sig_mut %>%
      select(pathway, NES_MUT = NES, padj_MUT = padj),
    by = "pathway"
  ) %>%
  mutate(
    same_direction = sign(NES_WT) == sign(NES_MUT)
  ) %>%
  arrange(desc(abs(NES_WT)))

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
write.csv(shared_df,
          file.path(out_dir, "shared_pathways.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# Plot: pathway concordance
# ------------------------------------------------------------
p <- ggplot(shared_df, aes(x = NES_WT, y = NES_MUT, label = pathway)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 30) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Shared biological programs",
    subtitle = paste0("n = ", length(shared)),
    x = "WT Young vs Old (NES)",
    y = "Mut Old vs WT Old (NES)"
  )

ggsave(file.path(out_dir, "shared_pathways.png"),
       p, width = 8, height = 6, dpi = 300)

print(p)

# ------------------------------------------------------------
# Bubble plot function
# ------------------------------------------------------------
plot_gsea_bubble <- function(gsea_df, title, filename) {
  
  gsea_df <- gsea_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      negLogFDR = -log10(padj),
      negLogFDR = pmin(negLogFDR, 10)   # cap values manually
    )
  
  top <- gsea_df %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(
      pathway = gsub("HALLMARK_", "", pathway),
      pathway = gsub("_", " ", pathway)
    )
  
  p <- ggplot(top, aes(x = NES, y = reorder(pathway, NES))) +
    geom_point(aes(size = negLogFDR, color = NES)) +
    
    # FIXED color scale (works with oob)
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      limits = c(-3, 3),
      oob = scales::squish
    ) +
    
    # size WITHOUT oob
    scale_size(
      range = c(3, 10),
      limits = c(0, 10)
    ) +
    
    theme_bw() +
    labs(
      title = title,
      x = "NES",
      y = "",
      size = "-log10(FDR)"
    )
  
  ggsave(file.path(out_dir, filename),
         p, width = 7, height = 5, dpi = 300)
  
  print(p)
}
# ------------------------------------------------------------
# Bubble plots (ALL comparisons)
# ------------------------------------------------------------
plot_gsea_bubble(
  gsea_wt,
  "GSEA Hallmark – WT Young vs Old",
  "bubble_wt_age.png"
)

plot_gsea_bubble(
  gsea_mut_old,
  "GSEA Hallmark – Mut Old vs WT Old",
  "bubble_mut_old.png"
)

plot_gsea_bubble(
  gsea_mut_young,
  "GSEA Hallmark – Mut Young vs WT Young",
  "bubble_mut_young.png"
)

plot_gsea_bubble(
  gsea_mut_age,
  "GSEA Hallmark – Mut Young vs Mut Old",
  "bubble_mut_age.png"
)

cat("\nScript 10 complete\n")



# ------------------------------------------------------------
# Prepare combined bubble dataset
# ------------------------------------------------------------
prep_bubble <- function(gsea_df, label) {
  
  gsea_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      negLogFDR = -log10(padj),
      negLogFDR = pmin(negLogFDR, 10),
      pathway = gsub("HALLMARK_", "", pathway),
      pathway = gsub("_", " ", pathway),
      contrast = label
    ) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 10)
}

bubble_all <- dplyr::bind_rows(
  prep_bubble(gsea_wt,        "WT Young vs Old"),
  prep_bubble(gsea_mut_old,   "Mut Old vs WT Old"),
  prep_bubble(gsea_mut_young, "Mut Young vs WT Young"),
  prep_bubble(gsea_mut_age,   "Mut Young vs Mut Old")
)

# ------------------------------------------------------------
# Combined bubble plot (FACET)
# ------------------------------------------------------------
p_combined <- ggplot(bubble_all, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = negLogFDR, color = NES)) +
  
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  
  scale_size(
    range = c(3, 10),
    limits = c(0, 10)
  ) +
  
  facet_wrap(~ contrast, scales = "free_y") +
  
  theme_bw() +
  labs(
    title = "GSEA Hallmark pathways across comparisons",
    x = "NES",
    y = "",
    size = "-log10(FDR)"
  )

ggsave(
  file.path(out_dir, "bubble_all_contrasts.png"),
  p_combined,
  width = 12,
  height = 8,
  dpi = 300
)

print(p_combined)
