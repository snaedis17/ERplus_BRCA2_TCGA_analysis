# ============================================================
#Compare_contrasts_overlap_concordance.R
# ============================================================
#Compare the aging signature (WT) with BRCA2 effects in both young and old tumors 
#to determine whether BRCA2 tumors imitate or differ from aging biology.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

results_dir <- "results/within_age_brca2"
out_dir     <- "results/contrast_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


filter  <- dplyr::filter
select  <- dplyr::select
arrange <- dplyr::arrange
mutate  <- dplyr::mutate
# ------------------------------------------------------------
# Load results
# ------------------------------------------------------------
res_wt_age  <- readRDS(file.path(results_dir, "res_wtyoung_vs_wtold.rds"))
res_mut_yng <- readRDS(file.path(results_dir, "res_mutyoung_vs_wtyoung.rds"))
res_mut_old <- readRDS(file.path(results_dir, "res_mutold_vs_wtold.rds"))

# ------------------------------------------------------------
# Prepare
# ------------------------------------------------------------
prep <- function(df, name){
  
  df %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::select(SYMBOL, log2FoldChange, padj) %>%
    dplyr::rename(
      log2FC = log2FoldChange,
      padj   = padj
    ) %>%
    dplyr::rename_with(~ paste0(., "_", name), -SYMBOL)
}

wt_age  <- prep(res_wt_age,  "WTage")
mut_yng <- prep(res_mut_yng, "MutY")
mut_old <- prep(res_mut_old, "MutO")

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------
collapse_symbol <- function(df, suffix){
  
  lfc_col  <- paste0("log2FC_", suffix)
  padj_col <- paste0("padj_", suffix)
  
  df %>%
    group_by(SYMBOL) %>%
    summarise(
      !!lfc_col := .data[[lfc_col]][which.max(abs(.data[[lfc_col]]))],
      !!padj_col := if(all(is.na(.data[[padj_col]]))) NA_real_ else min(.data[[padj_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

wt_age_clean  <- collapse_symbol(wt_age,  "WTage")
mut_yng_clean <- collapse_symbol(mut_yng, "MutY")
mut_old_clean <- collapse_symbol(mut_old, "MutO")

cmp_young <- wt_age_clean %>%
  inner_join(mut_yng_clean, by="SYMBOL")

cmp_old <- wt_age_clean %>%
  inner_join(mut_old_clean, by="SYMBOL")

# ------------------------------------------------------------
# CORRELATION
# ------------------------------------------------------------
cor_young <- cor.test(
  cmp_young$log2FC_WTage,
  cmp_young$log2FC_MutY,
  method="spearman"
)

cor_old <- cor.test(
  cmp_old$log2FC_WTage,
  cmp_old$log2FC_MutO,
  method="spearman"
)

print(cor_young)
#Spearman's rank correlation rho
#data:  cmp_young$log2FC_WTage and cmp_young$log2FC_MutY
#S = 6.9077e+11, p-value = 1.805e-09
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
 # rho 
#-0.04780204 

print(cor_old)

#Spearman's rank correlation rho
#data:  cmp_old$log2FC_WTage and cmp_old$log2FC_MutO
#S = 4.4704e+11, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
 # rho 
#0.3218999 


# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------
p1 <- ggplot(cmp_young, aes(log2FC_WTage, log2FC_MutY)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha=0.4) +
  geom_smooth(method="lm", color="red") +
  theme_bw() +
  labs(title="WT age vs Mut young")

p2 <- ggplot(cmp_old, aes(log2FC_WTage, log2FC_MutO)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha=0.4) +
  geom_smooth(method="lm", color="red") +
  theme_bw() +
  labs(title="WT age vs Mut old")



ggsave(file.path(out_dir, "concordance_young.png"), p1, width=6, height=5)
ggsave(file.path(out_dir, "concordance_old.png"), p2, width=6, height=5)

print(p1)
print(p2)

#summary from plots: The transcriptional programs associated with aging in ER+ wild-type tumors are largely absent or reversed in BRCA2-mutated tumors from young patients.

#BRCA2-mutated tumors in older patients partially recapitulate age-associated transcriptional changes observed in wild-type tumors.

#Age is not just a covariate — it fundamentally changes the biology of BRCA2-mutated tumors.
# ------------------------------------------------------------
# OVERLAP (significant genes)
# ------------------------------------------------------------
sig_wt  <- res_wt_age  %>% filter(!is.na(padj) & padj < 0.05)
sig_mut <- res_mut_yng %>% filter(!is.na(padj) & padj < 0.05)
sig_old <- res_mut_old %>% filter(!is.na(padj) & padj < 0.05)

overlap_young <- intersect(sig_wt$SYMBOL, sig_mut$SYMBOL)
overlap_old   <- intersect(sig_wt$SYMBOL, sig_old$SYMBOL)

cat("Overlap WT age vs Mut young:", length(overlap_young), "\n")
#Overlap WT age vs Mut young: 2 
cat("Overlap WT age vs Mut old:", length(overlap_old), "\n")
#Overlap WT age vs Mut old: 71

# save overlap
write.csv(data.frame(SYMBOL=overlap_young),
          file.path(out_dir, "overlap_young.csv"),
          row.names=FALSE)

write.csv(data.frame(SYMBOL=overlap_old),
          file.path(out_dir, "overlap_old.csv"),
          row.names=FALSE)

# ------------------------------------------------------------
# SAME DIRECTION 
# ------------------------------------------------------------
same_dir_young <- cmp_young %>%
  filter(
    padj_WTage < 0.05,
    padj_MutY  < 0.05,
    sign(log2FC_WTage) == sign(log2FC_MutY)
  )

write.csv(same_dir_young,
          file.path(out_dir, "same_direction_young.csv"),
          row.names=FALSE)

cat("Same-direction genes (young):", nrow(same_dir_young), "\n")
#Same-direction genes (young): 1 

cat("\nDone\n")
#Same-direction genes (young): 1 

# correlation summary table
cor_summary <- data.frame(
  comparison = c("WT vs MutYoung", "WT vs MutOld"),
  rho = c(cor_young$estimate, cor_old$estimate),
  pval = c(cor_young$p.value, cor_old$p.value)
)

write.csv(cor_summary,
          file.path(out_dir, "correlation_summary.csv"),
          row.names = FALSE)
#When I compared age-related gene expression in wild-type tumors to BRCA2-mutated tumors, I saw clear differences. In older patients, BRCA2 tumors showed some similarity to normal aging patterns. However, in young patients, there was almost no overlap, and some changes even went in the opposite direction. This suggests that BRCA2-mutated tumors in young individuals follow a different biology, rather than just being an accelerated form of aging..
