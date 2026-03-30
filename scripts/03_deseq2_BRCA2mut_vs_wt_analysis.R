rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
  library(AnnotationDbi)
  library(edgeR)
})

# shorthand
filter  <- dplyr::filter
select  <- dplyr::select
arrange <- dplyr::arrange

# ============================================================
# GOAL
# ============================================================
# Compare gene expression between:
#   BRCA2-mutated vs WT tumors (ER+ only)
#
# Using:
#   - DESeq2 for differential expression
#   - Shrinkage for stable effect sizes
#   - VST for visualization
# ============================================================


# ============================================================
# LOAD DATA
# ============================================================
# counts = raw RNA-seq counts (NOT normalized)
# annot  = sample metadata (patients, BRCA2 status, etc.)

counts_one <- readRDS("data/counts_ERpos_clean_final.rds")
annot_patient <- readRDS("data/annot_ERpos_clean_final.rds")


# ------------------------------------------------------------
# Make sure samples match between counts and metadata
# (very important sanity check)
# ------------------------------------------------------------
stopifnot(identical(colnames(counts_one), annot_patient$patient12))


# ============================================================
# BUILD DESIGN MATRIX
# ============================================================
# grp = comparison variable (WT vs Mut)

coldata <- annot_patient %>%
  transmute(
    patient12,
    grp = factor(BRCA2_Status, levels = c("WT", "Mut"))
  ) %>%
  as.data.frame()

rownames(coldata) <- coldata$patient12
coldata$patient12 <- NULL

stopifnot(identical(rownames(coldata), colnames(counts_one)))


# ============================================================
# FILTER LOW-EXPRESSION GENES
# ============================================================
# Remove genes that are too lowly expressed to be reliable
# (reduces noise, improves statistical power)

dge <- DGEList(counts = counts_one)

keep <- filterByExpr(
  y = dge,
  group = coldata$grp
)

counts_one <- counts_one[keep, , drop = FALSE]

# Now we only keep genes with meaningful expression


# ============================================================
# RUN DESeq2
# ============================================================
# This models count data using negative binomial distribution

dds <- DESeqDataSetFromMatrix(
  countData = counts_one,
  colData   = coldata,
  design    = ~ grp
)

dds <- DESeq(dds)

# coefficients (important for shrink later)
resultsNames(dds)


# ============================================================
# EXTRACT RESULTS
# ============================================================
# res_raw = original DESeq2 results
#   → contains log2FC (can be noisy for low counts)

res_raw <- results(dds, name = "grp_Mut_vs_WT")


# ============================================================
# SHRINKAGE (IMPORTANT)
# ============================================================
# This stabilizes log2 fold changes
# → prevents extreme values for low-count genes

res_shrunk <- lfcShrink(
  dds,
  coef = "grp_Mut_vs_WT",
  type = "apeglm"
)

# KEY IDEA:
# raw = statistical test
# shrunk = better effect size for interpretation


# ============================================================
# ADD GENE SYMBOLS
# ============================================================

add_symbols <- function(res_df) {
  res_df$ENSEMBL <- gsub("\\..*", "", rownames(res_df))
  res_df$SYMBOL <- mapIds(
    EnsDb.Hsapiens.v86,
    keys = res_df$ENSEMBL,
    keytype = "GENEID",
    column = "SYMBOL",
    multiVals = "first"
  )
  res_df %>% relocate(SYMBOL, ENSEMBL)
}

res_raw_df    <- add_symbols(as.data.frame(res_raw))
res_shrunk_df <- add_symbols(as.data.frame(res_shrunk))


# ============================================================
# COMBINE RAW + SHRUNK
# ============================================================
# This is what Adrian does manually
# → keep BOTH versions

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


# ============================================================
# DEFINE SIGNIFICANT GENES
# ============================================================
# Use:
#   - statistical significance (padj)
#   - biological effect size (shrunk LFC)

sig_results <- full_results %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FC_shrunk) >= 1
  ) %>%
  arrange(padj)


# ============================================================
# CLEAN VOLCANO (publication style)
# ============================================================

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

p_volcano <- ggplot(volcano_df, aes(x = log2FC_shrunk, y = negLogP)) +
  
  # points
  geom_point(aes(color = category), alpha = 0.7, size = 1.8) +
  
  # colors (clean)
  scale_color_manual(values = c(
    "Up" = "#D55E00",     # red/orange
    "Down" = "#0072B2",   # blue
    "NS" = "grey70"
  )) +
  
  # threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  
  # theme
  theme_classic(base_size = 14) +
  
  labs(
    title = "BRCA2-mutated vs WT ER+ tumors",
    x = "log2 fold change (shrunk)",
    y = "-log10 adjusted p-value",
    color = ""
  )

# ------------------------------------------------------------
# LABEL TOP GENES
# ------------------------------------------------------------
library(ggrepel)

top_genes <- volcano_df %>%
  filter(category != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 15)

p_volcano <- p_volcano +
  ggrepel::geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 20
  )

print(p_volcano)
# ============================================================
# VST TRANSFORMATION (VERY IMPORTANT)
# ============================================================
# This is used for:
#   - plots
#   - PCA
#   - ssGSEA
#
# It transforms counts → normalized expression values
# similar idea to log2(TPM + 1)

vst_mat <- assay(vst(dds, blind = TRUE))

saveRDS(vst_mat, "data/vst_expression_ERpos_cleanWT.rds")


# ============================================================
# QUALITY CHECK: raw vs shrunk
# ============================================================
# If shrink worked → points move toward 0

plot(res_raw$log2FoldChange, res_shrunk$log2FoldChange)
abline(0,1,col="red")
