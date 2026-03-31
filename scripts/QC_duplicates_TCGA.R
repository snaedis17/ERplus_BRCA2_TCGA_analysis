# ============================================================
# QC_duplicates_TCGA_precollapse.R
# QC of duplicate aliquots in FINAL analysis dataset
# (female + primary tumors + annotated, BEFORE collapsing)
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------
if (!file.exists("data/counts_precollapse_TCGA.rds")) {
  stop("Missing counts_precollapse_TCGA.rds. Run 01_prepare_TCGA_data.R first.")
}

if (!file.exists("data/annot_precollapse_TCGA.rds")) {
  stop("Missing annot_precollapse_TCGA.rds. Run 01_prepare_TCGA_data.R first.")
}

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
norm12 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 12)
norm16 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 16)

# ------------------------------------------------------------
# Load data (FINAL dataset BEFORE collapsing duplicates)
# ------------------------------------------------------------
counts_unique <- readRDS("data/counts_precollapse_TCGA.rds")
annot <- readRDS("data/annot_precollapse_TCGA.rds")

cat("Loaded matrix:\n")
print(dim(counts_unique))

cat("Loaded annotation:\n")
print(dim(annot))

# ------------------------------------------------------------
# Restrict to ER+ (same as DESeq analysis)
# ------------------------------------------------------------
annot <- annot %>%
  dplyr::filter(ER == "ER+")

# Keep only overlapping samples
common <- intersect(colnames(counts_unique), annot$aliquot)

# Subset BOTH
counts_unique <- counts_unique[, common, drop = FALSE]

annot <- annot %>%
  dplyr::filter(aliquot %in% common) %>%
  dplyr::distinct(aliquot, .keep_all = TRUE) %>%
  dplyr::arrange(match(aliquot, common))

# Check alignment
stopifnot(identical(colnames(counts_unique), annot$aliquot))


cat("After ER+ filtering:\n")
print(dim(counts_unique))

# ------------------------------------------------------------
# Identify duplicates
# ------------------------------------------------------------
aliquot2pat <- tibble(
  aliquot   = norm16(colnames(counts_unique)),
  patient12 = norm12(colnames(counts_unique))
)

duplicates <- aliquot2pat %>%
  dplyr::group_by(patient12) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(patient12)

cat("Patients with duplicates:\n")
print(length(unique(duplicates$patient12)))
#Patients with duplicates:6

cat("Total duplicate samples:\n")
print(nrow(duplicates))
#Total duplicate samples: 12

cat("Duplicate mapping:\n")
print(duplicates)

#Duplicate mapping:
> print(duplicates)

# A tibble: 12 × 2
# Groups:   patient12 [6]
#aliquot          patient12   
#<chr>            <chr>       
#  1 TCGA-A7-A0DB-01A TCGA-A7-A0DB
#2 TCGA-A7-A0DB-01C TCGA-A7-A0DB
#3 TCGA-A7-A0DC-01B TCGA-A7-A0DC
#4 TCGA-A7-A0DC-01A TCGA-A7-A0DC
#5 TCGA-A7-A13G-01B TCGA-A7-A13G
#6 TCGA-A7-A13G-01A TCGA-A7-A13G
#7 TCGA-A7-A26E-01A TCGA-A7-A26E
#8 TCGA-A7-A26E-01B TCGA-A7-A26E
#9 TCGA-A7-A26J-01A TCGA-A7-A26J
#10 TCGA-A7-A26J-01B TCGA-A7-A26J
#11 TCGA-AC-A3OD-01A TCGA-AC-A3OD
#12 TCGA-AC-A3OD-01B TCGA-AC-A3OD
# ------------------------------------------------------------
# Save duplicate mapping (for Adrian)
# ------------------------------------------------------------
write.csv(
  duplicates,
  "results/QC_duplicate_aliquots_ERpos.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Correlation between duplicate samples
# ------------------------------------------------------------
dup_groups <- duplicates %>%
  dplyr::group_by(patient12) %>%
  dplyr::summarise(samples = list(aliquot), .groups = "drop")

cor_results <- purrr::map_df(dup_groups$samples, function(samps) {
  
  if (length(samps) < 2) return(NULL)
  
  combs <- combn(samps, 2, simplify = FALSE)
  
  purrr::map_df(combs, function(pair) {
    data.frame(
      patient12   = norm12(pair[1]),
      sample1     = pair[1],
      sample2     = pair[2],
      correlation = cor(
        counts_unique[, pair[1]],
        counts_unique[, pair[2]],
        method = "pearson"
      )
    )
  })
})

# ------------------------------------------------------------
# Inspect correlations
# ------------------------------------------------------------
cat("Correlation summary:\n")
print(summary(cor_results$correlation))

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3025  0.7283  0.7531  0.6860  0.7735  0.8069 

cat("All duplicate correlations:\n")
print(cor_results)

#All duplicate correlations:
#> print(cor_results)
#patient12          sample1          sample2 correlation
#1 TCGA-A7-A0DB TCGA-A7-A0DB-01A TCGA-A7-A0DB-01C   0.8069017
#2 TCGA-A7-A0DC TCGA-A7-A0DC-01B TCGA-A7-A0DC-01A   0.7555540
#3 TCGA-A7-A13G TCGA-A7-A13G-01B TCGA-A7-A13G-01A   0.7208643
#4 TCGA-A7-A26E TCGA-A7-A26E-01A TCGA-A7-A26E-01B   0.7795402
#5 TCGA-A7-A26J TCGA-A7-A26J-01A TCGA-A7-A26J-01B   0.7506292
#6 TCGA-AC-A3OD TCGA-AC-A3OD-01A TCGA-AC-A3OD-01B   0.3025360
# ------------------------------------------------------------
# Flag potential problematic duplicates
# ------------------------------------------------------------
low_cor <- cor_results %>%
  dplyr::filter(correlation < 0.9)

cat("Potential problematic duplicates (cor < 0.9):\n")
print(nrow(low_cor))
print(low_cor)

# ------------------------------------------------------------
# Save QC outputs
# ------------------------------------------------------------
write.csv(
  cor_results,
  "results/QC_duplicate_correlations_ERpos.csv",
  row.names = FALSE
)

write.csv(
  low_cor,
  "results/QC_duplicate_lowcor_ERpos.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save expression matrix used for QC (optional, for sharing)
# ------------------------------------------------------------
saveRDS(
  counts_unique,
  "results/QC_counts_ERpos_with_duplicates.rds"
)

cat("Duplicate QC complete\n")
