# ------------------------------------------------------------
# OUTPUT:
# - counts_precollapse_TCGA.rds
# - annot_precollapse_TCGA.rds
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
})

filter   <- dplyr::filter
select   <- dplyr::select
distinct <- dplyr::distinct

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
norm12 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 12)
norm16 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 16)

# ------------------------------------------------------------
# 1) Download TCGA
# ------------------------------------------------------------
q <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(q)
se <- GDCprepare(q)

saveRDS(se, "data/se_TCGA_raw.rds")

#se    <- readRDS("data/se_TCGA_raw.rds")
counts <- assay(se)

# ------------------------------------------------------------
# Clean gene IDs
# ------------------------------------------------------------
ens <- sub("\\..*$", "", rownames(counts))
counts <- rowsum(as.matrix(counts), group = ens)
counts <- counts[rowSums(counts) > 0, ]
storage.mode(counts) <- "integer"
colnames(counts) <- norm16(colnames(counts))

# ------------------------------------------------------------
# Protein-coding filter
# ------------------------------------------------------------
# Protein-coding filter
gene_biotype <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = rownames(counts),
  keytype = "GENEID",
  column = "GENEBIOTYPE",
  multiVals = "first"
)

counts <- counts[gene_biotype == "protein_coding", ]
counts <- counts[!is.na(rownames(counts)), ]


# ------------------------------------------------------------
# Clinical
# ------------------------------------------------------------
clin <- as.data.frame(colData(se))

if ("patient" %in% names(clin)) {
  clin$patient12 <- norm12(clin$patient)
} else if ("bcr_patient_barcode" %in% names(clin)) {
  clin$patient12 <- norm12(clin$bcr_patient_barcode)
} else {
  clin$patient12 <- norm12(rownames(clin))
}

sex_col <- intersect(c("gender", "sex"), names(clin))
clin$sex <- tolower(trimws(as.character(clin[[sex_col[1]]])))

clin <- clin %>% filter(sex %in% c("female", "f"))

# ------------------------------------------------------------
# ER annotation
# Load PAM50 subtype data
# ------------------------------------------------------------
sub <- TCGAquery_subtype("BRCA")

sub$patient12 <- norm12(sub$patient)
clin_er <- sub %>%
  transmute(
    patient12,
    ER = case_when(
      BRCA_Subtype_PAM50 %in% c("LumA", "LumB") ~ "ER+",
      BRCA_Subtype_PAM50 == "Basal" ~ "ER-",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ER)) %>%   # removes HER2 + Normal
  distinct()


# ------------------------------------------------------------
# Merge ER into clinical
# ------------------------------------------------------------

clin <- clin %>%
  left_join(clin_er, by = "patient12") %>%
  mutate(ER = ifelse(patient12 == "TCGA-B6-A0I8", "ER+", ER)) %>%  # ← bæta við hér
  filter(ER == "ER+") %>%
  distinct(patient12, .keep_all = TRUE)

# ------------------------------------------------------------
# Add age (NO JOIN)
# ------------------------------------------------------------
clin <- clin %>%
  mutate(
    age_num = if ("age_at_diagnosis" %in% names(clin)) {
      as.numeric(age_at_diagnosis) / 365.25
    } else {
      abs(as.numeric(days_to_birth)) / 365.25
    }
  )

# ------------------------------------------------------------
# Build sample annotation
# ------------------------------------------------------------
sample_annot <- tibble(
  aliquot   = norm16(colnames(counts)),
  patient12 = norm12(colnames(counts))
) %>%
  left_join(clin, by = "patient12") %>%
  filter(!is.na(ER))   



# ------------------------------------------------------------
# Keep primary tumors
# ------------------------------------------------------------
sample_annot <- sample_annot %>%
  mutate(sample_type = substr(aliquot, 14, 15)) %>%
  filter(sample_type == "01")

# ------------------------------------------------------------
# Align counts
# ------------------------------------------------------------
counts <- counts[, sample_annot$aliquot, drop = FALSE]
stopifnot(all(colnames(counts) == sample_annot$aliquot))


cat("\nFinal dataset:\n")
print(dim(counts))
print(table(sample_annot$ER))


# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
saveRDS(counts, "data/counts_precollapse_TCGA.rds")
saveRDS(sample_annot, "data/annot_precollapse_TCGA.rds")

cat("DONE - precollapse dataset ready\n")
