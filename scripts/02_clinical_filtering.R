# ============================================================
# Clinical annotation and filtering
# TCGA BRCA ER+ BRCA2 analysis
#
# Steps
# 1. Load TCGA SummarizedExperiment
# 2. Extract clinical metadata
# 3. Restrict to female samples
# 4. Determine ER status
#    - clinical ER status if available
#    - PAM50 subtype fallback
# 5. Annotate BRCA2 mutation carriers using Maxwell file
# 6. Build patient-level clinical table
# 7. Build sample annotation
# 8. Identify duplicate aliquots
#
# Outputs
# data_processed/counts_unique.rds
# data_processed/sample_annot.rds
# data_processed/clin2_unique.rds
# data_processed/dup_info.rds
# ============================================================


suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(readr)
})

source("scripts/00_utils.R")

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------
se <- readRDS("data_raw/tcga_brca_se.rds")

# Maxwell file:

library(readxl)

Maxwell_BRCA2uppl <- read_excel("data/Maxwell_BRCA2uppl.xlsx")


# ------------------------------------------------------------
# Counts matrix
# ------------------------------------------------------------
count_matrix <- assay(se)

ens <- sub("\\..*$", "", rownames(count_matrix))
counts_unique <- rowsum(as.matrix(count_matrix), group = ens, reorder = FALSE)
counts_unique <- counts_unique[rowSums(counts_unique) > 0, , drop = FALSE]
storage.mode(counts_unique) <- "integer"
colnames(counts_unique) <- norm16(colnames(counts_unique))

message("#genes: ", nrow(counts_unique), "  #samples: ", ncol(counts_unique))

# ------------------------------------------------------------
# Clinical table
# ------------------------------------------------------------
clin <- as.data.frame(colData(se))

if ("patient" %in% names(clin)) {
  clin$patient12 <- norm12(clin$patient)
} else if ("bcr_patient_barcode" %in% names(clin)) {
  clin$patient12 <- norm12(clin$bcr_patient_barcode)
} else {
  clin$patient12 <- norm12(rownames(clin))
}

# ------------------------------------------------------------
# Female only
# ------------------------------------------------------------
sex_col <- intersect(c("gender", "sex"), names(clin))
stopifnot(length(sex_col) > 0)

clin$sex <- tolower(trimws(as.character(clin[[sex_col[1]]])))
clin <- clin[!is.na(clin$sex) & clin$sex %in% c("female", "f"), , drop = FALSE]

cat("Female N:", nrow(clin), "\n")

# ------------------------------------------------------------
# ER status from clinical
# ------------------------------------------------------------
cli_raw <- NULL
suppressWarnings(
  try({ cli_raw <- GDCquery_clinic("TCGA-BRCA", type = "clinical") }, silent = TRUE)
)

get_er_from_clinical <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(NULL)

  idcol <- NULL
  for (cand in c("bcr_patient_barcode", "case_submitter_id", "submitter_id")) {
    if (cand %in% names(df)) {
      idcol <- cand
      break
    }
  }
  if (is.null(idcol)) return(NULL)

  df$patient12 <- norm12(df[[idcol]])

  cand <- names(df)[
    grepl("(^|_)(er|estrogen).*status", names(df), ignore.case = TRUE) |
      grepl("^er[._]?status$", names(df), ignore.case = TRUE) |
      grepl("breast_carcinoma_estrogen_receptor_status", names(df), ignore.case = TRUE)
  ]

  if (!length(cand)) return(NULL)

  er_col <- cand[1]
  er_raw <- tolower(trimws(as.character(df[[er_col]])))

  ER <- ifelse(
    er_raw %in% c("positive", "+", "pos", "1", "true", "yes", "positive or strong"), "ER+",
    ifelse(
      er_raw %in% c("negative", "-", "neg", "0", "false", "no"), "ER-",
      ifelse(grepl("pos", er_raw), "ER+",
             ifelse(grepl("neg", er_raw), "ER-", NA))
    )
  )

  out <- df %>%
    transmute(patient12, ER = ER) %>%
    distinct(patient12, .keep_all = TRUE)

  attr(out, "er_source") <- paste0("clinical:", er_col)
  out
}

er_from_clin <- get_er_from_clinical(cli_raw)

# ------------------------------------------------------------
# PAM50 fallback
# ------------------------------------------------------------
sub <- NULL
suppressWarnings(try({ sub <- TCGAquery_subtype("BRCA") }, silent = TRUE))

if (!is.null(sub) && is.data.frame(sub)) {
  sub$patient12 <- norm12(sub$patient)
  pam_col <- intersect(c("BRCA_Subtype_PAM50", "PAM50", "PAM50.mRNA", "pam50"), names(sub))
  if (length(pam_col)) {
    sub$PAM50 <- toupper(trimws(as.character(sub[[pam_col[1]]])))
  } else {
    sub$PAM50 <- NA_character_
  }
} else {
  sub <- data.frame(patient12 = character(0), PAM50 = character(0))
}

er_from_pam50 <- sub %>%
  transmute(
    patient12,
    ER = case_when(
      PAM50 %in% c("LUMINAL A", "LUMINALA", "LUMA", "LUMINAL B", "LUMINALB", "LUMB") ~ "ER+",
      PAM50 %in% c("BASAL", "BASAL-LIKE", "BASALLIKE") ~ "ER-",
      TRUE ~ NA_character_
    ),
    PAM50 = PAM50
  ) %>%
  distinct(patient12, .keep_all = TRUE)

attr(er_from_pam50, "er_source") <- "PAM50 (proxy; HER2 set to NA)"

# ------------------------------------------------------------
# Merge ER calls
# ------------------------------------------------------------
clin_er <- clin[, c("patient12", "sex"), drop = FALSE]

if (!is.null(er_from_clin)) {
  clin_er <- clin_er %>%
    left_join(er_from_clin, by = "patient12") %>%
    rename(ER_clin = ER)
} else {
  clin_er$ER_clin <- NA_character_
}

clin_er <- clin_er %>%
  left_join(er_from_pam50[, c("patient12", "ER", "PAM50")], by = "patient12") %>%
  mutate(
    ER = ifelse(!is.na(ER_clin), ER_clin, ER),
    ER_source = ifelse(
      !is.na(ER_clin),
      if (!is.null(er_from_clin)) attr(er_from_clin, "er_source") else "clinical",
      attr(er_from_pam50, "er_source")
    )
  ) %>%
  dplyr::select(patient12, sex, ER, PAM50, ER_source)

cat("ER counts (female only):\n")
print(table(clin_er$ER, useNA = "ifany"))

# ------------------------------------------------------------
# Age
# ------------------------------------------------------------
age_years <- rep(NA_real_, nrow(clin))

if ("age_at_diagnosis" %in% names(clin)) {
  age_years <- suppressWarnings(as.numeric(clin$age_at_diagnosis) / 365.25)
} else if ("days_to_birth" %in% names(clin)) {
  age_years <- suppressWarnings(abs(as.numeric(clin$days_to_birth)) / 365.25)
}

clin_er$age_num <- age_years[match(clin_er$patient12, clin$patient12)]

# ------------------------------------------------------------
# Maxwell BRCA2 carriers
# ------------------------------------------------------------
mx <- Maxwell_BRCA2uppl
names(mx) <- make.names(names(mx))

stopifnot(all(c("Tumor.ID", "Mutation") %in% names(mx)))

maxwell_carriers <- mx %>%
  transmute(
    patient12 = norm12(Tumor.ID),
    BRCA2 = grepl("^\\s*BRCA2\\b", Mutation, ignore.case = TRUE)
  ) %>%
  filter(BRCA2) %>%
  distinct(patient12) %>%
  mutate(BRCA2_Status = factor("Mut", levels = c("WT", "Mut")))

clin2 <- clin_er %>%
  left_join(maxwell_carriers[, c("patient12", "BRCA2_Status")], by = "patient12") %>%
  mutate(
    BRCA2_Status = ifelse(is.na(BRCA2_Status), "WT", as.character(BRCA2_Status)),
    BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut"))
  )

cat("BRCA2 (female):\n")
print(table(clin2$BRCA2_Status, useNA = "ifany"))

# ------------------------------------------------------------
# Collapse duplicates to patient-level
# ------------------------------------------------------------
to_chr <- function(x) {
  if (is.list(x)) vapply(x, function(z) as.character(z)[1], character(1)) else as.character(x)
}

to_num <- function(x) {
  if (is.list(x)) vapply(x, function(z) suppressWarnings(as.numeric(z)[1]), numeric(1)) else suppressWarnings(as.numeric(x))
}

clin2_unique <- clin2 %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(
    patient12    = to_chr(patient12),
    ER           = to_chr(ER),
    BRCA2_Status = to_chr(BRCA2_Status),
    age_num      = to_num(age_num)
  ) %>%
  group_by(patient12) %>%
  summarise(
    ER = {
      vals <- unique(na.omit(ER))
      if (length(vals) == 1) vals else NA_character_
    },
    age_num = median(age_num, na.rm = TRUE),
    BRCA2_Status = if (any(BRCA2_Status == "Mut", na.rm = TRUE)) "Mut" else "WT",
    .groups = "drop"
  ) %>%
  mutate(BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut")))

cat("# collapsed patients:", nrow(clin2_unique), "\n")

# ------------------------------------------------------------
# Build sample annotation
# ------------------------------------------------------------
aliquots <- norm16(colnames(counts_unique))
patients <- norm12(colnames(counts_unique))

sample_annot <- tibble(
  aliquot   = aliquots,
  patient12 = patients
) %>%
  left_join(clin2_unique, by = "patient12") %>%
  mutate(sample_type = substr(aliquot, 14, 15)) %>%
  filter(sample_type == "01")

counts_unique <- counts_unique[, sample_annot$aliquot, drop = FALSE]

stopifnot(identical(colnames(counts_unique), sample_annot$aliquot))

cat("ER coverage:\n")
print(table(sample_annot$ER, useNA = "ifany"))

cat("BRCA2 status:\n")
print(table(sample_annot$BRCA2_Status, useNA = "ifany"))

# ------------------------------------------------------------
# Duplicate aliquot check after filtering
# ------------------------------------------------------------
aliquot2pat_filtered <- data.frame(
  aliquot   = colnames(counts_unique),
  patient12 = substr(colnames(counts_unique), 1, 12)
)

cat("Filtered matrix: aliquots per patient\n")
print(table(table(aliquot2pat_filtered$patient12)))

dup_patients <- aliquot2pat_filtered %>%
  as_tibble() %>%
  count(patient12, name = "n_aliquots") %>%
  filter(n_aliquots > 1)

dup_info <- dup_patients %>%
  left_join(
    sample_annot %>%
      distinct(patient12, .keep_all = TRUE) %>%
      select(patient12, BRCA2_Status, age_num, ER),
    by = "patient12"
  ) %>%
  mutate(
    AgeGroup = case_when(
      is.na(age_num) ~ NA_character_,
      age_num <= 40 ~ "Young",
      age_num > 40 ~ "Older"
    )
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
if (!dir.exists("data_processed")) dir.create("data_processed", recursive = TRUE)

saveRDS(counts_unique, "data_processed/counts_unique.rds")
saveRDS(clin2_unique, "data_processed/clin2_unique.rds")
saveRDS(sample_annot, "data_processed/sample_annot.rds")
saveRDS(dup_info, "data_processed/dup_info.rds")

write_csv(as_tibble(clin2_unique), "data_processed/clin2_unique.csv")
write_csv(as_tibble(sample_annot), "data_processed/sample_annot.csv")
write_csv(as_tibble(dup_info), "data_processed/dup_info.csv")

message("✓ Wrote:
- data_processed/counts_unique.rds
- data_processed/clin2_unique.rds
- data_processed/sample_annot.rds
- data_processed/dup_info.rds")
