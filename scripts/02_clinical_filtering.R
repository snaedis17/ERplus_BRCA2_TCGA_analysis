# ------------------------------------------------------------
# OUTPUT:
# - counts_ERpos_TCGA.rds
# - annot_ERpos_TCGA.rds
# - counts_precollapse_TCGA.rds
# - annot_precollapse_TCGA.rds
#
# DESCRIPTION:
# Prepares TCGA-BRCA RNA-seq data:
# - gene filtering (protein-coding + low count)
# - clinical + ER annotation
# - BRCA2 germline annotation
# - primary tumor selection
# - ER+ patient-level dataset
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(matrixStats)
})

filter   <- dplyr::filter
select   <- dplyr::select
arrange  <- dplyr::arrange
distinct <- dplyr::distinct

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
norm12 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 12)
norm16 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 16)

# as.character(x):ensures input is a text
#gsub("\\s+", "", x): removes all spaces 
# toupper(...): makes everything uppercase
#I standardize TCGA barcodes to ensure consistent matching between datasets.
#I use 12-character IDs for patient-level data like clinical variables and mutation status,and 16-character IDs for sample-level expression data.This prevents mismatches and ensures that each sample is correctly linked to the right patient

# ------------------------------------------------------------
# 1) Download TCGA-BRCA STAR counts
# ------------------------------------------------------------
message("Querying TCGA-BRCA STAR counts...")

q <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(q)
se <- GDCprepare(q)

# SAVE RAW SummarizedExperiment

saveRDS(se, "data/se_TCGA_raw.rds")

count_matrix <- assay(se)


# remove ENSEMBL version and collapse duplicate ENSG rows
ens <- sub("\\..*$", "", rownames(count_matrix))
counts_unique <- rowsum(as.matrix(count_matrix), group = ens, reorder = FALSE)
counts_unique <- counts_unique[rowSums(counts_unique) > 0, , drop = FALSE]
storage.mode(counts_unique) <- "integer"
colnames(counts_unique) <- norm16(colnames(counts_unique))

# SAVE CLEANED (no filters yet)
saveRDS(counts_unique, "data/counts_cleaned_TCGA.rds")

cat("Initial matrix dimensions:", dim(counts_unique), "\n")

#First, I clean the gene identifiers by removing ENSEMBL version numbers, since the same gene can appear multiple times with different versions. I then collapse duplicates by summing counts to obtain one value per gene. After that, I remove genes with zero expression across all samples and ensure the matrix is stored as integers, as required by DESeq2. Finally, I standardize sample IDs to maintain consistency across datasets.

# ------------------------------------------------------------
# 2) Protein-coding filter FIRST
# ------------------------------------------------------------
ens_ids <- rownames(counts_unique) #Get gene IDs

gene_biotype <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = ens_ids,
  keytype = "GENEID",
  column = "GENEBIOTYPE",
  multiVals = "first"
)
#Get gene type annotation

keep_pc <- !is.na(gene_biotype) & gene_biotype == "protein_coding" #Keep only protein-coding
counts_unique <- counts_unique[keep_pc, , drop = FALSE]
#	Filter matrix, remove: lncRNA, pseudogenes, miRNA, non-coding genes


cat("After protein-coding filter:", dim(counts_unique), "\n")



# ------------------------------------------------------------
# 3) Low-count filter FIRST
#    Rule: keep genes with >=10 counts in >=50 samples
# ------------------------------------------------------------
keep_lc <- rowSums(counts_unique >= 10) >= 50
counts_unique <- counts_unique[keep_lc, , drop = FALSE]

# SAVE FILTERED (ready for analysis)
saveRDS(counts_unique, "data/counts_filtered_TCGA.rds")

cat("After low-count filter:", dim(counts_unique), "\n")


#“I remove lowly expressed genes by requiring at least 10 counts in at least 50 samples. This reduces noise and improves statistical power, as genes with very low expression tend to have unstable variance and can lead to false positives in differential expression analysis.”
# ------------------------------------------------------------
# 4) Clinical annotation: extract, clean, and restrict cohort
# ------------------------------------------------------------

# Extract clinical/sample metadata from SummarizedExperiment
# colData(se) contains annotation for each sample (columns in count matrix)
clin <- as.data.frame(colData(se))

# ------------------------------------------------------------
# Create standardized patient ID (12-character TCGA barcode)
# TCGA data may store patient IDs under different column names,
# so we check multiple possibilities and normalize them
# ------------------------------------------------------------
if ("patient" %in% names(clin)) {
  clin$patient12 <- norm12(clin$patient)
} else if ("bcr_patient_barcode" %in% names(clin)) {
  clin$patient12 <- norm12(clin$bcr_patient_barcode)
} else {
  clin$patient12 <- norm12(rownames(clin))
}

# ------------------------------------------------------------
# Identify sex column (can be named "gender" or "sex")
# ------------------------------------------------------------
sex_col <- intersect(c("gender", "sex"), names(clin))

# Stop execution if no sex column is found (data integrity check)
stopifnot(length(sex_col) > 0)

# ------------------------------------------------------------
# Clean sex values:
# - convert to character
# - trim whitespace
# - convert to lowercase
# ------------------------------------------------------------
clin$sex <- tolower(trimws(as.character(clin[[sex_col[1]]])))

# ------------------------------------------------------------
# Restrict to female patients only
# This ensures a consistent cohort for breast cancer analysis
# and removes rare male breast cancer cases
# ------------------------------------------------------------
clin <- clin[
  !is.na(clin$sex) & clin$sex %in% c("female", "f"),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Print number of female samples retained
# ------------------------------------------------------------
cat("Female N:", nrow(clin), "\n")

# ------------------------------------------------------------
# 5) ER annotation: clinical extraction with PAM50 fallback
# ------------------------------------------------------------

# ------------------------------------------------------------
# Attempt to download detailed clinical data from TCGA
# This may fail depending on connection/API, so we wrap in try()
# ------------------------------------------------------------
cli_raw <- NULL
suppressWarnings(
  try({ cli_raw <- GDCquery_clinic("TCGA-BRCA", type = "clinical") }, silent = TRUE)
)

# ------------------------------------------------------------
# Function to extract ER status from clinical data
# Handles inconsistent column naming and value formats
# ------------------------------------------------------------
get_er_from_clinical <- function(df) {
  
  # Return NULL if input is missing or invalid
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  
  # ----------------------------------------------------------
  # Identify patient ID column (varies across TCGA tables)
  # ----------------------------------------------------------
  idcol <- NULL
  for (cand in c("bcr_patient_barcode", "case_submitter_id", "submitter_id")) {
    if (cand %in% names(df)) {
      idcol <- cand
      break
    }
  }
  if (is.null(idcol)) return(NULL)
  
  # Standardize patient ID to 12-character TCGA barcode
  df$patient12 <- norm12(df[[idcol]])
  
  # ----------------------------------------------------------
  # Identify ER status column using flexible pattern matching
  # ----------------------------------------------------------
  cand <- names(df)[
    grepl("(^|_)(er|estrogen).*status", names(df), ignore.case = TRUE) |
      grepl("^er[._]?status$", names(df), ignore.case = TRUE) |
      grepl("breast_carcinoma_estrogen_receptor_status", names(df), ignore.case = TRUE)
  ]
  if (!length(cand)) return(NULL)
  
  er_col <- cand[1]
  
  # ----------------------------------------------------------
  # Clean ER values (standardize text formatting)
  # ----------------------------------------------------------
  er_raw <- tolower(trimws(as.character(df[[er_col]])))
  
  # ----------------------------------------------------------
  # Convert raw values to binary ER status (ER+ / ER-)
  # Handles multiple possible encodings
  # ----------------------------------------------------------
  ER <- ifelse(
    er_raw %in% c("positive", "+", "pos", "1", "true", "yes", "positive or strong"), "ER+",
    ifelse(
      er_raw %in% c("negative", "-", "neg", "0", "false", "no"), "ER-",
      ifelse(
        grepl("pos", er_raw), "ER+",
        ifelse(grepl("neg", er_raw), "ER-", NA)
      )
    )
  )
  
  # ----------------------------------------------------------
  # Keep one ER annotation per patient
  # ----------------------------------------------------------
  out <- df %>%
    transmute(patient12, ER = ER) %>%
    distinct(patient12, .keep_all = TRUE)
  
  # Store metadata about source column
  attr(out, "er_source") <- paste0("clinical:", er_col)
  
  return(out)
}


# Extract ER from clinical data
er_from_clin <- get_er_from_clinical(cli_raw)

# ------------------------------------------------------------
# PAM50 fallback (proxy for ER status)
# Used when clinical ER is missing
# ------------------------------------------------------------
sub <- NULL
suppressWarnings(try({ sub <- TCGAquery_subtype("BRCA") }, silent = TRUE))

if (!is.null(sub) && is.data.frame(sub)) {
  
  sub$patient12 <- norm12(sub$patient)
  
  # Identify PAM50 subtype column
  pam_col <- intersect(c("BRCA_Subtype_PAM50", "PAM50", "PAM50.mRNA", "pam50"), names(sub))
  
  if (length(pam_col)) {
    sub$PAM50 <- toupper(trimws(as.character(sub[[pam_col[1]]])))
  } else {
    sub$PAM50 <- NA_character_
  }
  
} else {
  sub <- data.frame(patient12 = character(0), PAM50 = character(0))
}

# ------------------------------------------------------------
# Convert PAM50 subtype to ER status
# Luminal → ER+, Basal → ER-
# ------------------------------------------------------------
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

# Store source information
attr(er_from_pam50, "er_source") <- "PAM50 (proxy; HER2 set to NA)"

# ------------------------------------------------------------
# Merge ER annotations into clinical dataset
# Prioritize clinical ER over PAM50 proxy
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

# ------------------------------------------------------------
# Manual correction (known misannotation)
# ------------------------------------------------------------
clin_er$ER[clin_er$patient12 == "TCGA-EW-A1P7"] <- "ER+"

# ------------------------------------------------------------
# Print ER distribution
# ------------------------------------------------------------
cat("ER counts in female tumors:\n")
print(table(clin_er$ER, useNA = "ifany"))

# ER-  ER+ <NA> 
#215  870  132 

# ------------------------------------------------------------
# Extract age (convert days → years)
# ------------------------------------------------------------
age_years <- rep(NA_real_, nrow(clin))

if ("age_at_diagnosis" %in% names(clin)) {
  age_years <- suppressWarnings(as.numeric(clin$age_at_diagnosis) / 365.25)
} else if ("days_to_birth" %in% names(clin)) {
  age_years <- suppressWarnings(abs(as.numeric(clin$days_to_birth)) / 365.25)
}

clin_er$age_num <- age_years[match(clin_er$patient12, clin$patient12)]

# ------------------------------------------------------------
# Inspect age distribution
# ------------------------------------------------------------
summary(clin_er$age_num)
#> summary(clin_er$age_num)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  26.94   49.09   58.77   58.90   67.85   90.00      53 
# ------------------------------------------------------------
# 6) Annotate germline BRCA2 carriers and build patient-level cohort
# ------------------------------------------------------------

# Ensure that both external carrier tables are available in memory
# - Maxwell_BRCA2uppl: curated list of BRCA2 carriers from Maxwell et al.
# - yost_table: additional germline variant table used to capture carriers
stopifnot(exists("Maxwell_BRCA2uppl"))
stopifnot(exists("yost_table"))

# ------------------------------------------------------------
# Extract BRCA2 germline carriers from the Maxwell table
# ------------------------------------------------------------
# Standardize column names and confirm required variables exist
mx <- Maxwell_BRCA2uppl
names(mx) <- make.names(names(mx))
stopifnot(all(c("Tumor.ID", "Mutation") %in% names(mx)))

# Keep only entries annotated as BRCA2 mutation carriers
# Tumor.ID is converted to a standardized 12-character TCGA patient ID
maxwell_carriers <- mx %>%
  transmute(
    patient12 = norm12(Tumor.ID),
    BRCA2 = grepl("^\\s*BRCA2\\b", Mutation, ignore.case = TRUE)
  ) %>%
  filter(BRCA2) %>%
  distinct(patient12)

# ------------------------------------------------------------
# Extract pathogenic BRCA2 germline carriers from the Yost table
# ------------------------------------------------------------
# Restrict to BRCA2, pathogenic germline variants, and breast cancer cases
yost_table <- yost_table
names(yost_table) <- make.names(names(yost_table))

yost_carriers <- yost_table %>%
  filter(
    Gene == "BRCA2",
    GermlineVariantGroup == "PathogenicMutation",
    Cancer == "BR"
  ) %>%
  transmute(patient12 = norm12(ID)) %>%
  distinct(patient12)


# ------------------------------------------------------------
# Combine carrier lists from both sources
# ------------------------------------------------------------
# Any patient present in either source is labeled as BRCA2-mutated
all_brca2_carriers <- bind_rows(maxwell_carriers, yost_carriers) %>%
  distinct(patient12) %>%
  mutate(BRCA2_Status = factor("Mut", levels = c("WT", "Mut")))

# ------------------------------------------------------------
# Merge BRCA2 carrier status into clinical annotation
# ------------------------------------------------------------
# Patients not found in the carrier list are labeled as wild-type (WT)
clin2 <- clin_er %>%
  left_join(all_brca2_carriers, by = "patient12") %>%
  mutate(
    BRCA2_Status = ifelse(is.na(BRCA2_Status), "WT", as.character(BRCA2_Status)),
    BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut"))
  )

# Manual ER correction for one known sample
clin2 <- clin2 %>%
  mutate(
    ER = ifelse(patient12 == "TCGA-B6-A0I8", "ER+", ER)
  )

cat("BRCA2 status in female tumors:\n")
print(table(clin2$BRCA2_Status, useNA = "ifany"))

#  WT  Mut 
#1192   25 

# ------------------------------------------------------------
# Collapse to one row per patient
# ------------------------------------------------------------
# Some imported objects may contain list-columns or duplicated entries,
# so convert values to simple character/numeric vectors first
to_chr <- function(x) {
  if (is.list(x)) vapply(x, function(z) as.character(z)[1], character(1)) else as.character(x)
}
to_num <- function(x) {
  if (is.list(x)) vapply(x, function(z) suppressWarnings(as.numeric(z)[1]), numeric(1)) else suppressWarnings(as.numeric(x))
}

# Build a unique patient-level clinical table
# - ER is kept only if all non-missing values agree
# - age is summarized as the median
# - BRCA2 status is set to Mut if any row indicates mutation
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
  mutate(
    BRCA2_Status = factor(BRCA2_Status, levels = c("WT", "Mut"))
  )

cat("Collapsed patients:", nrow(clin2_unique), "\n")
##Collapsed patients: 1082 
# ------------------------------------------------------------
# 7) Build sample annotation aligned to the count matrix
# ------------------------------------------------------------
# Derive both sample-level (16-character aliquot) and patient-level (12-character) IDs
aliquots <- norm16(colnames(counts_unique))
patients <- norm12(colnames(counts_unique))

sample_annot <- tibble(
  aliquot   = aliquots,
  patient12 = patients
) %>%
  left_join(clin2_unique, by = "patient12")

# Restrict to primary tumor samples only
# In TCGA barcodes, sample type "01" corresponds to primary solid tumor
sample_annot <- sample_annot %>%
  mutate(sample_type = substr(aliquot, 14, 15)) %>%
  filter(sample_type == "01")

# Keep only columns corresponding to retained primary tumor samples
counts_unique <- counts_unique[, sample_annot$aliquot, drop = FALSE]

# SAVE PRIMARY TUMORS (still includes duplicates)
saveRDS(counts_unique, "data/counts_primary_tumor_TCGA.rds")


# Confirm that sample annotation and count matrix are aligned
stopifnot(identical(colnames(counts_unique), sample_annot$aliquot))

cat("ER coverage after sample_type filter:\n")
print(table(sample_annot$ER, useNA = "ifany"))

#ER-  ER+ <NA> 
#197  780  134 

cat("BRCA2 coverage after sample_type filter:\n")
print(table(sample_annot$BRCA2_Status, useNA = "ifany"))

#WT  Mut <NA> 
#1077 21   13 
# ------------------------------------------------------------
# SAVE dataset used for duplicate QC (AFTER filtering, BEFORE collapsing)
# ------------------------------------------------------------
counts_precollapse <- counts_unique
annot_precollapse  <- sample_annot

saveRDS(counts_precollapse, "data/counts_precollapse_TCGA.rds")
saveRDS(annot_precollapse,  "data/annot_precollapse_TCGA.rds")



#stop her





