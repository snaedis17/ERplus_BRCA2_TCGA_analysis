# ============================================================
 Build FINAL ER+ TCGA cohort
# ============================================================
#This script builds the final ER-positive TCGA breast cancer cohort used for all downstream analyses (DESeq2, GSEA, etc.).

#It performs three main steps:
#	1.	Cleans the dataset (removes unwanted samples)
#2.	Defines BRCA2 mutation status (Mut vs WT)
#3.	Constructs a patient-level expression dataset

#### Input data

##counts_precollapse_TCGA.rds

#Gene expression matrix (RNA-seq counts)
	#•	Rows = genes
	#•	Columns = aliquots (individual tumor samples)
	#•	Raw, unfiltered dataset

#Precollapse means:
	#•	multiple samples per patient may exist
 #•	no filtering has been applied yet
	#•	includes all tumor types before cleaning

##Sample annotation (metadata)

#Contains clinical and sample-level information:
	#•	aliquot → sample ID (matches columns in counts matrix)
	#•	patient12 → TCGA patient ID (first 12 characters of barcode)
	#•	ER → estrogen receptor status (ER+ / ER−)
	#•	age_num → age at diagnosis

#maf_TCGA_BRCA.rds = Somatic mutation data (MAF file)

#Contains mutation information per tumor sample:
	#•	Hugo_Symbol → gene name (e.g. BRCA1, BRCA2)
	#•	Variant_Classification → mutation type
	#•	Tumor_Sample_Barcode → sample ID

##Nik_germline
##maxwell_germline
#These contain germline mutation annotations from external sources.

# Used to identify: pathogenic BRCA1 carriers (to exclude) and	pathogenic BRCA2 carriers (to define the “Mut” group)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

filter   <- dplyr::filter
select   <- dplyr::select
arrange  <- dplyr::arrange
distinct <- dplyr::distinct
count <- dplyr::count
pull  <- dplyr::pull


norm12 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 12)
norm16 <- function(x) substr(toupper(gsub("\\s+", "", as.character(x))), 1, 16)

# ------------------------------------------------------------
# Load data (from Script 01)
# ------------------------------------------------------------
counts <- readRDS("data/counts_precollapse_TCGA.rds")
annot  <- readRDS("data/annot_precollapse_TCGA.rds")
maf    <- readRDS("data/maf_TCGA_BRCA.rds")


# external
Nik_germline <- yost_table_12
names(Nik_germline) <- make.names(names(Nik_germline))
names(maxwell_germline) <- make.names(names(maxwell_germline))

# ------------------------------------------------------------
# 1) Restrict to ER+
# ------------------------------------------------------------
annot <- annot %>% filter(ER == "ER+")
counts <- counts[, annot$aliquot, drop = FALSE]

# ------------------------------------------------------------
# 2) REMOVE BRCA1 (germline + somatic)
# ------------------------------------------------------------

# Germline (Nik)
brca1_nik <- Nik_germline %>%
  filter(toupper(Gene) == "BRCA1",
         GermlineVariantGroup == "PathogenicMutation") %>%
  transmute(patient12 = norm12(ID)) %>%
  distinct()

# Germline (Maxwell)
brca1_maxwell <- maxwell_germline %>%
  filter(grepl("^\\s*BRCA1\\b", Mutation, ignore.case = TRUE)) %>%
  transmute(patient12 = norm12(Tumor.ID)) %>%
  distinct()

# Somatic
brca1_lof <- maf %>%
  filter(Hugo_Symbol == "BRCA1",
         Variant_Classification %in% c(
           "Frame_Shift_Del","Frame_Shift_Ins",
           "Nonsense_Mutation","Splice_Site","Nonstop_Mutation"
         ))

somatic_brca1_ids <- unique(norm12(brca1_lof$Tumor_Sample_Barcode))

# Combine
brca1_all_ids <- unique(c(
  brca1_nik$patient12,
  brca1_maxwell$patient12,
  somatic_brca1_ids
))

# Filter
annot <- annot %>% filter(!patient12 %in% brca1_all_ids)
counts <- counts[, annot$aliquot, drop = FALSE]

# ------------------------------------------------------------
# 3) DEFINE BRCA2 GERMLINE (Mut group)
# ------------------------------------------------------------

brca2_nik <- Nik_germline %>%
  filter(toupper(Gene) == "BRCA2",
         GermlineVariantGroup == "PathogenicMutation") %>%
  transmute(patient12 = norm12(ID)) %>%
  distinct()

brca2_maxwell <- maxwell_germline %>%
  filter(grepl("^\\s*BRCA2\\b", Mutation, ignore.case = TRUE)) %>%
  transmute(patient12 = norm12(Tumor.ID)) %>%
  distinct()

brca2_all <- bind_rows(brca2_nik, brca2_maxwell) %>%
  distinct(patient12)

annot <- annot %>%
  mutate(
    BRCA2_Status = ifelse(
      patient12 %in% brca2_all$patient12,
      "Mut",
      "WT"
    )
  )

# ------------------------------------------------------------
# 4) REMOVE BRCA2 SOMATIC FROM WT
# ------------------------------------------------------------

brca2_lof <- maf %>%
  filter(Hugo_Symbol == "BRCA2",
         Variant_Classification %in% c(
           "Frame_Shift_Del","Frame_Shift_Ins",
           "Nonsense_Mutation","Splice_Site","Nonstop_Mutation"
         ))

somatic_brca2_ids <- unique(norm12(brca2_lof$Tumor_Sample_Barcode))

annot <- annot %>%
  filter(!(BRCA2_Status == "WT" & patient12 %in% somatic_brca2_ids))

counts <- counts[, annot$aliquot, drop = FALSE]

# ------------------------------------------------------------
# 5) REMOVE DUPLICATE PATIENTS (REMOVE ALL)
# ------------------------------------------------------------

dup_patients <- annot %>%
  count(patient12) %>%
  filter(n > 1) %>%
  pull(patient12)

annot <- annot %>% filter(!patient12 %in% dup_patients)
counts <- counts[, annot$aliquot, drop = FALSE]

cat("Removed duplicate patients:", length(dup_patients), "\n")
#Removed duplicate patients: 6 
# ------------------------------------------------------------
# 6) COLLAPSE TO PATIENT LEVEL
# ------------------------------------------------------------

# map aliquot -> patient
aliquot2pat <- tibble(
  aliquot   = norm16(colnames(counts)),
  patient12 = norm12(colnames(counts))
)

counts_one <- counts
colnames(counts_one) <- aliquot2pat$patient12

# patient-level annotation
annot_patient <- annot %>%
  distinct(patient12, .keep_all = TRUE) %>%
  select(patient12, ER, BRCA2_Status, age_num)

# ensure matching order
keep <- intersect(colnames(counts_one), annot_patient$patient12)

counts_one <- counts_one[, keep, drop = FALSE]

annot_patient <- annot_patient %>%
  filter(patient12 %in% keep) %>%
  arrange(match(patient12, keep)) %>%
  as.data.frame()

rownames(annot_patient) <- annot_patient$patient12

# ------------------------------------------------------------
# SAVE FINAL DATA
# ------------------------------------------------------------
saveRDS(counts_one,    "data/counts_ERpos_clean_final.rds")
saveRDS(annot_patient,"data/annot_ERpos_clean_final.rds")

cat("DONE\n")
print(table(annot_patient$BRCA2_Status))
print(dim(counts_one))





