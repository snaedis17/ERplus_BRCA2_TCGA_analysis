# ER+ BRCA2 Breast Cancer Transcriptomics

Transcriptomic analysis of ER-positive breast cancer in TCGA.

RNA-seq analysis of TCGA-BRCA tumors focusing on ER+ BRCA2 mutation carriers.

---

## Objectives

1. Compare ER+ BRCA2-mutated vs wild-type tumors
2. Compare young vs older ER+ tumors
3. Identify pathways associated with aggressive ER+ BRCA2 tumors

---

## Data

Source: TCGA-BRCA RNA-seq (STAR counts)

Downloaded using:
TCGAbiolinks

---

## Pipeline

1. Download TCGA RNA-seq data  
2. Clinical filtering (ER+, female)  
3. Patient-level count collapsing  
4. DESeq2 differential expression  
5. Hallmark pathway analysis (fgsea)  
6. ssGSEA pathway scoring  

---

## Additional metadata

BRCA2 mutation carriers were defined using curated mutation annotations.

The file `data/brca2_carriers_tcga.csv` contains the list of BRCA2 mutation carriers used in this analysis.

## Author

Snædís Ragnarsdóttir  
PhD student – University of Iceland
