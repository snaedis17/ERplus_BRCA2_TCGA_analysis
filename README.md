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

The files `yost_table_12.xlsx` and `maxwell_germline.xlsx` contain the list of BRCA2 mutation carriers used in this analysis.

---
## Results

### Differential expression: BRCA2-mutated vs WT (ER+)

Volcano plot showing differential gene expression between BRCA2-mutated and wild-type ER-positive breast tumors.

![Volcano](results/brca2/plots/BRCA2_volcano_cleanWT.png)

**Key observation:**  
Differential expression is driven by a relatively small subset of genes with strong effect sizes, rather than a global shift in the transcriptome. This is consistent with a structured but subtle transcriptional phenotype in BRCA2-mutated tumors.

---

### Biological interpretation

These findings align with a model where BRCA2-mutated ER+ tumors do not form a completely distinct transcriptional cluster, but instead exhibit targeted dysregulation of specific pathways. As shown in our analysis  [oai_citation:0‡GREIN 3_nyjast.docx](file-service://file-932ajgEjLdVrL1r2qXT4H3), this type of pattern reflects coordinated biological programs (e.g., proliferation, metabolic rewiring, EMT) rather than widespread gene-level changes.


### PCA – global transcriptomic structure

Samples do not clearly separate by BRCA2 status or age at a global level.

![PCA](results/age_interaction/plots/PCA_age_brca2.png)

---

### Interaction: Age × BRCA2

A subset of genes shows expression patterns that depend on both age and BRCA2 status.

![Interaction genes](results/age_interaction/plots/boxplot_interaction_genes.png)

---

### BRCA2 expression (control check)

BRCA2 expression is lower in mutated tumors, confirming expected biology.

![BRCA2 expression](results/age_interaction/plots/BRCA2_expression.png)
## Author



![test](results/age_interaction/plots/PCA_age_brca2.png) 
Snædís Ragnarsdóttir  
PhD student – University of Iceland
