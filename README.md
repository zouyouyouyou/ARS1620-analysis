# ARS-1620 Proteomics and Activity Analysis

This repository contains R code for analysis of ARS-1620 experimental data, including:

- one-pot SPROX analysis
- one-pot TPP analysis
- volcano plot generation
- ALDH1A3 enzymatic activity analysis

The workflow includes data import, filtering, normalization, Welch t tests, hit selection, UniProt annotation, and figure generation.

---

## Repository contents

This project includes analysis for:

1. UniProt annotation import and preprocessing
2. One-pot SPROX data analysis
3. One-pot TPP data analysis
4. Volcano plots using Z-scores
5. Volcano plots using log2 fold changes
6. ALDH1A3 enzymatic activity assay analysis

---

## Requirements

### Software

- R
- RStudio recommended

### Required R packages

```r
library(tidyverse)
library(tidyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(rentrez)
library(Biostrings)
library(openxlsx)
library(ggrepel)
library(readr)
library(stringr)
library(forcats)
library(RColorBrewer)
```

## UniProt annotation and one-pot SPROX analysis

### Read UniProt annotation file

The analysis begins by reading the UniProt human proteome annotation file and keeping the accession and gene name fields for downstream annotation.

```r
uniprot_proteome_Homo_Sapiens <- read_tsv("SPROX/Data/uniprotkb_proteome_UP000005640_2024_09_23.tsv")

names(uniprot_proteome_Homo_Sapiens)[1] <- "Master_Protein_Accession"
names(uniprot_proteome_Homo_Sapiens)[3] <- "Protein_names"
names(uniprot_proteome_Homo_Sapiens)[4] <- "Gene_Name"

uniprot_proteome <- uniprot_proteome_Homo_Sapiens %>%
  dplyr::select(Master_Protein_Accession, Gene_Name) %>%
  filter(!duplicated(.))
```

### Read one-pot SPROX data

```r
SPROX_Met_Enrich <- read_excel("SPROX/Data/ARS_1620_Enrich_non_scale.xlsx")
names(SPROX_Met_Enrich)[3] <- "Annotated_Sequence"
names(SPROX_Met_Enrich)[10] <- "Master_Protein_Accessions"
names(SPROX_Met_Enrich)[11] <- "Positions_in_Master_Proteins"
```


