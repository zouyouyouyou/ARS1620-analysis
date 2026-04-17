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
