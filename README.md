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

### SPROX analysis

```r
# Read one-pot SPROX data
SPROX_Met_Enrich <- read_excel("SPROX/Data/ARS_1620_Enrich_non_scale.xlsx")

# rename the columns
names(SPROX_Met_Enrich)[3] <- "Annotated_Sequence"
names(SPROX_Met_Enrich)[10] <- "Master_Protein_Accessions" 
names(SPROX_Met_Enrich)[11] <- "Positions_in_Master_Proteins" 
names(SPROX_Met_Enrich)[15] <- "Control_1"
names(SPROX_Met_Enrich)[16] <- "Control_2"
names(SPROX_Met_Enrich)[17] <- "Control_3"
names(SPROX_Met_Enrich)[18] <- "Control_4"
names(SPROX_Met_Enrich)[19] <- "ten_uM_1"
names(SPROX_Met_Enrich)[20] <- "ten_uM_2"
names(SPROX_Met_Enrich)[21] <- "ten_uM_3"
names(SPROX_Met_Enrich)[22] <- "hundred_uM_1"
names(SPROX_Met_Enrich)[23] <- "hundred_uM_2"
names(SPROX_Met_Enrich)[24] <- "hundred_uM_3"

# select the columns for further analysis
SPROX_Met_Enrich_fil <- SPROX_Met_Enrich %>% 
  mutate(sequence = sub("^\\[.+\\]\\.(.+)\\..+$", "\\1", `Annotated_Sequence`)) %>% 
  dplyr::select("Master_Protein_Accessions", "Annotated_Sequence", "sequence", "Positions_in_Master_Proteins", "Modifications", 15:24)

# remove missing value, filter for wt-Met
SPROX_Met_Enrich_fil_wt_Met <- SPROX_Met_Enrich_fil %>%
  filter(if_all(c(Control_1, Control_2, Control_3, Control_4,
                  ten_uM_1, ten_uM_2, ten_uM_3,
                  hundred_uM_1, hundred_uM_2, hundred_uM_3), ~ !is.na(.))) %>% 
  filter(grepl("M", sequence)) %>% 
  filter(!grepl("Oxidation", Modifications))

# calculate the normalization factor, This normalization step mitigates the TMT channel to channel errors (e.g., differential sample loss and/or isobaric mass tag labeling efficiency). 
SPROX_norm_1 = sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_2 = sum(SPROX_Met_Enrich_fil_wt_Met$Control_2, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_3 = sum(SPROX_Met_Enrich_fil_wt_Met$Control_3, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_4 = sum(SPROX_Met_Enrich_fil_wt_Met$Control_4, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_5 = sum(SPROX_Met_Enrich_fil_wt_Met$ten_uM_1, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_6 = sum(SPROX_Met_Enrich_fil_wt_Met$ten_uM_2, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_7 = sum(SPROX_Met_Enrich_fil_wt_Met$ten_uM_3, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_8 = sum(SPROX_Met_Enrich_fil_wt_Met$hundred_uM_1, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_9 = sum(SPROX_Met_Enrich_fil_wt_Met$hundred_uM_2, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)
SPROX_norm_10 = sum(SPROX_Met_Enrich_fil_wt_Met$hundred_uM_3, na.rm = TRUE)/sum(SPROX_Met_Enrich_fil_wt_Met$Control_1, na.rm = TRUE)

# print nromalization factors
SPROX_norm_1 
SPROX_norm_2 
SPROX_norm_3 
SPROX_norm_4 
SPROX_norm_5
SPROX_norm_6 
SPROX_norm_7
SPROX_norm_8 
SPROX_norm_9 
SPROX_norm_10

> SPROX_norm_1 
[1] 1
> SPROX_norm_2 
[1] 0.9998314
> SPROX_norm_3 
[1] 1.003346
> SPROX_norm_4 
[1] 1.000769
> SPROX_norm_5
[1] 1.001683
> SPROX_norm_6 
[1] 0.998944
> SPROX_norm_7
[1] 1.001597
> SPROX_norm_8 
[1] 0.9992734
> SPROX_norm_9 
[1] 1.003763
> SPROX_norm_10 
[1] 1.003005

# apply normalization factor to each column
SPROX_Met_Enrich_fil_wt_Met_norm <- SPROX_Met_Enrich_fil_wt_Met %>% 
  mutate(Control_1 = Control_1/SPROX_norm_1) %>% 
  mutate(Control_2 = Control_2/SPROX_norm_2) %>% 
  mutate(Control_3 = Control_3/SPROX_norm_3) %>% 
  mutate(Control_4 = Control_4/SPROX_norm_4) %>% 
  mutate(ten_uM_1 = ten_uM_1/SPROX_norm_5) %>% 
  mutate(ten_uM_2 = ten_uM_2/SPROX_norm_6) %>% 
  mutate(ten_uM_3 = ten_uM_3/SPROX_norm_7) %>% 
  mutate(hundred_uM_1 = hundred_uM_1/SPROX_norm_8) %>% 
  mutate(hundred_uM_2 = hundred_uM_2/SPROX_norm_9) %>% 
  mutate(hundred_uM_3 = hundred_uM_3/SPROX_norm_10) %>% 
  mutate(Master_Protein_Accession = str_extract(Master_Protein_Accessions, "^[^;-]+")) %>% 
  dplyr::select(Master_Protein_Accession, 2:15) 

# perform Welch's t-test for wt Met
SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test <- SPROX_Met_Enrich_fil_wt_Met_norm %>% 
  rowwise() %>%
  mutate(
    # =======================
    # 10 uM: 6:9 vs 10:12
    # =======================
    fold_change_ten_uM = {
      x <- as.numeric(c_across(6:9))
      y <- as.numeric(c_across(10:12))
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) > 0 && length(y) > 0) mean(y) / mean(x) else NA_real_
    },
    
    log2FC_ten_uM = if (!is.na(fold_change_ten_uM) && fold_change_ten_uM > 0) {
      log2(fold_change_ten_uM)
    } else {
      NA_real_
    },
    
    p_value_ten_uM = {
      x <- as.numeric(c_across(6:9))
      y <- as.numeric(c_across(10:12))
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) >= 2 && length(y) >= 2) t.test(x, y)$p.value else NA_real_
    },
    
    neg_log10_p_value_ten_uM = if (!is.na(p_value_ten_uM) && p_value_ten_uM > 0) {
      -log10(p_value_ten_uM)
    } else {
      NA_real_
    },
    
    # =======================
    # 100 uM: 6:9 vs 13:15
    # =======================
    fold_change_hundred_uM = {
      x <- as.numeric(c_across(6:9))
      y <- as.numeric(c_across(13:15))
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) > 0 && length(y) > 0) mean(y) / mean(x) else NA_real_
    },
    
    log2FC_hundred_uM = if (!is.na(fold_change_hundred_uM) && fold_change_hundred_uM > 0) {
      log2(fold_change_hundred_uM)
    } else {
      NA_real_
    },
    
    p_value_hundred_uM = {
      x <- as.numeric(c_across(6:9))
      y <- as.numeric(c_across(13:15))
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) >= 2 && length(y) >= 2) t.test(x, y)$p.value else NA_real_
    },
    
    neg_log10_p_value_hundred_uM = if (!is.na(p_value_hundred_uM) && p_value_hundred_uM > 0) {
      -log10(p_value_hundred_uM)
    } else {
      NA_real_
    }
  ) %>%
  ungroup() %>% 
  mutate(
    # BH-FDR adjusted p values
    adj_p_value_ten_uM = p.adjust(p_value_ten_uM, method = "BH"),
    adj_p_value_hundred_uM = p.adjust(p_value_hundred_uM, method = "BH"),
    
    # Z-scores
    log2FC_ten_uM_avg = mean(log2FC_ten_uM, na.rm = TRUE),
    log2FC_ten_uM_sd  = sd(log2FC_ten_uM, na.rm = TRUE),
    Z_Score_ten_uM    = (log2FC_ten_uM - log2FC_ten_uM_avg) / log2FC_ten_uM_sd,
    
    log2FC_hundred_uM_avg = mean(log2FC_hundred_uM, na.rm = TRUE),
    log2FC_hundred_uM_sd  = sd(log2FC_hundred_uM, na.rm = TRUE),
    Z_Score_hundred_uM    = (log2FC_hundred_uM - log2FC_hundred_uM_avg) / log2FC_hundred_uM_sd
  )

SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name <- SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test %>% 
  left_join(
    uniprot_proteome,
    by = "Master_Protein_Accession"
  )

# filter for hits
SPROX_Hits_ten_uM <- SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name %>%
  filter(p_value_ten_uM < 0.05) %>% 
  filter(abs(Z_Score_ten_uM) > 2) %>% 
  dplyr::select(Gene_Name) %>% 
  distinct()

SPROX_Hits_hundred_uM <- SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name %>%
  filter(p_value_hundred_uM < 0.05) %>% 
  filter(abs(Z_Score_hundred_uM) > 2) %>% 
  dplyr::select(Gene_Name) %>% 
  distinct()

SPROX_Hits_overlapped <- inner_join(
  SPROX_Hits_ten_uM, 
  SPROX_Hits_hundred_uM) %>% 
  distinct()

#plot the result

# plot the volcano plots with Z-Score.

SPROX_Hits_overlapped_plot <- SPROX_Hits_overlapped %>% 
  left_join(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name) %>% 
  filter(p_value_ten_uM < 0.05) %>% 
  filter(abs(Z_Score_ten_uM) > 2) %>% 
  filter(p_value_hundred_uM < 0.05) %>% 
  filter(abs(Z_Score_hundred_uM) > 2) %>% 
  filter(!Gene_Name == "ALDH1A3") %>% 
  filter(!Gene_Name == "KRAS")

# SPROX_10uM
SPROX_10uM <- ggplot() +
  # Background
  geom_point(data = SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM), color = "grey", size = 0.3, alpha = 0.6) +
  # other proteins
  geom_point(data = SPROX_Hits_overlapped_plot, aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM), color = "#4B9CD3", size = 1, alpha = 0.6) +
  # KRAS
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "KRAS" & sequence == "VKDSEDVPMVLVGNK"), 
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM), color = "red", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "KRAS" & sequence == "VKDSEDVPMVLVGNK"),
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = 0.3, nudge_y = 0, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # ALDH1A3
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ALDH1A3" &  sequence == "GLFIKPTVFSEVTDNMR"),
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM), color = "blue", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ALDH1A3" & sequence == "GLFIKPTVFSEVTDNMR"),
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = -0.2, nudge_y = 0.2, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # ADK
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ADK" &  sequence == "VAQWMIQQPHK"),
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM), color = "blue", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ADK" & sequence == "VAQWMIQQPHK"),
    aes(x = Z_Score_ten_uM, y = neg_log10_p_value_ten_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = -0.2, nudge_y = 0.2, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # Thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "33", color = "black", linewidth = 0.25) +
  geom_vline(xintercept = c(-2, 2), linetype = "33", color = "black", linewidth = 0.25) +
  labs(
    x = expression("Z-Score"),y = expression(-log[10]~"(p-value)")) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2), labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    plot.margin = margin(6, 10, 6, 10),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12))

SPROX_10uM

ggsave("volcano_SPROX_10uM.png", SPROX_10uM,
       device = ragg::agg_png, width = 3.45, height = 3.2, units = "in",
       dpi = 600)


# SPROX_100uM
SPROX_100uM <- ggplot() +
  # Background
  geom_point(data = SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM), color = "grey", size = 0.3, alpha = 0.6) +
  # other proteins
  geom_point(data = SPROX_Hits_overlapped_plot, aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM), color = "#4B9CD3", size = 1, alpha = 0.6) +
  # KRAS
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "KRAS" & sequence == "VKDSEDVPMVLVGNK"), 
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM), color = "red", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "KRAS" & sequence == "VKDSEDVPMVLVGNK"),
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = 0.2, nudge_y = 0, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # ALDH1A3
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ALDH1A3" & sequence == "GLFIKPTVFSEVTDNMR"),
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM), color = "blue", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ALDH1A3" & sequence == "GLFIKPTVFSEVTDNMR"),
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = -0.2, nudge_y = 0.1, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # ADK
  geom_point(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ADK" & sequence == "VAQWMIQQPHK"),
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM), color = "blue", size = 2.5, alpha = 0.95, shape = 16) +
  geom_text_repel(
    data = subset(SPROX_Met_Enrich_fil_wt_Met_norm_Welch_t_test_Gene_Name, Gene_Name == "ADK" & sequence == "VAQWMIQQPHK"),
    aes(x = Z_Score_hundred_uM, y = neg_log10_p_value_hundred_uM, label = Gene_Name), size = 5, box.padding = 0.25, point.padding = 0.25, nudge_x = -0.2, nudge_y = 0.1, color = "black", segment.color = "black", segment.size = 0.25, max.overlaps = Inf) +
  # Thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "33", color = "black", linewidth = 0.25) +
  geom_vline(xintercept = c(-2, 2), linetype = "33", color = "black", linewidth = 0.25) +
  labs(
    x = expression("Z-Score"),y = expression(-log[10]~"(p-value)")) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 2), labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    plot.margin = margin(6, 10, 6, 10),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12))

SPROX_100uM

ggsave("volcano_SPROX_100uM.png", SPROX_100uM,
       device = ragg::agg_png, width = 3.45, height = 3.2, units = "in",
       dpi = 600)
```


