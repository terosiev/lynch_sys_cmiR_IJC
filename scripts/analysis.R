# Introduction ----
# This script performs differential expression (DE) analysis of circulating microRNA (c-miR) expression
# in Lynch syndrome (LS) patients, healthy controls, and sporadic rectal cancer (SRME) patients.
# The goal is to identify differences in c-miR expression across groups and conditions.

# Script written by Tia-Marje Korhonen (Uni. of Jyväskylä) and Tero Sievänen (Uni. East. Finland)
# Last edited: 30.1.2025

# Load required packages ----
library(edgeR)       # For differential expression analysis
library(DESeq2)      # For DE analysis
library(tidyverse)   # For data manipulation and visualization
library(DT)          # For interactive tables
library(plotly)      # For interactive plots

# Data import ----
counts <- read.csv("rawCounts.tsv", header = TRUE, sep = "\t")  # Raw counts data
targets <- read.csv("phenodata_age.txt", sep = "\t", header = TRUE)  # Study design

# Add column names to counts
colnames(counts) <- targets$Type

# Filtering of raw c-miR counts ----
myDGEList <- DGEList(counts = counts)  # Create DGEList object
cpm <- cpm(myDGEList)  # Convert to counts per million (CPM)
keepers <- rowSums(cpm > 1) >= 108  # Keep c-miRs with >1 CPM in at least 70% of samples
myDGEList.filtered <- myDGEList[keepers, ]  # Filter DGEList
FilteredCounts <- counts[rownames(counts) %in% rownames(myDGEList.filtered$counts), ]  # Filter raw counts

# Save filtered counts to a file
write.table(FilteredCounts, "FilteredCounts.txt", sep = "\t")

# Function for DE analysis ----
perform_de_analysis <- function(filtered_counts, targets, condition_col, batch_col, alpha = 0.05) {
  # Subset counts and targets based on condition
  select <- targets[[condition_col]] %in% c("CTRL", "SRME", "LS")  # Adjust as needed
  Counts <- filtered_counts[, select]
  Targets <- targets[select, ]
  
  # Set up design matrix
  condition <- as.character(Targets[[condition_col]])
  batch <- as.character(Targets[[batch_col]])
  design <- data.frame(condition = as.factor(condition), batch = batch)
  rownames(design) <- colnames(Counts)
  
  # Perform DE analysis
  dds <- DESeqDataSetFromMatrix(countData = Counts, colData = design, design = ~ batch + condition)
  dds <- DESeq(dds)
  results(dds, alpha = alpha)
}

# Function to display results ----
display_results <- function(res, top_n = 20) {
  resOrdered <- res[order(res$padj), ]
  subset <- head(resOrdered, top_n) %>%
    as_tibble(rownames = "miR")
  
  # Interactive table
  datatable(subset, 
            extensions = c('KeyTable', "FixedHeader"), 
            caption = 'Differentially expressed c-miRs',
            options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 20, lengthMenu = c("5", "10", "15", "50"))) %>%
    formatRound(columns = c(2:7), digits = 3)
}

# Function to create a volcano plot ----
create_volcano_plot <- function(res, title = "Volcano plot", subtitle = "") {
  res.df <- as_tibble(res, rownames = "miR")
  ggplot(res.df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "longdash", colour = "grey", size = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "longdash", colour = c("#2C467A", "#BE684D"), size = 1) +
    labs(title = title, subtitle = subtitle, x = "Log2FC", y = "-log10(Padj)") +
    theme_bw()
}

# Perform DE analysis for different comparisons using filtered counts ----

# 1. Sporadic rectal cancer patients vs non-LS controls
res_srme_vs_ctrl <- perform_de_analysis(FilteredCounts, targets, "Type", "NGS")
display_results(res_srme_vs_ctrl)
create_volcano_plot(res_srme_vs_ctrl, title = "SRME vs CTRL")

# 2. Healthy LS vs diseased LS
res_ls_healthy_vs_diseased <- perform_de_analysis(FilteredCounts, targets, "Healthy_now", "NGS")
display_results(res_ls_healthy_vs_diseased)
create_volcano_plot(res_ls_healthy_vs_diseased, title = "Healthy LS vs Diseased LS")

# 3. Healthy LS vs non-LS controls
res_ls_healthy_vs_ctrl <- perform_de_analysis(FilteredCounts, targets, "Type", "NGS")
display_results(res_ls_healthy_vs_ctrl)
create_volcano_plot(res_ls_healthy_vs_ctrl, title = "Healthy LS vs CTRL")

# 4. Healthy LS vs SRME
res_ls_healthy_vs_srme <- perform_de_analysis(FilteredCounts, targets, "Type", "NGS")
display_results(res_ls_healthy_vs_srme)
create_volcano_plot(res_ls_healthy_vs_srme, title = "Healthy LS vs SRME")

# 5. Diseased LS vs SRME
res_ls_diseased_vs_srme <- perform_de_analysis(FilteredCounts, targets, "Type", "NGS")
display_results(res_ls_diseased_vs_srme)
create_volcano_plot(res_ls_diseased_vs_srme, title = "Diseased LS vs SRME")
