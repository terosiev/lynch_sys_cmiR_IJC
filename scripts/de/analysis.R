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

# Get helpers
source(utils.R)

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
