# This script contains helper functions for DE analysis

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
