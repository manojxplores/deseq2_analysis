# Mini-Project: Differential Gene Expression Analysis - GENERALIZED SCRIPT
# Comparison: Control (0hr) vs Dehydration (Specified Timepoint) in Soybean (GSE57252)

# --- 1. Load Libraries ---
# Ensure required packages are installed:
# install.packages(c("DESeq2", "ggplot2", "pheatmap", "EnhancedVolcano", "ggrepel"))
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)

# --- 2. Configuration ---

# !!! USER INPUT: SET THE DESIRED TIMEPOINT FOR ANALYSIS (1, 6, or 12) !!!
analysis_timepoint_hr <- 6  # <--- CHANGE THIS VALUE (e.g., 1, 6, or 12)

# -------------------------------------------------------------------------

# Validate timepoint input
valid_timepoints <- c(1, 6, 12)
if (!(analysis_timepoint_hr %in% valid_timepoints)) {
  stop("Invalid analysis_timepoint_hr. Please choose 1, 6, or 12.")
}
cat("--- Starting analysis for", analysis_timepoint_hr, "hr timepoint ---\n\n")

# Define input file path (adjust if necessary)
count_file <- "datasets\\GSE57252_RawRead_counts.txt" # Assumes file is in the working directory or provide full/relative path

# Define output directory dynamically based on timepoint
output_dir <- paste0("DGE_Analysis_Control_vs_Dehydration_", analysis_timepoint_hr, "hr")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Created output directory:", output_dir, "\n")
} else {
  cat("Output directory already exists:", output_dir, "\n")
}

# Define FDR threshold
fdr_cutoff <- 0.1

# Define Log2 Fold Change threshold for Volcano plot highlighting
lfc_cutoff <- 1.0

# --- 3. Load and Prepare Data ---

# Check if count file exists
if (!file.exists(count_file)) {
  stop("Count file not found: ", count_file)
}

# Read raw count data
cat("Reading raw count data from:", count_file, "\n")
raw_counts <- read.table(count_file, header = TRUE, row.names = 1, check.names = FALSE)
cat("Dimensions of the full count data matrix:", dim(raw_counts), "\n")

# Define samples for comparison dynamically
control_samples <- c("Co0hrR1", "Co0hrR2", "Co0hrR3")
# Construct dehydration sample names based on the chosen timepoint
dehydration_samples <- paste0("De", analysis_timepoint_hr, "hrR", 1:3)
selected_samples <- c(control_samples, dehydration_samples)

cat("\nSelected samples for comparison:\n")
print(selected_samples)

# Check if all selected samples are in the count data columns
if (!all(selected_samples %in% colnames(raw_counts))) {
    missing_samples <- selected_samples[!selected_samples %in% colnames(raw_counts)]
    stop("Some selected samples are not found in the columns of the count data: ", paste(missing_samples, collapse=", "))
}

# Subset the count data for selected samples
counts_subset <- raw_counts[, selected_samples]
cat("\nDimensions of the subsetted count data matrix:", dim(counts_subset), "\n")

# Create sample metadata (colData)
# The structure is always 3 control vs 3 dehydration for this comparison
condition <- factor(rep(c("Control", "Dehydration"), each = 3))
col_data <- data.frame(row.names = selected_samples, condition = condition)
cat("\nSample metadata (colData):\n")
print(col_data)

# --- 4. Create DESeqDataSet Object ---
cat("\nCreating DESeqDataSet object...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData = col_data,
  design = ~ condition
)

# --- 5. Pre-filtering ---
# Keep genes with at least 10 reads total across the selected samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Number of genes after filtering (>= 10 reads total):", nrow(dds), "\n")

# --- 6. Set Reference Level ---
# Ensure 'Control' is the reference level for fold change calculations
dds$condition <- relevel(dds$condition, ref = "Control")
cat("Reference level for comparison set to 'Control'.\n")

# --- 7. Run DESeq2 Analysis ---
cat("\nRunning DESeq2 analysis for", analysis_timepoint_hr, "hr (this may take a few minutes)...\n")
dds <- DESeq(dds)
cat("DESeq2 analysis complete.\n")

# --- 8. Get Results ---
cat("\nExtracting results (FDR <", fdr_cutoff, ")...\n")
# Contrast is always 'Dehydration' vs 'Control' based on the metadata factor levels
res <- results(dds, contrast = c("condition", "Dehydration", "Control"), alpha = fdr_cutoff)

cat("\nSummary of results for", analysis_timepoint_hr, "hr vs Control:\n")
summary(res)

# Sort results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# --- 9. Extract and Save Significant Genes ---
cat("\nIdentifying significant differentially expressed genes (DEGs)...\n")
sig_genes <- subset(res_ordered, padj < fdr_cutoff)
cat("Number of significant DEGs (FDR <", fdr_cutoff, "):", nrow(sig_genes), "\n")

# Prepare data frame for saving
if (nrow(sig_genes) > 0) {
    sig_genes_df <- as.data.frame(sig_genes)
    sig_genes_df$gene_id <- rownames(sig_genes_df)
    # Reorder columns for clarity
    sig_genes_df <- sig_genes_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
} else {
    # Create empty dataframe with correct columns if no significant genes found
    sig_genes_df <- data.frame(gene_id=character(), baseMean=numeric(), log2FoldChange=numeric(),
                               lfcSE=numeric(), stat=numeric(), pvalue=numeric(), padj=numeric())
    cat("No significant DEGs found at FDR <", fdr_cutoff, ".\n")
}

# Save the list of significant genes to a CSV file (dynamic filename)
output_csv <- file.path(output_dir, paste0("DEGs_Control_vs_Dehydration_", analysis_timepoint_hr, "hr_FDR", fdr_cutoff, ".csv"))
write.csv(sig_genes_df, output_csv, row.names = FALSE)
cat("Saved significant DEGs list to:", output_csv, "\n")

# --- 10. Generate Visualizations ---
cat("\nGenerating plots for", analysis_timepoint_hr, "hr...\n")

# Define file paths for plots dynamically
ma_plot_file <- file.path(output_dir, paste0("MA_plot_Control_vs_Dehydration_", analysis_timepoint_hr, "hr.png"))
pca_plot_file <- file.path(output_dir, paste0("PCA_plot_Control_vs_Dehydration_", analysis_timepoint_hr, "hr.png"))
volcano_plot_file <- file.path(output_dir, paste0("Volcano_plot_Control_vs_Dehydration_", analysis_timepoint_hr, "hr.png"))
heatmap_top30_file <- file.path(output_dir, paste0("Heatmap_Top30_DEGs_Control_vs_Dehydration_", analysis_timepoint_hr, "hr.png"))
count_plot_dir <- file.path(output_dir, paste0("Top_DEG_Count_Plots_", analysis_timepoint_hr, "hr")) # Subdirectory for count plots
if (!dir.exists(count_plot_dir)) {
  dir.create(count_plot_dir)
}

# 10.1 MA Plot
png(ma_plot_file, width = 800, height = 600, res = 100)
plot_title_ma <- paste("MA Plot: Control vs Dehydration", analysis_timepoint_hr, "hr")
plotMA(res, ylim = c(-5, 5), main = plot_title_ma, alpha = fdr_cutoff)
abline(h = c(-lfc_cutoff, lfc_cutoff), col = "dodgerblue", lty = 2) # Add LFC threshold lines
dev.off()
cat("Saved MA plot to:", ma_plot_file, "\n")

# 10.2 PCA Plot
# Use variance stabilizing transformation for PCA visualization
vsd <- vst(dds, blind = FALSE) # blind=FALSE uses the design for transformation

png(pca_plot_file, width = 800, height = 600, res = 100)
# Use ntop=500 for consistency, DESeq2 default
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pca_data, "percentVar"))
plot_title_pca <- paste("PCA Plot: Control vs Dehydration", analysis_timepoint_hr, "hr")
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(plot_title_pca) +
  scale_color_manual(values = c("Control" = "coral", "Dehydration" = "turquoise")) + # Consistent colors
  theme_bw()
dev.off()
cat("Saved PCA plot to:", pca_plot_file, "\n")

# 10.3 Volcano Plot
if (!is.null(res)) {
    png(volcano_plot_file, width = 1000, height = 800, res = 100)
    plot_title_volcano <- paste('Control vs Dehydration', analysis_timepoint_hr, 'hr')
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = plot_title_volcano,
                    pCutoff = fdr_cutoff,
                    FCcutoff = lfc_cutoff,
                    pointSize = 2.0,
                    labSize = 3.0,
                    max.overlaps = 20, # Increase if needed, but might get cluttered
                    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                    colAlpha = 0.6,
                    legendPosition = 'right',
                    legendLabSize = 12,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5) +
      labs(subtitle = paste(nrow(sig_genes), "DEGs found (FDR <", fdr_cutoff, ")"))
    dev.off()
    cat("Saved Volcano plot to:", volcano_plot_file, "\n")
} else {
    cat("Skipping Volcano plot generation as DESeq results are NULL.\n")
}


# 10.4 Heatmap of Top 30 Significant Genes
if (nrow(sig_genes) >= 1) {
  num_top_genes <- min(30, nrow(sig_genes))
  top_sig_genes <- head(rownames(sig_genes), num_top_genes)
  mat <- assay(vsd)[top_sig_genes, , drop = FALSE]

  if (nrow(mat) > 0 && ncol(mat) > 0) {
      mat <- mat - rowMeans(mat) # Center genes
      anno_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])
      anno_colors <- list(condition = c(Control = "coral", Dehydration = "turquoise"))
      plot_title_heatmap <- paste("Heatmap of Top", num_top_genes, "DEGs (Control vs Dehydration", analysis_timepoint_hr, "hr, FDR <", fdr_cutoff, ")")

      png(heatmap_top30_file, width = 800, height = 1000, res = 100)
      pheatmap(mat,
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               main = plot_title_heatmap,
               fontsize_row = 8,
               fontsize_col = 10,
               cluster_rows = TRUE,
               cluster_cols = TRUE, # Cluster cols is usually desired
               show_rownames = TRUE,
               scale = "none", # Already centered
               border_color = NA
               )
      dev.off()
      cat("Saved Heatmap of top DEGs to:", heatmap_top30_file, "\n")
   } else {
        cat("Skipping heatmap: Matrix for top significant genes is empty or invalid.\n")
   }
} else {
  cat("Skipping heatmap generation as no significant DEGs were found.\n")
}


# 10.5 Count plots for Top 6 Significant Genes
if (nrow(sig_genes) > 0) {
  cat("Generating count plots for top 6 DEGs...\n")
  num_count_plots <- min(6, nrow(sig_genes))
  top_genes_for_counts <- head(rownames(sig_genes), num_count_plots)
  normalized_counts <- counts(dds, normalized = TRUE)

  for (gene in top_genes_for_counts) {
    d <- data.frame(
      Sample = colnames(normalized_counts),
      Count = normalized_counts[gene, ],
      Condition = colData(dds)$condition
    )
    p <- ggplot(d, aes(x = Sample, y = Count, fill = Condition)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Control" = "coral", "Dehydration" = "turquoise")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste("Normalized counts for", gene)) +
      xlab("Sample") +
      ylab("Normalized Count")

    # Save plot with dynamic filename
    count_plot_filename <- file.path(count_plot_dir, paste0("Count_plot_", analysis_timepoint_hr, "hr_", gene, ".png"))
    ggsave(count_plot_filename, plot = p, width = 8, height = 6)
  }
  cat("Saved count plots to:", count_plot_dir, "\n")
} else {
    cat("Skipping count plot generation as no significant DEGs were found.\n")
}

# --- 11. Session Information ---
cat("\nAnalysis for", analysis_timepoint_hr, "hr complete. Session Information:\n")
sessionInfo()

# --- End of Script ---