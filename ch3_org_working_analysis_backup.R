## QC Analysis Script for Gene-Level Seurat Objects------
# Load essential packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(Matrix)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")
dir.create("./output_files/QC", recursive = TRUE, showWarnings = FALSE)

# Sample information
sample_ids <- c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C")
timepoints <- c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org")

# FLAMES output directory
FLAMES_OUTPUT_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"

## STEP 1: Create gene-level count matrices with gene symbols------

# Load isoform-gene dictionary (contains gene_id to gene_symbol mapping)
cat("=== Loading isoform-gene dictionary ===\n")
isoform_gene_dict <- read.csv("./output_files/ref_files/isoform_gene_dict.csv")

# Extract unique gene_id to gene_symbol mapping
gene_symbol_dict <- isoform_gene_dict %>%
  select(gene_id, gene_symbol) %>%
  distinct()

cat("Loaded", nrow(gene_symbol_dict), "unique gene mappings\n")

# Function to convert gene IDs to symbols in gene-level count matrices
convert_gene_counts_to_symbols <- function(gene_count_file, sample_id, gene_dict) {

  cat("\nProcessing:", sample_id, "\n")

  # Read gene count matrix
  cat("  Reading count matrix...\n")
  gene_counts <- read.csv(gene_count_file, row.names = 1, check.names = FALSE)

  # Convert NAs to zeros (FLAMES uses empty cells for zero counts)
  gene_counts[is.na(gene_counts)] <- 0
  cat("  Converted NA values to zeros\n")

  cat("  Original dimensions:", nrow(gene_counts), "genes x", ncol(gene_counts), "cells\n")

  # Create data frame with gene IDs
  gene_counts_df <- gene_counts %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id")

  # Remove version numbers from gene IDs (e.g., ENSG00000123.1 -> ENSG00000123)
  gene_counts_df$gene_id_clean <- str_replace(gene_counts_df$gene_id, "\\.\\d+$", "")

  # Also clean the dictionary gene IDs
  gene_dict_clean <- gene_dict %>%
    mutate(gene_id_clean = str_replace(gene_id, "\\.\\d+$", ""))

  # Merge with gene symbols
  cat("  Merging with gene symbols...\n")
  gene_counts_with_symbols <- gene_counts_df %>%
    left_join(gene_dict_clean, by = "gene_id_clean") %>%
    select(-gene_id_clean)

  # Check how many genes got symbols
  n_with_symbols <- sum(!is.na(gene_counts_with_symbols$gene_symbol))
  n_without_symbols <- sum(is.na(gene_counts_with_symbols$gene_symbol))
  cat("  Genes with symbols:", n_with_symbols, "\n")
  cat("  Genes without symbols:", n_without_symbols, "(will use gene_id)\n")

  # For genes without symbols, use the gene_id
  gene_counts_with_symbols <- gene_counts_with_symbols %>%
    mutate(gene_symbol = ifelse(is.na(gene_symbol), gene_id.x, gene_symbol))

  # Create combined gene_id_gene_symbol format for rownames (ensures uniqueness)
  gene_counts_with_symbols <- gene_counts_with_symbols %>%
    mutate(gene_name = paste0(gene_id.x, "_", gene_symbol))

  # Set combined names as rownames
  rownames(gene_counts_with_symbols) <- gene_counts_with_symbols$gene_name

  # Remove gene_id and gene_symbol columns, keeping only count data
  gene_counts_final <- gene_counts_with_symbols %>%
    select(-gene_id.x, -gene_id.y, -gene_symbol, -gene_name)

  cat("  Final dimensions:", nrow(gene_counts_final), "genes x", ncol(gene_counts_final), "cells\n")

  return(gene_counts_final)
}

# Process all samples and create gene symbol count matrices
cat("\n=== Converting gene IDs to symbols for all samples ===\n")

gene_symbol_count_matrices <- list()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  gene_count_file <- file.path(FLAMES_OUTPUT_DIR, paste0(sample_id, "_gene_count.csv"))

  gene_symbol_count_matrices[[sample_id]] <- convert_gene_counts_to_symbols(
    gene_count_file = gene_count_file,
    sample_id = sample_id,
    gene_dict = gene_symbol_dict
  )
}

cat("\n=== All count matrices converted successfully ===\n")

## STEP 2: Create Seurat objects from gene symbol count matrices------

# Function to create Seurat objects (adapted from ch3_qc_flames_v220_p2.R)
create_seurat_from_counts <- function(count_matrix_df, sample_id,
                                      assay_type = "RNA",
                                      MIN.CELL = 3,
                                      MIN.FEATURES = 20) {

  cat("\n=== Creating Seurat object for:", sample_id, "===\n")

  # Convert to matrix
  count_matrix <- as.matrix(count_matrix_df)

  # Check for and handle duplicate rownames
  if(any(duplicated(rownames(count_matrix)))) {
    cat("Warning: Found", sum(duplicated(rownames(count_matrix))), "duplicate gene names in", sample_id, "\n")
    cat("Making rownames unique...\n")
    rownames(count_matrix) <- make.unique(rownames(count_matrix))
  }

  # PRE-FILTER before converting to sparse matrix
  cat("Filtering genes: keeping genes detected in at least", MIN.CELL, "cells\n")
  genes_keep <- rowSums(count_matrix > 0, na.rm = TRUE) >= MIN.CELL
  count_matrix <- count_matrix[genes_keep, ]

  cat("Filtering cells: keeping cells with at least", MIN.FEATURES, "features\n")
  cells_keep <- colSums(count_matrix > 0, na.rm = TRUE) >= MIN.FEATURES
  count_matrix <- count_matrix[, cells_keep]

  cat("Filtered dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "cells\n")

  # Convert to sparse matrix
  count_sparse <- Matrix::Matrix(count_matrix, sparse = TRUE)

  # Create Seurat object WITHOUT additional filtering (already filtered)
  seu_obj <- CreateSeuratObject(
    counts = count_sparse,
    project = sample_id,
    assay = assay_type
  )

  # Add sample metadata
  seu_obj$sample_id <- sample_id

  cat("Seurat object created successfully!\n")

  return(seu_obj)
}

# Create and save Seurat objects for all samples
cat("\n=== Creating Seurat objects for all samples ===\n")

seurat_objects <- list()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]

  # Create Seurat object
  seurat_objects[[sample_id]] <- create_seurat_from_counts(
    count_matrix_df = gene_symbol_count_matrices[[sample_id]],
    sample_id = sample_id,
    assay_type = "RNA",
    MIN.CELL = 3,
    MIN.FEATURES = 20
  )

  # Add timepoint metadata
  seurat_objects[[sample_id]]$timepoint <- timepoints[i]

  # Save immediately to disk
  output_file <- paste0("./output_files/seu_objects/", sample_id, "_gene_seurat_with_symbols.rds")
  saveRDS(seurat_objects[[sample_id]], file = output_file)
  cat("Saved:", output_file, "\n")

  # Print summary
  print(seurat_objects[[sample_id]])
}

cat("\n========================================\n")
cat("=== All Seurat objects created and saved successfully! ===\n")
cat("========================================\n")

## STEP 3: QC Analysis on Gene-Level Seurat Objects------

cat("\n=== Starting QC Analysis ===\n")

# Calculate mitochondrial and ribosomal percentages for all samples
for(sample_id in sample_ids) {
  cat("\nCalculating QC metrics for:", sample_id, "\n")

  # Mitochondrial genes (pattern matches MT- anywhere in feature name)
  # Feature names are formatted as: ENSG00000198888.2-MT-ND1
  seurat_objects[[sample_id]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_objects[[sample_id]],
    pattern = "-MT-"
  )

  # Ribosomal genes (pattern matches RP[SL] anywhere in feature name)
  seurat_objects[[sample_id]][["percent.ribo"]] <- PercentageFeatureSet(
    seurat_objects[[sample_id]],
    pattern = "-RP[SL]"
  )

  cat("  Median percent.mt:", median(seurat_objects[[sample_id]]$percent.mt), "\n")
  cat("  Median percent.ribo:", median(seurat_objects[[sample_id]]$percent.ribo), "\n")
}

# Create QC violin plots for all samples
cat("\n=== Creating QC violin plots ===\n")

for(sample_id in sample_ids) {
  cat("Plotting:", sample_id, "\n")

  # Use VlnPlot with default labels
  # Default labels are clear: nFeature_RNA, nCount_RNA, percent.mt, percent.ribo
  p <- VlnPlot(
    seurat_objects[[sample_id]],
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    ncol = 4,
    pt.size = 0.1
  )

  # Use png() device instead of ggsave() to avoid S4SXP error with Seurat v5
  output_file <- paste0("./output_files/QC/", sample_id, "_QC_violin.png")
  png(output_file, width = 16*300, height = 6*300, res = 300)
  print(p)
  dev.off()
  cat("  Saved:", output_file, "\n")
}

# Create scatter plots showing relationships between QC metrics
cat("\n=== Creating QC scatter plots ===\n")

for(sample_id in sample_ids) {
  cat("Plotting:", sample_id, "\n")

  # nCount vs nFeature - with intuitive labels
  p1 <- FeatureScatter(
    seurat_objects[[sample_id]],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ) +
    xlab("Total UMI Counts (Library Size)") +
    ylab("Genes Detected per Cell") +
    ggtitle("Library Complexity")

  # nCount vs percent.mt - with intuitive labels
  p2 <- FeatureScatter(
    seurat_objects[[sample_id]],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  ) +
    xlab("Total UMI Counts (Library Size)") +
    ylab("Mitochondrial Content (%)") +
    ggtitle("Cell Viability")

  # Combine plots
  p <- p1 + p2 + plot_annotation(title = paste("QC Metrics -", sample_id))

  # Save plot
  output_file <- paste0("./output_files/QC/", sample_id, "_QC_scatter.png")
  ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
  cat("  Saved:", output_file, "\n")
}

# Generate summary statistics table
cat("\n=== Generating summary statistics ===\n")

qc_summary <- data.frame()

for(i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  seu <- seurat_objects[[sample_id]]

  qc_summary <- rbind(qc_summary, data.frame(
    sample_id = sample_id,
    timepoint = timepoints[i],
    n_cells = ncol(seu),
    n_genes = nrow(seu),
    median_nCount = median(seu$nCount_RNA),
    median_nFeature = median(seu$nFeature_RNA),
    median_percent_mt = median(seu$percent.mt),
    median_percent_ribo = median(seu$percent.ribo),
    mean_nCount = mean(seu$nCount_RNA),
    mean_nFeature = mean(seu$nFeature_RNA),
    mean_percent_mt = mean(seu$percent.mt),
    mean_percent_ribo = mean(seu$percent.ribo)
  ))
}


# Print summary
cat("\n=== QC Summary Statistics ===\n")
print(qc_summary)

# Save summary table
write.csv(qc_summary, file = "./output_files/QC/qc_summary_statistics.csv", row.names = FALSE)
cat("\nSaved: ./output_files/QC/qc_summary_statistics.csv\n")

cat("\n========================================\n")
cat("=== QC Analysis Complete! ===\n")
cat("========================================\n")
cat("\nSeurat objects with QC metrics are stored in 'seurat_objects' list\n")
cat("QC plots saved in: ./output_files/QC/\n")

