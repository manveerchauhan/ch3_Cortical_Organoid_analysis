#!/usr/bin/env Rscript
# =============================================================================
# FLAMES Cortical Organoid Analysis Pipeline - Main Entry Point
# =============================================================================
# Author: Manveer Chauhan and Sefi Prawer
# Description: Main entry point for the modular FLAMES analysis pipeline
#              Runs complete analysis from data loading to final reports
# =============================================================================

# Clear environment and set up
rm(list = ls())
cat("\n")
cat("================================================================================\n")
cat("      FLAMES Cortical Organoid Differentiation Analysis Pipeline v1.0.0       \n")
cat("================================================================================\n")
cat("\n")

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# Check if running in correct directory
if (!file.exists("analysis_pipeline/functions/config.R")) {
  stop("Please run this script from the main project directory containing analysis_pipeline/")
}

# Load configuration first
cat("[INIT] Loading configuration...\n")
source("analysis_pipeline/functions/config.R")

# Create output directories
create_output_dirs()

# Initialize logging
start_time <- Sys.time()
log_message(sprintf("Starting pipeline at %s", start_time), "INFO")
log_message(sprintf("Pipeline version: %s", PIPELINE_VERSION), "INFO")

# Load all required libraries
cat("[INIT] Loading required libraries...\n")
required_packages <- c(
  "Seurat", "Matrix", "dplyr", "ggplot2", "grid", "gridExtra",
  "cowplot", "qs", "data.table", "viridis", "RColorBrewer",
  "ggrepel", "scales", "stringr", "patchwork"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(sprintf("Missing required packages: %s\nPlease install with: install.packages(c(%s))",
               paste(missing_packages, collapse = ", "),
               paste0("'", paste(missing_packages, collapse = "', '"), "'")))
}

# Load libraries quietly
suppressPackageStartupMessages({
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }
})

# Load all functional modules
cat("[INIT] Loading functional modules...\n")
source("analysis_pipeline/functions/flames_io.R")
source("analysis_pipeline/functions/plotting_utils.R")
source("analysis_pipeline/functions/id_mapping_utils.R")
source("analysis_pipeline/functions/seurat_qc_modules.R")

# =============================================================================
# PIPELINE CONTROL FUNCTIONS
# =============================================================================

#' Run pipeline step with error handling and progress tracking
#'
#' @param step_name Name of the analysis step
#' @param step_function Function to execute
#' @param dependencies List of required previous steps
#' @param resume_on_error Should pipeline continue if this step fails?
run_pipeline_step <- function(step_name, step_function, dependencies = NULL, resume_on_error = FALSE) {
  log_message(sprintf("Starting Step: %s", step_name), "INFO")
  step_start_time <- Sys.time()

  # Check dependencies
  if (!is.null(dependencies)) {
    for (dep in dependencies) {
      dep_file <- file.path(OUTPUT_BASE_DIR, "data", "metadata", paste0(dep, "_completed.rds"))
      if (!file.exists(dep_file)) {
        stop(sprintf("Dependency not met: %s must be completed before %s", dep, step_name))
      }
    }
  }

  # Progress indicator
  if (ENABLE_PROGRESS_BARS) {
    cat(sprintf("\n=== %s ===\n", toupper(step_name)))
  }

  tryCatch({
    # Execute step function
    result <- step_function()

    # Mark step as completed
    step_end_time <- Sys.time()
    step_duration <- as.numeric(difftime(step_end_time, step_start_time, units = "mins"))

    completion_info <- list(
      step_name = step_name,
      start_time = step_start_time,
      end_time = step_end_time,
      duration_minutes = step_duration,
      status = "SUCCESS"
    )

    # Save completion marker
    saveRDS(completion_info, file.path(OUTPUT_BASE_DIR, "data", "metadata", paste0(gsub(" ", "_", tolower(step_name)), "_completed.rds")))

    log_message(sprintf("Completed Step: %s (%.2f minutes)", step_name, step_duration), "INFO")

    return(result)

  }, error = function(e) {
    step_end_time <- Sys.time()
    step_duration <- as.numeric(difftime(step_end_time, step_start_time, units = "mins"))

    error_info <- list(
      step_name = step_name,
      start_time = step_start_time,
      end_time = step_end_time,
      duration_minutes = step_duration,
      status = "FAILED",
      error_message = e$message
    )

    # Save error info
    saveRDS(error_info, file.path(OUTPUT_BASE_DIR, "data", "metadata", paste0(gsub(" ", "_", tolower(step_name)), "_error.rds")))

    log_message(sprintf("FAILED Step: %s - %s", step_name, e$message), "ERROR")

    if (resume_on_error) {
      log_message(sprintf("Continuing pipeline despite error in %s", step_name), "WARN")
      return(NULL)
    } else {
      stop(sprintf("Pipeline stopped due to error in %s: %s", step_name, e$message))
    }
  })
}

# =============================================================================
# ANALYSIS STEP FUNCTIONS
# =============================================================================

#' Step 1: Data Loading and Initial QC
step_data_loading <- function() {
  log_message("Loading FLAMES data and performing initial QC...", "INFO")

  # Load all FLAMES data
  flames_data <- load_all_flames_data(FLAMES_OUTPUT_DIR, save_raw = TRUE)

  # Extract components
  sample_metadata <- flames_data$sample_metadata
  count_matrices <- flames_data$count_matrices
  gene_count_matrices <- flames_data$gene_count_matrices
  summary_metrics <- flames_data$summary_metrics

  # Add gene-level statistics to metadata
  for (i in seq_len(nrow(sample_metadata))) {
    sample_id <- sample_metadata$sample_id[i]
    if (sample_id %in% names(gene_count_matrices)) {
      gene_counts <- gene_count_matrices[[sample_id]]
      sample_metadata$n_genes[i] <- nrow(gene_counts)
      sample_metadata$total_gene_counts[i] <- sum(gene_counts)
    }
  }

  # Save updated metadata
  write.csv(sample_metadata, file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.csv"), row.names = FALSE)
  saveRDS(sample_metadata, file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.rds"))

  # Create initial QC plots
  qc_plots <- list()

  # Plot 1: Sample overview (cells per sample)
  p1 <- ggplot(sample_metadata, aes(x = sample_id, y = n_cells, fill = timepoint)) +
    geom_bar(stat = "identity") +
    labs(title = "Detected Cells per Sample (Unfiltered)",
         x = "Sample", y = "Number of Cells") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Total cells: %d", sum(sample_metadata$n_cells)),
             hjust = 1.1, vjust = 1.5,
             size = 3, fontface = "bold", color = "darkblue")

  # Plot 2: Isoforms (transcripts) vs Cells
  p2 <- ggplot(sample_metadata, aes(x = n_cells, y = n_features, color = timepoint, label = sample_id)) +
    geom_point(size = 3) +
    geom_text_repel() +
    labs(title = "Isoforms (Transcripts) vs Cells by Sample (Unfiltered)",
         x = "Number of Cells",
         y = "Number of Isoforms Detected") +
    theme_publication() +
    scale_color_manual(values = TIMEPOINT_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Avg isoforms: %.0f", mean(sample_metadata$n_features)),
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold", color = "darkblue")

  # Plot 3: Genes vs Cells
  p3 <- ggplot(sample_metadata, aes(x = n_cells, y = n_genes, color = timepoint, label = sample_id)) +
    geom_point(size = 3) +
    geom_text_repel() +
    labs(title = "Genes vs Cells by Sample (Unfiltered)",
         x = "Number of Cells",
         y = "Number of Genes Detected") +
    theme_publication() +
    scale_color_manual(values = TIMEPOINT_COLORS) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Avg genes: %.0f", mean(sample_metadata$n_genes)),
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold", color = "darkblue")

  # Plot 4: Isoforms per gene ratio
  sample_metadata$isoforms_per_gene <- sample_metadata$n_features / sample_metadata$n_genes
  p4 <- ggplot(sample_metadata, aes(x = timepoint, y = isoforms_per_gene, fill = timepoint)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 3) +
    labs(title = "Isoforms per Gene Ratio (Unfiltered)",
         x = "Timepoint",
         y = "Isoforms / Genes") +
    theme_publication() +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    theme(legend.position = "none") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Overall avg: %.2f", mean(sample_metadata$isoforms_per_gene)),
             hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold", color = "darkblue")

  # Plot 5: Total isoform counts by timepoint
  p5 <- ggplot(sample_metadata, aes(x = timepoint, y = total_counts, fill = timepoint)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 3) +
    labs(title = "Total Isoform Counts by Timepoint (Unfiltered)",
         x = "Timepoint",
         y = "Total Isoform Counts") +
    theme_publication() +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    theme(legend.position = "none")

  # Plot 6: Total gene counts by timepoint
  p6 <- ggplot(sample_metadata, aes(x = timepoint, y = total_gene_counts, fill = timepoint)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 3) +
    labs(title = "Total Gene Counts by Timepoint (Unfiltered)",
         x = "Timepoint",
         y = "Total Gene Counts") +
    theme_publication() +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    theme(legend.position = "none")

  qc_plots <- list(p1, p2, p3, p4, p5, p6)

  # Save individual plots
  save_figure(p1, "01_cells_per_sample", "single", "both", "qc")
  save_figure(p2, "02_isoforms_vs_cells", "single", "both", "qc")
  save_figure(p3, "03_genes_vs_cells", "single", "both", "qc")
  save_figure(p4, "04_isoforms_per_gene_ratio", "single", "both", "qc")
  save_figure(p5, "05_total_isoform_counts", "single", "both", "qc")
  save_figure(p6, "06_total_gene_counts", "single", "both", "qc")

  # Create PDF report
  create_multipanel_figure(
    qc_plots,
    "01_Data_Loading_QC",
    ncol = 2,
    title = "Data Loading and Initial Quality Control - Isoform and Gene Level",
    subdir = "../reports"
  )

  log_message("Data loading and initial QC completed", "INFO")
  log_message(sprintf("Average isoforms per gene: %.2f", mean(sample_metadata$isoforms_per_gene)), "INFO")

  return(list(
    sample_metadata = sample_metadata,
    count_matrices = count_matrices,
    gene_count_matrices = gene_count_matrices,
    summary_metrics = summary_metrics
  ))
}

#' Step 2: Save Per-Cell Metadata (No QC Plots - Done in Step 3)
step_save_metadata <- function() {
  log_message("Saving per-cell metadata for reference", "INFO")

  # Load data from previous step
  isoform_count_matrices <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "raw_counts", "transcript_count_matrices.rds"))
  gene_count_matrices <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "raw_counts", "gene_count_matrices.rds"))
  sample_metadata <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.rds"))

  detailed_qc_isoform <- list()
  detailed_qc_gene <- list()

  # Calculate basic per-cell metrics for each sample
  for (sample_id in sample_metadata$sample_id) {
    log_message(sprintf("Processing metadata for sample: %s", sample_id), "DEBUG")

    # ========== ISOFORM-LEVEL BASIC METRICS ==========
    isoform_counts <- isoform_count_matrices[[sample_id]]

    # Calculate per-cell isoform metrics (NO MT% - not relevant for isoform assay)
    cells_total_isoform_counts <- Matrix::colSums(isoform_counts)
    cells_n_isoforms <- Matrix::colSums(isoform_counts > 0)

    # Store basic isoform metadata
    qc_data_isoform <- data.frame(
      sample_id = sample_id,
      cell_barcode = colnames(isoform_counts),
      total_isoform_counts = cells_total_isoform_counts,
      n_isoforms = cells_n_isoforms,
      timepoint = sample_metadata$timepoint[sample_metadata$sample_id == sample_id],
      stringsAsFactors = FALSE
    )
    detailed_qc_isoform[[sample_id]] <- qc_data_isoform

    # ========== GENE-LEVEL BASIC METRICS ==========
    gene_counts <- gene_count_matrices[[sample_id]]

    # Calculate per-cell gene metrics (NO MT% yet - Seurat will handle this in Step 3)
    cells_total_gene_counts <- Matrix::colSums(gene_counts)
    cells_n_genes <- Matrix::colSums(gene_counts > 0)

    # Store basic gene metadata
    qc_data_gene <- data.frame(
      sample_id = sample_id,
      cell_barcode = colnames(gene_counts),
      total_gene_counts = cells_total_gene_counts,
      n_genes = cells_n_genes,
      timepoint = sample_metadata$timepoint[sample_metadata$sample_id == sample_id],
      stringsAsFactors = FALSE
    )
    detailed_qc_gene[[sample_id]] <- qc_data_gene
  }

  # Combine all metadata
  all_qc_isoform <- do.call(rbind, detailed_qc_isoform)
  all_qc_gene <- do.call(rbind, detailed_qc_gene)

  # Save metadata CSVs for reference (no plots - Seurat will handle QC visualization)
  write.csv(all_qc_isoform, file.path(OUTPUT_BASE_DIR, "data", "metadata", "per_cell_isoform_metadata.csv"), row.names = FALSE)
  write.csv(all_qc_gene, file.path(OUTPUT_BASE_DIR, "data", "metadata", "per_cell_gene_metadata.csv"), row.names = FALSE)

  saveRDS(all_qc_isoform, file.path(OUTPUT_BASE_DIR, "data", "metadata", "per_cell_isoform_metadata.rds"))
  saveRDS(all_qc_gene, file.path(OUTPUT_BASE_DIR, "data", "metadata", "per_cell_gene_metadata.rds"))

  # Log summary statistics
  log_message("Per-cell metadata saved", "INFO")
  log_message(sprintf("Isoform assay - Total cells: %d, Median isoforms/cell: %.0f, Median counts/cell: %.0f",
                      nrow(all_qc_isoform),
                      median(all_qc_isoform$n_isoforms),
                      median(all_qc_isoform$total_isoform_counts)), "INFO")
  log_message(sprintf("Gene assay - Total cells: %d, Median genes/cell: %.0f, Median counts/cell: %.0f",
                      nrow(all_qc_gene),
                      median(all_qc_gene$n_genes),
                      median(all_qc_gene$total_gene_counts)), "INFO")
  log_message("All QC visualization will be performed in Step 3 using Seurat", "INFO")

  return(list(
    metadata_isoform = all_qc_isoform,
    metadata_gene = all_qc_gene
  ))
}

# =============================================================================
# STEP 3: SEURAT QC AND FILTERING
# =============================================================================

step_seurat_qc <- function() {
  log_message("Starting Seurat-based QC workflow", "INFO")

  # Load sample metadata from Step 1
  sample_metadata <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.rds"))

  # === Step 3.1: Create ID Mappings (Once) ===
  log_message("Creating ID mappings from FLAMES GTF", "INFO")

  id_maps <- create_id_mappings_from_flames(
    flames_dir = FLAMES_OUTPUT_DIR,
    reference_gtf = REFERENCE_GTF,
    cache_file = file.path(OUTPUT_BASE_DIR, "data", "id_mappings_cache.rds")
  )

  log_message(sprintf("ID mappings created: %d genes, %d transcripts",
                      id_maps$metadata$n_genes,
                      id_maps$metadata$n_transcripts), "INFO")

  # === Step 3.2: Process Each Sample ===
  qc_summary_list <- list()
  seurat_objects_list <- list()

  for (i in 1:nrow(sample_metadata)) {
    sample_id <- sample_metadata$sample_id[i]
    timepoint <- sample_metadata$timepoint[i]

    log_message(sprintf("Processing sample %d/%d: %s (%s)",
                        i, nrow(sample_metadata), sample_id, timepoint), "INFO")

    # --- 3.2.1: Load Raw Counts ---
    log_message(sprintf("Loading raw counts for %s", sample_id), "DEBUG")

    # Load gene counts
    gene_count_file <- file.path(FLAMES_OUTPUT_DIR, paste0(sample_id, "_gene_count.csv"))
    gene_counts <- read.csv(gene_count_file, row.names = 1)

    # Load isoform counts
    iso_counts <- load_flames_10x(FLAMES_OUTPUT_DIR, sample_id)

    # --- 3.2.2: Initialize Seurat Object ---
    seurat_obj <- initializeSeuratObjs(
      sample_id = sample_id,
      gene_counts = gene_counts,
      iso_counts = iso_counts,
      id_maps = id_maps,
      min_cells = MIN_CELLS,
      min_features = MIN_FEATURES
    )

    # Store timepoint metadata
    seurat_obj$timepoint <- timepoint

    # --- 3.2.3: QC Filtering ---
    seurat_obj <- runSeuratPreliminaryFiltering(
      seurat_obj = seurat_obj,
      timepoint = timepoint,
      qc_thresholds = QC_THRESHOLDS
    )

    # --- 3.2.4: Normalization and Cell Cycle Scoring ---
    seurat_obj <- NormaliseScaleAndElbow(
      seurat_obj = seurat_obj,
      cell_cycle_genes = CELL_CYCLE_GENES,
      n_pcs_test = N_PCS_TEST
    )

    optimal_pcs <- seurat_obj@misc$optimal_pcs

    # --- 3.2.5: Clustering Optimization ---
    optimal_res <- find.optimal.cluster.res(
      seurat_obj = seurat_obj,
      resolutions = CLUSTERING_RESOLUTIONS
    )

    # --- 3.2.6: UMAP and Doublet Removal ---
    seurat_obj <- runUMAP.removeDoublets(
      seurat_obj = seurat_obj,
      optimal_pcs = optimal_pcs,
      optimal_res = optimal_res,
      doublet_rate = DOUBLET_RATE
    )

    # --- 3.2.7: Save Individual Seurat Object ---
    output_file <- file.path(OUTPUT_BASE_DIR, "data", "seurat_objects", paste0(sample_id, "_qc.rds"))
    log_message(sprintf("Saving Seurat object: %s", basename(output_file)), "DEBUG")
    saveRDS(seurat_obj, output_file)

    # --- 3.2.8: Collect QC Summary Metrics ---
    qc_summary <- data.frame(
      sample_id = sample_id,
      timepoint = timepoint,
      n_cells_initial = seurat_obj@misc$qc_stats$n_cells_before,
      n_cells_filtered = seurat_obj@misc$qc_stats$n_cells_after,
      pct_retained_filtering = seurat_obj@misc$qc_stats$pct_retained,
      n_doublets = seurat_obj@misc$doublet_stats$n_doublets,
      pct_doublets = seurat_obj@misc$doublet_stats$pct_doublets,
      n_cells_final = ncol(seurat_obj),
      n_genes = nrow(seurat_obj[["RNA"]]),
      n_isoforms = nrow(seurat_obj[["iso"]]),
      optimal_pcs = optimal_pcs,
      optimal_resolution = optimal_res,
      n_clusters = length(unique(seurat_obj$seurat_clusters)),
      median_nFeature = median(seurat_obj$nFeature_RNA),
      median_nCount = median(seurat_obj$nCount_RNA),
      median_pct_mt = median(seurat_obj$percent.mt),
      stringsAsFactors = FALSE
    )

    qc_summary_list[[sample_id]] <- qc_summary
    seurat_objects_list[[sample_id]] <- seurat_obj

    # --- 3.2.9: Generate Per-Sample QC Report ---
    generate_sample_qc_report(seurat_obj, sample_id)

    log_message(sprintf("Completed QC for %s: %d cells final, %d clusters",
                        sample_id, ncol(seurat_obj), qc_summary$n_clusters), "INFO")
  }

  # === Step 3.3: Combine QC Summary Table ===
  qc_summary_df <- do.call(rbind, qc_summary_list)

  # Save summary table
  write.csv(qc_summary_df,
            file.path(OUTPUT_BASE_DIR, "data", "seurat_objects", "qc_summary_table.csv"),
            row.names = FALSE)

  saveRDS(qc_summary_df,
          file.path(OUTPUT_BASE_DIR, "data", "seurat_objects", "qc_summary_table.rds"))

  log_message("QC summary table saved", "INFO")

  # === Step 3.4: Generate Summary Comparison Report ===
  generate_qc_summary_report(qc_summary_df, seurat_objects_list)

  log_message("Seurat QC workflow completed for all samples", "INFO")

  return(list(
    qc_summary = qc_summary_df,
    seurat_objects = seurat_objects_list
  ))
}

# =============================================================================
# STEP 3 HELPER FUNCTIONS
# =============================================================================

#' Generate per-sample QC report PDF
#'
#' @param seurat_obj Seurat object with QC results
#' @param sample_id Sample identifier
generate_sample_qc_report <- function(seurat_obj, sample_id) {

  pdf_file <- file.path(OUTPUT_BASE_DIR, "reports", paste0(sample_id, "_QC_Report.pdf"))

  pdf(pdf_file, width = 11, height = 14)

  # Page 1: QC Filtering Comparison
  qc_plots <- seurat_obj@misc$qc_plots
  grid.arrange(
    qc_plots$nFeature,
    qc_plots$nCount,
    qc_plots$percent_mt,
    ncol = 2,
    top = textGrob(paste(sample_id, "- QC Filtering Results"),
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )

  # Page 2: Elbow Plot
  print(seurat_obj@misc$elbow_plot)

  # Page 3: Clustree
  print(seurat_obj@misc$clustree_plot)

  # Page 4: Final UMAP
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    labs(title = paste(sample_id, "- Final UMAP (Singlets Only)")) +
    theme_publication()

  print(umap_plot)

  # Page 5: Cell Cycle UMAP
  cc_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") +
    labs(title = paste(sample_id, "- Cell Cycle Phase")) +
    theme_publication()

  print(cc_plot)

  dev.off()

  log_message(sprintf("Generated QC report: %s", basename(pdf_file)), "DEBUG")

  # Also save individual PNG figures
  ggsave(file.path(OUTPUT_BASE_DIR, "figures", "qc", paste0(sample_id, "_qc_before_after.png")),
         grid.arrange(qc_plots$nFeature, qc_plots$nCount, qc_plots$percent_mt, ncol = 2),
         width = 10, height = 8, dpi = 150)

  ggsave(file.path(OUTPUT_BASE_DIR, "figures", "qc", paste0(sample_id, "_elbow.png")),
         seurat_obj@misc$elbow_plot, width = 8, height = 6, dpi = 150)

  ggsave(file.path(OUTPUT_BASE_DIR, "figures", "qc", paste0(sample_id, "_clustree.png")),
         seurat_obj@misc$clustree_plot, width = 10, height = 8, dpi = 150)

  ggsave(file.path(OUTPUT_BASE_DIR, "figures", "qc", paste0(sample_id, "_umap_final.png")),
         umap_plot, width = 8, height = 6, dpi = 150)
}

#' Generate summary QC comparison report across all samples
#'
#' @param qc_summary_df QC summary data frame
#' @param seurat_objects_list List of Seurat objects
generate_qc_summary_report <- function(qc_summary_df, seurat_objects_list) {

  pdf_file <- file.path(OUTPUT_BASE_DIR, "reports", "03_Seurat_QC_Summary.pdf")

  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Cell retention across samples
  p1 <- ggplot(qc_summary_df, aes(x = sample_id)) +
    geom_col(aes(y = n_cells_initial), fill = "lightgray", alpha = 0.5) +
    geom_col(aes(y = n_cells_final), fill = "steelblue") +
    geom_text(aes(y = n_cells_final, label = n_cells_final), vjust = -0.5, size = 3) +
    facet_wrap(~ timepoint, scales = "free_x") +
    labs(title = "Cell Retention Across QC Pipeline",
         subtitle = "Gray: Initial | Blue: Final (after filtering + doublet removal)",
         x = "Sample", y = "Number of Cells") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p1)

  # Page 2: QC metrics by timepoint
  p2 <- ggplot(qc_summary_df, aes(x = timepoint, y = median_nFeature, fill = timepoint)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(title = "Median Features per Cell by Timepoint", y = "Median nFeature") +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    theme_publication()

  p3 <- ggplot(qc_summary_df, aes(x = timepoint, y = median_pct_mt, fill = timepoint)) +
    geom_boxplot() +
    geom_point(size = 3) +
    labs(title = "Median MT% by Timepoint", y = "Median MT%") +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    theme_publication()

  print(grid.arrange(p2, p3, ncol = 2))

  # Page 3: Optimal PCs and clustering
  p4 <- ggplot(qc_summary_df, aes(x = sample_id, y = optimal_pcs, fill = timepoint)) +
    geom_col() +
    geom_text(aes(label = optimal_pcs), vjust = -0.5, size = 3) +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    labs(title = "Optimal PCs Selected per Sample", y = "Optimal PCs") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p5 <- ggplot(qc_summary_df, aes(x = sample_id, y = n_clusters, fill = timepoint)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.1f", optimal_resolution)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = TIMEPOINT_COLORS) +
    labs(title = "Number of Clusters (Optimal Resolution Shown)", y = "N Clusters") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(grid.arrange(p4, p5, ncol = 2))

  # Page 4: All samples UMAP grid
  umap_plots <- lapply(names(seurat_objects_list), function(sid) {
    DimPlot(seurat_objects_list[[sid]], reduction = "umap", group.by = "seurat_clusters") +
      labs(title = sid) +
      theme_publication() +
      theme(legend.position = "none")
  })

  do.call(grid.arrange, c(umap_plots, ncol = 4))

  dev.off()

  log_message(sprintf("Generated summary QC report: %s", basename(pdf_file)), "INFO")
}

# =============================================================================
# MAIN PIPELINE EXECUTION
# =============================================================================

main <- function() {
  cat("\n[PIPELINE] Starting full analysis pipeline...\n")

  # Ensure output directories exist before saving metadata
  create_output_dirs()

  # Save analysis metadata
  save_analysis_metadata()

  # Step 1: Data Loading and Initial QC
  step1_result <- run_pipeline_step(
    "Data Loading",
    step_data_loading,
    dependencies = NULL,
    resume_on_error = FALSE
  )

  # Step 2: Save Per-Cell Metadata (No Plots)
  step2_result <- run_pipeline_step(
    "Save Metadata",
    step_save_metadata,
    dependencies = c("data_loading"),
    resume_on_error = FALSE
  )

  # Step 3: Seurat QC and Filtering (All QC Visualization Here)
  step3_result <- run_pipeline_step(
    "Seurat QC",
    step_seurat_qc,
    dependencies = c("save_metadata"),
    resume_on_error = FALSE
  )

  # Additional steps would be added here following the same pattern...
  # step4_result <- run_pipeline_step("Integration", step_integration, dependencies = c("seurat_qc"))
  # step5_result <- run_pipeline_step("Cell Annotation", step_cell_annotation, dependencies = c("integration"))
  # etc.

  # Pipeline completion
  end_time <- Sys.time()
  total_duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

  log_message(sprintf("Pipeline completed successfully in %.2f minutes", total_duration), "INFO")

  # Create final summary
  pipeline_summary <- list(
    start_time = start_time,
    end_time = end_time,
    total_duration_minutes = total_duration,
    completed_steps = c("Data Loading", "Save Metadata", "Seurat QC"),
    status = "SUCCESS"
  )

  saveRDS(pipeline_summary, file.path(OUTPUT_BASE_DIR, "summary", "pipeline_summary.rds"))

  cat("\n")
  cat("================================================================================\n")
  cat("                         PIPELINE COMPLETED SUCCESSFULLY                       \n")
  cat("================================================================================\n")
  cat(sprintf("Total runtime: %.2f minutes\n", total_duration))
  cat(sprintf("Output directory: %s\n", OUTPUT_BASE_DIR))
  cat(sprintf("Reports available in: %s\n", file.path(OUTPUT_BASE_DIR, "reports")))
  cat("\n")

  return(pipeline_summary)
}

# =============================================================================
# EXECUTION
# =============================================================================

# Check if script is being run directly (not sourced)
if (sys.nframe() == 0) {
  # Run main pipeline
  result <- main()
} else {
  log_message("Main script sourced - use main() to run pipeline", "INFO")
}