# =============================================================================
# FLAMES Data I/O and Auto-Detection Functions
# =============================================================================
# Author: Manveer Chauhan and Sefi Prawer
# Description: Functions for auto-detecting and loading FLAMES output data
#              in 10X format (MTX + barcodes + features)
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(dplyr)
})

# =============================================================================
# FLAMES AUTO-DETECTION FUNCTIONS
# =============================================================================

#' Auto-detect FLAMES samples from directory structure
#'
#' @param flames_dir Path to FLAMES output directory
#' @return Data frame with sample metadata
auto_detect_flames_samples <- function(flames_dir) {
  log_message("Starting FLAMES sample auto-detection...", "INFO")

  if (!dir.exists(flames_dir)) {
    stop(sprintf("FLAMES directory does not exist: %s", flames_dir))
  }

  # List all files in FLAMES directory
  all_files <- list.files(flames_dir, full.names = FALSE)

  # Extract unique sample IDs based on file patterns
  sample_ids <- c()

  # Look for count matrix files (.count.mtx)
  mtx_files <- grep(FLAMES_FILE_PATTERNS$count_matrix, all_files, value = TRUE)
  if (length(mtx_files) > 0) {
    sample_ids <- c(sample_ids, gsub("\\.count\\.mtx$", "", mtx_files))
  }

  # Also check for gene count files as backup
  gene_count_files <- grep(FLAMES_FILE_PATTERNS$gene_counts, all_files, value = TRUE)
  if (length(gene_count_files) > 0) {
    backup_samples <- gsub("_gene_count\\.csv$", "", gene_count_files)
    sample_ids <- unique(c(sample_ids, backup_samples))
  }

  if (length(sample_ids) == 0) {
    stop("No FLAMES samples detected. Check directory structure and file naming.")
  }

  # Filter to expected sample pattern (org_1A, org_3B, etc.)
  valid_samples <- sample_ids[grepl(SAMPLE_PATTERN, sample_ids)]

  if (length(valid_samples) == 0) {
    warning("No samples matching expected pattern found. Using all detected samples.")
    valid_samples <- sample_ids
  }

  # Extract timepoint and replicate information
  sample_metadata <- data.frame(
    sample_id = valid_samples,
    stringsAsFactors = FALSE
  )

  # Parse timepoint and replicate from sample names (assuming org_XY format)
  sample_metadata$timepoint_num <- as.numeric(gsub("org_([0-9]+)[A-C]", "\\1", sample_metadata$sample_id))
  sample_metadata$replicate <- gsub("org_[0-9]+([A-C])", "\\1", sample_metadata$sample_id)
  sample_metadata$timepoint <- paste0(sample_metadata$timepoint_num, "M_Org")

  # Convert timepoint to days (assuming months)
  sample_metadata$timepoint_days <- sample_metadata$timepoint_num * 30

  # Validate that required files exist for each sample
  for (sample in sample_metadata$sample_id) {
    required_files <- c(
      paste0(sample, ".count.mtx"),
      paste0(sample, ".barcodes.txt"),
      paste0(sample, ".features.txt")
    )

    missing_files <- required_files[!file.exists(file.path(flames_dir, required_files))]

    if (length(missing_files) > 0) {
      stop(sprintf("Missing required files for sample %s: %s",
                   sample, paste(missing_files, collapse = ", ")))
    }
  }

  log_message(sprintf("Auto-detected %d valid FLAMES samples", nrow(sample_metadata)), "INFO")
  log_message(sprintf("Samples: %s", paste(sample_metadata$sample_id, collapse = ", ")), "INFO")

  return(sample_metadata)
}

# =============================================================================
# FLAMES DATA LOADING FUNCTIONS
# =============================================================================

#' Load FLAMES 10X format data for a single sample
#'
#' @param flames_dir Path to FLAMES output directory
#' @param sample_id Sample identifier
#' @return Sparse count matrix with transcript IDs as rows, cells as columns
load_flames_10x <- function(flames_dir, sample_id) {
  log_message(sprintf("Loading FLAMES data for sample: %s", sample_id), "INFO")

  # File paths
  mtx_file <- file.path(flames_dir, paste0(sample_id, ".count.mtx"))
  barcodes_file <- file.path(flames_dir, paste0(sample_id, ".barcodes.txt"))
  features_file <- file.path(flames_dir, paste0(sample_id, ".features.txt"))

  # Check file existence
  files_to_check <- c(mtx_file, barcodes_file, features_file)
  missing_files <- files_to_check[!file.exists(files_to_check)]

  if (length(missing_files) > 0) {
    stop(sprintf("Missing files for sample %s: %s",
                 sample_id, paste(basename(missing_files), collapse = ", ")))
  }

  # Load count matrix
  log_message(sprintf("Reading count matrix: %s", basename(mtx_file)), "DEBUG")
  counts <- readMM(mtx_file)

  # Load barcodes (cell IDs)
  log_message(sprintf("Reading barcodes: %s", basename(barcodes_file)), "DEBUG")
  barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1

  # Load features (transcript IDs)
  log_message(sprintf("Reading features: %s", basename(features_file)), "DEBUG")
  features <- read.table(features_file, header = FALSE, stringsAsFactors = FALSE)$V1

  # Validate dimensions and transpose if necessary
  if (nrow(counts) == length(barcodes) && ncol(counts) == length(features)) {
    # Matrix is transposed - features should be rows, cells should be columns
    log_message(sprintf("Transposing matrix for %s (was %d x %d)", sample_id, nrow(counts), ncol(counts)), "DEBUG")
    counts <- Matrix::t(counts)
  }

  # Validate final dimensions
  if (nrow(counts) != length(features)) {
    stop(sprintf("Matrix rows (%d) don't match features (%d)", nrow(counts), length(features)))
  }

  if (ncol(counts) != length(barcodes)) {
    stop(sprintf("Matrix columns (%d) don't match barcodes (%d)", ncol(counts), length(barcodes)))
  }

  # Set row and column names
  rownames(counts) <- features
  colnames(counts) <- barcodes

  log_message(sprintf("Loaded %s: %d features × %d cells",
                      sample_id, nrow(counts), ncol(counts)), "INFO")

  return(counts)
}

#' Load gene-level counts from FLAMES CSV output
#'
#' @param flames_dir Path to FLAMES output directory
#' @param sample_id Sample identifier
#' @return Matrix with gene counts
load_flames_gene_counts <- function(flames_dir, sample_id) {
  gene_count_file <- file.path(flames_dir, paste0(sample_id, "_gene_count.csv"))

  if (!file.exists(gene_count_file)) {
    warning(sprintf("Gene count file not found for sample %s", sample_id))
    return(NULL)
  }

  log_message(sprintf("Loading gene counts for sample: %s", sample_id), "INFO")

  # Read gene count matrix
  gene_counts <- read.csv(gene_count_file, row.names = 1, check.names = FALSE)

  log_message(sprintf("Loaded gene counts for %s: %d genes × %d cells",
                      sample_id, nrow(gene_counts), ncol(gene_counts)), "INFO")

  return(as.matrix(gene_counts))
}

#' Load FLAMES summary statistics
#'
#' @param flames_dir Path to FLAMES output directory
#' @param sample_id Sample identifier
#' @return List with summary metrics
load_flames_summary <- function(flames_dir, sample_id) {
  summary_file <- file.path(flames_dir, paste0(sample_id, "_summary.txt"))

  if (!file.exists(summary_file)) {
    warning(sprintf("Summary file not found for sample %s", sample_id))
    return(NULL)
  }

  # Read summary file
  summary_lines <- readLines(summary_file)

  # Extract key metrics
  total_reads_line <- grep("Total reads:", summary_lines, value = TRUE)
  cells_line <- grep("Identified # of cells:", summary_lines, value = TRUE)
  reads_in_cells_line <- grep("Total reads in cells:", summary_lines, value = TRUE)

  summary_metrics <- list(
    sample_id = sample_id,
    total_reads = ifelse(length(total_reads_line) > 0,
                         as.numeric(gsub("[^0-9]", "", total_reads_line)), NA),
    identified_cells = ifelse(length(cells_line) > 0,
                              as.numeric(gsub(".*: ", "", cells_line)), NA),
    reads_in_cells_count = ifelse(length(reads_in_cells_line) > 0,
                                  as.numeric(gsub(".*: ([0-9,]+) .*", "\\1",
                                                  gsub(",", "", reads_in_cells_line))), NA),
    reads_in_cells_percent = ifelse(length(reads_in_cells_line) > 0,
                                    as.numeric(gsub(".*\\(([0-9.]+)%\\).*", "\\1", reads_in_cells_line)), NA)
  )

  return(summary_metrics)
}

# =============================================================================
# BATCH LOADING FUNCTIONS
# =============================================================================

#' Load all FLAMES samples and create sample metadata
#'
#' @param flames_dir Path to FLAMES output directory
#' @param save_raw Should raw count matrices be saved to output directory?
#' @return List containing sample metadata and count matrices
load_all_flames_data <- function(flames_dir, save_raw = TRUE) {
  log_message("Starting batch loading of all FLAMES data...", "INFO")

  # Auto-detect samples
  sample_metadata <- auto_detect_flames_samples(flames_dir)

  # Initialize storage
  count_matrices <- list()
  gene_count_matrices <- list()
  summary_metrics <- list()

  # Load data for each sample
  for (i in seq_len(nrow(sample_metadata))) {
    sample_id <- sample_metadata$sample_id[i]

    tryCatch({
      # Load transcript-level counts (10X format)
      counts <- load_flames_10x(flames_dir, sample_id)
      count_matrices[[sample_id]] <- counts

      # Load gene-level counts (CSV format)
      gene_counts <- load_flames_gene_counts(flames_dir, sample_id)
      if (!is.null(gene_counts)) {
        gene_count_matrices[[sample_id]] <- gene_counts
      }

      # Load summary statistics
      summary <- load_flames_summary(flames_dir, sample_id)
      if (!is.null(summary)) {
        summary_metrics[[sample_id]] <- summary
      }

      # Add cell and feature counts to metadata
      sample_metadata$n_cells[i] <- ncol(counts)
      sample_metadata$n_features[i] <- nrow(counts)
      sample_metadata$total_counts[i] <- sum(counts)

    }, error = function(e) {
      stop(sprintf("Failed to load sample %s: %s", sample_id, e$message))
    })
  }

  # Save raw data if requested
  if (save_raw) {
    log_message("Saving raw count matrices...", "INFO")
    saveRDS(count_matrices, file.path(OUTPUT_BASE_DIR, "data", "raw_counts", "transcript_count_matrices.rds"))

    if (length(gene_count_matrices) > 0) {
      saveRDS(gene_count_matrices, file.path(OUTPUT_BASE_DIR, "data", "raw_counts", "gene_count_matrices.rds"))
    }

    # Save sample metadata
    write.csv(sample_metadata, file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.csv"), row.names = FALSE)
    saveRDS(sample_metadata, file.path(OUTPUT_BASE_DIR, "data", "metadata", "sample_metadata.rds"))

    # Save summary metrics
    if (length(summary_metrics) > 0) {
      summary_df <- do.call(rbind, lapply(summary_metrics, as.data.frame))
      write.csv(summary_df, file.path(OUTPUT_BASE_DIR, "data", "metadata", "flames_summary_metrics.csv"), row.names = FALSE)
    }
  }

  log_message(sprintf("Successfully loaded %d samples", length(count_matrices)), "INFO")

  return(list(
    sample_metadata = sample_metadata,
    count_matrices = count_matrices,
    gene_count_matrices = gene_count_matrices,
    summary_metrics = summary_metrics
  ))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Get basic statistics for a count matrix
#'
#' @param counts Sparse count matrix
#' @return Data frame with basic statistics
get_basic_stats <- function(counts) {
  stats <- data.frame(
    n_features = nrow(counts),
    n_cells = ncol(counts),
    total_counts = sum(counts),
    median_counts_per_cell = median(Matrix::colSums(counts)),
    median_features_per_cell = median(Matrix::colSums(counts > 0)),
    features_detected = sum(Matrix::rowSums(counts) > 0),
    stringsAsFactors = FALSE
  )
  return(stats)
}

log_message("FLAMES I/O functions loaded successfully", "INFO")