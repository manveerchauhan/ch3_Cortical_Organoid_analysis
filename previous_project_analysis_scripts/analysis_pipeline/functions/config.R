# =============================================================================
# Global Configuration and Figure Parameters
# =============================================================================
# Author: Manveer Chauhan and Sefi Prawer
# Description: Global settings, figure parameters, and publication presets
#              for the cortical organoid FLAMES analysis pipeline
# =============================================================================

# =============================================================================
# PROJECT CONFIGURATION
# =============================================================================

# Project metadata
PROJECT_NAME <- "Cortical_Organoid_Differentiation"
ANALYSIS_DATE <- Sys.Date()
PIPELINE_VERSION <- "1.0.0"

# Data paths - MODIFY THESE FOR YOUR PROJECT
FLAMES_OUTPUT_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"
REFERENCE_GTF <- "/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf"
REFERENCE_GENOME <- "/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_genome.fa"

# Output directory (relative to main.R location)
OUTPUT_BASE_DIR <- file.path(getwd(), "analysis_pipeline", "output")

# =============================================================================
# FIGURE PARAMETERS
# =============================================================================

# Default figure dimensions (inches)
DEFAULT_FIGURE_WIDTH <- 8
DEFAULT_FIGURE_HEIGHT <- 6
MULTIPANEL_WIDTH <- 11
MULTIPANEL_HEIGHT <- 8.5

# DPI settings
DEFAULT_DPI <- 150  # Development/preview
PUBLICATION_DPI <- 300  # Publication quality

# Color palettes
TIMEPOINT_COLORS <- c("1M_Org" = "#E31A1C", "3M_Org" = "#1F78B4", "6M_Org" = "#33A02C")
QC_COLORS <- c("Pass" = "#2E8B57", "Fail" = "#DC143C", "Borderline" = "#FF8C00")

# Font settings
DEFAULT_FONT_SIZE <- 12
TITLE_FONT_SIZE <- 14
AXIS_FONT_SIZE <- 10

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Quality Control Thresholds
QC_PARAMS <- list(
  min_features_per_cell = 200,
  max_features_per_cell = 8000,
  min_counts_per_cell = 500,
  max_counts_per_cell = 50000,
  max_mt_percent = 20,
  min_cells_per_feature = 3
)

# Integration Parameters
INTEGRATION_PARAMS <- list(
  dims_to_use = 1:30,
  resolution = c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0),
  k_neighbors = 20
)

# Differential Analysis Parameters
DE_PARAMS <- list(
  min_pct = 0.1,
  logfc_threshold = 0.25,
  test_use = "wilcox",
  p_val_cutoff = 0.05
)

# DTU Analysis Parameters
DTU_PARAMS <- list(
  dif_cutoff = 0.25,
  alpha = 0.01,
  min_gene_expression = 30,
  min_isoform_expression = 20
)

# =============================================================================
# SEURAT QC PARAMETERS (Step 3)
# =============================================================================

# Timepoint-specific QC thresholds (fine-tunable)
QC_THRESHOLDS <- list(
  "1M_Org" = list(
    mt_percent = 10,          # Fixed MT% threshold
    sd_multiplier = 1.5       # Dynamic threshold: mean ± 1.5 SD
  ),
  "3M_Org" = list(
    mt_percent = 10,
    sd_multiplier = 1.5
  ),
  "6M_Org" = list(
    mt_percent = 10,
    sd_multiplier = 1.5
  )
)

# Seurat object creation parameters
MIN_CELLS <- 5                # Minimum cells required per feature
MIN_FEATURES <- 500           # Minimum features required per cell

# Clustering parameters
CLUSTERING_RESOLUTIONS <- c(0.1, 0.3, 0.5, 0.8, 1.0)  # Resolutions to test
USE_SILHOUETTE_OPTIMIZATION <- TRUE                    # Auto-select optimal resolution

# Doublet removal parameters
DOUBLET_RATE <- 0.016         # Expected doublet rate (1.6%)

# Cell cycle and normalization
CELL_CYCLE_GENES <- "/data/gpfs/projects/punim0646/manveer/cycle.rda"
USE_SCTRANSFORM <- TRUE       # Use SCTransform for normalization

# PCA parameters
N_PCS_TEST <- 50              # Number of PCs to test for elbow detection
PCA_METHOD <- "elbow"         # Method for optimal PC selection ("elbow" or "variance")

# =============================================================================
# PUBLICATION PRESETS
# =============================================================================

# Publication mode toggle
PUBLICATION_MODE <- FALSE

# Function to set publication parameters
set_publication_mode <- function(enabled = TRUE) {
  if (enabled) {
    assign("PUBLICATION_MODE", TRUE, envir = .GlobalEnv)
    assign("DEFAULT_DPI", PUBLICATION_DPI, envir = .GlobalEnv)
    assign("DEFAULT_FONT_SIZE", 10, envir = .GlobalEnv)
    message("✓ Publication mode enabled: High DPI, optimized fonts")
  } else {
    assign("PUBLICATION_MODE", FALSE, envir = .GlobalEnv)
    assign("DEFAULT_DPI", 150, envir = .GlobalEnv)
    assign("DEFAULT_FONT_SIZE", 12, envir = .GlobalEnv)
    message("✓ Development mode enabled: Standard DPI, larger fonts")
  }
}

# Function to get figure parameters
get_figure_params <- function(type = "single") {
  if (type == "single") {
    return(list(width = DEFAULT_FIGURE_WIDTH, height = DEFAULT_FIGURE_HEIGHT, dpi = DEFAULT_DPI))
  } else if (type == "multipanel") {
    return(list(width = MULTIPANEL_WIDTH, height = MULTIPANEL_HEIGHT, dpi = DEFAULT_DPI))
  }
}

# =============================================================================
# SAMPLE AUTO-DETECTION CONFIGURATION
# =============================================================================

# Expected FLAMES output file patterns
FLAMES_FILE_PATTERNS <- list(
  count_matrix = "\\.count\\.mtx$",
  barcodes = "\\.barcodes\\.txt$",
  features = "\\.features\\.txt$",
  gene_counts = "_gene_count\\.csv$",
  summary = "_summary\\.txt$"
)

# Sample naming pattern (for auto-detection)
SAMPLE_PATTERN <- "^org_[0-9][A-C]$"

# =============================================================================
# LOGGING AND PROGRESS CONFIGURATION
# =============================================================================

# Progress tracking
ENABLE_PROGRESS_BARS <- TRUE
ENABLE_TIMESTAMPS <- TRUE
LOG_LEVEL <- "INFO"  # DEBUG, INFO, WARN, ERROR

# Create log directory
LOG_DIR <- file.path(OUTPUT_BASE_DIR, "logs")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Function to create timestamped messages
log_message <- function(message, level = "INFO") {
  if (ENABLE_TIMESTAMPS) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, message))
  } else {
    cat(sprintf("%s: %s\n", level, message))
  }
}

# Function to create output directories
create_output_dirs <- function() {
  dirs <- c(
    file.path(OUTPUT_BASE_DIR, "data", "raw_counts"),
    file.path(OUTPUT_BASE_DIR, "data", "seurat_objects"),
    file.path(OUTPUT_BASE_DIR, "data", "processed_data"),
    file.path(OUTPUT_BASE_DIR, "data", "metadata"),
    file.path(OUTPUT_BASE_DIR, "figures", "individual"),
    file.path(OUTPUT_BASE_DIR, "figures", "qc"),
    file.path(OUTPUT_BASE_DIR, "figures", "integration"),
    file.path(OUTPUT_BASE_DIR, "figures", "cell_types"),
    file.path(OUTPUT_BASE_DIR, "figures", "isoform_analysis"),
    file.path(OUTPUT_BASE_DIR, "reports"),
    file.path(OUTPUT_BASE_DIR, "summary"),
    LOG_DIR
  )

  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      log_message(sprintf("Created directory: %s", dir))
    }
  }
}

# Function to save analysis metadata
save_analysis_metadata <- function() {
  metadata <- list(
    project_name = PROJECT_NAME,
    analysis_date = ANALYSIS_DATE,
    pipeline_version = PIPELINE_VERSION,
    flames_output_dir = FLAMES_OUTPUT_DIR,
    publication_mode = PUBLICATION_MODE,
    qc_params = QC_PARAMS,
    integration_params = INTEGRATION_PARAMS,
    de_params = DE_PARAMS,
    dtu_params = DTU_PARAMS,
    r_version = R.version.string,
    session_info = capture.output(sessionInfo())
  )

  # Fix: Use correct path structure
  metadata_file <- file.path(OUTPUT_BASE_DIR, "data", "metadata", "analysis_metadata.rds")
  log_message(sprintf("Saving metadata to: %s", metadata_file), "DEBUG")

  # Check if directory exists
  metadata_dir <- dirname(metadata_file)
  if (!dir.exists(metadata_dir)) {
    log_message(sprintf("Creating metadata directory: %s", metadata_dir), "INFO")
    dir.create(metadata_dir, recursive = TRUE)
  }

  saveRDS(metadata, metadata_file)

  # Also save as readable text file for thesis methods
  methods_text <- sprintf("
METHODS SUMMARY - %s
Generated: %s
Pipeline Version: %s

=== QUALITY CONTROL PARAMETERS ===
- Minimum features per cell: %d
- Maximum features per cell: %d
- Minimum counts per cell: %d
- Maximum counts per cell: %d
- Maximum mitochondrial percentage: %d%%
- Minimum cells per feature: %d

=== INTEGRATION PARAMETERS ===
- Dimensions used: %s
- Clustering resolutions tested: %s
- K-neighbors: %d

=== DIFFERENTIAL EXPRESSION PARAMETERS ===
- Minimum percent expression: %s
- Log fold-change threshold: %s
- Statistical test: %s
- P-value cutoff: %s

=== DTU ANALYSIS PARAMETERS ===
- dIF cutoff: %s
- Alpha level: %s
- Minimum gene expression: %d
- Minimum isoform expression: %d

=== SOFTWARE VERSIONS ===
R Version: %s
Analysis Date: %s
FLAMES Output Directory: %s
",
    PROJECT_NAME, ANALYSIS_DATE, PIPELINE_VERSION,
    QC_PARAMS$min_features_per_cell, QC_PARAMS$max_features_per_cell,
    QC_PARAMS$min_counts_per_cell, QC_PARAMS$max_counts_per_cell,
    QC_PARAMS$max_mt_percent, QC_PARAMS$min_cells_per_feature,
    paste(INTEGRATION_PARAMS$dims_to_use[1], "-", tail(INTEGRATION_PARAMS$dims_to_use, 1)),
    paste(INTEGRATION_PARAMS$resolution, collapse = ", "),
    INTEGRATION_PARAMS$k_neighbors,
    DE_PARAMS$min_pct, DE_PARAMS$logfc_threshold, DE_PARAMS$test_use, DE_PARAMS$p_val_cutoff,
    DTU_PARAMS$dif_cutoff, DTU_PARAMS$alpha, DTU_PARAMS$min_gene_expression, DTU_PARAMS$min_isoform_expression,
    R.version.string, ANALYSIS_DATE, FLAMES_OUTPUT_DIR
  )

  writeLines(methods_text, file.path(OUTPUT_BASE_DIR, "summary", "Methods_Summary.txt"))
  log_message("Analysis metadata and methods summary saved")
}

# Initialize configuration
log_message("Configuration loaded successfully", "INFO")
log_message(sprintf("Project: %s", PROJECT_NAME), "INFO")
log_message(sprintf("FLAMES directory: %s", FLAMES_OUTPUT_DIR), "INFO")