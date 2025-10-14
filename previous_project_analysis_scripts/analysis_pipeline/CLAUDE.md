# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a modular R pipeline for analyzing FLAMES (Full-Length Analysis of Matched End Sequencing) isoform data from cortical brain organoids across developmental timepoints (1, 3, and 6 months). The pipeline is designed for thesis writing and reproducible research, generating PDF reports at each analysis step using grid/gridExtra.

**Key Design Philosophy:**
- Functional module architecture with reusable components
- Auto-detection of FLAMES samples from directory structure
- Granular PDF reports per analysis step (no workflowr dependency)
- Complete reproducibility with automatic intermediate file saving
- Publication-ready figure generation with configurable parameters

## Running the Pipeline

### Quick Start

The pipeline is designed to be run from the parent directory:

```r
# From: /data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/previous_project_analysis_scripts/
source("analysis_pipeline/main.R")

# Or run main function directly after sourcing
main()
```

**Important:** The pipeline expects to be run from the parent directory containing `analysis_pipeline/`. It will fail with a clear error if run from the wrong location.

### Individual Analysis Steps

To run individual steps without the full pipeline:

```r
# Load configuration and modules first
source("analysis_pipeline/functions/config.R")
source("analysis_pipeline/functions/flames_io.R")
source("analysis_pipeline/functions/plotting_utils.R")

# Then run specific step functions
step1_result <- step_data_loading()
step2_result <- step_quality_control()
```

### Configuration

All project-specific settings are in `functions/config.R`:

```r
# Modify these for your project:
FLAMES_OUTPUT_DIR <- "/path/to/your/flames/output"
REFERENCE_GTF <- "/path/to/reference.gtf"
REFERENCE_GENOME <- "/path/to/genome.fa"

# Toggle publication mode for high-DPI figures
set_publication_mode(TRUE)  # 300 DPI for publication
set_publication_mode(FALSE) # 150 DPI for development (default)
```

## Pipeline Architecture

### Entry Point: `main.R`

The main script orchestrates the entire analysis workflow:
- **Dependency checking:** Validates all required R packages before starting
- **Progress tracking:** Timestamped logging and error handling for each step
- **Automatic resumption:** Tracks completed steps to enable resuming after failures
- **Metadata saving:** Captures all parameters and session info for reproducibility

### Functional Modules

#### `functions/config.R`
Global configuration including:
- Project metadata and data paths
- Figure parameters (dimensions, DPI, fonts)
- Analysis parameters (QC thresholds, integration settings, DTU cutoffs)
- Publication presets for high-quality outputs
- Logging and progress tracking settings

**Key parameters:**
```r
QC_PARAMS$min_features_per_cell = 200
QC_PARAMS$max_mt_percent = 20
INTEGRATION_PARAMS$dims_to_use = 1:30
DTU_PARAMS$dif_cutoff = 0.25
```

#### `functions/flames_io.R`
FLAMES data loading and auto-detection:
- `auto_detect_flames_samples()`: Automatically finds all org_XY samples in FLAMES directory
- `load_flames_10x()`: Loads transcript-level count matrices in 10X MTX format
- `load_flames_gene_counts()`: Loads gene-level aggregated counts from CSV
- `load_all_flames_data()`: Batch loads all samples with metadata generation

**Expected FLAMES file structure:**
```
flame_outs_standard_ref/
├── org_1A.count.mtx          # Transcript counts (sparse matrix)
├── org_1A.barcodes.txt        # Cell barcodes
├── org_1A.features.txt        # Transcript IDs (ENST format)
├── org_1A_gene_count.csv      # Gene-level counts
├── org_1A_summary.txt         # FLAMES summary statistics
└── [repeat for org_1B, org_3A, org_3B, org_3C, org_6A, org_6B, org_6C]
```

**Sample naming convention:** `org_[timepoint][replicate]` (e.g., org_1A, org_3B, org_6C)
- Timepoint: 1, 3, or 6 (months)
- Replicate: A, B, or C

#### `functions/plotting_utils.R`
Publication-ready plotting with grid/gridExtra:
- `theme_publication()`: Consistent ggplot theme for all figures
- `save_figure()`: Saves plots as PNG, PDF, or both with proper parameters
- `create_multipanel_figure()`: Combines multiple plots into multi-panel PDFs
- `plot_qc_histogram()`, `plot_qc_scatter()`, `plot_qc_violin()`: QC-specific plots
- `plot_umap()`, `plot_feature()`: UMAP visualization functions
- `plot_volcano()`: Differential expression volcano plots
- `plot_expression_heatmap()`: Gene/isoform expression heatmaps

**Color schemes:**
```r
TIMEPOINT_COLORS: 1M_Org (red), 3M_Org (blue), 6M_Org (green)
QC_COLORS: Pass (green), Fail (red), Borderline (orange)
```

### Pipeline Steps (Current Implementation)

#### Step 1: Data Loading and Initial QC (`step_data_loading()`)
- Auto-detects and loads all 8 organoid samples
- Creates sample metadata with timepoint and replicate info
- Generates initial QC plots (cells per sample, features vs cells, total counts)
- **Output:** `reports/01_Data_Loading_QC.pdf`

#### Step 2: Detailed Quality Control (`step_quality_control()`)
- Calculates per-cell QC metrics (total counts, features, MT%)
- Creates comprehensive QC visualizations stratified by timepoint
- Saves detailed QC metrics table
- **Output:** `reports/02_Quality_Control.pdf`, `data/metadata/detailed_qc_metrics.csv`

### Planned Steps (Not Yet Implemented)

#### Step 3: Sample Integration and Batch Correction
**Implementation needed:** Integration strategy questions must be answered first (see Implementation Status section in parent CLAUDE.md)

Expected workflow:
```r
# Create Seurat objects for each sample
# Apply SCTransform normalization
# Harmony batch correction across timepoints
# UMAP dimensional reduction
# Generate integration QC plots
```

#### Step 4: Cell Type Annotation
Clustering and cell type identification with reference-based annotation

#### Step 5: Isoform-Level Analysis
DTU analysis, trajectory analysis, and isoform switching

#### Step 6: Specialized Analysis
Single-cell entropy, microexon detection, SQANTI quality assessment

## Output Directory Structure

The pipeline generates organized outputs:

```
analysis_pipeline/output/
├── data/
│   ├── raw_counts/
│   │   ├── transcript_count_matrices.rds  # Auto-loaded FLAMES data
│   │   └── gene_count_matrices.rds
│   ├── seurat_objects/                    # Processed Seurat objects (.qs format)
│   ├── processed_data/                    # Filtered matrices
│   └── metadata/
│       ├── sample_metadata.csv            # Sample info with timepoints
│       ├── detailed_qc_metrics.csv        # Per-cell QC metrics
│       ├── analysis_metadata.rds          # Complete pipeline parameters
│       └── [step]_completed.rds           # Step completion markers
├── figures/
│   ├── individual/                        # Individual PNG/PDF figures
│   ├── qc/                               # QC-specific plots
│   ├── integration/                      # Integration analysis (planned)
│   ├── cell_types/                       # Cell annotation (planned)
│   └── isoform_analysis/                 # DTU analysis (planned)
├── reports/
│   ├── 01_Data_Loading_QC.pdf           ✅ IMPLEMENTED
│   ├── 02_Quality_Control.pdf           ✅ IMPLEMENTED
│   ├── 03_Integration.pdf               🚧 PLANNED
│   ├── 04_Cell_Annotation.pdf           🚧 PLANNED
│   ├── 05_Isoform_Analysis.pdf          🚧 PLANNED
│   └── 06_Specialized_Analysis.pdf      🚧 PLANNED
├── summary/
│   ├── Methods_Summary.txt              ✅ Auto-generated for thesis
│   ├── pipeline_summary.rds
│   └── Complete_Analysis.pdf            🚧 PLANNED
└── logs/
    └── [timestamped logs]
```

## Required R Packages

### Core Packages (Required for Current Implementation)
```r
required_packages <- c(
  "Seurat",      # Single-cell analysis framework
  "Matrix",      # Sparse matrix operations for FLAMES data
  "dplyr",       # Data manipulation
  "ggplot2",     # Plotting
  "grid",        # Grid graphics for multi-panel figures
  "gridExtra",   # Arranging multiple plots
  "cowplot",     # Publication-ready plot themes
  "qs",          # Fast serialization for Seurat objects
  "data.table",  # Fast data operations
  "viridis",     # Color palettes
  "RColorBrewer",# Color palettes
  "ggrepel",     # Label repulsion in plots
  "scales",      # Scale functions for ggplot
  "stringr",     # String manipulation
  "patchwork"    # Combining plots
)
```

### Additional Packages (For Future Steps)
- **Integration:** `harmony`, `SeuratWrappers`
- **Cell annotation:** `clustifyr`, `celldex`, `SingleR`
- **Clustering:** `clustree`, `igraph`
- **DTU analysis:** `IsoformSwitchAnalyzeR`, `DRIMSeq`, `DEXSeq`, `stageR`
- **Genomics:** `rtracklayer`, `GenomicFeatures`, `biomaRt`, `BSgenome.Hsapiens.UCSC.hg38`
- **Doublet removal:** `DoubletFinder`

## Key Functions Reference

### Data Loading
```r
# Auto-detect all FLAMES samples
metadata <- auto_detect_flames_samples(FLAMES_OUTPUT_DIR)

# Load single sample (10X MTX format)
counts <- load_flames_10x(FLAMES_OUTPUT_DIR, "org_1A")

# Batch load all samples
all_data <- load_all_flames_data(FLAMES_OUTPUT_DIR, save_raw = TRUE)
```

### Figure Generation
```r
# Save single figure
save_figure(plot, "figure_name", type = "single", format = "both", subdir = "qc")

# Create multi-panel PDF
create_multipanel_figure(
  plot_list = list(p1, p2, p3, p4),
  filename = "01_QC_Report",
  ncol = 2,
  title = "Quality Control Analysis",
  subdir = "../reports"
)
```

### Pipeline Control
```r
# Run step with dependency checking and error handling
run_pipeline_step(
  "Data Loading",
  step_data_loading,
  dependencies = NULL,
  resume_on_error = FALSE
)
```

## Error Handling and Resumption

The pipeline tracks completed steps using marker files:
```r
# Check if step completed
step_marker <- file.path(OUTPUT_BASE_DIR, "data", "metadata", "data_loading_completed.rds")
if (file.exists(step_marker)) {
  # Step already completed, can skip
}

# Error information saved for failed steps
error_info <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "metadata", "quality_control_error.rds"))
```

## Development Workflow

### Adding a New Analysis Step

1. **Create step function in `main.R`:**
```r
step_new_analysis <- function() {
  log_message("Starting new analysis...", "INFO")

  # Load required data from previous step
  data <- readRDS(file.path(OUTPUT_BASE_DIR, "data", "..."))

  # Perform analysis
  results <- analyze_data(data)

  # Create plots
  plots <- list(p1, p2, p3)

  # Save results
  saveRDS(results, file.path(OUTPUT_BASE_DIR, "data", "results.rds"))

  # Generate report
  create_multipanel_figure(plots, "XX_New_Analysis", ncol = 2,
                          title = "New Analysis", subdir = "../reports")

  return(results)
}
```

2. **Add step to main pipeline:**
```r
stepX_result <- run_pipeline_step(
  "New Analysis",
  step_new_analysis,
  dependencies = c("previous_step"),
  resume_on_error = FALSE
)
```

### Adding New Plotting Functions

Add custom plot functions to `functions/plotting_utils.R`:
```r
plot_custom <- function(data, params) {
  p <- ggplot(data, aes(...)) +
    geom_...() +
    theme_publication()
  return(p)
}
```

## File Path Conventions

- **All paths in config.R are absolute:** Full paths to FLAMES data and reference files
- **Output paths are relative:** Relative to `OUTPUT_BASE_DIR`
- **FLAMES auto-detection:** Expects samples matching pattern `^org_[0-9][A-C]$`

## Current Data Location

**Project-specific paths** (modify in config.R for different projects):
```r
FLAMES_OUTPUT_DIR: /data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref
REFERENCE_GTF: /data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf
REFERENCE_GENOME: /data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_genome.fa
```

## Logging and Progress Tracking

The pipeline provides detailed logging:
```r
[2025-09-27 08:32:15] INFO: Loading configuration...
[2025-09-27 08:32:16] INFO: Auto-detected 8 valid FLAMES samples
[2025-09-27 08:32:17] INFO: Starting Step: Data Loading
[2025-09-27 08:35:42] INFO: Completed Step: Data Loading (3.42 minutes)
```

**Log levels:** DEBUG, INFO, WARN, ERROR (configurable in config.R)

## Relationship to Legacy Code

The parent directory contains legacy `.Rmd` files from a previous workflowr-based approach:
- `setup_adapted.Rmd`, `geneQC_adapted.Rmd`: Previous analysis notebooks
- This modular pipeline **replaces** the workflowr approach
- Legacy code serves as reference for implementing remaining steps
- Key improvements: No website generation, direct PDF output, better modularity

## Testing and Validation

To test the pipeline on a subset of samples:

```r
# Temporarily modify config.R to use only specific samples
# Or manually filter in auto_detect_flames_samples()

# Run pipeline in development mode (lower DPI, faster)
set_publication_mode(FALSE)
source("analysis_pipeline/main.R")
```

## Common Issues and Solutions

**Issue:** Pipeline fails with "Please run this script from the main project directory"
**Solution:** Run from `/previous_project_analysis_scripts/`, not from `analysis_pipeline/`

**Issue:** "Missing required files for sample org_XY"
**Solution:** Check FLAMES directory contains `.count.mtx`, `.barcodes.txt`, and `.features.txt` for all samples

**Issue:** "Matrix rows don't match features"
**Solution:** The pipeline auto-transposes matrices if needed - check FLAMES output format

## Implementation Status

See parent `CLAUDE.md` for detailed implementation status and planning questions for Steps 3-6. The current pipeline (Steps 1-2) is **fully functional and tested**.
