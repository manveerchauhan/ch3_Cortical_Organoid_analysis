# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a neurodevelopmental models analysis project focused on organoid timeseries analysis using single-cell RNA sequencing (scRNA-seq) data. The codebase implements a comprehensive Seurat-based pipeline for quality control, processing, and analysis of multiple organoid samples across different timepoints.

## Main Components

### Core Analysis Script
- `seurat_batch_QC-script1.R`: Main batch processing script for multiple organoid samples
- `Seurat_QC.Rmd`: R Markdown file (currently empty) for QC documentation

### Key Functions

The main R script contains several modular functions:

1. **Data Import Functions**:
   - `countMatrixConvertIDstoSymbols()`: Converts ENSG IDs to gene symbols using a dictionary
   - `readInCountMatrices()`: Reads and formats gene and isoform count matrices
   - `initializeSeuratObjs()`: Creates initial Seurat objects from count matrices

2. **Quality Control Functions**:
   - `runSeuratPreliminaryFiltering()`: Performs initial cell filtering based on gene counts, mitochondrial content, and other metrics
   - `NormaliseScaleAndElbow()`: Handles normalization, scaling, PCA, and elbow plot generation

3. **Clustering Functions**:
   - `find.optimal.cluster.res()`: Determines optimal clustering resolution using silhouette analysis
   - `runUMAP.removeDoublets()`: Runs UMAP dimensionality reduction and removes doublets using DoubletFinder

4. **Batch Processing**:
   - `process_single_sample()`: Wrapper function that processes a single sample through the entire pipeline
   - Main execution loop for batch processing multiple samples

### Sample Configuration

The script uses a configuration-driven approach with `sample_config` data frame containing:
- Sample IDs (org1A, org1B, org3A, etc.)
- Timepoints (1month, 3month, 6month)
- File paths for gene and isoform matrices
- Sample-specific parameters (max_counts, mt_threshold, optimal_pcs, cluster_resolution)

## Dependencies

### R Packages Required
- Seurat: Main single-cell analysis framework
- tidyverse: Data manipulation and visualization
- reticulate: Python integration
- DoubletFinder: Doublet detection and removal
- gridExtra, grid: Plot arrangement
- clustree: Cluster tree visualization
- RColorBrewer: Color palettes
- ggmin: Minimal ggplot theme

### External Resources
- Cell cycle genes (`cycle.rda` from `/data/gpfs/projects/punim0646/manveer/`)
- Silhouette analysis functions (`silhouette.R` from `/data/gpfs/projects/punim2251/scripts/`)
- Gene ID conversion dictionary (`v41_ENSG_ID_GENEsymbol.csv`)

## File Structure

```
/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/
├── seurat_batch_QC-script1.R    # Main analysis script
├── Seurat_QC.Rmd               # R Markdown documentation (empty)
└── [Generated output files]     # Created during analysis runs
```

## Running the Analysis

1. **Prerequisites**: Ensure all required R packages are installed and external resource files are accessible
2. **Configuration**: Update file paths in `sample_config` data frame to point to actual count matrix files
3. **Execution**: Run the entire script or individual functions as needed
4. **Output**: Results are saved in a timestamped directory (`batch_processed_YYYY-MM-DD`) containing:
   - Individual processed Seurat objects (.rds files)
   - QC plots and summaries (PDF files)
   - Processing summary (.csv file)

## Key Parameters

- **Filtering thresholds**: Customizable per sample (max_counts, mt_threshold)
- **PCA dimensions**: Sample-specific optimal PC counts
- **Clustering resolution**: Tunable clustering parameters
- **Doublet detection**: Uses 1.6% expected doublet rate

## Output Files

- `processed_[sample_id].rds`: Final processed Seurat objects
- `[sample_id]-filter-figs-QC.pdf`: QC filtering visualizations
- `[sample_id]-cellCycle-PCs.pdf`: Cell cycle and PCA analysis plots
- `[sample_id]Clustree.pdf`: Clustering resolution comparison
- `processing_summary.csv`: Batch processing results summary