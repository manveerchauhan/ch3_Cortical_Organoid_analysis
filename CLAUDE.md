# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a comprehensive neurodevelopmental models analysis project that integrates TWAS (Transcriptome-Wide Association Study) data with single-cell RNA sequencing (scRNA-seq) analysis of brain organoids. The project analyzes organoid data across three developmental timepoints (1, 3, and 6 months) to identify developmental patterns and disease-associated genes, with particular focus on psychiatric traits including schizophrenia (SCZ), ADHD, Alzheimer's disease, autism spectrum disorder (ASD), and bipolar disorder.

## Repository Structure

### Main Analysis Scripts
- `DE_script.R`: Original differential expression analysis comparing organoid data with TWAS results
- `DE_updated_script.R`: Updated analysis using GWAS-filtered TWAS databases with enhanced logging and error handling
- `cellTypeProportion_script.R`: Cell type proportion analysis across timepoints using speckle package
- `face_validity_check.R`: Comprehensive isoform-level expression analysis with co-expression and pseudotime components

### Data Processing Pipelines
- `Org_Timeseries_Analysis/`: Complete Seurat-based scRNA-seq processing pipeline with batch processing capabilities
  - `seurat_batch_QC-script1.R`: Main batch processing script for multiple organoid samples
  - Contains modular functions for QC, normalization, clustering, and doublet removal

### Analysis Databases
- `Updated_Mitch_Filtered_Database/`: GWAS-filtered TWAS results
  - Contains trait-specific filtered CSV files (e.g., `SCZ_filtered_5e-8_Sig_isoTWAS_developmental.csv`)
  - `isotwas_filtering_analysis.R`: Analysis and reporting script for filtering validation

### Output Directories
- `Updated_GWAS_Filtered_Heatmaps/`: GWAS-filtered trait-specific heatmaps
- `geneTWAS_Developmental_Heatmaps/`: Gene-level TWAS association heatmaps
- `isoTWAS_Developmental_Heatmaps/`: Isoform-level TWAS association heatmaps

## Core Analysis Workflow

### 1. Data Import and Setup
```r
# Standard working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis")

# Load processed Seurat objects
one_month_org <- readRDS("one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("six_month_seurat.intergrated_harm.isofrom.rds")

# Set default identity
Idents(seurat_object) <- "cluster_annotations"
```

### 2. TWAS-Differential Expression Integration
The main analysis function `returnTWAS_DE_Features()` identifies overlapping differentially expressed transcripts:
- Runs `FindAllMarkers()` on isoform assays (`assay = "iso"`)
- Filters for adjusted p-value ≤ 0.05
- Matches with TWAS transcript databases
- Returns trait-stratified results

### 3. Visualization Pipeline
- Generates pseudobulked expression heatmaps using `pheatmap`
- Creates trait-specific visualizations for each timepoint
- Exports publication-ready PNG files with consistent formatting

## Required R Packages

### Core Analysis
- `Seurat`: Single-cell analysis framework
- `tidyverse`: Data manipulation and visualization
- `presto`: Fast differential expression testing

### Specialized Packages
- `DoubletFinder`: Doublet detection and removal
- `speckle`: Cell type proportion analysis
- `clustree`: Clustering resolution visualization
- `pheatmap`: Heatmap generation
- `gridExtra`, `grid`, `patchwork`: Plot arrangement

### Data Import
- `readxl`: Excel file reading for TWAS data
- `reticulate`: Python integration (for Seurat pipeline)

## Key Data File Locations

### Seurat Objects (Main Analysis)
- `/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/[timepoint]_seurat.intergrated_harm.isofrom.rds`

### TWAS Reference Data
- Original: `/data/gpfs/projects/punim2251/Mitch/FILTERED_Bhattacharya2023*`
- GWAS-filtered: `Updated_Mitch_Filtered_Database/[trait]_filtered_5e-8_Sig_isoTWAS_developmental.csv`

### Resource Files
- Gene ID dictionary: `/data/gpfs/projects/punim2251/resoucres/v41_ENSG_ID_GENEsymbol.csv`
- Cell cycle genes: `/data/gpfs/projects/punim0646/manveer/cycle.rda`
- Silhouette analysis functions: `/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R`

## Seurat Object Structure

### Assays
- `RNA`: Gene expression data
- `iso`: Isoform expression data (primary assay for TWAS analysis)

### Metadata Columns
- `cluster_annotations`: Cell type identities (primary grouping variable)
- `orig.ident`: Original sample identifiers
- `timepoint`: Developmental timepoint (1M_Org, 3M_Org, 6M_Org)

### Dimensional Reductions
- `umap.harm`: Harmony-corrected UMAP coordinates

## Running Analyses

### Standard Analysis Parameters
- Identity column: `cluster_annotations`
- Default assay: `iso` for TWAS integration, `RNA` for gene-level analysis
- DE test method: `wilcox`
- Significance threshold: `p_val_adj <= 0.05`

### Batch Processing Configuration
The Seurat pipeline uses sample-specific configuration with parameters for:
- Maximum counts and mitochondrial thresholds
- Optimal PCA dimensions
- Clustering resolutions
- Doublet detection rate (typically 1.6%)

### Output Generation
```r
# Typical heatmap generation workflow
createPDFSummary(isoLvl = TRUE, output_dir = "Updated_GWAS_Filtered_Heatmaps")

# Cell type proportion analysis
MakeGlobalCellTypePropPlot(seurat_object)
MakeReplicateCellTypePropPlot(seurat_object)
```

## Analysis Features

### Filtering and Quality Control
- Removes Bambu-generated transcripts from analysis
- Applies GWAS significance filtering (typically 5e-8)
- Sample-specific QC thresholds for cell filtering

### Visualization Standards
- Heatmaps: Row-scaled pseudobulked expression with hierarchical clustering
- Cell type proportions: Global and replicate-specific views
- Co-expression analysis: Blended feature plots with standardized color schemes

### Key Gene Sets
- Splicing factors: PTBP1, PTBP2, RBFOX1-3, NOVA1-2
- Synaptic genes: DLG4, NRXN1-3, NLGN1-4
- Focus on neurodevelopmental trajectories and psychiatric disease associations

## Development Environment

This project runs on high-performance computing infrastructure with access to:
- Large memory requirements for Seurat object processing
- Parallel processing capabilities for batch analyses
- Shared storage systems for large RDS files and TWAS databases

The analysis pipeline is designed for R environments with specific package versions compatible with large-scale single-cell datasets and integrative genomics workflows.