# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a neurodevelopmental models analysis project focused on integrating TWAS (Transcriptome-Wide Association Study) data with single-cell RNA sequencing analysis of brain organoids. The project analyzes organoid data across three developmental timepoints (1, 3, and 6 months) to identify developmental patterns and disease-associated genes.

## Repository Structure

The current directory (`Updated_Mitch_Filtered_Database/`) contains filtered TWAS results:
- `Sig_isoTWAS_developmental.csv`: Significant isoform-level TWAS results for developmental traits
- `Sig_TWAS_developmental.csv`: Significant gene-level TWAS results for developmental traits

The parent directory contains the broader analysis framework:
- `DE_script.R`: Main differential expression analysis comparing organoid data with TWAS results
- `cellTypeProportion_script.R`: Cell type proportion analysis across timepoints
- `Org_Timeseries_Analysis/`: Complete Seurat-based scRNA-seq processing pipeline
- Processed Seurat objects (`.rds` files) for 1, 3, and 6-month organoids
- TWAS data files from Bhattacharya2023 (both original and filtered versions)

## Core Analysis Components

### 1. TWAS Data Integration (`DE_script.R`)
- Integrates isoform-level TWAS data with single-cell differential expression
- Main function: `returnTWAS_DE_Features()` - identifies overlapping differentially expressed transcripts
- Uses Seurat's `FindAllMarkers()` on isoform assays
- Processes multiple psychiatric traits (SCZ, ADHD, ALZ, etc.)

### 2. Cell Type Analysis (`cellTypeProportion_script.R`)
- Uses the `speckle` package for cell type proportion analysis
- Main function: `MakeGlobalCellTypePropPlot()` - creates proportion visualizations
- Compares cell type distributions across developmental timepoints
- Generates both global and replicate-specific proportion plots

### 3. Single-cell Processing Pipeline (`Org_Timeseries_Analysis/`)
Contains a comprehensive Seurat-based pipeline with modular functions:
- Data import and ID conversion
- Quality control and filtering
- Normalization, scaling, and PCA
- Optimal clustering resolution finding
- UMAP and doublet removal
- Batch processing capabilities

## Key Data Files

### TWAS Results Format
Both CSV files contain standardized TWAS results with columns:
- Gene/Transcript identifiers (ENSG IDs, HGNC symbols, Ensembl transcript IDs)
- Statistical measures (Z-scores, P-values, FDR)
- Genomic coordinates and biotype information
- GWAS SNP information

### Seurat Objects
Processed organoid data contains:
- Gene expression (`RNA` assay) and isoform expression (`iso` assay)
- Cell type annotations in `cluster_annotations` metadata
- Timepoint information for developmental comparisons

## Running Analyses

### R Environment Setup
All scripts require these core packages:
- `Seurat`: Single-cell analysis framework
- `tidyverse`: Data manipulation and visualization
- `readxl`: Excel file reading for TWAS data
- `presto`: Fast differential expression
- `speckle`: Cell type proportion analysis

### Typical Analysis Workflow
1. Load processed Seurat objects from parent directory
2. Set active identity to `cluster_annotations` 
3. Run differential expression analysis with TWAS integration
4. Generate visualizations (heatmaps stored in trait-specific directories)
5. Export results and plots

## File Paths and Dependencies

### Standard Data Locations
- Seurat objects: `/data/gpfs/projects/punim2251/Mitch/[timepoint]_seurat.intergrated_harm.isofrom.rds`
- TWAS data: `/data/gpfs/projects/punim2251/Mitch/FILTERED_Bhattacharya2023*`
- Output directories: `geneTWAS_Developmental_Heatmaps/` and `isoTWAS_Developmental_Heatmaps/`

### Key Resource Files
- Gene ID dictionary: `/data/gpfs/projects/punim2251/resoucres/v41_ENSG_ID_GENEsymbol.csv`
- Cell cycle genes: `/data/gpfs/projects/punim0646/manveer/cycle.rda`
- Silhouette analysis functions: `/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R`

## Analysis Parameters

### Quality Control Thresholds (from pipeline)
- Sample-specific maximum counts and mitochondrial thresholds
- Doublet detection rate: 1.6%
- PCA dimensions and clustering resolutions optimized per sample

### TWAS Integration
- Uses transcript-level analysis for isoform TWAS
- Filters for significant associations (typically FDR < 0.05)
- Focuses on developmental traits (SCZ, ADHD, ALZ, etc.)

## Output Visualization

The analysis generates trait-specific heatmaps showing:
- Gene-level TWAS associations (`geneTWAS_Developmental_Heatmaps/`)
- Isoform-level TWAS associations (`isoTWAS_Developmental_Heatmaps/`)
- Cell type proportion plots (PDF format in parent directory)
- UMAP visualizations of organoid data