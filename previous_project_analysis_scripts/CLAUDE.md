# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a single-cell RNA sequencing (scRNA-seq) analysis project focused on FLAMES (Full-Length Analysis of Matched End Sequencing) isoform data from brain organoids. The project uses a modular R script pipeline that outputs PDF reports and figures using grid/gridExtra, replacing the previous workflowr-dependent approach.

**Current Project**: Cortical Organoid Differentiation Analysis
- **Data Location**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref`
- **Timepoints**: 1, 3, and 6 months of cortical organoid differentiation
- **Samples**: 8 samples total (org_1A, org_1B, org_3A, org_3B, org_3C, org_6A, org_6B, org_6C)
- **Expected cells**: 600-2000 cells per sample

## Analysis Pipeline Architecture

The project uses a functional module approach optimized for thesis writing and lab reproducibility:
- **Entry Point**: `main.R` script that runs the entire pipeline in logical order
- **Functional Modules**: Reusable functions grouped by analysis type for easy maintenance
- **Auto-Detection**: Automatically detects FLAMES samples from directory structure
- **Granular Output**: Individual PDFs per analysis step with summary statistics and plots
- **Publication-Ready**: Configurable figure dimensions and publication-quality outputs
- **Complete Reproducibility**: All intermediate results saved automatically between steps
- **R-Contained**: All configuration and settings handled within R for simplicity

### Design Philosophy
- **Current Focus**: Optimized for cortical organoid differentiation analysis
- **Future Extension**: Designed for easy sharing and extension to full modularity
- **Co-Authorship Ready**: Clean, reproducible codebase suitable for FLAMES paper collaboration

## Repository Structure

### Functional Module Pipeline Structure
```
analysis_pipeline/
├── main.R                           # Entry point - runs full pipeline with progress tracking
├── functions/
│   ├── config.R                     # Global configuration and figure parameters
│   ├── flames_io.R                  # FLAMES auto-detection and data loading
│   ├── qc_modules.R                 # Quality control analysis and plotting
│   ├── integration_modules.R        # Sample integration and batch correction
│   ├── cell_annotation_modules.R    # Cell type identification and clustering
│   ├── isoform_modules.R            # DTU, trajectory, entropy analysis
│   ├── plotting_utils.R             # Publication-ready plotting with grid/gridExtra
│   └── report_generators.R          # PDF report generation per analysis step
└── output/
    ├── data/
    │   ├── raw_counts/              # Auto-loaded FLAMES count matrices
    │   ├── seurat_objects/          # Intermediate Seurat objects (.qs)
    │   ├── processed_data/          # Filtered and processed matrices
    │   └── metadata/                # Sample metadata and QC tables
    ├── figures/
    │   ├── individual/              # Individual PNG figures (configurable dimensions)
    │   ├── qc/                      # QC-specific figures
    │   ├── integration/             # Integration analysis figures
    │   ├── cell_types/              # Cell annotation figures
    │   └── isoform_analysis/        # DTU and specialized analysis figures
    ├── reports/
    │   ├── 01_Data_Loading_QC.pdf   # FLAMES loading + initial QC
    │   ├── 02_Quality_Control.pdf   # Detailed QC metrics and filtering
    │   ├── 03_Integration.pdf       # Sample integration and batch correction
    │   ├── 04_Cell_Annotation.pdf   # Cell type identification and clustering
    │   ├── 05_Isoform_Analysis.pdf  # DTU and isoform-level analysis
    │   └── 06_Specialized_Analysis.pdf # Trajectory, entropy, microexons
    └── summary/
        ├── Complete_Analysis.pdf    # Comprehensive pipeline summary
        ├── Methods_Summary.pdf      # Detailed methods for thesis/publication
        └── QC_Summary_Table.csv     # Tabular summary of all QC metrics
```

### Legacy Files (Previous Workflowr Approach)
- `setup_adapted.Rmd`: Previous setup approach (now replaced by modular scripts)
- `geneQC_adapted.Rmd`: Previous QC approach (now replaced by modular scripts)
- `_site.yml`: Workflowr configuration (deprecated)
- Various `.Rmd` files: Legacy analysis notebooks

## Key Analysis Pipeline

### 1. Data Setup and Preprocessing
The analysis begins with FLAMES output processing:
- GTF file processing to create isoform-gene dictionaries
- Gene count matrix conversion from ENSG IDs to gene symbols
- Quality filtering to remove Bambu-generated transcripts

### 2. Quality Control and Integration
- Sample-specific QC metrics and filtering
- Seurat-based integration using Harmony for batch correction
- Doublet detection and removal with DoubletFinder

### 3. Cell Type Analysis
- Clustering at multiple resolutions using clustree visualization
- Cell type annotation using clustifyr and reference datasets
- UMAP visualization with harmony-corrected embeddings

### 4. Isoform-Level Analysis
- Addition of isoform counts to Seurat objects as separate assay ("iso")
- Trajectory analysis for both gene and isoform expression
- Differential transcript usage (DTU) analysis comparing conditions/timepoints
- Single-cell entropy calculations for isoform diversity

### 5. Specialized Analyses
- SQANTI-based transcript quality assessment
- Microexon detection and functional analysis
- Gene-specific isoform switching analysis (e.g., PKM)

## Required R Packages

### Core Single-Cell Analysis
- `Seurat`: Primary framework for scRNA-seq analysis
- `qs2`/`qs`: Fast serialization for large R objects
- `workflowr`: Reproducible research workflow framework

### Data Processing and Quality Control
- `rtracklayer`: GTF/GFF file handling
- `DropletUtils`: Single-cell data utilities
- `DoubletFinder`: Doublet detection
- `celda`: Decontamination and clustering

### Visualization and Plotting
- `ggplot2`, `cowplot`, `patchwork`: Advanced plotting
- `ComplexHeatmap`, `pheatmap`: Heatmap generation
- `clustifyr`: Cell type classification plots
- `clustree`: Clustering resolution visualization

### Specialized Isoform Analysis
- `IsoformSwitchAnalyzeR`: Differential transcript usage analysis
- `DRIMSeq`, `DEXSeq`, `stageR`: DTU statistical testing
- `ORFik`: Open reading frame analysis
- `GenomicFeatures`, `Biostrings`: Genomic annotation handling

### Genomics and Annotation
- `biomaRt`: Gene annotation and ID conversion
- `BSgenome.Hsapiens.UCSC.hg38`: Human genome sequences
- `Gviz`: Genomic visualization

## Key File Locations and Naming Conventions

### Current Project FLAMES Data
- **FLAMES Output Directory**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref`
- **Count Matrices**: `*.count.mtx` (10X format sparse matrices)
- **Cell Barcodes**: `*.barcodes.txt` (one barcode per line)
- **Feature IDs**: `*.features.txt` (transcript IDs in ENST format)
- **Gene Counts**: `*_gene_count.csv` (gene-level aggregated counts)
- **FLAMES GTF**: `isoform_annotated.gtf` (filtered) and `isoform_annotated_unfiltered.gtf`
- **Transcript Assembly**: `transcript_assembly.fa`

### Reference Files for Current Project
- **Reference GTF**: `/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf`
- **Reference Genome**: `/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_genome.fa`
- **FLAMES Config**: `/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/config_file_org_resume.json`

### Seurat Objects
- Main integrated objects stored as `.qs` files in `output/seurat_objects/`
- Standard naming: `seurat_filt_harm_annotated_iso.qs`
- Objects contain both `RNA` and `iso` assays

### Analysis Outputs
- QC results: `output/qc_results/` directory
- DTU results: `output/DTU/` directory with CSV summaries and plots
- Cell type plots and annotations: `figure/` directory
- Cached intermediate results: Various `*_cache/` directories
- Processed count data: `data/gene_counts/` directory

## Analysis Standards and Conventions

### Seurat Object Structure
- **Assays**: `RNA` (gene-level), `iso` (isoform-level)
- **Reductions**: `umap.harm` (harmony-corrected UMAP)
- **Identity**: Typically set to cell type annotations
- **Metadata**: Includes sample IDs, timepoints, cell type classifications

### Common Analysis Parameters
- **Clustering**: Multiple resolutions (0.1 to 1.0 by 0.1)
- **Integration**: Harmony-based batch correction
- **DTU significance**: q-value < 0.01, dIF cutoff > 0.25
- **Quality filtering**: Sample-specific thresholds for UMI counts and mitochondrial percentage

### Modular Pipeline Workflow
The analysis is run through the main entry point:
```r
# Run the complete pipeline
source("main.R")

# Or run individual analysis steps
source("scripts/02_quality_control.R")

# Generate specific reports
source("functions/plotting_utils.R")
generate_qc_report()
```

### File Reading Patterns
- **FLAMES 10X format**: Use `Matrix::readMM()` for `.mtx` files
- **Seurat objects**: `qs2::qs_read()` for `.qs` files
- **Count matrices**: `read.csv()` with `row.names = 1` for gene counts
- **GTF files**: `rtracklayer::import()`
- **FLAMES loading function**:
```r
load_flames_10x <- function(flames_dir, sample_id) {
  mtx_file <- file.path(flames_dir, paste0(sample_id, ".count.mtx"))
  barcodes_file <- file.path(flames_dir, paste0(sample_id, ".barcodes.txt"))
  features_file <- file.path(flames_dir, paste0(sample_id, ".features.txt"))

  counts <- readMM(mtx_file)
  barcodes <- read.table(barcodes_file, header = FALSE)$V1
  features <- read.table(features_file, header = FALSE)$V1

  rownames(counts) <- features
  colnames(counts) <- barcodes
  return(counts)
}
```

## Development Environment

This project is designed for high-performance computing environments with:
- Large memory requirements for single-cell datasets (>500GB for FLAMES processing)
- Access to FLAMES pipeline outputs (>2TB of data)
- R environment with Bioconductor packages
- Workflowr for reproducible research documentation

### Current Project Environment
- **HPC Cluster**: Spartan (University of Melbourne)
- **Partition**: sapphire (high-memory nodes)
- **FLAMES Environment**: conda environment FLAMES_v2.2
- **Data Storage**: `/data/gpfs/projects/punim2251/` (shared project space)

### Sample Information for Current Project
```r
sample_metadata <- data.frame(
  sample_id = c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C"),
  timepoint = c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org"),
  timepoint_days = c(30, 30, 90, 90, 90, 180, 180, 180),
  replicate = c("A", "B", "A", "B", "C", "A", "B", "C"),
  expected_cells = c(900, 1600, 1000, 1300, 2000, 600, 1200, 1100)
)
```

## Getting Started with Current Project

### Quick Start Commands
```bash
# Navigate to analysis directory
cd /data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/previous_project_analysis_scripts

# Run the complete modular pipeline
R
> source("main.R")

# Or run specific analysis steps
> source("scripts/01_data_loading.R")
> source("scripts/02_quality_control.R")
```

### Pipeline Configuration
Before running, configure your project in `config/project_config.R`:
```r
# Project-specific settings
PROJECT_NAME <- "Cortical_Organoid_Differentiation"
FLAMES_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref"
OUTPUT_DIR <- file.path(getwd(), "output")

# Sample metadata
SAMPLE_METADATA <- data.frame(
  sample_id = c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C"),
  timepoint = c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org"),
  expected_cells = c(900, 1600, 1000, 1300, 2000, 600, 1200, 1100)
)
```

### Output Structure
The pipeline generates organized outputs:
```
output/
├── data/
│   ├── seurat_objects/              # Processed Seurat objects (.qs files)
│   ├── count_matrices/              # Processed count matrices
│   └── metadata/                    # Sample metadata and QC metrics
├── figures/
│   ├── qc/                         # Quality control plots
│   ├── integration/                # Integration and batch correction plots
│   ├── cell_types/                 # Cell type identification plots
│   └── isoform_analysis/           # DTU and isoform analysis plots
├── reports/
│   ├── 01_QC_Report.pdf            # Quality control summary
│   ├── 02_Integration_Report.pdf   # Integration analysis
│   ├── 03_CellTypes_Report.pdf     # Cell type identification
│   └── 04_Isoform_Report.pdf       # Isoform analysis results
└── summary/
    ├── pipeline_summary.pdf        # Complete analysis overview
    └── qc_metrics_table.csv        # Tabular QC metrics
```

### PDF Report Generation
Each analysis step generates publication-ready PDF reports using grid/gridExtra:
```r
# Example: Generate QC report
library(grid)
library(gridExtra)

# Create multi-panel figure
qc_plots <- list(plot1, plot2, plot3, plot4)
pdf("output/reports/01_QC_Report.pdf", width = 11, height = 8.5)
do.call(grid.arrange, c(qc_plots, ncol = 2))
dev.off()
```

## Implementation Status

### ✅ Completed Components (Ready for Testing)

**Core Infrastructure:**
- `config.R`: Complete global configuration with publication presets
- `flames_io.R`: Auto-detection and loading of FLAMES 10X format data
- `plotting_utils.R`: Publication-ready plotting utilities with grid/gridExtra
- `main.R`: Entry point with progress tracking and error handling

**Implemented Analysis Steps:**
1. **Data Loading (Step 1)**:
   - Auto-detects 8 organoid samples from FLAMES directory
   - Loads transcript-level count matrices (.mtx format)
   - Creates initial QC overview plots
   - Generates `01_Data_Loading_QC.pdf` report

2. **Quality Control (Step 2)**:
   - Detailed per-cell QC metrics (counts, features, MT%)
   - Comprehensive QC visualizations (histograms, scatter, violin plots)
   - Timepoint-stratified analysis
   - Generates `02_Quality_Control.pdf` report

### 🚧 Planned Implementation (Next Steps)

**Step 3: Sample Integration and Batch Correction**
```r
step_integration <- function() {
  # Create Seurat objects for each sample
  # Apply SCTransform normalization
  # Harmony batch correction across timepoints
  # UMAP dimensional reduction
  # Generate integration QC plots
  # Output: 03_Integration.pdf + integrated Seurat object
}
```

**📋 Step 3 Planning Questions (To Be Addressed Before Implementation):**

*Integration Strategy 1: All Samples Together*
1. **Batch correction variables**: Should Harmony correct for `timepoint` only, or also `replicate` (A/B/C), or both?
2. **Normalization**: Should I use SCTransform (recommended for integration) or standard LogNormalize?
3. **Features for integration**: How many variable features should I use? (default 2000-3000?)

*Integration Strategy 2: Timepoint-Specific Integration*
4. **Within-timepoint integration**: For timepoints with multiple replicates (3M has 3 reps, 6M has 3 reps), should I integrate by `replicate` or just combine without integration since they're biological replicates?
5. **Separate object naming**: How should I name the timepoint-specific objects? (e.g., `seurat_1M`, `seurat_3M`, `seurat_6M`?)

*Analysis Parameters*
6. **Clustering resolutions**: Should I test the same resolutions for both strategies (0.1, 0.3, 0.5, 0.8, 1.0) or focus on a specific range?
7. **UMAP parameters**: Any specific UMAP settings (n.neighbors, min.dist) or use Seurat defaults?
8. **Variable features**: Should I use the same variable features across both integration strategies for comparison?

*Output Organization*
9. **Comparison plots**: Do you want side-by-side UMAP comparisons between the two integration strategies in the same PDF?
10. **Seurat object storage**: Should both integration strategies be saved as separate .qs files, or do you want to choose the best one after seeing results?

*Quality Control Integration*
11. **Integration metrics**: Should I include integration-specific QC like mixing metrics, silhouette scores, or kBET statistics?
12. **Batch effect assessment**: Do you want before/after integration comparisons to show batch effect removal?

*Assay Strategy*
13. **Primary assay**: Should this integration step focus on the `RNA` assay (gene-level), with isoform-level integration (`iso` assay) happening later in Step 5?
14. **Dual assay approach**: Or should I integrate both `RNA` and `iso` assays in this step?

**Step 4: Cell Type Annotation**
```r
step_cell_annotation <- function() {
  # Clustering at multiple resolutions
  # Find cluster biomarkers
  # Cell type annotation (manual + reference-based)
  # Generate cell type plots (UMAP, dotplots, heatmaps)
  # Output: 04_Cell_Annotation.pdf + annotated Seurat object
}
```

**Step 5: Isoform-Level Analysis**
```r
step_isoform_analysis <- function() {
  # Add isoform assay to Seurat objects
  # Differential transcript usage (DTU) analysis
  # Trajectory analysis (pseudotime)
  # Isoform switching analysis
  # Output: 05_Isoform_Analysis.pdf + DTU results
}
```

**Step 6: Specialized Analysis**
```r
step_specialized_analysis <- function() {
  # Single-cell entropy analysis
  # Microexon detection
  # Gene-specific isoform analysis (e.g., PKM)
  # SQANTI transcript quality assessment
  # Output: 06_Specialized_Analysis.pdf + specialized results
}
```

### 🎯 Implementation Priorities

**Next Development Session:**
1. **Test current implementation** (Steps 1-2)
2. **Implement Step 3**: Integration module with Seurat workflow
3. **Add Seurat helper functions** in new `seurat_helpers.R`
4. **Expand plotting utilities** for integration-specific plots

**Required Additional Functions:**
- `seurat_helpers.R`: Seurat object creation, normalization, integration
- `integration_modules.R`: Harmony integration, clustering, UMAP
- `cell_annotation_modules.R`: Cell type identification and annotation
- `isoform_modules.R`: DTU analysis, trajectory analysis, isoform switching
- `report_generators.R`: Automated PDF report generation per step

### 🔧 Technical Notes for Future Implementation

**Seurat Integration Workflow:**
```r
# Standard workflow to implement
seurat_obj <- CreateSeuratObject(counts, assay = "RNA")
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "timepoint")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.3, 0.5))
```

**Isoform Assay Addition:**
```r
# Add transcript-level data as separate assay
seurat_obj[["iso"]] <- CreateAssayObject(counts = transcript_counts)
DefaultAssay(seurat_obj) <- "iso"
# Run isoform-specific analysis
```

**DTU Analysis Framework:**
```r
# Use existing DTU code from legacy scripts
# Implement IsoformSwitchAnalyzeR workflow
# Generate switching plots and summary tables
```

### 📊 Expected Final Output Structure

**Complete Pipeline Outputs:**
```
analysis_pipeline/output/
├── reports/
│   ├── 01_Data_Loading_QC.pdf      ✅ IMPLEMENTED
│   ├── 02_Quality_Control.pdf      ✅ IMPLEMENTED
│   ├── 03_Integration.pdf          🚧 PLANNED
│   ├── 04_Cell_Annotation.pdf      🚧 PLANNED
│   ├── 05_Isoform_Analysis.pdf     🚧 PLANNED
│   └── 06_Specialized_Analysis.pdf 🚧 PLANNED
├── summary/
│   ├── Complete_Analysis.pdf       🚧 PLANNED
│   ├── Methods_Summary.txt         ✅ IMPLEMENTED
│   └── QC_Summary_Table.csv        ✅ IMPLEMENTED
└── data/seurat_objects/
    ├── integrated_seurat_object.qs 🚧 PLANNED
    └── final_annotated_object.qs   🚧 PLANNED
```

The analysis pipeline emphasizes isoform-level analysis and differential transcript usage, making it particularly suitable for studying alternative splicing and transcript diversity in single-cell contexts during cortical organoid differentiation.