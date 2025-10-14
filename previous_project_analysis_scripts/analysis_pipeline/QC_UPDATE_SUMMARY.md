# QC Update Summary - Isoform and Gene-Level Separation

## Overview
Updated the analysis pipeline to clearly distinguish between **isoform-level (transcript)** and **gene-level** data in all QC steps, plots, and debug statements.

## Key Changes

### 1. Step 1: Data Loading and Initial QC

**New plots generated (6 total):**
1. `01_cells_per_sample.png/pdf` - Cells detected per sample
2. `02_isoforms_vs_cells.png/pdf` - **Isoforms** (transcripts) vs cells by sample
3. `03_genes_vs_cells.png/pdf` - **Genes** vs cells by sample
4. `04_isoforms_per_gene_ratio.png/pdf` - Isoforms per gene ratio by timepoint
5. `05_total_isoform_counts.png/pdf` - Total **isoform** counts by timepoint
6. `06_total_gene_counts.png/pdf` - Total **gene** counts by timepoint

**New metadata columns:**
- `n_features` → Number of isoforms (transcripts) detected
- `n_genes` → Number of genes detected
- `total_counts` → Total isoform counts
- `total_gene_counts` → Total gene counts
- `isoforms_per_gene` → Ratio of isoforms to genes

**Debug output:**
```
[INFO] Average isoforms per gene: X.XX
```

### 2. Step 2: Detailed Quality Control

**New QC metrics separated by data type:**

#### Isoform-level metrics:
- `total_isoform_counts` - Total UMI counts from isoforms per cell
- `n_isoforms` - Number of isoforms detected per cell
- `mt_isoform_percent` - Mitochondrial percentage from isoform counts

#### Gene-level metrics:
- `total_gene_counts` - Total UMI counts from genes per cell
- `n_genes` - Number of genes detected per cell
- `mt_gene_percent` - Mitochondrial percentage from gene counts

#### Combined metrics:
- `isoforms_per_gene` - Per-cell ratio of isoforms to genes

**New plots generated (12 total):**

| # | Filename | Description | Data Type |
|---|----------|-------------|-----------|
| 1 | `07_total_isoform_counts_histogram` | Total isoform counts per cell distribution | Isoform |
| 2 | `08_isoforms_detected_histogram` | Number of isoforms detected per cell | Isoform |
| 3 | `09_total_gene_counts_histogram` | Total gene counts per cell distribution | Gene |
| 4 | `10_genes_detected_histogram` | Number of genes detected per cell | Gene |
| 5 | `11_isoform_counts_vs_isoforms_scatter` | Isoform counts vs isoforms detected | Isoform |
| 6 | `12_gene_counts_vs_genes_scatter` | Gene counts vs genes detected | Gene |
| 7 | `13_isoforms_per_gene_histogram` | Isoforms per gene ratio per cell | Combined |
| 8 | `14_isoform_counts_by_timepoint_violin` | Isoform counts by timepoint | Isoform |
| 9 | `15_gene_counts_by_timepoint_violin` | Gene counts by timepoint | Gene |
| 10 | `16_mt_isoform_percent_violin` | MT% from isoform data | Isoform |
| 11 | `17_mt_gene_percent_violin` | MT% from gene data | Gene |
| 12 | `18_isoforms_per_gene_ratio_violin` | Isoforms/gene ratio by timepoint | Combined |

**New output files:**
- `detailed_qc_isoform.csv/rds` - Per-cell isoform QC metrics
- `detailed_qc_gene.csv/rds` - Per-cell gene QC metrics
- `detailed_qc_combined.csv/rds` - Merged isoform + gene metrics

**Debug output:**
```
[DEBUG] Processing QC for sample: org_XY
[INFO] Average isoforms per cell (median): XXXX
[INFO] Average genes per cell (median): XXXX
[INFO] Average isoforms per gene per cell: X.XX
```

## Nomenclature Standards

### Variable Naming Convention
- **Isoform data**: `*_isoform_*`, `n_isoforms`, `isoform_counts`
- **Gene data**: `*_gene_*`, `n_genes`, `gene_counts`
- **Combined metrics**: `*_combined`, `isoforms_per_gene`

### Plot Titles
- Clear indication of data type in every title
- Examples:
  - "Total **Isoform** Counts per Cell"
  - "Total **Gene** Counts per Cell"
  - "**Isoforms** (Transcripts) vs Cells by Sample"
  - "**Genes** vs Cells by Sample"

### File Naming
- Isoform plots: `*_isoform_*`
- Gene plots: `*_gene_*`
- Combined plots: `*_ratio_*` or descriptive names

## Expected Data Dimensions

Based on FLAMES output:
- **Isoforms (transcripts)**: ~385,659 features per sample (ENST IDs)
- **Genes**: ~52,622 features per sample (ENSG IDs)
- **Average ratio**: ~7.3 isoforms per gene

## PDF Reports

### Report 1: `01_Data_Loading_QC.pdf`
- 6 plots in 2×3 grid
- Title: "Data Loading and Initial Quality Control - Isoform and Gene Level"

### Report 2: `02_Quality_Control.pdf`
- 12 plots in 2×6 grid
- Title: "Detailed Quality Control Analysis - Isoform and Gene Level"

## Benefits

1. **Clarity**: Immediate visual distinction between isoform and gene-level data
2. **Comparability**: Side-by-side comparison of gene vs isoform QC metrics
3. **Interpretability**: Clear understanding of alternative splicing complexity (isoforms/gene ratio)
4. **Debugging**: Explicit debug statements for both data types
5. **Reproducibility**: Separate CSV files for isoform and gene QC metrics

## Next Steps

When implementing **Step 3 (Integration)**, ensure:
- Both `RNA` (gene) and `iso` (isoform) assays are clearly labeled in Seurat objects
- QC plots for integrated data maintain this nomenclature
- Debug statements specify which assay is being processed
