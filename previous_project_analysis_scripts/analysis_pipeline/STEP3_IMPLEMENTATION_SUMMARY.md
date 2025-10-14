# Step 3 Implementation Summary - Seurat QC Workflow

## Implementation Complete ✅

**Date:** 2025-10-05
**Author:** Manveer Chauhan

## What Was Implemented

### 1. New File: `functions/id_mapping_utils.R`
**GTF-based hash map for fast ID conversion**

**Main Function:**
- `create_id_mappings_from_flames()` - Creates gene and transcript ID mappings from FLAMES GTF
  - Reads `isoform_annotated.gtf` from FLAMES output
  - Creates gene_map: ENSG → Gene Symbol
  - Creates transcript_map: ENST → GeneName_TxID
  - Returns hash maps (named vectors) for O(1) lookup
  - Caches to RDS for instant reload
  - Merges with reference GTF to fill missing annotations

**Helper Functions:**
- `apply_gene_symbol_conversion()` - Converts gene count matrix rownames (ENSG → Symbol)
- `apply_transcript_label_conversion()` - Converts isoform count matrix rownames (ENST → GeneName_TxID)

**Key Features:**
- ✅ No biomaRt dependencies (uses GTF directly)
- ✅ Hash map speed (O(1) lookup)
- ✅ Automatic caching to RDS
- ✅ Handles Bambu transcripts gracefully
- ✅ Handles duplicate gene symbols (sum counts by default)

---

### 2. New File: `functions/seurat_qc_modules.R`
**Five main QC functions for dual-assay Seurat workflow**

#### Function 1: `initializeSeuratObjs()`
Creates Seurat object with both RNA (gene) and iso (isoform) assays

**Steps:**
1. Convert gene counts: ENSG → Gene Symbol using gene_hash
2. Convert isoform counts: ENST → GeneName_TxID using tx_hash
3. Handle duplicate gene symbols (sum counts)
4. Create Seurat object with RNA assay
5. Add isoform assay (match cells between assays)
6. Calculate MT% for both assays
7. Apply basic filtering: min.cells=5, min.features=500

**Output:** Dual-assay Seurat object with unfiltered data

---

#### Function 2: `runSeuratPreliminaryFiltering()`
Applies timepoint-specific dynamic QC thresholds

**Steps:**
1. Calculate dynamic thresholds: mean ± 1.5 SD (timepoint-specific)
2. Apply fixed MT% threshold: 10%
3. Filter gene assay (RNA) based on thresholds
4. Apply same cell list to isoform assay (iso)
5. Generate before/after QC comparison plots
6. Store QC statistics in object metadata

**Thresholds:**
- nFeature_RNA: [mean - 1.5×SD, mean + 1.5×SD]
- nCount_RNA: [mean - 1.5×SD, mean + 1.5×SD]
- percent.mt: < 10%

**Output:** Filtered Seurat object + QC comparison plots

---

#### Function 3: `NormaliseScaleAndElbow()`
SCTransform normalization and cell cycle scoring

**Steps:**
1. SCTransform normalization on RNA assay
2. Run PCA (test up to 50 PCs)
3. Determine optimal PCs using elbow method
4. Load cell cycle genes from `/data/gpfs/projects/punim0646/manveer/cycle.rda`
5. Score cell cycle phases (G1/S/G2M)
6. **Do NOT regress** cell cycle effects (label only)
7. Generate elbow plot with optimal PC marked

**Output:** Normalized Seurat object with optimal PC selection

---

#### Function 4: `find.optimal.cluster.res()`
Silhouette-based automatic clustering resolution selection

**Steps:**
1. Build SNN graph using optimal PCs
2. Test all resolutions: 0.1, 0.3, 0.5, 0.8, 1.0
3. Calculate silhouette score for each resolution
4. Automatically select resolution with highest silhouette score
5. Re-run clustering with optimal resolution
6. Generate clustree visualization

**Output:** Optimal resolution value + clustree plot

---

#### Function 5: `runUMAP.removeDoublets()`
UMAP dimensionality reduction and DoubletFinder

**Steps:**
1. Run UMAP with optimal PCs
2. Estimate DoubletFinder pK parameter via BCmvn
3. Calculate expected doublets (1.6% doublet rate)
4. Adjust for homotypic doublets
5. Run DoubletFinder
6. Remove doublets (keep only singlets)
7. Re-run UMAP on clean data

**Output:** Final QC'd Seurat object (singlets only)

---

### 3. Updated: `functions/config.R`
**Added Seurat QC parameters**

```r
# Timepoint-specific QC thresholds (fine-tunable)
QC_THRESHOLDS <- list(
  "1M_Org" = list(mt_percent = 10, sd_multiplier = 1.5),
  "3M_Org" = list(mt_percent = 10, sd_multiplier = 1.5),
  "6M_Org" = list(mt_percent = 10, sd_multiplier = 1.5)
)

# Seurat parameters
MIN_CELLS <- 5
MIN_FEATURES <- 500
CLUSTERING_RESOLUTIONS <- c(0.1, 0.3, 0.5, 0.8, 1.0)
DOUBLET_RATE <- 0.016  # 1.6%
CELL_CYCLE_GENES <- "/data/gpfs/projects/punim0646/manveer/cycle.rda"
N_PCS_TEST <- 50
```

**Key Features:**
- ✅ Timepoint-specific thresholds (easily tunable)
- ✅ Fixed 10% MT threshold across all timepoints
- ✅ Dynamic nFeature/nCount thresholds (mean ± 1.5 SD)

---

### 4. Updated: `main.R`
**Added Step 3: Seurat QC workflow**

**Main Step 3 Function: `step_seurat_qc()`**

**Workflow:**
1. Load sample metadata from Step 1
2. Create ID mappings from FLAMES GTF (cached)
3. **For each of 8 samples:**
   - Load raw gene and isoform counts
   - Initialize dual-assay Seurat object
   - Apply QC filtering (timepoint-specific)
   - Normalize and score cell cycle
   - Optimize clustering resolution (silhouette)
   - Run UMAP and remove doublets
   - Save individual Seurat RDS file
   - Generate per-sample QC PDF report
4. Combine QC summary table across all samples
5. Generate summary comparison PDF report

**Helper Functions:**
- `generate_sample_qc_report()` - Creates per-sample PDF with 5 pages:
  - Page 1: QC filtering comparison (before/after)
  - Page 2: Elbow plot (optimal PCs)
  - Page 3: Clustree (resolution selection)
  - Page 4: Final UMAP (clusters)
  - Page 5: Cell cycle UMAP

- `generate_qc_summary_report()` - Creates summary PDF with 4 pages:
  - Page 1: Cell retention bar plots (by timepoint)
  - Page 2: QC metrics boxplots (nFeature, MT%)
  - Page 3: Optimal PCs and cluster counts
  - Page 4: All samples UMAP grid (4×2 layout)

---

## Expected Output Structure

```
analysis_pipeline/output/
├── data/
│   ├── id_mappings_cache.rds              # Cached GTF mappings (reusable)
│   └── seurat_objects/
│       ├── org_1A_qc.rds                  # 8 individual Seurat objects
│       ├── org_1B_qc.rds
│       ├── org_3A_qc.rds
│       ├── org_3B_qc.rds
│       ├── org_3C_qc.rds
│       ├── org_6A_qc.rds
│       ├── org_6B_qc.rds
│       ├── org_6C_qc.rds
│       └── qc_summary_table.csv/rds       # Summary statistics table
│
├── figures/qc/
│   ├── org_1A_qc_before_after.png         # Per-sample QC comparison
│   ├── org_1A_elbow.png                   # PCA elbow plots
│   ├── org_1A_clustree.png                # Clustering resolution trees
│   ├── org_1A_umap_final.png              # Final UMAP
│   └── ... (repeated for all 8 samples)
│
└── reports/
    ├── org_1A_QC_Report.pdf               # Per-sample detailed reports (8 files)
    ├── org_1B_QC_Report.pdf
    ├── ... (6 more)
    └── 03_Seurat_QC_Summary.pdf           # Combined summary report
```

---

## Seurat Object Structure

**Each `*_qc.rds` file contains:**

### Assays:
- `RNA` - Gene-level expression (gene symbols as rownames)
- `SCT` - SCTransform-normalized gene expression
- `iso` - Isoform-level expression (GeneName_TxID as rownames)

### Reductions:
- `pca` - Principal component analysis
- `umap` - UMAP coordinates

### Metadata Columns:
- `sample_id` - Sample identifier (e.g., "org_1A")
- `timepoint` - Timepoint (e.g., "1M_Org")
- `nFeature_RNA` - Number of genes detected
- `nCount_RNA` - Total gene counts
- `percent.mt` - Mitochondrial percentage (gene-level)
- `percent.mt.iso` - Mitochondrial percentage (isoform-level)
- `Phase` - Cell cycle phase (G1/S/G2M)
- `S.Score` - S phase score
- `G2M.Score` - G2M phase score
- `seurat_clusters` - Cluster assignments (optimal resolution)
- `doublet_classification` - Singlet/Doublet (all are Singlet in final object)
- `qc_status` - "filtered"

### Misc (stored metadata):
- `optimal_pcs` - Optimal number of PCs selected
- `optimal_resolution` - Optimal clustering resolution
- `qc_thresholds` - Applied QC thresholds
- `qc_stats` - QC filtering statistics
- `doublet_stats` - Doublet removal statistics
- `qc_plots` - Before/after QC comparison plots
- `elbow_plot` - PCA elbow plot
- `clustree_plot` - Clustering resolution tree
- `silhouette_scores` - Silhouette scores for all tested resolutions

---

## QC Summary Table Columns

**File:** `output/data/seurat_objects/qc_summary_table.csv`

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier |
| `timepoint` | Timepoint (1M_Org, 3M_Org, 6M_Org) |
| `n_cells_initial` | Cells before filtering |
| `n_cells_filtered` | Cells after QC filtering |
| `pct_retained_filtering` | % cells retained after filtering |
| `n_doublets` | Number of doublets detected |
| `pct_doublets` | % doublets detected |
| `n_cells_final` | Final cell count (singlets only) |
| `n_genes` | Number of genes in RNA assay |
| `n_isoforms` | Number of isoforms in iso assay |
| `optimal_pcs` | Optimal PCs selected |
| `optimal_resolution` | Optimal clustering resolution |
| `n_clusters` | Number of clusters |
| `median_nFeature` | Median features per cell |
| `median_nCount` | Median counts per cell |
| `median_pct_mt` | Median MT% |

---

## Key Implementation Features

### 1. GTF-Based ID Mapping (No biomaRt)
- ✅ Fast: O(1) hash map lookup
- ✅ Reliable: No external API dependencies
- ✅ Cached: Instant reload on re-runs
- ✅ Robust: Handles Bambu transcripts and missing annotations

### 2. Dynamic Timepoint-Specific Thresholds
- ✅ Calculated per timepoint: mean ± 1.5 SD
- ✅ Easily tunable in `config.R`
- ✅ Fixed MT% threshold: 10%

### 3. Silhouette-Based Clustering
- ✅ Automatic optimal resolution selection
- ✅ Tests 5 resolutions: 0.1, 0.3, 0.5, 0.8, 1.0
- ✅ User can override with clustree visualization

### 4. Comprehensive QC Reporting
- ✅ Per-sample detailed PDFs (5 pages each)
- ✅ Summary comparison PDF (4 pages)
- ✅ Individual PNG figures for each plot
- ✅ CSV summary table for easy inspection

### 5. Dual-Assay Support
- ✅ RNA assay: Gene symbols (standard Seurat workflows)
- ✅ iso assay: GeneName_TxID format (preserves isoform identity)
- ✅ Matched cells between assays
- ✅ Separate MT% calculations

---

## How to Run

### Option 1: Run Full Pipeline (Steps 1-3)
```bash
cd /data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/previous_project_analysis_scripts

R
source("analysis_pipeline/main.R")
# Pipeline will auto-run Steps 1, 2, and 3
```

### Option 2: Run Only Step 3 (if Steps 1-2 already complete)
```r
source("analysis_pipeline/main.R")

# Run only Step 3
step3_result <- step_seurat_qc()
```

### Option 3: Load Individual Seurat Objects
```r
library(Seurat)

# Load a specific sample
seurat_obj <- readRDS("analysis_pipeline/output/data/seurat_objects/org_1A_qc.rds")

# Inspect object
seurat_obj  # Shows assays, reductions, metadata

# View QC metrics
seurat_obj@misc$qc_stats
seurat_obj@misc$optimal_pcs
seurat_obj@misc$optimal_resolution

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase")
```

---

## Next Steps (Future Implementation)

### Step 4: Integration (Harmony-based)
- Integrate samples per timepoint using Harmony
- Create timepoint-specific integrated objects
- Compare within-timepoint vs across-timepoint integration strategies

### Step 5: Cell Type Annotation
- Cell type identification using reference-based annotation
- Manual curation of cluster identities
- Generate cell type marker plots

### Step 6: Isoform-Level Analysis
- Differential transcript usage (DTU) analysis
- Trajectory analysis (pseudotime)
- Isoform switching analysis

---

## Dependencies

**R Packages Required:**
```r
# Core single-cell
library(Seurat)
library(DoubletFinder)
library(clustree)

# Data manipulation
library(dplyr)
library(Matrix)
library(data.table)

# Visualization
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(scales)

# Genomics
library(rtracklayer)

# Utilities
library(stringr)
```

**External Files Required:**
- FLAMES GTF: `{FLAMES_DIR}/isoform_annotated.gtf`
- Reference GTF: `/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf`
- Cell cycle genes: `/data/gpfs/projects/punim0646/manveer/cycle.rda`
- Gene count CSVs: `{FLAMES_DIR}/*_gene_count.csv`
- Isoform count MTX: `{FLAMES_DIR}/*.count.mtx`

---

## Testing Recommendations

1. **Test on single sample first:**
   ```r
   # Load only org_1A for testing
   sample_metadata <- sample_metadata[1, ]  # First sample only
   step3_result <- step_seurat_qc()
   ```

2. **Check ID mapping cache:**
   ```r
   id_maps <- readRDS("analysis_pipeline/output/data/id_mappings_cache.rds")
   head(id_maps$gene_map)
   head(id_maps$transcript_map)
   ```

3. **Inspect QC summary:**
   ```r
   qc_summary <- read.csv("analysis_pipeline/output/data/seurat_objects/qc_summary_table.csv")
   print(qc_summary)
   ```

4. **View generated PDFs:**
   - Check `output/reports/org_1A_QC_Report.pdf` for detailed QC
   - Check `output/reports/03_Seurat_QC_Summary.pdf` for overview

---

## Known Considerations

1. **Memory usage:** Each sample requires ~5-10GB RAM (385K isoforms × cells)
2. **Runtime:** ~10-15 minutes per sample (total ~2 hours for 8 samples)
3. **DoubletFinder:** May take 5-10 minutes per sample for pK estimation
4. **Cell cycle genes:** Ensure `cycle.rda` file exists and is accessible
5. **Bambu transcripts:** Novel transcripts retain ENSG_ENST format (by design)

---

## Summary

✅ **Implementation Complete**
- 3 new/updated files
- 5 main QC functions
- 2 helper report generation functions
- GTF-based ID mapping system
- Dual-assay Seurat workflow
- Comprehensive QC reporting

**Ready for testing on full 8-sample dataset.**
