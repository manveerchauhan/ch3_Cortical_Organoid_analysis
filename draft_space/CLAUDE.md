# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the draft workspace for Chapter 3 of the neurodevelopmental organoid analysis project. It implements a complete single-cell RNA-seq analysis pipeline from FLAMES v2.2.0 isoform quantification outputs through to validated, annotated cell type clusters. The workspace handles 8 organoid samples across three developmental timepoints (1, 3, and 6 months) with both gene-level and isoform-level quantification.

## Pipeline Architecture

The analysis follows a modular 6-stage pipeline:

1. **Data Import** → Raw FLAMES outputs to Seurat objects
2. **QC & Processing** → Cell/gene filtering, normalization, doublet removal
3. **Integration** → Harmony-based batch correction by timepoint
4. **Validation** → Compare new clusters to ground truth annotations
5. **Marker Analysis** → Validate cell type marker expression patterns
6. **Annotation** → Automated cell type assignment via scType

## Core Analysis Scripts

### Stage 1: Data Import

#### `ch3_qc_flames_v220.R`
Primary pipeline for processing FLAMES outputs into Seurat-compatible count matrices:
- Creates isoform-to-gene symbol dictionary from FLAMES GTF
- Converts Oarfish quantification outputs (`.count.mtx`, `.barcodes.txt`, `.features.txt`) into CSV count matrices
- Generates sample metadata tibble with embedded count matrices
- Creates separate gene-level and isoform-level Seurat objects for each sample
- **Memory-optimized**: Processes samples sequentially to avoid memory overload

#### `ch3_qc_flames_v220_p2.R`
Lightweight version for creating Seurat objects from pre-processed count matrices:
- Loads pre-saved `org_flames_v2.2.0_experiment_data.rds`
- Uses configurable filtering parameters (`MIN.CELL`, `MIN.FEATURES`)
- Processes samples one at a time with immediate RDS saving

#### `ch3_org_working_analysis.R`
Converts gene IDs to gene symbols and creates final Seurat objects:
- Uses combined `gene_id_gene_symbol` format for uniqueness (e.g., `ENSG00000002587.10_CD99`)
- Handles NA values in FLAMES gene count CSVs
- Adds timepoint metadata

### Stage 2: QC & Processing

#### `qc_functions.R`
Comprehensive QC function library implementing 6-step processing:
1. `normalize_and_find_features()` - Normalization + variable feature selection
2. `cell_cycle_and_scale()` - Cell cycle scoring + data scaling
3. `quantitative_elbow()` - Automated optimal PC determination (with +5 PC adjustment)
4. `optimize_clustering()` - Silhouette-based resolution optimization
5. `umap_and_doublets()` - UMAP generation + DoubletFinder
6. Helper: `convert_cc_genes_to_seurat_format()` - Converts cell cycle gene symbols to ENSG-SYMBOL format

**Key features**:
- Quantitative elbow uses two metrics: cumulative variance >90% & individual <5%, and consecutive change <0.1%
- Silhouette analysis tests resolutions 0.1-1.2 (step 0.1)
- DoubletFinder with homotypic doublet correction
- All functions return list with `seurat_obj` and `plots` components

### Stage 3: Integration

#### `ch3_org_integration_harmony.R`
End-to-end Harmony integration pipeline for timepoint-grouped samples:
- Merges samples by timepoint (1M_Org: 2 samples, 3M_Org: 3 samples, 6M_Org: 3 samples)
- Runs pre-integration processing with quantitative elbow
- Applies Harmony batch correction on `sample_id`
- Creates before/after integration UMAPs
- Optimizes clustering with silhouette analysis
- **Configurable resolution override** via `clustering_params` list (NA = use calculated optimal)
- Generates comprehensive PDF reports with elbow plots, integration comparison, clustree, and QC metrics
- Saves integrated objects to `./output_files/integrated_objects/[timepoint]_integrated_harmony.rds`

**Critical configuration**:
```r
clustering_params <- list(
  "1M_Org" = 0.3,   # Custom override
  "3M_Org" = NA,    # Use calculated optimal
  "6M_Org" = 0.5    # Custom override
)
```

### Stage 4: Validation

#### `ch3_cluster_validation.R`
Validates new integrated clusters against ground truth manual annotations:
- Compares to published organoid objects: `[one/three/six]_month_seurat.intergrated_harm.isofrom.rds`
- **Barcode normalization**: Strips leading underscores from new objects, uses `Barcode` metadata column from ground truth
- Creates confusion matrices (raw counts + row-normalized percentages)
- Generates distribution summaries showing top 3 cluster assignments per annotation
- Produces heatmap visualizations with cluster counts and percentages
- Outputs CSVs and PDFs to `./output_files/cluster_validation/`

### Stage 5: Marker Validation

#### `ch3_marker_validation.R`
Validates 14 cell type marker gene sets across integrated timepoints:
- **Marker sets**: oRG, vRG, SST interneurons, parvalbumin interneurons, intermediate progenitors, stem cells, glutamatergic neurons, microglia, oligodendrocytes, pan-neuronal, astrocytes, Cajal-Retzius
- Resolves gene symbols to ENSG-SYMBOL format via regex matching
- Calculates module scores using `AddModuleScore()`
- Identifies top 5 most expressed genes per marker set
- Creates comprehensive PDFs with:
  - Summary page (dataset info + missing genes)
  - Module score UMAP per marker set
  - Individual gene expression UMAPs (top 5 genes)
- Outputs to `./output_files/marker_validation/[timepoint]_marker_validation.pdf`

**Important**: Uses `JoinLayers()` for Seurat v5 compatibility before visualization

### Stage 6: Annotation

#### `cellannotation-script3.R`
Automated cell type annotation using scType:
- Downloads scType gene sets from GitHub or uses custom Excel references
- Calculates scType scores on `scale.data` layer
- Assigns cell types to clusters (low-confidence → "Unknown" if score < ncells/4)
- Creates UMAP comparisons (labeled vs unsupervised clusters)
- Generates confidence pie charts for each cluster
- Exports PDFs with labels and confidence visualizations

## Sample Information

```r
sample_ids <- c("org_1A", "org_1B", "org_3A", "org_3B", "org_3C", "org_6A", "org_6B", "org_6C")
timepoints <- c("1M_Org", "1M_Org", "3M_Org", "3M_Org", "3M_Org", "6M_Org", "6M_Org", "6M_Org")
```

## Key Data Locations

### Input Data
- **FLAMES output directory**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models/Cortical_Org_Diff/flame_outs_standard_ref`
  - Contains per-sample Oarfish outputs: `[sample_id].count.mtx`, `[sample_id].barcodes.txt`, `[sample_id].features.txt`
  - Gene-level counts: `[sample_id]_gene_count.csv`
  - FLAMES GTF: `isoform_annotated.gtf`
- **Reference GTF**: `/data/gpfs/projects/punim2251/Organoid_NDR_FLAMES_ref/modified_gencode_v47_primary_annotation.gtf`

### Output Structure
```
output_files/
├── ref_files/                # isoform_gene_dict.csv - transcript ID to gene symbol mapping
├── counts/                   # gene_symbol_[sample_id]_counts.csv - formatted isoform count matrices
├── seu_objects/              # Individual Seurat objects
│   ├── [sample_id]_gene_seurat.rds
│   ├── [sample_id]_isoform_seurat.rds
│   └── [sample_id]_QCed_final.rds       # Post-QC objects (input for integration)
├── integrated_objects/       # Timepoint-integrated Seurat objects
│   ├── 1M_Org_integrated_harmony.rds
│   ├── 3M_Org_integrated_harmony.rds
│   └── 6M_Org_integrated_harmony.rds
├── integration_QC/           # Integration pipeline PDFs
│   └── [timepoint]_integration_report.pdf
├── cluster_validation/       # Validation vs ground truth
│   ├── [timepoint]_overlap_stats.csv
│   ├── [timepoint]_confusion_matrix.csv
│   ├── [timepoint]_confusion_matrix_pct.csv
│   └── [timepoint]_confusion_plots.pdf
├── marker_validation/        # Cell type marker expression
│   └── [timepoint]_marker_validation.pdf
└── QC/                       # Individual sample QC reports
```

### Intermediate Data
- `org_flames_v2.2.0_experiment_data.rds`: Sample metadata tibble with embedded count matrices (509 MB)
- Speeds up re-running analysis without re-reading CSVs

## Complete Pipeline Workflow

### Stage 1: Create Initial Seurat Objects

**Option A: From scratch (ch3_qc_flames_v220.R)**
```r
# 1. Create isoform-gene dictionary (takes 5-10 minutes)
isoform_gene_dict <- make_isoform_gene_symbol_dict(FLAMES_gtf_file, reference_gtf_file, output_file)

# 2. Process Oarfish outputs to formatted count matrices
process_oarfish_files_to_counts_matrix(flames_output_folder, sample_name, output_dir)

# 3. Build sample metadata tibble and save
sample_metadata <- tibble(...)
saveRDS(sample_metadata, file = "org_flames_v2.2.0_experiment_data.rds")

# 4. Create Seurat objects (sequential, memory-safe)
for(i in 1:nrow(sample_metadata)) {
  process_single_sample(i, sample_metadata)
}
```

**Option B: From pre-processed data (ch3_org_working_analysis.R - RECOMMENDED)**
```r
# Converts FLAMES gene counts to gene symbols and creates Seurat objects
# Handles NA values and duplicate gene names automatically
# Outputs: [sample_id]_gene_seurat_with_symbols.rds
source("ch3_org_working_analysis.R")
```

### Stage 2: QC Processing

**Apply comprehensive QC to each sample using qc_functions.R:**
```r
source("qc_functions.R")
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")  # Load s_genes, g2m_genes
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")

# For each sample:
seu <- readRDS("./output_files/seu_objects/org_1A_gene_seurat_with_symbols.rds")

# Step 1: Normalization
result1 <- normalize_and_find_features(seu, "org_1A")
seu <- result1$seurat_obj

# Step 2: Cell cycle scoring
result2 <- cell_cycle_and_scale(seu, "org_1A")
seu <- result2$seurat_obj

# Step 3: Determine optimal PCs
result3 <- quantitative_elbow(seu, "org_1A", max_pcs = 50, pc_adjustment = 5)
optimal_pcs <- result3$optimal_pcs_used

# Step 4: Optimize clustering
result4 <- optimize_clustering(seu, "org_1A", optimal_pcs)
seu <- result4$seurat_obj

# Step 5: UMAP and doublet removal
result5 <- umap_and_doublets(seu, "org_1A", optimal_pcs, doublet_rate = 0.016)
seu <- result5$seurat_obj

# Save QC'd object
saveRDS(seu, "./output_files/seu_objects/org_1A_QCed_final.rds")
```

### Stage 3: Integration

**Run Harmony integration by timepoint:**
```r
source("ch3_org_integration_harmony.R")

# The script will:
# 1. Load all QCed samples (*_QCed_final.rds)
# 2. Merge by timepoint (1M_Org, 3M_Org, 6M_Org)
# 3. Run Harmony integration on sample_id
# 4. Optimize clustering (with optional resolution override)
# 5. Generate comprehensive PDF reports
# 6. Save integrated objects

# Outputs:
# - ./output_files/integrated_objects/[timepoint]_integrated_harmony.rds
# - ./output_files/integration_QC/[timepoint]_integration_report.pdf
```

### Stage 4: Validation

**Validate clusters against ground truth:**
```r
source("ch3_cluster_validation.R")

# Compares integrated objects to published organoid annotations
# Generates confusion matrices, overlap statistics, heatmaps
# Outputs to ./output_files/cluster_validation/
```

### Stage 5: Marker Validation

**Validate cell type marker expression:**
```r
source("ch3_marker_validation.R")

# Tests 14 marker gene sets across all integrated timepoints
# Calculates module scores and top gene expression
# Generates multi-page PDFs with module scores and individual gene UMAPs
# Outputs to ./output_files/marker_validation/
```

### Stage 6: Annotation (Optional)

**Automated cell type assignment:**
```r
source("cellannotation-script3.R")

# Uses scType for automated annotation
# Generates pie charts showing confidence scores
# Exports labeled UMAPs and confidence visualizations
```

## Key Function Reference

### Data Import Functions (ch3_qc_flames_v220.R)

#### `make_isoform_gene_symbol_dict(FLAMES_gtf, reference_gtf, output_file)`
Merges FLAMES transcript IDs with GENCODE gene symbols. Returns data frame with:
- `transcript_id`: FLAMES transcript ID
- `gene_id`: ENSEMBL gene ID
- `gene_symbol`: Human-readable gene name

#### `process_oarfish_files_to_counts_matrix(flames_output_folder, sample_name, output_dir)`
Converts Oarfish Matrix Market format to CSV with gene symbols:
- Reads `.mtx`, `.barcodes.txt`, `.features.txt`
- Transposes matrix (cells = columns, features = rows)
- Appends gene symbols: `[transcript_id]_[gene_symbol]`
- Saves as `gene_symbol_[sample_name]_counts.csv`

#### `create_seurat_from_counts(count_matrix_df, sample_id, assay_type, MIN.CELL, MIN.FEATURES)`
Creates Seurat objects with pre-filtering:
- Handles duplicate feature names via `make.unique()`
- Pre-filters genes and cells before sparse matrix conversion (saves memory)
- Returns Seurat object with `sample_id` in metadata

#### `process_single_sample(sample_idx, metadata)`
Memory-safe sequential processing:
- Creates both gene-level and isoform-level Seurat objects
- Saves immediately to disk as `.rds`
- Clears from memory with `rm()` and `gc()`

### QC Functions (qc_functions.R)

#### `normalize_and_find_features(seurat_obj, sample_id)`
Returns: `list(seurat_obj, plots = list(var_feature_plot), top10)`
- `NormalizeData()` + `FindVariableFeatures()` (2000 features)
- Creates labeled variable feature plot

#### `cell_cycle_and_scale(seurat_obj, sample_id, gene_dict_path)`
Returns: `list(seurat_obj, plots = list(cc_pc1_pc2, cc_pc2_pc3))`
- Converts cell cycle genes to ENSG-SYMBOL format
- `CellCycleScoring()` + `ScaleData()` + `RunPCA()`
- Creates cell cycle phase plots on PC1/PC2 and PC2/PC3

#### `quantitative_elbow(seurat_obj, sample_id, max_pcs = 50, pc_adjustment = 5)`
Returns: `list(optimal_pcs, optimal_pcs_used, plots = list(elbow_plot, quant_plot), metrics)`
- Calculates two metrics: cumulative >90% & individual <5%, consecutive change <0.1%
- Returns both calculated optimal and adjusted (+5) values
- Creates annotated elbow plot with both values marked

#### `optimize_clustering(seurat_obj, sample_id, optimal_pcs)`
Returns: `list(seurat_obj, optimal_res, sil_results, plots = list(clustree))`
- `FindNeighbors()` on specified PCs
- Silhouette analysis on resolutions 0.1-1.2 (step 0.1)
- Returns top silhouette resolution
- Generates clustree visualization

#### `umap_and_doublets(seurat_obj, sample_id, optimal_pcs, doublet_rate)`
Returns: `list(seurat_obj, doublet_stats, plots = list(doublet_umap, final_umap, nfeature_umap, ncount_umap, cellcycle_umap))`
- `RunUMAP()` + DoubletFinder parameter sweep
- Removes doublets via `subset()`
- Returns singlet-only Seurat object with QC plots

### Integration Functions (ch3_org_integration_harmony.R)

#### `merge_samples_by_timepoint(sample_ids, seurat_list, timepoint_name)`
Merges multiple Seurat objects with `sample_id` and `timepoint` metadata

#### `pre_integration_processing(seurat_obj, timepoint_name)`
Returns: `list(seurat_obj, optimal_pcs, elbow_plots)`
- Normalization + variable features + scaling + PCA
- Quantitative elbow analysis

#### `harmony_integration(seurat_obj, timepoint_name, optimal_pcs)`
Returns: `list(seurat_obj, plots = list(before, after))`
- Creates unintegrated UMAP for comparison
- `IntegrateLayers()` with HarmonyIntegration
- Creates harmony-corrected UMAP

#### `post_integration_clustering(seurat_obj, timepoint_name, optimal_pcs, custom_resolution = NA)`
Returns: `list(seurat_obj, optimal_res, final_res, resolution_method, sil_results, plots)`
- `FindNeighbors()` on harmony reduction
- Silhouette optimization
- Allows custom resolution override (NA = use calculated)
- Generates clustree and final cluster UMAP

### Validation Functions (ch3_cluster_validation.R)

#### `analyze_barcode_overlap(new_obj, gt_obj, timepoint_name)`
Returns: `list(overlap_cells_new, overlap_cells_gt, overlap_barcodes_clean, stats)`
- Strips leading underscores from new object barcodes
- Uses `Barcode` metadata column from ground truth
- Returns mapping between normalized and original cell names

#### `create_confusion_matrix(new_obj, gt_obj, overlap_result, timepoint_name)`
Returns: `list(confusion_mat, confusion_df, confusion_pct, new_clusters, gt_annotations)`
- Creates cluster × annotation contingency table
- Row-normalizes to percentages

#### `visualize_confusion_matrix(confusion_result, timepoint_name)`
Returns: `list(count_plot, pct_plot)`
- Heatmaps with cluster counts and percentages

### Marker Validation Functions (ch3_marker_validation.R)

#### `resolve_gene_name(gene_symbol, seurat_obj)`
Finds full ENSG-SYMBOL feature name via regex pattern matching

#### `validate_gene_set(gene_symbols, seurat_obj, set_name)`
Returns: `data.frame(gene_symbol, full_name, found)`
- Resolves all genes in a marker set
- Reports matching statistics

#### `get_top_genes(gene_full_names, seurat_obj, n = 5)`
Returns top n genes by mean expression across all cells

#### `create_marker_plot(seurat_obj, marker_set_name, validation_results, module_score_name, timepoint_name)`
Returns: `list(title_info, module_plot, cluster_plot, gene_plots)`
- Creates module score UMAP
- Creates cluster UMAP
- Creates individual gene UMAPs (top 5 genes)

## Important Implementation Notes

### Memory Management
- **Do NOT** load all 8 samples into memory simultaneously
- **Do NOT** create all Seurat objects in parallel
- **Always** use `process_single_sample()` sequential approach
- Expected Seurat object sizes: ~50-200 MB per sample

### Rowname Handling
- FLAMES may produce duplicate transcript IDs
- Use `make.unique()` to append `.1`, `.2`, etc. to duplicates
- **Gene-level naming**: `ENSG00000002587.10_CD99` format (gene_id_gene_symbol)
  - Multiple ENSG IDs can map to same gene symbol (paralogs, sex chromosomes)
  - Combined format ensures uniqueness while preserving interpretability
- **Isoform-level naming**: `ENST00000123456_GENE1` format (transcript_id_gene_symbol)

### Assay Naming Convention
- Gene-level: `assay = "RNA"` (standard Seurat convention)
- Isoform-level: `assay = "isoform"` or `assay = "iso"` (varies by script)
- **Note**: Parent directory uses `assay = "iso"`, this workspace uses `"isoform"`

### Filtering Philosophy
- Gene-level: More stringent (≥3 cells, ≥20 features) for standard QC
- Isoform-level: Permissive (≥1 cell, ≥1 feature) to retain rare isoforms

## Required R Packages

### Core Pipeline
- `Seurat` (v5): Single-cell object framework with layer support
- `Matrix`: Sparse matrix operations
- `tidyverse`: Data manipulation (tibble, dplyr, ggplot2, stringr)

### QC & Processing
- `DoubletFinder`: Doublet detection and removal
- `clustree`: Clustering resolution visualization

### Integration
- `harmony`: Batch correction via HarmonyIntegration

### Visualization
- `cowplot`, `patchwork`, `gridExtra`, `grid`: Multi-plot layouts

### Annotation
- `openxlsx`: Reading scType marker databases
- `HGNChelper`: Gene symbol validation

### Optional (Data Import)
- `rtracklayer`: GTF import
- `GenomicFeatures`: GTF manipulation

## Running on HPC (SLURM)

### `org_create_seu_objs.slurm`
SLURM job script for Seurat object creation:
```bash
#!/bin/bash
#SBATCH --partition="sapphire"
#SBATCH --account="punim2251"
#SBATCH --cpus-per-task=8
#SBATCH --mem=300gb
#SBATCH --time=0-03:30:00

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.4.0

Rscript ch3_qc_flames_v220_p2.R
```

**Resource recommendations**:
- **Data import**: 300 GB RAM, 3.5 hours (processes all 8 samples sequentially)
- **QC per sample**: 64-128 GB RAM, 2-4 hours (DoubletFinder is memory-intensive)
- **Integration**: 128-256 GB RAM, 4-6 hours (depends on merged cell count)
- **Validation/Marker analysis**: 64 GB RAM, 1-2 hours

**Module loading on Spartan**:
```bash
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.4.0
```

## Design Patterns and Best Practices

### Function Return Pattern
All major functions follow a consistent return pattern:
```r
list(
  seurat_obj = modified_object,
  plots = list(plot1, plot2, ...),
  [stage-specific results]
)
```
This allows chaining operations: `result <- func1(seu); seu <- result$seurat_obj`

### Metadata Conventions
- `sample_id`: Original sample name (org_1A, org_1B, etc.)
- `timepoint`: Developmental stage (1M_Org, 3M_Org, 6M_Org)
- `Phase`: Cell cycle phase (G1, S, G2M) from `CellCycleScoring()`
- `seurat_clusters`: Final cluster assignments
- `cluster_annotations`: Manual annotations (ground truth objects only)
- `DF.classifications`: DoubletFinder results (Singlet/Doublet)
- `Barcode`: Clean barcode without sample prefix (ground truth objects only)

### Reduction Naming
- `pca`: Standard PCA reduction (pre-integration)
- `harmony`: Harmony-corrected embedding (post-integration)
- `umap.unintegrated`: Pre-integration UMAP (for comparison)
- `umap.harmony`: Post-integration UMAP (primary visualization)
- `umap`: Standard UMAP (individual samples)

### Gene Naming Standards
- **Gene-level (RNA assay)**: `ENSG00000002587.10_CD99` (gene_id_gene_symbol)
  - Ensures uniqueness for paralogs/sex chromosomes
  - Preserves interpretability with gene symbol
- **Isoform-level (isoform/iso assay)**: `ENST00000123456_GENE1` (transcript_id_gene_symbol)
  - FLAMES may generate duplicate transcript IDs → use `make.unique()`
- **Cell cycle genes**: Must be converted to ENSG-SYMBOL format via `convert_cc_genes_to_seurat_format()`

### Resolution Configuration Pattern
```r
# Define custom overrides (NA = use calculated optimal)
clustering_params <- list(
  "1M_Org" = 0.3,   # Force 0.3
  "3M_Org" = NA,    # Use silhouette optimal
  "6M_Org" = 0.5    # Force 0.5
)
```

### Seurat v5 Compatibility
- Use `JoinLayers()` before plotting with FeaturePlot/DimPlot if working with split/integrated objects
- Access data with `layer` parameter: `GetAssayData(obj, assay = "RNA", layer = "data")`
- `IntegrateLayers()` replaces legacy `RunHarmony()` workflow

### Validation Architecture
The validation stage compares new analysis to ground truth via:
1. **Barcode normalization** → Handle different cell naming conventions
2. **Overlap analysis** → Identify shared cells between objects
3. **Confusion matrices** → Quantify cluster-annotation correspondence
4. **Visualization** → Heatmaps showing assignment patterns

### Marker Validation Strategy
1. **Define biological marker sets** → 14 cell type gene lists
2. **Resolve gene names** → Convert symbols to ENSG-SYMBOL format
3. **Calculate module scores** → `AddModuleScore()` for gene set expression
4. **Identify top genes** → Rank by mean expression
5. **Visualize patterns** → Module score + individual gene UMAPs

## Common Issues and Solutions

### Issue: "duplicate rownames not allowed" when converting gene IDs to symbols
**Problem**: Multiple ENSG IDs map to the same gene symbol (e.g., CD99, ASMTL on X/Y chromosomes, pseudogenes)

**Solution**: Use combined `gene_id_gene_symbol` format instead of gene symbol alone
```r
# In convert_gene_counts_to_symbols():
gene_counts_with_symbols <- gene_counts_with_symbols %>%
  mutate(gene_name = paste0(gene_id, "_", gene_symbol))
rownames(gene_counts_with_symbols) <- gene_counts_with_symbols$gene_name
```
This creates rownames like `ENSG00000002587.10_CD99`, ensuring uniqueness while preserving interpretability.

### Issue: "duplicate rownames not allowed" in isoform data
**Solution**: Use `make.unique()` in `create_seurat_from_counts()` (already implemented)

### Issue: Out of memory errors
**Solution**: Use `ch3_qc_flames_v220_p2.R` sequential processing instead of parallel approaches

### Issue: Seurat objects have wrong assay names
**Solution**: Check `assay_type` parameter - use `"RNA"` for genes, `"isoform"` for transcripts

### Issue: Count matrices have unnamed first column
**Solution**: CSV export includes rownames as first column - handled by `count_matrix_df[[1]]` indexing

### Issue: Gene count matrices full of NA values
**Problem**: FLAMES gene count CSV files use empty cells (NAs) to represent zero counts, causing Seurat object creation to fail and QC metrics to be all NA

**Solution**: Convert NAs to zeros immediately after reading CSV
```r
gene_counts <- read.csv(gene_count_file, row.names = 1, check.names = FALSE)
gene_counts[is.na(gene_counts)] <- 0  # Convert empty cells to zeros
```
This is standard practice for sparse matrix CSV formats where zeros are omitted to save space.

## Typical Analysis Scenarios

### Starting a new analysis from FLAMES outputs
1. Run `ch3_org_working_analysis.R` to create Seurat objects with gene symbols
2. Create QC script for each sample using functions from `qc_functions.R`
3. Run `ch3_org_integration_harmony.R` to integrate by timepoint
4. Run `ch3_cluster_validation.R` and `ch3_marker_validation.R` to assess quality

### Re-running integration with different resolution
1. Edit `clustering_params` in `ch3_org_integration_harmony.R`
2. Source the script to regenerate integrated objects
3. Re-run validation scripts to compare results

### Adding new samples to existing analysis
1. Process new samples through Stage 1 and 2 (data import + QC)
2. Update `sample_groups` in `ch3_org_integration_harmony.R`
3. Re-run integration pipeline
4. Compare new integrated objects to previous versions

### Investigating specific marker genes
1. Load integrated object: `seu <- readRDS("./output_files/integrated_objects/[timepoint]_integrated_harmony.rds")`
2. Use `resolve_gene_name("GENE_SYMBOL", seu)` to find full feature name
3. Use `FeaturePlot()` or `DotPlot()` for visualization
4. Consider adding to marker gene sets in `ch3_marker_validation.R` for systematic tracking

### Troubleshooting cluster assignments
1. Check confusion matrices in `./output_files/cluster_validation/`
2. Review module scores in marker validation PDFs
3. Try different clustering resolutions via `clustering_params` override
4. Examine silhouette analysis results in integration reports

## Related Resources

- **Parent project documentation**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/CLAUDE.md`
- **Cell cycle markers**: `/data/gpfs/projects/punim0646/manveer/cycle.rda` (requires `s_genes`, `g2m_genes`)
- **Silhouette functions**: `/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R` (provides `optimize_silhouette()`)
- **FLAMES documentation**: [GitHub - LongGF/FLAMES](https://github.com/LongGF/FLAMES)
- **Ground truth objects**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/[one/three/six]_month_seurat.intergrated_harm.isofrom.rds`

## Development Environment

- **HPC system**: Spartan (University of Melbourne)
- **Storage**: GPFS shared filesystem (`punim2251`, `punim0646`)
- **R version**: 4.4.0 (GCC/11.3.0, OpenMPI/4.1.4)
- **Typical memory**: 64-300 GB depending on pipeline stage
- **Working directory**: `/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space`
