# Ready for Testing - QC Pipeline Updates

## Summary of Changes

All QC plots and debug statements now clearly distinguish between **isoform-level** and **gene-level** data with proper nomenclature and statistics.

## Files Modified

1. **`main.R`** - Steps 1 and 2 updated with isoform/gene separation
2. **`functions/plotting_utils.R`** - Updated QC plotting functions to include statistics
3. **`CLAUDE.md`** - Created comprehensive documentation
4. **`QC_UPDATE_SUMMARY.md`** - Detailed change log

## What to Expect When Running

### Console Output (DEBUG/INFO statements)

```r
[INFO] Loading configuration...
[INFO] Project: Cortical_Organoid_Differentiation
[INFO] FLAMES directory: /data/gpfs/projects/punim2251/Neurodevelopmental_Models/...
[INFO] Starting pipeline at 2025-XX-XX XX:XX:XX
[INFO] Loading FLAMES data and performing initial QC...
[INFO] Auto-detected 8 valid FLAMES samples
[INFO] Samples: org_1A, org_1B, org_3A, org_3B, org_3C, org_6A, org_6B, org_6C
[INFO] Loading FLAMES data for sample: org_1A
[INFO] Loaded org_1A: 385659 features × XXXX cells
[INFO] Loading gene counts for sample: org_1A
[INFO] Loaded gene counts for org_1A: 52622 genes × XXXX cells
[... repeated for all samples ...]
[INFO] Data loading and initial QC completed
[INFO] Average isoforms per gene: X.XX

[INFO] Starting Step: Quality Control
[DEBUG] Processing QC for sample: org_1A
[... repeated for all samples ...]
[INFO] Detailed quality control completed
[INFO] Average isoforms per cell (median): XXXX
[INFO] Average genes per cell (median): XXXX
[INFO] Average isoforms per gene per cell: X.XX
```

### Output Files Generated

#### Step 1: Data Loading QC
**PDF Report:** `output/reports/01_Data_Loading_QC.pdf` (6 plots in 2×3 grid)

**Individual Plots:**
1. `output/figures/qc/01_cells_per_sample.png/pdf` - Total cells: X (Unfiltered)
2. `output/figures/qc/02_isoforms_vs_cells.png/pdf` - Avg isoforms: XXXX (Unfiltered)
3. `output/figures/qc/03_genes_vs_cells.png/pdf` - Avg genes: XXXX (Unfiltered)
4. `output/figures/qc/04_isoforms_per_gene_ratio.png/pdf` - Overall avg: X.XX (Unfiltered)
5. `output/figures/qc/05_total_isoform_counts.png/pdf` - (Unfiltered)
6. `output/figures/qc/06_total_gene_counts.png/pdf` - (Unfiltered)

**Data Files:**
- `output/data/metadata/sample_metadata.csv/rds` - With new columns: `n_genes`, `total_gene_counts`, `isoforms_per_gene`

#### Step 2: Detailed Quality Control
**PDF Report:** `output/reports/02_Quality_Control.pdf` (10 plots in 2×5 grid)

**Individual Plots:** (numbered 07-16)
- Isoform-level: 07, 09, 11, 13, 15 (odd numbers)
- Gene-level: 08, 10, 12, 14, 16 (even numbers)

All plots include:
- **"(Unfiltered)"** label in title
- **Statistics box** in top-right corner with:
  - n = number of cells
  - Median/Mean values
  - Q25-Q75 quartiles (histograms)
  - Correlation (scatter plots)
  - Per-group statistics (violin plots)

**Data Files:**
- `output/data/metadata/detailed_qc_isoform.csv/rds`
- `output/data/metadata/detailed_qc_gene.csv/rds`
- `output/data/metadata/detailed_qc_combined.csv/rds`

## Expected QC Metrics

Based on FLAMES output structure:

| Metric | Expected Range |
|--------|---------------|
| **Isoforms per sample** | ~385,000 transcripts (ENST IDs) |
| **Genes per sample** | ~52,000 genes (ENSG IDs) |
| **Isoforms per gene** | ~7-8 on average |
| **Cells per sample** | 600-2000 (from sample metadata) |
| **Isoforms per cell** | Variable (QC will show distribution) |
| **Genes per cell** | Variable (QC will show distribution) |

## Plot Features

### All Histograms Include:
- Blue bars with distribution
- Red dashed line at median
- Statistics box (top-right):
  - n = total cells
  - Median
  - Mean
  - Q25-Q75 quartiles

### All Scatter Plots Include:
- Points colored by timepoint
- Sample labels (repelled text)
- Statistics box (top-right):
  - n = total cells
  - Correlation coefficient

### All Violin Plots Include:
- Violin + boxplot overlay
- Colored by timepoint
- Statistics box (top-right):
  - Per-timepoint: n, median

## How to Run

```bash
cd /data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/previous_project_analysis_scripts

# Open R
R

# Run the pipeline
source("analysis_pipeline/main.R")
```

Or if the script is being sourced:
```r
main()
```

## Potential Issues to Watch For

1. **Missing gene count files**: If any `org_XY_gene_count.csv` files are missing, that sample will have NA for gene metrics
2. **Memory usage**: Loading all 8 samples × 385K isoforms may require significant RAM
3. **MT gene detection**:
   - Isoforms: Pattern `^.*-MT-` (transcript names like ENST-MT-XXX)
   - Genes: Pattern `^MT-` (gene names like MT-CO1)
4. **Plot rendering**: 12-plot PDF may take time to render

## Verification Checklist

After running, verify:
- [ ] All 8 samples detected (org_1A through org_6C)
- [ ] Both isoform and gene counts loaded for all samples
- [ ] Isoforms per gene ratio is ~7-8
- [ ] All plots have "(Unfiltered)" label
- [ ] All plots have statistics annotations
- [ ] Both PDF reports generated (01 and 02)
- [ ] All individual PNG/PDF files generated
- [ ] CSV files contain expected columns with no NAs

## Next Steps After Testing

Once Steps 1-2 are verified:
1. Implement **Step 3: Integration** with both RNA and iso assays
2. Ensure filtered data gets `filtered = TRUE` parameter in plots
3. Create comparison plots (filtered vs unfiltered)
