# Script used to perform QC for LR scRNA-seq analysis

## Import seurat object 'atlases' from silvia's group----
one_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/six_month_seurat.intergrated_harm.isofrom.rds")

Idents(one_month_org) <- "cluster_annotations"
Idents(three_month_org) <- "cluster_annotations"
Idents(six_month_org) <- "cluster_annotations"

### Extract raw count matrices (to use as proxy for FLAMES v2.2.0 data)----

## Extract raw count matrices from one_month_org
org1A_counts <- LayerData(one_month_org, assay = "joined", layer = "counts.org_1A")
org1B_counts <- LayerData(one_month_org, assay = "joined", layer = "counts.org_1B")

## Extract raw count matrices from three_month_org
org3A_counts <- LayerData(three_month_org, assay = "joined", layer = "counts.org_3A")
org3B_counts <- LayerData(three_month_org, assay = "joined", layer = "counts.org_3B")
org3C_counts <- LayerData(three_month_org, assay = "joined", layer = "counts.org_3C")

## Extract raw count matrices from six_month_org
org6A_counts <- LayerData(six_month_org, assay = "joined", layer = "counts.org_6A")
org6B_counts <- LayerData(six_month_org, assay = "joined", layer = "counts.org_6B")
org6C_counts <- LayerData(six_month_org, assay = "joined", layer = "counts.org_6C")

## Convert raw count matrices into dataframes to inspect
org1A_counts_df <- as.data.frame(as.matrix(org1A_counts))
org1B_counts_df <- as.data.frame(as.matrix(org1B_counts))

org3A_counts_df <- as.data.frame(as.matrix(org3A_counts))
org3B_counts_df <- as.data.frame(as.matrix(org3B_counts))
org3C_counts_df <- as.data.frame(as.matrix(org3C_counts))

org6A_counts_df <- as.data.frame(as.matrix(org6A_counts))
org6B_counts_df <- as.data.frame(as.matrix(org6B_counts))
org6C_counts_df <- as.data.frame(as.matrix(org6C_counts))

### Extract isoform count matrices----

## Extract isoform counts and filter by sample-specific cells
# Get cell barcodes for each sample from gene count matrices
org1A_cells <- colnames(org1A_counts)
org1B_cells <- colnames(org1B_counts)
org3A_cells <- colnames(org3A_counts)  
org3B_cells <- colnames(org3B_counts)
org3C_cells <- colnames(org3C_counts)
org6A_cells <- colnames(org6A_counts)
org6B_cells <- colnames(org6B_counts)
org6C_cells <- colnames(org6C_counts)

# Extract full isoform matrices from each timepoint object
org1_iso_counts <- LayerData(one_month_org, assay = "iso", layer = "counts")
org3_iso_counts <- LayerData(three_month_org, assay = "iso", layer = "counts") 
org6_iso_counts <- LayerData(six_month_org, assay = "iso", layer = "counts")

# Filter isoform matrices by sample-specific cells
org1A_iso_counts <- org1_iso_counts[, org1A_cells]
org1B_iso_counts <- org1_iso_counts[, org1B_cells]
org3A_iso_counts <- org3_iso_counts[, org3A_cells]
org3B_iso_counts <- org3_iso_counts[, org3B_cells] 
org3C_iso_counts <- org3_iso_counts[, org3C_cells]
org6A_iso_counts <- org6_iso_counts[, org6A_cells]
org6B_iso_counts <- org6_iso_counts[, org6B_cells]
org6C_iso_counts <- org6_iso_counts[, org6C_cells]

# Convert isoform matrices to dataframes
org1A_iso_counts_df <- as.data.frame(as.matrix(org1A_iso_counts))
org1B_iso_counts_df <- as.data.frame(as.matrix(org1B_iso_counts))
org3A_iso_counts_df <- as.data.frame(as.matrix(org3A_iso_counts))
org3B_iso_counts_df <- as.data.frame(as.matrix(org3B_iso_counts))
org3C_iso_counts_df <- as.data.frame(as.matrix(org3C_iso_counts))
org6A_iso_counts_df <- as.data.frame(as.matrix(org6A_iso_counts))
org6B_iso_counts_df <- as.data.frame(as.matrix(org6B_iso_counts))
org6C_iso_counts_df <- as.data.frame(as.matrix(org6C_iso_counts))

## Export count matrices to organized directory structure----
# Create main output directory and subdirectories
output_dir <- "previous_raw_data"
gene_dir <- file.path(output_dir, "gene_counts")
isoform_dir <- file.path(output_dir, "isoform_counts")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
if (!dir.exists(gene_dir)) {
  dir.create(gene_dir)
}
if (!dir.exists(isoform_dir)) {
  dir.create(isoform_dir)
}

# Gene count matrices
gene_matrices <- list(
  "org1A_gene_counts.csv" = org1A_counts_df,
  "org1B_gene_counts.csv" = org1B_counts_df,
  "org3A_gene_counts.csv" = org3A_counts_df,
  "org3B_gene_counts.csv" = org3B_counts_df,
  "org3C_gene_counts.csv" = org3C_counts_df,
  "org6A_gene_counts.csv" = org6A_counts_df,
  "org6B_gene_counts.csv" = org6B_counts_df,
  "org6C_gene_counts.csv" = org6C_counts_df
)

# Isoform count matrices
isoform_matrices <- list(
  "org1A_transcript_counts.csv" = org1A_iso_counts_df,
  "org1B_transcript_counts.csv" = org1B_iso_counts_df,
  "org3A_transcript_counts.csv" = org3A_iso_counts_df,
  "org3B_transcript_counts.csv" = org3B_iso_counts_df,
  "org3C_transcript_counts.csv" = org3C_iso_counts_df,
  "org6A_transcript_counts.csv" = org6A_iso_counts_df,
  "org6B_transcript_counts.csv" = org6B_iso_counts_df,
  "org6C_transcript_counts.csv" = org6C_iso_counts_df
)

# Export gene count matrices
for (filename in names(gene_matrices)) {
  filepath <- file.path(gene_dir, filename)
  write.csv(gene_matrices[[filename]], filepath, row.names = TRUE)
  message(paste("Exported gene counts:", filepath))
}

# Export isoform count matrices
for (filename in names(isoform_matrices)) {
  filepath <- file.path(isoform_dir, filename)
  write.csv(isoform_matrices[[filename]], filepath, row.names = TRUE)
  message(paste("Exported isoform counts:", filepath))
}

message(paste("All count matrices exported to:", output_dir))
message(paste("Gene matrices in:", gene_dir))
message(paste("Isoform matrices in:", isoform_dir))