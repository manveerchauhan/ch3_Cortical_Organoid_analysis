#!/usr/bin/env Rscript
# QC Visualization Script for Filtered Seurat Objects
# Author: Manveer Chauhan
# Description: Generates comprehensive QC metrics and visualizations for presentation

library(Seurat)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(scales)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis")

# Create output directory
output_dir <- "QC_Summary_Visualizations"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ===== Load Seurat Objects =====
cat("Loading Seurat objects...\n")
one_month_org <- readRDS("one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("six_month_seurat.intergrated_harm.isofrom.rds")

Idents(one_month_org) <- "cluster_annotations"
Idents(three_month_org) <- "cluster_annotations"
Idents(six_month_org) <- "cluster_annotations"

cat("Successfully loaded all Seurat objects\n\n")

# ===== Function: Extract comprehensive QC metrics =====
extract_QC_metrics <- function(seu.obj, timepoint_label) {
  cat("Extracting QC metrics for", timepoint_label, "...\n")

  metadata <- seu.obj@meta.data

  # Basic cell metrics
  n_cells <- ncol(seu.obj)
  n_genes_detected <- nrow(seu.obj[["RNA"]])
  n_isoforms_detected <- nrow(seu.obj[["iso"]])

  # Cell type distribution
  celltype_counts <- table(Idents(seu.obj))
  n_celltypes <- length(unique(Idents(seu.obj)))

  # Sample/batch distribution
  sample_counts <- table(metadata$orig.ident)
  n_samples <- length(unique(metadata$orig.ident))

  # QC metric ranges - Gene level (RNA assay)
  median_nCount_RNA = median(metadata$nCount_RNA, na.rm = TRUE)
  mean_nCount_RNA = mean(metadata$nCount_RNA, na.rm = TRUE)
  median_nFeature_RNA = median(metadata$nFeature_RNA, na.rm = TRUE)
  mean_nFeature_RNA = mean(metadata$nFeature_RNA, na.rm = TRUE)

  # QC metric ranges - Isoform level (iso assay)
  median_nCount_iso = median(metadata$nCount_iso, na.rm = TRUE)
  mean_nCount_iso = mean(metadata$nCount_iso, na.rm = TRUE)
  median_nFeature_iso = median(metadata$nFeature_iso, na.rm = TRUE)
  mean_nFeature_iso = mean(metadata$nFeature_iso, na.rm = TRUE)

  # Mitochondrial percentage (same for both)
  median_percent_mt = median(metadata$percent.mt, na.rm = TRUE)
  mean_percent_mt = mean(metadata$percent.mt, na.rm = TRUE)

  metrics <- list(
    timepoint = timepoint_label,
    n_cells = n_cells,
    n_genes = n_genes_detected,
    n_isoforms = n_isoforms_detected,
    n_celltypes = n_celltypes,
    n_samples = n_samples,

    # Gene-level metrics
    median_nCount_RNA = median_nCount_RNA,
    mean_nCount_RNA = mean_nCount_RNA,
    median_nFeature_RNA = median_nFeature_RNA,
    mean_nFeature_RNA = mean_nFeature_RNA,

    # Isoform-level metrics
    median_nCount_iso = median_nCount_iso,
    mean_nCount_iso = mean_nCount_iso,
    median_nFeature_iso = median_nFeature_iso,
    mean_nFeature_iso = mean_nFeature_iso,

    # Mitochondrial percentage
    median_percent_mt = median_percent_mt,
    mean_percent_mt = mean_percent_mt,

    # Store tables
    celltype_distribution = celltype_counts,
    sample_distribution = sample_counts
  )

  return(metrics)
}

# ===== Extract metrics for all timepoints =====
qc_1M <- extract_QC_metrics(one_month_org, "1M")
qc_3M <- extract_QC_metrics(three_month_org, "3M")
qc_6M <- extract_QC_metrics(six_month_org, "6M")

# ===== Generate Summary Statistics Table =====
cat("\n=== GENERATING SUMMARY STATISTICS TABLE ===\n")

summary_stats <- data.frame(
  Metric = c("Total Cells", "Detected Genes (RNA)", "Detected Isoforms",
             "Cell Types", "Samples/Replicates",
             "Median Gene UMI Count", "Mean Gene UMI Count",
             "Median Gene Features/Cell", "Mean Gene Features/Cell",
             "Median Isoform UMI Count", "Mean Isoform UMI Count",
             "Median Isoform Features/Cell", "Mean Isoform Features/Cell",
             "Median % Mitochondrial", "Mean % Mitochondrial"),
  `1_Month` = c(
    format(qc_1M$n_cells, big.mark = ","),
    format(qc_1M$n_genes, big.mark = ","),
    format(qc_1M$n_isoforms, big.mark = ","),
    qc_1M$n_celltypes,
    qc_1M$n_samples,
    format(round(qc_1M$median_nCount_RNA, 0), big.mark = ","),
    format(round(qc_1M$mean_nCount_RNA, 0), big.mark = ","),
    format(round(qc_1M$median_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_1M$mean_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_1M$median_nCount_iso, 0), big.mark = ","),
    format(round(qc_1M$mean_nCount_iso, 0), big.mark = ","),
    format(round(qc_1M$median_nFeature_iso, 0), big.mark = ","),
    format(round(qc_1M$mean_nFeature_iso, 0), big.mark = ","),
    sprintf("%.2f%%", qc_1M$median_percent_mt),
    sprintf("%.2f%%", qc_1M$mean_percent_mt)
  ),
  `3_Month` = c(
    format(qc_3M$n_cells, big.mark = ","),
    format(qc_3M$n_genes, big.mark = ","),
    format(qc_3M$n_isoforms, big.mark = ","),
    qc_3M$n_celltypes,
    qc_3M$n_samples,
    format(round(qc_3M$median_nCount_RNA, 0), big.mark = ","),
    format(round(qc_3M$mean_nCount_RNA, 0), big.mark = ","),
    format(round(qc_3M$median_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_3M$mean_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_3M$median_nCount_iso, 0), big.mark = ","),
    format(round(qc_3M$mean_nCount_iso, 0), big.mark = ","),
    format(round(qc_3M$median_nFeature_iso, 0), big.mark = ","),
    format(round(qc_3M$mean_nFeature_iso, 0), big.mark = ","),
    sprintf("%.2f%%", qc_3M$median_percent_mt),
    sprintf("%.2f%%", qc_3M$mean_percent_mt)
  ),
  `6_Month` = c(
    format(qc_6M$n_cells, big.mark = ","),
    format(qc_6M$n_genes, big.mark = ","),
    format(qc_6M$n_isoforms, big.mark = ","),
    qc_6M$n_celltypes,
    qc_6M$n_samples,
    format(round(qc_6M$median_nCount_RNA, 0), big.mark = ","),
    format(round(qc_6M$mean_nCount_RNA, 0), big.mark = ","),
    format(round(qc_6M$median_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_6M$mean_nFeature_RNA, 0), big.mark = ","),
    format(round(qc_6M$median_nCount_iso, 0), big.mark = ","),
    format(round(qc_6M$mean_nCount_iso, 0), big.mark = ","),
    format(round(qc_6M$median_nFeature_iso, 0), big.mark = ","),
    format(round(qc_6M$mean_nFeature_iso, 0), big.mark = ","),
    sprintf("%.2f%%", qc_6M$median_percent_mt),
    sprintf("%.2f%%", qc_6M$mean_percent_mt)
  ),
  check.names = FALSE
)

# Save as CSV
write.csv(summary_stats,
          file.path(output_dir, "Summary_Statistics_Table.csv"),
          row.names = FALSE)

# Print to console
cat("\nSummary Statistics:\n")
print(summary_stats)

# ===== Cell Type Distribution Plots =====
cat("\n=== GENERATING CELL TYPE DISTRIBUTION PLOTS ===\n")

create_celltype_barplot <- function(seu.obj, timepoint_label, fill_color = "#4A90E2") {
  df <- data.frame(CellType = Idents(seu.obj)) %>%
    group_by(CellType) %>%
    summarise(Count = n()) %>%
    arrange(desc(Count))

  p <- ggplot(df, aes(x = reorder(CellType, Count), y = Count)) +
    geom_bar(stat = "identity", fill = fill_color, alpha = 0.8) +
    geom_text(aes(label = format(Count, big.mark = ",")),
              hjust = -0.2, size = 3) +
    coord_flip() +
    labs(title = paste0(timepoint_label, " Cell Type Distribution"),
         x = "Cell Type",
         y = "Number of Cells") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank())

  return(p)
}

# Generate individual barplots
p_celltype_1M <- create_celltype_barplot(one_month_org, "1 Month", "#E74C3C")
p_celltype_3M <- create_celltype_barplot(three_month_org, "3 Month", "#F39C12")
p_celltype_6M <- create_celltype_barplot(six_month_org, "6 Month", "#27AE60")

# Combined plot
combined_celltype <- p_celltype_1M / p_celltype_3M / p_celltype_6M
ggsave(file.path(output_dir, "CellType_Distribution_AllTimepoints.png"),
       combined_celltype, width = 12, height = 14, dpi = 300)

cat("Saved: CellType_Distribution_AllTimepoints.png\n")

# ===== Stacked cell type comparison =====
celltype_comparison_df <- bind_rows(
  data.frame(CellType = Idents(one_month_org), Timepoint = "1M"),
  data.frame(CellType = Idents(three_month_org), Timepoint = "3M"),
  data.frame(CellType = Idents(six_month_org), Timepoint = "6M")
) %>%
  group_by(Timepoint, CellType) %>%
  summarise(Count = n(), .groups = "drop")

p_stacked <- ggplot(celltype_comparison_df,
                    aes(x = Timepoint, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = format(Count, big.mark = ",")),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Cell Type Composition Across Developmental Timepoints",
       x = "Timepoint",
       y = "Number of Cells",
       fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(output_dir, "CellType_Stacked_Comparison.png"),
       p_stacked, width = 14, height = 8, dpi = 300)

cat("Saved: CellType_Stacked_Comparison.png\n")

# ===== QC Metric Violin Plots =====
cat("\n=== GENERATING QC VIOLIN PLOTS ===\n")

create_QC_violin <- function(seu.obj, timepoint_label, metric, y_label) {
  df <- seu.obj@meta.data %>%
    mutate(CellType = Idents(seu.obj))

  p <- ggplot(df, aes(x = CellType, y = !!sym(metric), fill = CellType)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    labs(title = paste0(timepoint_label, ": ", y_label),
         x = "Cell Type",
         y = y_label) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")

  return(p)
}

# === GENE-LEVEL VIOLINS (RNA assay) ===

# Gene UMI count violins
v_umi_gene_1M <- create_QC_violin(one_month_org, "1M", "nCount_RNA", "Gene UMI Counts")
v_umi_gene_3M <- create_QC_violin(three_month_org, "3M", "nCount_RNA", "Gene UMI Counts")
v_umi_gene_6M <- create_QC_violin(six_month_org, "6M", "nCount_RNA", "Gene UMI Counts")

combined_umi_gene <- v_umi_gene_1M / v_umi_gene_3M / v_umi_gene_6M
ggsave(file.path(output_dir, "QC_Gene_UMI_Counts_Violins.png"),
       combined_umi_gene, width = 14, height = 16, dpi = 300)

cat("Saved: QC_Gene_UMI_Counts_Violins.png\n")

# Gene feature count violins
v_feat_gene_1M <- create_QC_violin(one_month_org, "1M", "nFeature_RNA", "Gene Features Detected")
v_feat_gene_3M <- create_QC_violin(three_month_org, "3M", "nFeature_RNA", "Gene Features Detected")
v_feat_gene_6M <- create_QC_violin(six_month_org, "6M", "nFeature_RNA", "Gene Features Detected")

combined_feat_gene <- v_feat_gene_1M / v_feat_gene_3M / v_feat_gene_6M
ggsave(file.path(output_dir, "QC_Gene_Feature_Counts_Violins.png"),
       combined_feat_gene, width = 14, height = 16, dpi = 300)

cat("Saved: QC_Gene_Feature_Counts_Violins.png\n")

# === ISOFORM-LEVEL VIOLINS (iso assay) ===

# Isoform UMI count violins
v_umi_iso_1M <- create_QC_violin(one_month_org, "1M", "nCount_iso", "Isoform UMI Counts")
v_umi_iso_3M <- create_QC_violin(three_month_org, "3M", "nCount_iso", "Isoform UMI Counts")
v_umi_iso_6M <- create_QC_violin(six_month_org, "6M", "nCount_iso", "Isoform UMI Counts")

combined_umi_iso <- v_umi_iso_1M / v_umi_iso_3M / v_umi_iso_6M
ggsave(file.path(output_dir, "QC_Isoform_UMI_Counts_Violins.png"),
       combined_umi_iso, width = 14, height = 16, dpi = 300)

cat("Saved: QC_Isoform_UMI_Counts_Violins.png\n")

# Isoform feature count violins
v_feat_iso_1M <- create_QC_violin(one_month_org, "1M", "nFeature_iso", "Isoform Features Detected")
v_feat_iso_3M <- create_QC_violin(three_month_org, "3M", "nFeature_iso", "Isoform Features Detected")
v_feat_iso_6M <- create_QC_violin(six_month_org, "6M", "nFeature_iso", "Isoform Features Detected")

combined_feat_iso <- v_feat_iso_1M / v_feat_iso_3M / v_feat_iso_6M
ggsave(file.path(output_dir, "QC_Isoform_Feature_Counts_Violins.png"),
       combined_feat_iso, width = 14, height = 16, dpi = 300)

cat("Saved: QC_Isoform_Feature_Counts_Violins.png\n")

# === MITOCHONDRIAL PERCENTAGE (same for both) ===

# Mitochondrial percentage violins
v_mt_1M <- create_QC_violin(one_month_org, "1M", "percent.mt", "% Mitochondrial")
v_mt_3M <- create_QC_violin(three_month_org, "3M", "percent.mt", "% Mitochondrial")
v_mt_6M <- create_QC_violin(six_month_org, "6M", "percent.mt", "% Mitochondrial")

combined_mt <- v_mt_1M / v_mt_3M / v_mt_6M
ggsave(file.path(output_dir, "QC_Mitochondrial_Violins.png"),
       combined_mt, width = 14, height = 16, dpi = 300)

cat("Saved: QC_Mitochondrial_Violins.png\n")

# ===== UMAP Visualizations =====
cat("\n=== GENERATING UMAP VISUALIZATIONS ===\n")

# Function to create UMAP plots
create_umap_plot <- function(seu.obj, timepoint_label) {
  p <- DimPlot(seu.obj,
               reduction = "umap.harm",
               label = TRUE,
               label.size = 3.5,
               repel = TRUE) +
    labs(title = paste0(timepoint_label, " Organoids - UMAP")) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  return(p)
}

umap_1M <- create_umap_plot(one_month_org, "1 Month")
umap_3M <- create_umap_plot(three_month_org, "3 Month")
umap_6M <- create_umap_plot(six_month_org, "6 Month")

# Save individual UMAPs
ggsave(file.path(output_dir, "UMAP_1M_Organoids.png"),
       umap_1M, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "UMAP_3M_Organoids.png"),
       umap_3M, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "UMAP_6M_Organoids.png"),
       umap_6M, width = 10, height = 8, dpi = 300)

# Combined UMAP panel
combined_umap <- umap_1M + umap_3M + umap_6M +
  plot_layout(ncol = 3)

ggsave(file.path(output_dir, "UMAP_Combined_AllTimepoints.png"),
       combined_umap, width = 24, height = 8, dpi = 300)

cat("Saved: UMAP visualizations\n")

# ===== Sample Distribution Analysis =====
cat("\n=== GENERATING SAMPLE DISTRIBUTION PLOTS ===\n")

create_sample_distribution <- function(seu.obj, timepoint_label) {
  df <- seu.obj@meta.data %>%
    group_by(orig.ident, cluster_annotations) %>%
    summarise(Count = n(), .groups = "drop")

  p <- ggplot(df, aes(x = orig.ident, y = Count, fill = cluster_annotations)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = paste0(timepoint_label, " - Cell Type Proportions by Sample"),
         x = "Sample ID",
         y = "Proportion",
         fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5))

  return(p)
}

sample_1M <- create_sample_distribution(one_month_org, "1M")
sample_3M <- create_sample_distribution(three_month_org, "3M")
sample_6M <- create_sample_distribution(six_month_org, "6M")

combined_samples <- sample_1M / sample_3M / sample_6M
ggsave(file.path(output_dir, "Sample_CellType_Proportions.png"),
       combined_samples, width = 14, height = 16, dpi = 300)

cat("Saved: Sample_CellType_Proportions.png\n")

# ===== Export detailed cell type count table =====
celltype_detailed <- bind_rows(
  data.frame(Timepoint = "1M",
             CellType = names(table(Idents(one_month_org))),
             Count = as.numeric(table(Idents(one_month_org)))),
  data.frame(Timepoint = "3M",
             CellType = names(table(Idents(three_month_org))),
             Count = as.numeric(table(Idents(three_month_org)))),
  data.frame(Timepoint = "6M",
             CellType = names(table(Idents(six_month_org))),
             Count = as.numeric(table(Idents(six_month_org))))
) %>%
  pivot_wider(names_from = Timepoint, values_from = Count, values_fill = 0)

write.csv(celltype_detailed,
          file.path(output_dir, "Detailed_CellType_Counts.csv"),
          row.names = FALSE)

cat("Saved: Detailed_CellType_Counts.csv\n")

# ===== Final Summary Report =====
cat("\n========================================\n")
cat("QC VISUALIZATION SCRIPT COMPLETE\n")
cat("========================================\n\n")

cat("Output Directory:", output_dir, "\n\n")

cat("Generated Files:\n")
cat("  1. Summary_Statistics_Table.csv - Main summary metrics (gene + isoform)\n")
cat("  2. Detailed_CellType_Counts.csv - Cell type counts per timepoint\n")
cat("  3. CellType_Distribution_AllTimepoints.png - Bar plots of cell types\n")
cat("  4. CellType_Stacked_Comparison.png - Stacked comparison\n")
cat("  5. QC_Gene_UMI_Counts_Violins.png - Gene UMI distribution by cell type\n")
cat("  6. QC_Gene_Feature_Counts_Violins.png - Gene feature distribution by cell type\n")
cat("  7. QC_Isoform_UMI_Counts_Violins.png - Isoform UMI distribution by cell type\n")
cat("  8. QC_Isoform_Feature_Counts_Violins.png - Isoform feature distribution by cell type\n")
cat("  9. QC_Mitochondrial_Violins.png - MT% distribution by cell type\n")
cat(" 10. UMAP_Combined_AllTimepoints.png - Combined UMAP visualization\n")
cat(" 11. UMAP_[1M/3M/6M]_Organoids.png - Individual UMAP plots\n")
cat(" 12. Sample_CellType_Proportions.png - Sample composition analysis\n\n")

cat("Summary of Datasets:\n")
cat("  1M Organoids:", format(ncol(one_month_org), big.mark = ","),
    "cells across", length(unique(Idents(one_month_org))), "cell types\n")
cat("  3M Organoids:", format(ncol(three_month_org), big.mark = ","),
    "cells across", length(unique(Idents(three_month_org))), "cell types\n")
cat("  6M Organoids:", format(ncol(six_month_org), big.mark = ","),
    "cells across", length(unique(Idents(six_month_org))), "cell types\n")

cat("\nAll visualizations ready for presentation!\n")
cat("\nNote: For isoform-level expression plotting, use 'isoform_expression_plotter.R'\n")
cat("Example: source('isoform_expression_plotter.R')\n")
