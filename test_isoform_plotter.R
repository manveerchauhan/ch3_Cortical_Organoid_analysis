#!/usr/bin/env Rscript
# Test Script for Isoform Expression Plotter
# Author: Manveer Chauhan
# Description: Example usage of isoform expression plotting functions

library(Seurat)
library(tidyverse)
library(patchwork)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis")

# Source the isoform plotting functions
source("isoform_expression_plotter.R")

cat("\n=== Loading Seurat Objects ===\n")

# Load Seurat objects
one_month_org <- readRDS("one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("six_month_seurat.intergrated_harm.isofrom.rds")

# Set identities
Idents(one_month_org) <- "cluster_annotations"
Idents(three_month_org) <- "cluster_annotations"
Idents(six_month_org) <- "cluster_annotations"

cat("Seurat objects loaded successfully\n")

# Create named list of Seurat objects
seurat_timepoints <- list(
  "1M" = one_month_org,
  "3M" = three_month_org,
  "6M" = six_month_org
)

# ===== Example 1: Generate NRXN1 isoform report =====
cat("\n=== Example 1: Generating NRXN1 Isoform Report ===\n")

nrxn1_report <- generate_isoform_report(
  seurat_list = seurat_timepoints,
  gene_symbol = "NRXN1",
  output_filename = "QC_Summary_Visualizations/NRXN1_Isoform_Expression_Report.pdf"
)

cat("\nNRXN1 report generated:", nrxn1_report, "\n")

# ===== Example 2: Generate reports for synaptic adhesion genes =====
cat("\n=== Example 2: Generating Reports for Synaptic Adhesion Genes ===\n")

# PTBP family
generate_isoform_report(seurat_timepoints, "PTBP1", "QC_Summary_Visualizations/PTBP1_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "PTBP2", "QC_Summary_Visualizations/PTBP2_Isoform_Expression_Report.pdf")

# NRXN family
generate_isoform_report(seurat_timepoints, "NRXN2", "QC_Summary_Visualizations/NRXN2_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "NRXN3", "QC_Summary_Visualizations/NRXN3_Isoform_Expression_Report.pdf")

# NLGN family
generate_isoform_report(seurat_timepoints, "NLGN1", "QC_Summary_Visualizations/NLGN1_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN2", "QC_Summary_Visualizations/NLGN2_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN3", "QC_Summary_Visualizations/NLGN3_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN4X", "QC_Summary_Visualizations/NLGN4X_Isoform_Expression_Report.pdf")
generate_isoform_report(seurat_timepoints, "NLGN4Y", "QC_Summary_Visualizations/NLGN4Y_Isoform_Expression_Report.pdf")


# ===== Example 3: Using individual plotting functions =====
cat("\n=== Example 3: Using Individual Plotting Functions ===\n")

# Get all NRXN1 isoforms from 6M organoids
nrxn1_isoforms_6M <- get_isoforms_for_gene(six_month_org, "NRXN1")
cat("Found", length(nrxn1_isoforms_6M), "NRXN1 isoforms in 6M organoids:\n")
cat(paste(nrxn1_isoforms_6M, collapse = "\n"), "\n\n")

# Create individual plots for the first isoform (if available)
if(length(nrxn1_isoforms_6M) > 0) {
  first_isoform <- nrxn1_isoforms_6M[1]
  cat("Creating example UMAP plot for:", first_isoform, "\n")

  # UMAP plot with custom gradient
  p_umap <- plot_isoform_umap(six_month_org, first_isoform, "6M")

  # Save plot
  ggsave("QC_Summary_Visualizations/Example_Individual_Isoform_UMAP.png",
         p_umap, width = 8, height = 6, dpi = 300)

  cat("Saved example UMAP plot to: QC_Summary_Visualizations/Example_Individual_Isoform_UMAP.png\n")
}

cat("\n=== Test Script Complete ===\n")
cat("All isoform reports have been generated successfully!\n")
