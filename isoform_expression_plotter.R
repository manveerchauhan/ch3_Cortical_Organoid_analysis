#!/usr/bin/env Rscript
# Isoform Expression Plotting Functions
# Author: Manveer Chauhan
# Description: Generalized functions for plotting isoform-level expression patterns
#              for genes across developmental timepoints

# Required libraries
library(Seurat)
library(tidyverse)
library(patchwork)

# ===== Core Functions =====

#' Extract all isoforms for a given parent gene symbol
#'
#' @param seu.obj Seurat object with 'iso' assay
#' @param gene_symbol Parent gene symbol (e.g., "NRXN1")
#' @return Character vector of isoform IDs matching the gene symbol
get_isoforms_for_gene <- function(seu.obj, gene_symbol) {
  # Get all features from the isoform assay
  all_isoforms <- rownames(seu.obj[["iso"]])

  # Filter for isoforms containing the gene symbol
  # Pattern: gene symbol appears after a hyphen (e.g., "ENST00000335137.4-NRXN1")
  pattern <- paste0("-", gene_symbol, "$")
  matching_isoforms <- grep(pattern, all_isoforms, value = TRUE)

  return(matching_isoforms)
}

#' Create violin plot for a single isoform across cell types
#'
#' @param seu.obj Seurat object with 'iso' assay
#' @param isoform_id Full isoform ID (e.g., "ENST00000335137.4-NRXN1")
#' @param timepoint_label Label for timepoint (e.g., "1M", "3M", "6M")
#' @return ggplot object
plot_isoform_expression <- function(seu.obj, isoform_id, timepoint_label) {
  # Extract expression data
  expr_data <- FetchData(seu.obj,
                         vars = c(isoform_id, "cluster_annotations"),
                         slot = "data",
                         assay = "iso")

  colnames(expr_data) <- c("Expression", "CellType")

  # Create violin plot
  p <- ggplot(expr_data, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(trim = FALSE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    labs(title = paste0(timepoint_label, ": ", isoform_id),
         x = "Cell Type",
         y = "Normalized Expression (log)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "none",
          panel.grid.minor = element_blank())

  return(p)
}

#' Create UMAP feature plot for a single isoform
#'
#' @param seu.obj Seurat object with 'iso' assay and 'umap.harm' reduction
#' @param isoform_id Full isoform ID
#' @param timepoint_label Label for timepoint
#' @return ggplot object
plot_isoform_umap <- function(seu.obj, isoform_id, timepoint_label) {
  # Get expression data to calculate max
  expr_data <- FetchData(seu.obj, vars = isoform_id, slot = "data", assay = "iso")
  maxExpr <- max(expr_data, na.rm = TRUE)

  p <- FeaturePlot(seu.obj,
                   features = isoform_id,
                   reduction = "umap.harm",
                   order = TRUE,
                   pt.size = 0.5) +
    scale_colour_gradient(low = "lightgrey",
                         high = "#e31837",
                         limits = c(0, maxExpr),
                         na.value = "#e31837") +
    labs(title = paste0(timepoint_label, ": ", isoform_id)) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))

  return(p)
}

#' Create dot plot showing isoform expression across cell types
#'
#' @param seu.obj Seurat object with 'iso' assay
#' @param isoforms Vector of isoform IDs
#' @param gene_symbol Parent gene symbol
#' @param timepoint_label Label for timepoint
#' @return ggplot object or NULL if no isoforms
plot_isoform_dotplot <- function(seu.obj, isoforms, gene_symbol, timepoint_label) {
  if(length(isoforms) == 0) return(NULL)

  p <- DotPlot(seu.obj,
               features = isoforms,
               assay = "iso",
               dot.scale = 8) +
    coord_flip() +
    labs(title = paste0(timepoint_label, ": ", gene_symbol, " Isoforms"),
         x = "Isoform ID",
         y = "Cell Type") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position = "right")

  return(p)
}

#' Generate comprehensive isoform expression report for a gene across timepoints
#'
#' @param seurat_list Named list of Seurat objects (e.g., list("1M" = obj1, "3M" = obj2, "6M" = obj3))
#' @param gene_symbol Parent gene symbol to analyze
#' @param output_filename Path to output PDF (optional, auto-generated if NULL)
#' @return Path to generated PDF report
generate_isoform_report <- function(seurat_list,
                                   gene_symbol,
                                   output_filename = NULL) {

  cat("\n=== Generating isoform expression report for", gene_symbol, "===\n")

  if(is.null(output_filename)) {
    output_filename <- paste0(gene_symbol, "_Isoform_Expression_Report.pdf")
  }

  # Get timepoint labels
  timepoint_labels <- names(seurat_list)

  # Get isoforms for each timepoint
  isoforms_by_timepoint <- list()
  for(tp in timepoint_labels) {
    isoforms_by_timepoint[[tp]] <- get_isoforms_for_gene(seurat_list[[tp]], gene_symbol)
    cat("Found", length(isoforms_by_timepoint[[tp]]), "isoforms in", tp, "\n")
  }

  # Get unique isoforms across all timepoints
  all_isoforms <- unique(unlist(isoforms_by_timepoint))

  if(length(all_isoforms) == 0) {
    cat("WARNING: No isoforms found for", gene_symbol, "\n")
    cat("Checking alternative patterns...\n")

    # Try alternative pattern (gene symbol anywhere in the name)
    alt_pattern <- gene_symbol
    for(tp in timepoint_labels) {
      all_features <- rownames(seurat_list[[tp]][["iso"]])
      isoforms_by_timepoint[[tp]] <- grep(alt_pattern, all_features, value = TRUE)
    }
    all_isoforms <- unique(unlist(isoforms_by_timepoint))

    cat("Found", length(all_isoforms), "isoforms with alternative pattern\n")

    if(length(all_isoforms) == 0) {
      cat("ERROR: Still no isoforms found. Aborting.\n")
      return(NULL)
    }
  }

  cat("Total unique isoforms:", length(all_isoforms), "\n")
  cat("Isoforms:", paste(all_isoforms, collapse = ", "), "\n")

  # Open PDF device
  pdf(output_filename, width = 16, height = 10)

  # Generate UMAP plots for each isoform
  cat("Creating UMAP plots for isoforms...\n")

  for(isoform in all_isoforms) {
    cat("  Processing", isoform, "...\n")

    plot_list <- list()

    # Check presence in each timepoint and create UMAP plots
    for(tp in timepoint_labels) {
      if(isoform %in% isoforms_by_timepoint[[tp]]) {
        umap_plot <- plot_isoform_umap(seurat_list[[tp]], isoform, tp)
        plot_list <- c(plot_list, list(umap_plot))
      }
    }

    # Arrange plots on page (up to 3 timepoints per isoform)
    if(length(plot_list) > 0) {
      combined <- wrap_plots(plot_list, ncol = 3) +
        plot_annotation(title = paste0(gene_symbol, ": ", isoform),
                       theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
      print(combined)
    }
  }

  # Close PDF
  dev.off()

  cat("\nSaved isoform report to:", output_filename, "\n")
  cat("Report contains:\n")
  cat("  - UMAP plots for", length(all_isoforms), "isoforms across all timepoints\n")
  cat("  - Custom color gradient: lightgrey -> #e31837\n")

  return(output_filename)
}

cat("Isoform expression plotting functions loaded successfully.\n")
cat("Available functions:\n")
cat("  - get_isoforms_for_gene(seu.obj, gene_symbol)\n")
cat("  - plot_isoform_expression(seu.obj, isoform_id, timepoint_label)\n")
cat("  - plot_isoform_umap(seu.obj, isoform_id, timepoint_label)\n")
cat("  - plot_isoform_dotplot(seu.obj, isoforms, gene_symbol, timepoint_label)\n")
cat("  - generate_isoform_report(seurat_list, gene_symbol, output_filename = NULL)\n")
