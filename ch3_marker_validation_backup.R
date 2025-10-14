## Marker Gene Validation Script
# Author: Manveer Chauhan
# This script validates cell type marker genes across integrated organoid timepoints
# using module scores and top individual gene expression patterns

# Load essential packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)

# Set working directory
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/draft_space")

# Create output directory
dir.create("./output_files/marker_validation", recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("MARKER GENE VALIDATION ANALYSIS\n")
cat("========================================\n\n")

## DEFINE MARKER GENE SETS------
# 14 marker gene sets from original ch3-marker-check-script.R

marker_sets <- list(
  "Outer Radial Glia (oRG)" = c("HOPX", "TNC", "PTPRZ1", "FAM107A", "LIFR", "MOXD1"),

  "Ventricular Radial Glia (vRG)" = c("VIM", "NES", "PAX6", "HES1", "SLC1A3", "CDH2", "SOX2"),

  "SST Interneurons" = c("SST", "VIP", "CALB2", "CALB1", "NPY"),

  "Parvalbumin Interneurons" = c("PVALB", "SLIT2"),

  "Other Interneurons" = c("NOS1", "CCK", "SLC32A1"),

  "Intermediate Progenitors (IP)" = c("NEUROD6", "PPP1R17", "SOX4", "EOMES", "NEUROG2",
                                       "ASCL1", "DLL1", "DCX", "NES", "PAX6", "SOX2",
                                       "SOX9", "HES1"),

  "Neuronal Intermediate Progenitors" = c("EOMES", "PPP1R17", "NEUROG1"),

  "Stem Cells" = c("SOX2", "NES", "MYC", "PROM1", "POU5F1", "OLIG2", "PAX6", "SOX1"),

  "Glutamatergic Neurons" = c("SLC17A6", "SLC17A7", "SLC17A8"),

  "Microglia" = c("CX3CR1", "OBIF", "SP1999"),

  "Oligodendrocytes" = c("OLIG2", "OLIG1", "MBP", "SPG2", "CNP", "MOG"),

  "Pan-Neuronal Markers" = c("MAPT", "GABBR2", "KLC1", "CACNA1C", "DOC2A",
                              "RBFOX1", "RBFOX2", "CLCN3", "REST", "SRRM4"),

  "Astrocytes" = c("GFAP", "S100B", "ALDH1L1", "GLAST", "SLC1A2", "SLC1A3"),

  "Cajal-Retzius/Subplate" = c("RELN", "SSTR1", "NPY", "CCK", "CR", "CRH", "LHX6", "SLC6A1")
)

cat("Defined", length(marker_sets), "marker gene sets\n")
total_genes <- sum(sapply(marker_sets, length))
cat("Total genes to validate:", total_genes, "\n")

## FUNCTION 1: Resolve gene symbol to full feature name------
resolve_gene_name <- function(gene_symbol, seurat_obj) {
  # Search for gene in format: ENSG...-GENESYMBOL
  pattern <- paste0('-', gene_symbol, '$')
  matches <- grep(pattern, rownames(seurat_obj), value = TRUE)

  if(length(matches) == 0) {
    return(NA)
  } else if(length(matches) == 1) {
    return(matches[1])
  } else {
    # Multiple matches - take first (shouldn't happen often)
    warning("Multiple matches for ", gene_symbol, ": ", paste(matches, collapse = ", "))
    return(matches[1])
  }
}

## FUNCTION 2: Validate gene set against Seurat object------
validate_gene_set <- function(gene_symbols, seurat_obj, set_name) {

  cat("\n  Validating", set_name, ":", length(gene_symbols), "genes\n")

  validation_results <- data.frame(
    gene_symbol = gene_symbols,
    full_name = NA_character_,
    found = FALSE,
    stringsAsFactors = FALSE
  )

  for(i in seq_along(gene_symbols)) {
    full_name <- resolve_gene_name(gene_symbols[i], seurat_obj)
    validation_results$full_name[i] <- full_name
    validation_results$found[i] <- !is.na(full_name)
  }

  n_found <- sum(validation_results$found)
  n_missing <- sum(!validation_results$found)
  pct_found <- round(100 * n_found / length(gene_symbols), 1)

  cat("    Found:", n_found, "/", length(gene_symbols), "(", pct_found, "%)\n")

  if(n_missing > 0) {
    missing_genes <- validation_results$gene_symbol[!validation_results$found]
    cat("    Missing:", paste(missing_genes, collapse = ", "), "\n")
  }

  return(validation_results)
}

## FUNCTION 3: Calculate top expressed genes from a set------
get_top_genes <- function(gene_full_names, seurat_obj, n = 5) {

  # Remove NAs
  valid_genes <- gene_full_names[!is.na(gene_full_names)]

  if(length(valid_genes) == 0) {
    return(character(0))
  }

  # Calculate mean expression for each gene (Seurat v5 compatible)
  mean_expr <- sapply(valid_genes, function(gene) {
    mean(expm1(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[gene, ]))
  })

  # Sort and take top n
  top_genes <- names(sort(mean_expr, decreasing = TRUE))[1:min(n, length(mean_expr))]

  return(top_genes)
}

## FUNCTION 4: Create marker set validation plot------
create_marker_plot <- function(seurat_obj, marker_set_name, validation_results,
                               module_score_name, timepoint_name) {

  cat("\n  Creating plots for", marker_set_name, "\n")

  # Get valid genes
  valid_genes <- validation_results$full_name[validation_results$found]
  missing_genes <- validation_results$gene_symbol[!validation_results$found]

  # Check if we have enough genes
  if(length(valid_genes) == 0) {
    cat("    WARNING: No valid genes found - creating empty plot\n")
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(marker_set_name, "\n\nNo valid genes found"),
               size = 6, color = "gray50") +
      theme_void()

    # Return consistent structure even for empty case
    gene_list_text <- paste(validation_results$gene_symbol, collapse = ", ")
    return(list(
      title_info = list(
        marker_set_name = marker_set_name,
        gene_list = gene_list_text,
        missing_genes = missing_genes
      ),
      module_plot = empty_plot,
      cluster_plot = empty_plot,
      gene_plots = list()
    ))
  }

  if(length(valid_genes) < length(validation_results$gene_symbol) * 0.5) {
    cat("    WARNING: <50% genes found - proceeding with caution\n")
  }

  # Get top 5 genes by expression
  top_genes <- get_top_genes(valid_genes, seurat_obj, n = 5)
  cat("    Top genes:", paste(gsub(".*-", "", top_genes), collapse = ", "), "\n")

  # Determine max expression for consistent color scaling (Seurat v5 compatible)
  max_expr <- max(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[valid_genes, ], na.rm = TRUE)

  # Plot 1: Module score
  module_plot <- FeaturePlot(seurat_obj,
                              features = module_score_name,
                              reduction = "umap.harmony") +
    scale_colour_gradient(low = "lightgrey", high = "#e31837",
                          name = "Module\nScore") +
    ggtitle(paste0(marker_set_name, " - Module Score")) +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Plot 2: DimPlot with clusters
  cluster_plot <- DimPlot(seurat_obj,
                          reduction = "umap.harmony",
                          label = TRUE,
                          label.size = 3) +
    ggtitle("Clusters") +
    theme(plot.title = element_text(size = 10, face = "bold"))

  # Plots 3-7: Top individual genes
  gene_plots <- list()
  for(i in seq_along(top_genes)) {
    gene_symbol <- gsub(".*-", "", top_genes[i])

    gene_plots[[i]] <- FeaturePlot(seurat_obj,
                                    features = top_genes[i],
                                    reduction = "umap.harmony",
                                    min.cutoff = 'q10') +
      scale_colour_gradient(low = "lightgrey", high = "#e31837",
                            limits = c(0, max_expr),
                            na.value = "#e31837",
                            name = "Expression") +
      ggtitle(gene_symbol) +
      theme(plot.title = element_text(size = 10, face = "bold"))
  }

  # Prepare title information
  gene_list_text <- paste(validation_results$gene_symbol, collapse = ", ")

  # Return individual plot components instead of combined plot
  # This avoids patchwork layout conflicts in PDF generation
  return(list(
    title_info = list(
      marker_set_name = marker_set_name,
      gene_list = gene_list_text,
      missing_genes = missing_genes
    ),
    module_plot = module_plot,
    cluster_plot = cluster_plot,
    gene_plots = gene_plots
  ))
}

## FUNCTION 5: Create summary page------
create_summary_page <- function(timepoint_name, seurat_obj, all_validation_results) {

  cat("\n=== Creating summary page for", timepoint_name, "===\n")

  # Collect all missing genes across all sets
  missing_summary <- data.frame(
    marker_set = character(),
    gene_symbol = character(),
    stringsAsFactors = FALSE
  )

  for(set_name in names(all_validation_results)) {
    validation <- all_validation_results[[set_name]]
    missing_genes <- validation$gene_symbol[!validation$found]

    if(length(missing_genes) > 0) {
      missing_summary <- rbind(
        missing_summary,
        data.frame(
          marker_set = rep(set_name, length(missing_genes)),
          gene_symbol = missing_genes,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Create summary statistics
  total_markers <- length(marker_sets)
  total_genes_tested <- sum(sapply(all_validation_results, nrow))
  total_genes_found <- sum(sapply(all_validation_results, function(x) sum(x$found)))
  total_genes_missing <- total_genes_tested - total_genes_found
  pct_found <- round(100 * total_genes_found / total_genes_tested, 1)

  # Create summary text plot
  summary_text <- paste0(
    "MARKER GENE VALIDATION SUMMARY\n",
    "Timepoint: ", timepoint_name, "\n\n",
    "Dataset Information:\n",
    "  Total Cells: ", ncol(seurat_obj), "\n",
    "  Number of Clusters: ", length(unique(seurat_obj$seurat_clusters)), "\n",
    "  Optimal Resolution: ", unique(seurat_obj$seurat_clusters) %>% length(), " clusters\n\n",
    "Validation Statistics:\n",
    "  Marker Sets Tested: ", total_markers, "\n",
    "  Total Genes Tested: ", total_genes_tested, "\n",
    "  Genes Found: ", total_genes_found, " (", pct_found, "%)\n",
    "  Genes Missing: ", total_genes_missing, "\n"
  )

  summary_plot <- ggplot() +
    annotate("text", x = 0, y = 1,
             label = summary_text,
             hjust = 0, vjust = 1,
             size = 5, family = "mono") +
    theme_void() +
    xlim(0, 1) + ylim(0, 1)

  # Create missing genes table plot
  if(nrow(missing_summary) > 0) {
    missing_table <- gridExtra::tableGrob(missing_summary, rows = NULL)

    combined <- grid.arrange(
      summary_plot,
      missing_table,
      nrow = 2,
      heights = c(1, 1.5),
      top = textGrob(paste(timepoint_name, "- Marker Validation Summary"),
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )
  } else {
    combined <- grid.arrange(
      summary_plot,
      nrow = 1,
      top = textGrob(paste(timepoint_name, "- Marker Validation Summary"),
                     gp = gpar(fontsize = 16, fontface = "bold"))
    )
  }

  return(combined)
}

## FUNCTION 6: Process one timepoint------
process_timepoint <- function(timepoint_name, file_path) {

  cat("\n\n####################################################\n")
  cat("PROCESSING TIMEPOINT:", timepoint_name, "\n")
  cat("####################################################\n")

  # Load Seurat object
  cat("\nLoading object from:", file_path, "\n")
  seurat_obj <- readRDS(file_path)
  cat("  Loaded:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "features\n")

  # Join layers for Seurat v5 compatibility
  cat("  Joining layers for v5 compatibility...\n")
  seurat_obj <- JoinLayers(seurat_obj)
  cat("  Layers joined successfully\n")

  # Validate all marker sets
  cat("\n=== Validating marker gene sets ===\n")
  all_validation_results <- list()

  for(set_name in names(marker_sets)) {
    all_validation_results[[set_name]] <- validate_gene_set(
      marker_sets[[set_name]],
      seurat_obj,
      set_name
    )
  }

  # Add module scores for each marker set
  cat("\n=== Adding module scores ===\n")
  module_score_names <- character()

  for(i in seq_along(marker_sets)) {
    set_name <- names(marker_sets)[i]
    validation <- all_validation_results[[set_name]]
    valid_genes <- validation$full_name[validation$found]

    if(length(valid_genes) > 0) {
      module_name <- paste0("MarkerSet_", i)
      cat("  Adding module score for", set_name, "\n")
      seurat_obj <- AddModuleScore(
        seurat_obj,
        features = list(valid_genes),
        name = module_name
      )
      module_score_names[set_name] <- paste0(module_name, "1")  # Seurat adds "1" suffix
    } else {
      cat("  Skipping module score for", set_name, "(no valid genes)\n")
      module_score_names[set_name] <- NA
    }
  }

  # Create plots for each marker set
  cat("\n=== Creating plots ===\n")
  all_plots <- list()

  for(set_name in names(marker_sets)) {
    all_plots[[set_name]] <- create_marker_plot(
      seurat_obj,
      set_name,
      all_validation_results[[set_name]],
      module_score_names[set_name],
      timepoint_name
    )
  }

  # Generate PDF
  cat("\n=== Generating PDF ===\n")
  pdf_file <- paste0("./output_files/marker_validation/",
                     timepoint_name, "_marker_validation.pdf")
  pdf(pdf_file, width = 14, height = 10)

  # Page 1: Summary
  create_summary_page(timepoint_name, seurat_obj, all_validation_results)

  # Pages 2-15: Marker set plots
  for(set_name in names(all_plots)) {
    cat("  Plotting", set_name, "\n")
    plot_data <- all_plots[[set_name]]

    # Create title text
    title_info <- plot_data$title_info
    if(length(title_info$missing_genes) > 0) {
      title_text <- paste0(
        title_info$marker_set_name, "\n",
        "Genes: ", title_info$gene_list, "\n",
        "[Missing: ", paste(title_info$missing_genes, collapse = ", "), "]"
      )
    } else {
      title_text <- paste0(
        title_info$marker_set_name, "\n",
        "Genes: ", title_info$gene_list
      )
    }

    # Create title grob
    title_grob <- textGrob(
      title_text,
      gp = gpar(fontsize = 10, fontface = "bold")
    )

    # Page 1 for this marker set: Module score only
    grid.arrange(
      plot_data$module_plot,
      ncol = 1,
      nrow = 1,
      top = title_grob
    )

    # Page 2 for this marker set: Cluster + individual genes
    if(length(plot_data$gene_plots) > 0) {
      # Layout: cluster plot on top, gene plots below
      grid.arrange(
        plot_data$cluster_plot,
        arrangeGrob(grobs = plot_data$gene_plots, ncol = 3),
        ncol = 1,
        nrow = 2,
        heights = c(1, ceiling(length(plot_data$gene_plots) / 3) * 0.8),
        top = textGrob(
          paste0(title_info$marker_set_name, " - Individual Genes"),
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      )
    } else {
      # Only cluster plot if no individual genes
      grid.arrange(
        plot_data$cluster_plot,
        ncol = 1,
        nrow = 1,
        top = textGrob(
          paste0(title_info$marker_set_name, " - Clusters"),
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      )
    }
  }

  dev.off()
  cat("\nPDF saved:", pdf_file, "\n")

  cat("\n====================================\n")
  cat("VALIDATION COMPLETE:", timepoint_name, "\n")
  cat("====================================\n")

  return(list(
    validation_results = all_validation_results,
    plots = all_plots
  ))
}

## MAIN EXECUTION: Process all timepoints------
cat("\n\n####################################################\n")
cat("STARTING MARKER VALIDATION FOR ALL TIMEPOINTS\n")
cat("####################################################\n\n")

# Define timepoint file paths
timepoint_files <- list(
  "1M_Org" = "./output_files/integrated_objects/1M_Org_integrated_harmony.rds",
  "3M_Org" = "./output_files/integrated_objects/3M_Org_integrated_harmony.rds",
  "6M_Org" = "./output_files/integrated_objects/6M_Org_integrated_harmony.rds"
)

# Process each timepoint
results <- list()

for(timepoint_name in names(timepoint_files)) {
  results[[timepoint_name]] <- process_timepoint(
    timepoint_name,
    timepoint_files[[timepoint_name]]
  )

  # Force garbage collection
  gc()
}

cat("\n\n####################################################\n")
cat("ALL MARKER VALIDATIONS COMPLETE\n")
cat("####################################################\n")
cat("\nPDFs saved in: ./output_files/marker_validation/\n")
cat("\nGenerated", length(results), "validation reports\n")
