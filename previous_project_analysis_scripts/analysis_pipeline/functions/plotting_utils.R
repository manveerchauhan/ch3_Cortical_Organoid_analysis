# =============================================================================
# Publication-Ready Plotting Utilities with Grid/GridExtra
# =============================================================================
# Author: Manveer Chauhan and Sefi Prawer
# Description: Comprehensive plotting functions using grid and gridExtra
#              for publication-ready figures with flexible parameters
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(cowplot)
  library(RColorBrewer)
  library(viridis)
  library(ggrepel)
  library(scales)
})

# =============================================================================
# THEME AND STYLING FUNCTIONS
# =============================================================================

#' Create publication-ready ggplot theme
#'
#' @param base_size Base font size
#' @return ggplot theme object
theme_publication <- function(base_size = DEFAULT_FONT_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(size = TITLE_FONT_SIZE, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = AXIS_FONT_SIZE),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = AXIS_FONT_SIZE),
      strip.text = element_text(size = AXIS_FONT_SIZE, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.8),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

#' Get color palette for different data types
#'
#' @param type Type of data ("timepoint", "qc", "continuous", "discrete")
#' @param n Number of colors needed
#' @return Vector of colors
get_color_palette <- function(type = "timepoint", n = NULL) {
  switch(type,
    "timepoint" = TIMEPOINT_COLORS,
    "qc" = QC_COLORS,
    "continuous" = viridis::viridis(n %||% 10),
    "discrete" = RColorBrewer::brewer.pal(min(n %||% 8, 8), "Set2"),
    "diverging" = RColorBrewer::brewer.pal(min(n %||% 11, 11), "RdBu")
  )
}

# =============================================================================
# FIGURE SAVING FUNCTIONS
# =============================================================================

#' Save individual figure with flexible parameters
#'
#' @param plot ggplot object or grob
#' @param filename Filename (without extension)
#' @param type Figure type ("single", "multipanel")
#' @param format File format ("png", "pdf", "both")
#' @param subdir Subdirectory within figures folder
save_figure <- function(plot, filename, type = "single", format = "both", subdir = "individual") {
  # Get figure parameters
  params <- get_figure_params(type)

  # Create output paths
  fig_dir <- file.path(OUTPUT_BASE_DIR, "figures", subdir)
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  # Save PNG
  if (format %in% c("png", "both")) {
    png_path <- file.path(fig_dir, paste0(filename, ".png"))
    ggsave(png_path, plot,
           width = params$width, height = params$height,
           dpi = params$dpi, units = "in")
    log_message(sprintf("Saved PNG: %s", basename(png_path)), "DEBUG")
  }

  # Save PDF
  if (format %in% c("pdf", "both")) {
    pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
    ggsave(pdf_path, plot,
           width = params$width, height = params$height,
           units = "in", device = "pdf")
    log_message(sprintf("Saved PDF: %s", basename(pdf_path)), "DEBUG")
  }

  log_message(sprintf("Figure saved: %s", filename), "INFO")
}

#' Create and save multi-panel figure using gridExtra
#'
#' @param plot_list List of ggplot objects
#' @param filename Output filename
#' @param ncol Number of columns
#' @param nrow Number of rows (auto-calculated if NULL)
#' @param title Overall title for the figure
#' @param subdir Subdirectory within figures folder
create_multipanel_figure <- function(plot_list, filename, ncol = 2, nrow = NULL,
                                     title = NULL, subdir = "individual") {
  # Calculate rows if not specified
  if (is.null(nrow)) {
    nrow <- ceiling(length(plot_list) / ncol)
  }

  # Get figure parameters for multipanel
  params <- get_figure_params("multipanel")

  # Create output path
  fig_dir <- file.path(OUTPUT_BASE_DIR, "figures", subdir)
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  # Create multipanel figure
  pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))

  pdf(pdf_path, width = params$width, height = params$height)

  if (!is.null(title)) {
    # Create title grob
    title_grob <- textGrob(title, gp = gpar(fontsize = TITLE_FONT_SIZE, fontface = "bold"))

    # Arrange plots with title
    grid.arrange(
      title_grob,
      arrangeGrob(grobs = plot_list, ncol = ncol, nrow = nrow),
      heights = c(0.1, 0.9)
    )
  } else {
    # Arrange plots without title
    do.call(grid.arrange, c(plot_list, ncol = ncol, nrow = nrow))
  }

  dev.off()

  log_message(sprintf("Multipanel figure saved: %s", filename), "INFO")
  return(pdf_path)
}

# =============================================================================
# QC PLOTTING FUNCTIONS
# =============================================================================

#' Create QC metrics histogram with statistics
#'
#' @param data Vector of values
#' @param title Plot title
#' @param xlabel X-axis label
#' @param bins Number of histogram bins
#' @param filtered Is this filtered data? (default: FALSE)
#' @return ggplot object
plot_qc_histogram <- function(data, title, xlabel, bins = 50, filtered = FALSE) {
  # Calculate statistics
  data_median <- median(data, na.rm = TRUE)
  data_mean <- mean(data, na.rm = TRUE)
  data_q25 <- quantile(data, 0.25, na.rm = TRUE)
  data_q75 <- quantile(data, 0.75, na.rm = TRUE)
  n_cells <- length(data)

  # Create statistics label
  stats_label <- sprintf(
    "n = %d\nMedian = %.0f\nMean = %.0f\nQ25-Q75 = %.0f-%.0f",
    n_cells, data_median, data_mean, data_q25, data_q75
  )

  # Add filtered/unfiltered label to title
  full_title <- if (filtered) paste(title, "(Filtered)") else paste(title, "(Unfiltered)")

  p <- ggplot(data.frame(values = data), aes(x = values)) +
    geom_histogram(bins = bins, fill = "lightblue", alpha = 0.7, color = "black") +
    labs(title = full_title, x = xlabel, y = "Number of Cells") +
    theme_publication() +
    geom_vline(xintercept = data_median,
               color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = Inf, y = Inf,
             label = stats_label,
             hjust = 1.1, vjust = 1.5,
             size = 3, fontface = "bold",
             color = "darkblue")

  return(p)
}

#' Create QC scatter plot with statistics
#'
#' @param data Data frame with x and y columns
#' @param x_col Column name for x-axis
#' @param y_col Column name for y-axis
#' @param color_col Column name for coloring (optional)
#' @param title Plot title
#' @param filtered Is this filtered data? (default: FALSE)
#' @return ggplot object
plot_qc_scatter <- function(data, x_col, y_col, color_col = NULL, title, filtered = FALSE) {
  # Calculate statistics
  n_cells <- nrow(data)
  cor_val <- cor(data[[x_col]], data[[y_col]], use = "complete.obs")

  # Create statistics label
  stats_label <- sprintf("n = %d\nCorr = %.2f", n_cells, cor_val)

  # Add filtered/unfiltered label to title
  full_title <- if (filtered) paste(title, "(Filtered)") else paste(title, "(Unfiltered)")

  p <- ggplot(data, aes_string(x = x_col, y = y_col))

  if (!is.null(color_col)) {
    p <- p + geom_point(aes_string(color = color_col), alpha = 0.6)
  } else {
    p <- p + geom_point(alpha = 0.6, color = "steelblue")
  }

  p <- p +
    labs(title = full_title, x = gsub("_", " ", str_to_title(x_col)),
         y = gsub("_", " ", str_to_title(y_col))) +
    theme_publication() +
    annotate("text", x = Inf, y = Inf,
             label = stats_label,
             hjust = 1.1, vjust = 1.5,
             size = 3, fontface = "bold",
             color = "darkblue")

  return(p)
}

#' Create violin plot for QC metrics by group with statistics
#'
#' @param data Data frame
#' @param x_col Grouping column
#' @param y_col Value column
#' @param title Plot title
#' @param filtered Is this filtered data? (default: FALSE)
#' @return ggplot object
plot_qc_violin <- function(data, x_col, y_col, title, filtered = FALSE) {
  # Calculate statistics by group
  stats_by_group <- data %>%
    group_by(.data[[x_col]]) %>%
    summarise(
      n = n(),
      median = median(.data[[y_col]], na.rm = TRUE),
      mean = mean(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # Create overall statistics label
  stats_label <- paste(
    sapply(1:nrow(stats_by_group), function(i) {
      sprintf("%s: n=%d, med=%.0f",
              stats_by_group[[x_col]][i],
              stats_by_group$n[i],
              stats_by_group$median[i])
    }),
    collapse = "\n"
  )

  # Add filtered/unfiltered label to title
  full_title <- if (filtered) paste(title, "(Filtered)") else paste(title, "(Unfiltered)")

  p <- ggplot(data, aes_string(x = x_col, y = y_col, fill = x_col)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
    labs(title = full_title, x = gsub("_", " ", str_to_title(x_col)),
         y = gsub("_", " ", str_to_title(y_col))) +
    theme_publication() +
    theme(legend.position = "none") +
    scale_fill_manual(values = get_color_palette("discrete", length(unique(data[[x_col]])))) +
    annotate("text", x = Inf, y = Inf,
             label = stats_label,
             hjust = 1.1, vjust = 1.2,
             size = 2.5, fontface = "bold",
             color = "darkblue")

  return(p)
}

# =============================================================================
# INTEGRATION PLOTTING FUNCTIONS
# =============================================================================

#' Create UMAP plot with flexible coloring
#'
#' @param data Data frame with UMAP coordinates
#' @param color_by Column to color by
#' @param title Plot title
#' @param point_size Point size
#' @return ggplot object
plot_umap <- function(data, color_by, title, point_size = 0.5) {
  p <- ggplot(data, aes_string(x = "UMAP_1", y = "UMAP_2", color = color_by)) +
    geom_point(size = point_size, alpha = 0.7) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_publication() +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # Use appropriate color scale
  if (is.numeric(data[[color_by]])) {
    p <- p + scale_color_viridis_c()
  } else if (color_by == "timepoint") {
    p <- p + scale_color_manual(values = TIMEPOINT_COLORS)
  } else {
    n_groups <- length(unique(data[[color_by]]))
    p <- p + scale_color_manual(values = get_color_palette("discrete", n_groups))
  }

  return(p)
}

#' Create feature plot for gene/transcript expression
#'
#' @param data Data frame with coordinates and expression
#' @param feature_col Column with expression values
#' @param title Plot title
#' @return ggplot object
plot_feature <- function(data, feature_col, title) {
  p <- ggplot(data, aes_string(x = "UMAP_1", y = "UMAP_2", color = feature_col)) +
    geom_point(size = 0.5, alpha = 0.7) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2",
         color = "Expression") +
    theme_publication() +
    scale_color_viridis_c()

  return(p)
}

# =============================================================================
# STATISTICAL PLOTTING FUNCTIONS
# =============================================================================

#' Create volcano plot for differential expression
#'
#' @param de_results Data frame with DE results
#' @param title Plot title
#' @param pval_cutoff P-value cutoff for significance
#' @param fc_cutoff Fold change cutoff
#' @return ggplot object
plot_volcano <- function(de_results, title, pval_cutoff = 0.05, fc_cutoff = 0.25) {
  # Add significance categories
  de_results$significance <- "Not Significant"
  de_results$significance[de_results$p_val_adj < pval_cutoff & abs(de_results$avg_log2FC) > fc_cutoff] <- "Significant"
  de_results$significance[de_results$p_val_adj < pval_cutoff & de_results$avg_log2FC > fc_cutoff] <- "Up-regulated"
  de_results$significance[de_results$p_val_adj < pval_cutoff & de_results$avg_log2FC < -fc_cutoff] <- "Down-regulated"

  p <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.6) +
    labs(title = title, x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    theme_publication() +
    scale_color_manual(values = c("Not Significant" = "grey",
                                  "Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "Significant" = "darkgreen")) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black")

  return(p)
}

# =============================================================================
# HEATMAP FUNCTIONS
# =============================================================================

#' Create expression heatmap
#'
#' @param expr_matrix Expression matrix (genes x samples)
#' @param title Plot title
#' @param cluster_rows Should rows be clustered?
#' @param cluster_cols Should columns be clustered?
#' @return ggplot object
plot_expression_heatmap <- function(expr_matrix, title, cluster_rows = TRUE, cluster_cols = TRUE) {
  # Convert matrix to long format for ggplot
  expr_df <- reshape2::melt(as.matrix(expr_matrix))
  colnames(expr_df) <- c("Gene", "Sample", "Expression")

  p <- ggplot(expr_df, aes(x = Sample, y = Gene, fill = Expression)) +
    geom_tile() +
    labs(title = title, x = "Sample", y = "Gene") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6)) +
    scale_fill_viridis_c(name = "Expression")

  return(p)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Add statistical annotations to plots
#'
#' @param plot ggplot object
#' @param data Data for statistical test
#' @param x_col X variable
#' @param y_col Y variable
#' @param test Statistical test to use
#' @return ggplot object with annotations
add_statistical_annotations <- function(plot, data, x_col, y_col, test = "wilcox") {
  # This is a placeholder for statistical annotation functionality
  # Could be implemented with ggpubr::stat_compare_means() if needed
  return(plot)
}

#' Create publication-ready figure caption
#'
#' @param title Main title
#' @param description Figure description
#' @param methods Brief methods description
#' @return Formatted caption string
create_figure_caption <- function(title, description, methods = NULL) {
  caption <- paste0("Figure: ", title, ". ", description)
  if (!is.null(methods)) {
    caption <- paste0(caption, " Methods: ", methods)
  }
  return(caption)
}

log_message("Plotting utilities loaded successfully", "INFO")