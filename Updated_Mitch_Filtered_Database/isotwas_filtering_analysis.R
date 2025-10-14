#!/usr/bin/env Rscript
# isoTWAS Filtering Analysis and Report Generation
# Focused analysis of Sig_isoTWAS_developmental.csv
# Based on Mitchell's filtering criteria with emphasis on adjusted p-values
# 
# This script analyzes filtering parameters applied to isoTWAS data and generates
# a comprehensive PDF report with trait-stratified visualizations and summary statistics

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(knitr)
library(rmarkdown)

# Set working directory and theme
setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Updated_Mitch_Filtered_Database")
theme_set(theme_minimal())

# ==============================================================================
# DATA LOADING AND VALIDATION
# ==============================================================================

#' Load and validate isoTWAS data
#' @return Loaded dataframe with validation messages
load_isotwas_data <- function() {
  file_path <- "Sig_isoTWAS_developmental.csv"
  cat("Loading isoTWAS data from:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Load data
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Print basic info
  cat("Loaded", nrow(data), "rows and", ncol(data), "columns\n")
  cat("Column names:", paste(names(data), collapse = ", "), "\n")
  
  # Expected columns for isoTWAS
  expected_cols <- c("Gene", "Dataset", "Trait", "Method", "HGNC", "Transcript", 
                    "Chromosome", "Start", "End", "Biotype", "Z", "P", 
                    "Permutation.P", "R2", "Top.GWAS.SNP", "Top.GWAS.P", 
                    "Screen.P", "Screening.Adjusted.P", "Confirmation.P")
  
  # Check for required columns
  missing_cols <- setdiff(expected_cols, names(data))
  if (length(missing_cols) > 0) {
    cat("Warning: Missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
  }
  
  return(data)
}

#' Get unique traits from the dataset
get_unique_traits <- function(data) {
  traits <- unique(data$Trait)
  traits <- traits[!is.na(traits) & traits != "Trait"]
  return(sort(traits))
}

# ==============================================================================
# FILTERING VERIFICATION FUNCTIONS  
# ==============================================================================

#' Verify isoTWAS filtering status based on Mitchell's criteria
#' IMPORTANT: Uses Screening Adjusted P as primary p-value metric (accounts for multiple testing)
#' 
#' Filtering criteria:
#' - Screening Adjusted P < 0.05 (PRIMARY - accounts for multiple hypothesis testing)
#' - Permutation P < 0.05
#' - R2 >= 0.01  
#' - Confirmation P < 0.05
#' - Optional: GWAS SNP P analysis with 5e-8 threshold
#'
#' @param data isoTWAS dataframe
#' @param trait_filter Optional trait to filter by
verify_isotwas_filtering <- function(data, trait_filter = NULL) {
  
  # Filter by trait if specified
  if (!is.null(trait_filter)) {
    data <- data[data$Trait == trait_filter, ]
  }
  
  results <- list()
  results$total_records <- nrow(data)
  results$trait <- trait_filter
  
  # PRIMARY FILTERING CRITERIA - Using Adjusted P-values
  # This is the key change - we use Screening Adjusted P instead of raw P
  results$screening_adj_p_violations <- sum(data$Screening.Adjusted.P >= 0.05, na.rm = TRUE)
  results$permutation_p_violations <- sum(data$Permutation.P >= 0.05, na.rm = TRUE)  
  results$r2_violations <- sum(data$R2 < 0.01, na.rm = TRUE)
  results$confirmation_p_violations <- sum(data$Confirmation.P >= 0.05, na.rm = TRUE)
  
  # SECONDARY METRICS - Raw P-values for reference
  results$raw_p_violations <- sum(data$P >= 0.05, na.rm = TRUE)
  results$screen_p_violations <- sum(data$Screen.P >= 0.05, na.rm = TRUE)
  
  # GWAS SNP p-value analysis (5e-8 threshold)
  results$gwas_p_genome_wide_significant <- sum(data$Top.GWAS.P < 5e-8, na.rm = TRUE)
  results$gwas_p_not_significant <- sum(data$Top.GWAS.P >= 5e-8, na.rm = TRUE)
  
  # SUMMARY FLAGS
  results$any_primary_violations <- (results$screening_adj_p_violations > 0 | 
                                   results$permutation_p_violations > 0 | 
                                   results$r2_violations > 0 | 
                                   results$confirmation_p_violations > 0)
  
  return(results)
}

#' Generate summary statistics table for isoTWAS filtering results
generate_isotwas_summary <- function(filtering_results) {
  
  stats <- data.frame(
    Metric = character(),
    Value = character(),
    Priority = character(),
    stringsAsFactors = FALSE
  )
  
  trait_label <- ifelse(is.null(filtering_results$trait), "All Traits", filtering_results$trait)
  stats <- rbind(stats, data.frame(Metric = "Trait", Value = trait_label, Priority = "Info"))
  stats <- rbind(stats, data.frame(Metric = "Total Records", Value = filtering_results$total_records, Priority = "Info"))
  
  # PRIMARY FILTERING VIOLATIONS (based on adjusted p-values)
  stats <- rbind(stats, data.frame(Metric = "=== PRIMARY VIOLATIONS ===", Value = "", Priority = ""))
  stats <- rbind(stats, data.frame(Metric = "Screening Adj P >= 0.05", Value = filtering_results$screening_adj_p_violations, Priority = "HIGH"))
  stats <- rbind(stats, data.frame(Metric = "Permutation P >= 0.05", Value = filtering_results$permutation_p_violations, Priority = "HIGH"))
  stats <- rbind(stats, data.frame(Metric = "R2 < 0.01", Value = filtering_results$r2_violations, Priority = "HIGH"))
  stats <- rbind(stats, data.frame(Metric = "Confirmation P >= 0.05", Value = filtering_results$confirmation_p_violations, Priority = "HIGH"))
  
  # SECONDARY METRICS (reference only)
  stats <- rbind(stats, data.frame(Metric = "=== REFERENCE METRICS ===", Value = "", Priority = ""))
  stats <- rbind(stats, data.frame(Metric = "Raw P >= 0.05", Value = filtering_results$raw_p_violations, Priority = "LOW"))
  stats <- rbind(stats, data.frame(Metric = "Screen P >= 0.05", Value = filtering_results$screen_p_violations, Priority = "LOW"))
  
  # GWAS ANALYSIS
  stats <- rbind(stats, data.frame(Metric = "=== GWAS ANALYSIS ===", Value = "", Priority = ""))
  stats <- rbind(stats, data.frame(Metric = "GWAS P < 5e-8 (genome-wide)", Value = filtering_results$gwas_p_genome_wide_significant, Priority = "MED"))
  stats <- rbind(stats, data.frame(Metric = "GWAS P >= 5e-8", Value = filtering_results$gwas_p_not_significant, Priority = "MED"))
  
  return(stats)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Create histogram for p-values with log scale
#' @param data isoTWAS dataframe
#' @param column_name Name of the column to plot
#' @param title Plot title
#' @param trait_filter Optional trait to filter by
#' @param log_scale Whether to use log scale (default TRUE)
create_isotwas_histogram <- function(data, column_name, title, trait_filter = NULL, log_scale = TRUE) {
  
  # Filter by trait if specified
  if (!is.null(trait_filter)) {
    data <- data[data$Trait == trait_filter, ]
    title <- paste(title, "-", trait_filter)
  }
  
  # Remove NA values and zeros for log scale
  clean_data <- data[[column_name]]
  clean_data <- clean_data[!is.na(clean_data) & clean_data > 0]
  
  if (length(clean_data) == 0) {
    # Return empty plot if no data
    return(ggplot() + 
           geom_text(aes(x = 0.5, y = 0.5, label = "No data available"), size = 5) +
           labs(title = title) + 
           theme_void())
  }
  
  p <- ggplot(data.frame(values = clean_data), aes(x = values)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
    labs(title = title, x = column_name, y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  if (log_scale && length(clean_data) > 0) {
    p <- p + scale_x_log10() +
      labs(x = paste(column_name, "(log10 scale)"))
  }
  
  # Add significance thresholds
  if (grepl("GWAS", column_name)) {
    p <- p + geom_vline(xintercept = 5e-8, color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = 5e-8, y = Inf, label = "5e-8", 
               vjust = 1.5, hjust = -0.1, color = "red", size = 3)
  } else {
    p <- p + geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = 0.05, y = Inf, label = "0.05", 
               vjust = 1.5, hjust = 1.1, color = "red", size = 3)
  }
  
  return(p)
}

#' Create all histograms for isoTWAS data, optionally filtered by trait
#' @param data isoTWAS dataframe
#' @param trait_filter Optional trait to filter by
create_all_isotwas_histograms <- function(data, trait_filter = NULL) {
  
  # Filter data by trait if specified
  if (!is.null(trait_filter)) {
    data <- data[data$Trait == trait_filter, ]
  }
  
  plots <- list()
  trait_suffix <- ifelse(is.null(trait_filter), "", paste("-", trait_filter))
  
  # PRIMARY P-VALUE DISTRIBUTIONS (adjusted for multiple testing)
  plots[["Screening_Adj_P"]] <- create_isotwas_histogram(data, "Screening.Adjusted.P", 
                                                        paste("Screening Adjusted P", trait_suffix), trait_filter)
  plots[["Permutation_P"]] <- create_isotwas_histogram(data, "Permutation.P", 
                                                      paste("Permutation P", trait_suffix), trait_filter)
  plots[["Confirmation_P"]] <- create_isotwas_histogram(data, "Confirmation.P", 
                                                       paste("Confirmation P", trait_suffix), trait_filter)
  
  # SECONDARY P-VALUE DISTRIBUTIONS (reference only)  
  plots[["Raw_P"]] <- create_isotwas_histogram(data, "P", 
                                              paste("Raw P (reference)", trait_suffix), trait_filter)
  plots[["Screen_P"]] <- create_isotwas_histogram(data, "Screen.P", 
                                                 paste("Screen P (unadjusted)", trait_suffix), trait_filter)
  
  # GWAS SNP P-VALUES
  plots[["GWAS_P"]] <- create_isotwas_histogram(data, "Top.GWAS.P", 
                                               paste("GWAS SNP P", trait_suffix), trait_filter)
  
  # R2 DISTRIBUTION (not log scale)
  r2_data <- data$R2
  r2_data <- r2_data[!is.na(r2_data)]
  
  if (length(r2_data) > 0) {
    plots[["R2"]] <- ggplot(data.frame(values = r2_data), aes(x = values)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
      labs(title = paste("R2 Distribution", trait_suffix), x = "R2", y = "Count") +
      geom_vline(xintercept = 0.01, color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = 0.01, y = Inf, label = "0.01", 
               vjust = 1.5, hjust = -0.1, color = "red", size = 3) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
  } else {
    plots[["R2"]] <- ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "No R2 data available"), size = 5) +
      labs(title = paste("R2 Distribution", trait_suffix)) + 
      theme_void()
  }
  
  return(plots)
}

#' Create trait distribution plot for isoTWAS data
create_isotwas_trait_distribution <- function(data) {
  
  trait_counts <- data %>%
    count(Trait, sort = TRUE) %>%
    filter(Trait != "Trait")  # Remove header artifact if present
  
  ggplot(trait_counts, aes(x = reorder(Trait, n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "isoTWAS Records by Disease/Trait", 
         x = "Trait", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14))
}

#' Create biotype distribution plot for isoTWAS data
#' @param data isoTWAS dataframe
#' @param title_prefix Prefix for plot title
create_isotwas_biotype_distribution <- function(data, title_prefix = "") {
  
  # Handle potential missing or empty biotype values
  data_clean <- data %>%
    filter(!is.na(Biotype) & Biotype != "" & Biotype != "Biotype") %>%
    count(Biotype, sort = TRUE)
  
  if (nrow(data_clean) == 0) {
    return(ggplot() + 
           geom_text(aes(x = 0.5, y = 0.5, label = "No biotype data available"), size = 5) +
           labs(title = paste(title_prefix, "Biotype Distribution")) + 
           theme_void())
  }
  
  # Create horizontal bar plot
  p <- ggplot(data_clean, aes(x = reorder(Biotype, n), y = n)) +
    geom_bar(stat = "identity", fill = "darkorange", alpha = 0.7) +
    coord_flip() +
    labs(title = paste(title_prefix, "Biotype Distribution"), 
         x = "Biotype", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12),
          axis.text.y = element_text(size = 9))
  
  # Add count labels on bars
  p <- p + geom_text(aes(label = n), hjust = -0.1, size = 3)
  
  return(p)
}

#' Create biotype distribution comparison (before/after GWAS filtering)
#' @param original_data Original isoTWAS data
#' @param gwas_filtered_data GWAS filtered data
create_biotype_comparison <- function(original_data, gwas_filtered_data) {
  
  # Prepare data for comparison
  original_counts <- original_data %>%
    filter(!is.na(Biotype) & Biotype != "" & Biotype != "Biotype") %>%
    count(Biotype, sort = TRUE) %>%
    mutate(Dataset = "Original")
  
  gwas_counts <- gwas_filtered_data %>%
    filter(!is.na(Biotype) & Biotype != "" & Biotype != "Biotype") %>%
    count(Biotype, sort = TRUE) %>%
    mutate(Dataset = "GWAS Filtered (p<5e-8)")
  
  # Combine datasets
  combined_data <- bind_rows(original_counts, gwas_counts) %>%
    replace_na(list(n = 0))
  
  # Create grouped bar plot
  ggplot(combined_data, aes(x = reorder(Biotype, n), y = n, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    coord_flip() +
    labs(title = "Biotype Distribution: Before vs After GWAS Filtering", 
         x = "Biotype", y = "Count",
         fill = "Dataset") +
    scale_fill_manual(values = c("Original" = "steelblue", "GWAS Filtered (p<5e-8)" = "darkorange")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 12),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom")
}

# ==============================================================================
# PDF REPORT GENERATION
# ==============================================================================

#' Create a trait-specific page in the PDF report
#' @param data isoTWAS dataset
#' @param trait_name Name of the trait
create_isotwas_trait_page <- function(data, trait_name) {
  
  # Get filtering results for this trait
  trait_filtering <- verify_isotwas_filtering(data, trait_name)
  trait_stats <- generate_isotwas_summary(trait_filtering)
  
  # Create histograms for this trait
  trait_plots <- create_all_isotwas_histograms(data, trait_name)
  
  # Create new page
  grid.newpage()
  
  # Title for the trait page
  grid.text(paste("isoTWAS Analysis -", trait_name), 
            x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontweight = "bold"))
  
  # Summary statistics table (left side)
  pushViewport(viewport(x = 0.2, y = 0.7, width = 0.35, height = 0.45))
  grid.text("Filtering Summary", y = 0.95, gp = gpar(fontsize = 12, fontweight = "bold"))
  
  # Create a simplified summary table for display
  display_stats <- trait_stats[trait_stats$Priority %in% c("Info", "HIGH"), c("Metric", "Value")]
  grid.table(display_stats, rows = NULL)
  popViewport()
  
  # Arrange plots in a structured grid (4x2 layout)
  plot_positions <- list(
    list(x = 0.65, y = 0.85, width = 0.25, height = 0.2),   # Screening Adj P (most important)
    list(x = 0.65, y = 0.6, width = 0.25, height = 0.2),    # Permutation P
    list(x = 0.65, y = 0.35, width = 0.25, height = 0.2),   # Confirmation P
    list(x = 0.65, y = 0.1, width = 0.25, height = 0.2),    # GWAS P
    list(x = 0.2, y = 0.35, width = 0.25, height = 0.2),    # R2
    list(x = 0.2, y = 0.1, width = 0.25, height = 0.2),     # Raw P (reference)
    list(x = 0.42, y = 0.1, width = 0.18, height = 0.15)    # Screen P (reference, smaller)
  )
  
  # Plot order prioritizes adjusted p-values
  plot_order <- c("Screening_Adj_P", "Permutation_P", "Confirmation_P", "GWAS_P", "R2", "Raw_P", "Screen_P")
  
  # Draw each plot
  for (i in seq_along(plot_order)) {
    plot_name <- plot_order[i]
    if (plot_name %in% names(trait_plots) && i <= length(plot_positions)) {
      pos <- plot_positions[[i]]
      pushViewport(viewport(x = pos$x, y = pos$y, width = pos$width, height = pos$height))
      grid.draw(ggplotGrob(trait_plots[[plot_name]]))
      popViewport()
    }
  }
}

#' Generate comprehensive PDF report for isoTWAS filtering analysis
generate_isotwas_report <- function() {
  
  # Load data
  isotwas_data <- load_isotwas_data()
  
  # Get unique traits
  traits <- get_unique_traits(isotwas_data)
  cat("Traits found:", paste(traits, collapse = ", "), "\n")
  
  # Overall filtering analysis
  overall_filtering <- verify_isotwas_filtering(isotwas_data) 
  overall_stats <- generate_isotwas_summary(overall_filtering)
  
  # Create GWAS-filtered dataset (p < 5e-8)
  cat("Creating GWAS-filtered dataset (p < 5e-8)...\n")
  gwas_filtered_data <- isotwas_data %>%
    filter(Top.GWAS.P < 5e-8)
  
  cat("Original records:", nrow(isotwas_data), "\n")
  cat("GWAS-filtered records:", nrow(gwas_filtered_data), "\n")
  cat("Filtering reduction:", round((1 - nrow(gwas_filtered_data)/nrow(isotwas_data)) * 100, 1), "%\n")
  
  # Save GWAS-filtered dataset as CSV
  output_filename <- "GWAS_filtered_5e-8_Sig_isoTWAS_developmental.csv"
  write.csv(gwas_filtered_data, output_filename, row.names = FALSE)
  cat("GWAS-filtered data saved as:", output_filename, "\n")
  
  # Generate PDF
  pdf("isoTWAS_Filtering_Analysis_By_Trait.pdf", width = 16, height = 12)
  
  # Title page
  grid.newpage()
  grid.text("isoTWAS Filtering Analysis Report", 
            x = 0.5, y = 0.8, gp = gpar(fontsize = 28, fontweight = "bold"))
  grid.text("Analysis of Mitchell's Filtering Parameters", 
            x = 0.5, y = 0.7, gp = gpar(fontsize = 20))
  grid.text("Focus: Sig_isoTWAS_developmental.csv", 
            x = 0.5, y = 0.65, gp = gpar(fontsize = 18))
  grid.text("PRIMARY METRIC: Screening Adjusted P-values", 
            x = 0.5, y = 0.6, gp = gpar(fontsize = 16, col = "red"))
  grid.text("(Accounts for Multiple Hypothesis Testing)", 
            x = 0.5, y = 0.55, gp = gpar(fontsize = 14, col = "red"))
  grid.text(paste("Generated on:", Sys.Date()), 
            x = 0.5, y = 0.4, gp = gpar(fontsize = 14))
  
  # Overall summary page
  grid.newpage()
  grid.text("Overall isoTWAS Summary Statistics", 
            x = 0.5, y = 0.9, gp = gpar(fontsize = 18, fontweight = "bold"))
  grid.table(overall_stats[, c("Metric", "Value")], rows = NULL)
  
  # Trait distribution plot
  trait_plot <- create_isotwas_trait_distribution(isotwas_data)
  print(trait_plot)
  
  # Biotype distribution plots
  biotype_original <- create_isotwas_biotype_distribution(isotwas_data, "Original Dataset:")
  print(biotype_original)
  
  biotype_gwas_filtered <- create_isotwas_biotype_distribution(gwas_filtered_data, "GWAS Filtered (p<5e-8):")
  print(biotype_gwas_filtered)
  
  # Biotype comparison plot
  biotype_comparison <- create_biotype_comparison(isotwas_data, gwas_filtered_data)
  print(biotype_comparison)
  
  # Create separate pages for each trait
  for (trait in traits) {
    cat("Creating page for trait:", trait, "\n")
    create_isotwas_trait_page(isotwas_data, trait)
  }
  
  dev.off()
  
  cat("PDF report generated: isoTWAS_Filtering_Analysis_By_Trait.pdf\n")
  
  # Return results for console output
  return(list(
    overall_stats = overall_stats,
    overall_filtering = overall_filtering,
    traits = traits,
    total_records = nrow(isotwas_data),
    gwas_filtered_records = nrow(gwas_filtered_data),
    gwas_filtered_filename = output_filename
  ))
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main <- function() {
  cat("Starting isoTWAS Filtering Analysis...\n")
  cat("=====================================\n")
  cat("Focus: Sig_isoTWAS_developmental.csv\n")
  cat("Primary Metric: Screening Adjusted P-values (multiple testing corrected)\n\n")
  
  # Generate comprehensive report
  results <- generate_isotwas_report()
  
  cat("\n=====================================\n")
  cat("Analysis Complete!\n")
  cat("Key Findings:\n")
  
  cat("\nOverall isoTWAS Analysis Summary:\n")
  cat("- Total records:", results$total_records, "\n")
  cat("- Traits analyzed:", length(results$traits), "\n")
  cat("- Traits:", paste(results$traits, collapse = ", "), "\n")
  
  cat("\nPRIMARY FILTERING VIOLATIONS (Multiple Testing Corrected):\n")
  cat("- Screening Adjusted P >= 0.05:", results$overall_filtering$screening_adj_p_violations, "\n")
  cat("- Permutation P >= 0.05:", results$overall_filtering$permutation_p_violations, "\n") 
  cat("- R2 < 0.01:", results$overall_filtering$r2_violations, "\n")
  cat("- Confirmation P >= 0.05:", results$overall_filtering$confirmation_p_violations, "\n")
  
  cat("\nREFERENCE METRICS (Raw P-values):\n")
  cat("- Raw P >= 0.05:", results$overall_filtering$raw_p_violations, "\n")
  cat("- Screen P >= 0.05:", results$overall_filtering$screen_p_violations, "\n")
  
  cat("\nGWAS SNP ANALYSIS:\n")
  cat("- GWAS P < 5e-8 (genome-wide sig):", results$overall_filtering$gwas_p_genome_wide_significant, "\n")
  cat("- GWAS P >= 5e-8:", results$overall_filtering$gwas_p_not_significant, "\n")
  
  if (results$overall_filtering$screening_adj_p_violations > 0) {
    cat("\n*** IMPORTANT: Found Screening Adjusted P violations! ***\n")
    cat("These represent failures in the multiple testing corrected analysis.\n")
    cat("Mitchell should investigate these ", results$overall_filtering$screening_adj_p_violations, " records.\n")
  } else {
    cat("\n✓ All records pass Screening Adjusted P < 0.05 threshold!\n")
  }
  
  cat("\nGWAS FILTERING RESULTS:\n")
  cat("- Original records:", results$total_records, "\n")
  cat("- GWAS-filtered records (p<5e-8):", results$gwas_filtered_records, "\n")
  cat("- Reduction:", round((1 - results$gwas_filtered_records/results$total_records) * 100, 1), "% of records removed\n")
  cat("- Filtered data saved as:", results$gwas_filtered_filename, "\n")
  
  cat("\nThe PDF report 'isoTWAS_Filtering_Analysis_By_Trait.pdf' contains:\n")
  cat("- Overall summary with emphasis on adjusted p-values\n")
  cat("- Trait distribution analysis\n")
  cat("- Biotype distribution analysis (before/after GWAS filtering)\n")
  cat("- Biotype comparison plots\n") 
  cat("- Separate detailed pages for each of", length(results$traits), "traits\n")
  cat("- Prioritized histograms (adjusted p-values first)\n")
  cat("- Comprehensive filtering violation summaries\n")
}

#' Generate trait-specific filtered datasets and summaries
#' @param data Original isoTWAS dataset
#' @param trait_name Name of trait to filter (SCZ, MDD, BP)
generate_trait_filtered_data <- function(data, trait_name) {
  
  cat("=====================================\n")
  cat("Generating trait-specific analysis for:", trait_name, "\n")
  cat("=====================================\n")
  
  # Filter data for specific trait AND GWAS significance (p < 5e-8)
  trait_data <- data %>%
    filter(Trait == trait_name & Top.GWAS.P < 5e-8)
  
  cat("Records for", trait_name, ":", nrow(trait_data), "\n")
  
  if (nrow(trait_data) == 0) {
    cat("No records found for trait:", trait_name, "\n")
    return(NULL)
  }
  
  # Generate filtering analysis for this trait
  trait_filtering <- verify_isotwas_filtering(trait_data)
  trait_stats <- generate_isotwas_summary(trait_filtering)
  
  # Print summary statistics
  cat("\nFiltering Summary for", trait_name, ":\n")
  cat("- Total records:", trait_filtering$total_records, "\n")
  cat("- Screening Adjusted P >= 0.05:", trait_filtering$screening_adj_p_violations, "\n")
  cat("- Permutation P >= 0.05:", trait_filtering$permutation_p_violations, "\n")
  cat("- R2 < 0.01:", trait_filtering$r2_violations, "\n")
  cat("- Confirmation P >= 0.05:", trait_filtering$confirmation_p_violations, "\n")
  cat("- GWAS P < 5e-8 (genome-wide sig):", trait_filtering$gwas_p_genome_wide_significant, "\n")
  cat("- GWAS P >= 5e-8:", trait_filtering$gwas_p_not_significant, "\n")
  
  # Save trait-specific CSV file
  output_filename <- paste0(trait_name, "_filtered_5e-8_Sig_isoTWAS_developmental.csv")
  write.csv(trait_data, output_filename, row.names = FALSE)
  cat("Trait-specific data saved as:", output_filename, "\n")
  
  return(list(
    trait_name = trait_name,
    data = trait_data,
    filtering_results = trait_filtering,
    summary_stats = trait_stats,
    output_filename = output_filename
  ))
}

# Run main function
main()

# ==============================================================================
# TRAIT-SPECIFIC SUBFILTERING (SCZ, MDD, BP)
# ==============================================================================

cat("\n\n=====================================\n")
cat("TRAIT-SPECIFIC SUBFILTERING ANALYSIS\n")
cat("=====================================\n")
cat("Generating trait-specific filtered datasets for SCZ, MDD, BP, and ASD...\n")

# Load the original data again for trait-specific analysis
isotwas_data_for_traits <- load_isotwas_data()

# Define the traits of interest
target_traits <- c("SCZ", "MDD", "BP", "ASD")

# Generate trait-specific analyses
trait_results <- list()
for (trait in target_traits) {
  trait_result <- generate_trait_filtered_data(isotwas_data_for_traits, trait)
  if (!is.null(trait_result)) {
    trait_results[[trait]] <- trait_result
  }
}

# Summary of trait-specific results
cat("\n=====================================\n")
cat("TRAIT-SPECIFIC SUBFILTERING SUMMARY\n")
cat("=====================================\n")

for (trait in names(trait_results)) {
  result <- trait_results[[trait]]
  cat("\n", trait, "Summary:\n")
  cat("- Records:", result$filtering_results$total_records, "\n")
  cat("- Primary violations (Screening Adj P >= 0.05):", result$filtering_results$screening_adj_p_violations, "\n")
  cat("- GWAS genome-wide significant (p < 5e-8):", result$filtering_results$gwas_p_genome_wide_significant, "\n")
  cat("- Output file:", result$output_filename, "\n")
}

cat("\nTrait-specific filtering complete!\n")
cat("Generated", length(trait_results), "trait-specific CSV files with the naming structure:\n")
cat("- [TRAIT]_filtered_5e-8_Sig_isoTWAS_developmental.csv\n")
