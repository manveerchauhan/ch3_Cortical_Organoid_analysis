# ============================================================================
# Single-Cell Isoform Expression Analysis
# Neurodevelopmental Models Analysis
# ============================================================================
# Author: [Manveer Chauhan]
# Date: [20.06.25]
# Description: Analysis of isoform-level expression patterns in brain organoids
#              across developmental timepoints (1, 3, and 6 months)
# ============================================================================

# REQUIRED LIBRARIES --------------------------------------------------------
library(Seurat)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)  # For plot_annotation function

# GLOBAL PARAMETERS ----------------------------------------------------------
# Directory paths
BASE_DIR <- "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis"
setwd(BASE_DIR)

# Input file paths
ONE_MONTH_FILE <- file.path(BASE_DIR, "one_month_seurat.intergrated_harm.isofrom.rds")
THREE_MONTH_FILE <- file.path(BASE_DIR, "three_month_seurat.intergrated_harm.isofrom.rds")
SIX_MONTH_FILE <- file.path(BASE_DIR, "six_month_seurat.intergrated_harm.isofrom.rds")

# Analysis parameters
IDENTITY_COLUMN <- "cluster_annotations"
DEFAULT_ASSAY <- "iso"
DEFAULT_REDUCTION <- "umap.harm"
DEFAULT_SLOT <- "counts"
DE_TEST_METHOD <- "wilcox"

# Features of interest for splicing analysis
SPLICING_FACTORS <- c("PTBP1", "PTBP2", "RBFOX1", "RBFOX2", "RBFOX3", "NOVA1", "NOVA2")
SYNAPTIC_GENES <- c("DLG4", "NRXN1", "NRXN2", "NRXN3", 
                    "NLGN1", "NLGN2", "NLGN3", "NLGN4X",
                    "NLGN4Y")
ALL_FEATURES <- c(SPLICING_FACTORS, SYNAPTIC_GENES)

# Plotting parameters
DEFAULT_NCOL <- 2
DEFAULT_POINT_SIZE <- 0.3
CHUNK_SIZE <- 4
PT_SIZE <- 0.7
COEXPRESSION_COLOURS <- c("grey90", "#0072B2", "orange")
#COEXPRESSION_COLOURS = c("grey90", "#7570b3", "#1b9e77")

# ============================================================================
# FUNCTION DEFINITIONS
# ============================================================================

#' Check available metadata in Seurat objects
#'
#' @param seu_obj Seurat object
#' @param obj_name Name of the object for reporting
check_metadata <- function(seu_obj, obj_name) {
  cat("Metadata columns in", obj_name, ":\n")
  cat(paste(colnames(seu_obj@meta.data), collapse = ", "), "\n\n")
  
  # Check for common pseudotime-related columns
  pseudotime_cols <- grep("pseudotime|monocle|trajectory|time", 
                          colnames(seu_obj@meta.data), 
                          value = TRUE, ignore.case = TRUE)
  
  if (length(pseudotime_cols) > 0) {
    cat("Found potential pseudotime columns:", paste(pseudotime_cols, collapse = ", "), "\n")
  } else {
    cat("No pseudotime-related columns found\n")
  }
  cat("---\n")
}

#' Load and prepare Seurat objects
#'
#' @description Loads Seurat objects and sets the default identity
#' @return List containing the three Seurat objects
load_seurat_objects <- function() {
  cat("Loading Seurat objects...\n")
  
  # Load objects
  one_month_org <- readRDS(ONE_MONTH_FILE)
  three_month_org <- readRDS(THREE_MONTH_FILE)
  six_month_org <- readRDS(SIX_MONTH_FILE)
  
  # Set identities
  Idents(one_month_org) <- IDENTITY_COLUMN
  Idents(three_month_org) <- IDENTITY_COLUMN
  Idents(six_month_org) <- IDENTITY_COLUMN
  
  cat("Successfully loaded all objects\n")
  
  return(list(
    one_month = one_month_org,
    three_month = three_month_org,
    six_month = six_month_org
  ))
}

#' Find differentially expressed transcripts/isoforms
#'
#' @param seu.obj Seurat object
#' @param assay.type Assay to use for analysis (default: "iso")
#' @return DataFrame with differential expression results
returnTWAS_DE_Features <- function(seu.obj, assay.type = DEFAULT_ASSAY) {
  cat("Running FindAllMarkers for isoform-level analysis...\n")
  
  cellTypeMarkers_sc <- FindAllMarkers(
    seu.obj,
    assay = assay.type,
    slot = DEFAULT_SLOT,
    test.use = DE_TEST_METHOD
  )
  
  cat("Found", nrow(cellTypeMarkers_sc), "total markers\n")
  
  return(cellTypeMarkers_sc)
}

#' Create co-expression plot between two features
#'
#' @param seu.obj Seurat object
#' @param feature1 First gene/feature name
#' @param feature2 Second gene/feature name
#' @param reduction Dimensional reduction to use
#' @param title_prefix Prefix for plot title (e.g., "1 Month")
#' @return ggplot object
makeCoExpressionPlot <- function(seu.obj, 
                                 feature1 = "PTBP1", 
                                 feature2 = "PTBP2",
                                 reduction = DEFAULT_REDUCTION,
                                 title_prefix = "") {
  
  plt.output <- FeaturePlot(
    seu.obj,
    features = c(feature1, feature2), 
    blend = TRUE,
    reduction = reduction,
    pt.size = PT_SIZE,
    cols = COEXPRESSION_COLOURS,
    order = TRUE
  )
  
  # Add title if prefix provided
  if (title_prefix != "") {
    plot_title <- paste(title_prefix, ":", feature1, "vs", feature2, "Co-expression")
    plt.output <- plt.output + 
      plot_annotation(title = plot_title, theme = theme(plot.title = element_text(hjust = 0.5)))
  }
  
  return(plt.output)
}

#' Plot features in chunks using violin plots
#'
#' @param seurat_obj Seurat object
#' @param features Vector of feature names
#' @param chunk_size Number of features per chunk
#' @param ncol Number of columns in plot
#' @param title_prefix Prefix for plot titles
plot_feature_chunks <- function(seurat_obj, 
                                features, 
                                chunk_size = CHUNK_SIZE, 
                                ncol = DEFAULT_NCOL, 
                                title_prefix = "") {
  
  # Split features into chunks
  feature_chunks <- split(features, ceiling(seq_along(features) / chunk_size))
  
  # Loop over chunks and plot
  for (i in seq_along(feature_chunks)) {
    print(
      VlnPlot(seurat_obj, features = feature_chunks[[i]], ncol = ncol)
    )
  }
}

#' Create pseudotime comparison plots
#'
#' @param seu_obj Seurat object
#' @param title_prefix Prefix for plot titles
#' @param reduction Dimensional reduction to use
#' @return Combined ggplot object or NULL if pseudotime not available
plot_pseudotime_comparison <- function(seu_obj, 
                                       title_prefix = "", 
                                       reduction = DEFAULT_REDUCTION) {
  
  p1 <- DimPlot(
    seu_obj, 
    reduction = reduction, 
    group.by = IDENTITY_COLUMN,
    pt.size = DEFAULT_POINT_SIZE
  ) +
    ggtitle(paste(title_prefix, "- Cell types"))
  
  # Check if pseudotime exists in the metadata
  if ("pseudotime" %in% colnames(seu_obj@meta.data)) {
    p2 <- FeaturePlot(
      seu_obj, 
      features = "pseudotime", 
      reduction = reduction,
      pt.size = DEFAULT_POINT_SIZE
    ) +
      scale_color_viridis_c(name = "Pseudotime") +
      ggtitle(paste(title_prefix, "- Pseudotime"))
    
    return(p1 | p2)
  } else {
    cat("Warning: pseudotime not found in", title_prefix, "metadata. Skipping pseudotime plot.\n")
    return(p1)
  }
}

#' Create plots with and without cell type annotation labels shown
#'
#' @param seu_obj Seurat object
#' @param title_prefix Prefix for plot titles
#' @param reduction Dimensional reduction to use
#' @return Combined ggplot object or NULL if pseudotime not available
plot_annotation_comparison <- function(seu_obj, 
                                       title_prefix = "", 
                                       reduction = DEFAULT_REDUCTION) {
  
  p1 <- DimPlot(
    seu_obj, 
    reduction = reduction, 
    group.by = IDENTITY_COLUMN,
    pt.size = DEFAULT_POINT_SIZE
  ) +
    ggtitle(paste(title_prefix, "- Cell types"))
  
  # Check if pseudotime exists in the metadata
  p2 <- DimPlot(
    seu_obj, 
    reduction = reduction, 
    group.by = IDENTITY_COLUMN,
    pt.size = DEFAULT_POINT_SIZE,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste(title_prefix, "- Cell types")) +
    NoLegend()
  
  return(p1 | p2)
}

#' Generate ridge plots for multiple features
#'
#' @param seu_obj Seurat object
#' @param features Vector of feature names
#' @param ncol Number of columns in plot
#' @return ggplot object
plot_ridge_features <- function(seu_obj, 
                                features = ALL_FEATURES, 
                                ncol = DEFAULT_NCOL) {
  
  theme_set(theme_minimal())
  return(RidgePlot(seu_obj, features = features, ncol = ncol))
}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Load data ------------------------------------------------------------------
cat("=== LOADING DATA ===\n")
seurat_objects <- load_seurat_objects()
one_month_org <- seurat_objects$one_month
three_month_org <- seurat_objects$three_month
six_month_org <- seurat_objects$six_month

# Display cell type annotations  ---------------------------------------------
cat("\n=== PLOT CELL TYPE ANNOTATIONS ===\n")
plot_annotation_comparison(one_month_org, "1 month")
plot_annotation_comparison(three_month_org, "3 month") 
plot_annotation_comparison(six_month_org, "6 month")

# Check metadata for debugging -----------------------------------------------
cat("\n=== METADATA CHECK ===\n")
check_metadata(one_month_org, "1 Month")
check_metadata(three_month_org, "3 Month") 
check_metadata(six_month_org, "6 Month")

# Co-expression analysis -----------------------------------------------------
cat("\n=== CO-EXPRESSION ANALYSIS ===\n")

trace("FeaturePlot", quote(brewer.gran <- 2), at = 22L)

# PTBP1 vs PTBP2 across timepoints
cat("Analyzing PTBP1 vs PTBP2 co-expression...\n")
makeCoExpressionPlot(one_month_org, "PTBP1", "PTBP2", title_prefix = "1 Month")
makeCoExpressionPlot(three_month_org, "PTBP1", "PTBP2", title_prefix = "3 Month")
makeCoExpressionPlot(six_month_org, "PTBP1", "PTBP2", title_prefix = "6 Month")

# PTBP2 vs DLG4 across timepoints
cat("Analyzing PTBP2 vs DLG4 co-expression...\n")
makeCoExpressionPlot(one_month_org, "PTBP2", "DLG4", title_prefix = "1 Month")
makeCoExpressionPlot(three_month_org, "PTBP2", "DLG4", title_prefix = "3 Month")
makeCoExpressionPlot(six_month_org, "PTBP2", "DLG4", title_prefix = "6 Month")

# PTBP1 vs DLG4 across timepoints
cat("Analyzing PTBP1 vs DLG4 co-expression...\n")
makeCoExpressionPlot(one_month_org, "PTBP1", "DLG4", title_prefix = "1 Month")
makeCoExpressionPlot(three_month_org, "PTBP1", "DLG4", title_prefix = "3 Month")
makeCoExpressionPlot(six_month_org, "PTBP1", "DLG4", title_prefix = "6 Month")

# RBFOX family analysis
cat("Analyzing RBFOX family co-expression...\n")

# RBFOX1 vs RBFOX2
cat("RBFOX1 vs RBFOX2:\n")
makeCoExpressionPlot(one_month_org, "RBFOX1", "RBFOX2", title_prefix = "1 Month")  # 1 group of colocalization
makeCoExpressionPlot(three_month_org, "RBFOX1", "RBFOX2", title_prefix = "3 Month")  # 2 groups of colocalization
makeCoExpressionPlot(six_month_org, "RBFOX1", "RBFOX2", title_prefix = "6 Month")  # 2 groups of colocalization

# RBFOX2 vs RBFOX3
cat("RBFOX2 vs RBFOX3:\n")
# Note: Clear cell-type colocalization for RBFOX2 and RBFOX3 in neuronal subtypes
makeCoExpressionPlot(one_month_org, "RBFOX2", "RBFOX3", title_prefix = "1 Month")
makeCoExpressionPlot(three_month_org, "RBFOX2", "RBFOX3", title_prefix = "3 Month")  # Strong colocalization in 2 neuronal subtypes
makeCoExpressionPlot(six_month_org, "RBFOX2", "RBFOX3", title_prefix = "6 Month")  # Colocalization in 2 neuronal subtypes (RBFOX2 dominant)

# RBFOX1 vs RBFOX3
cat("RBFOX1 vs RBFOX3:\n")
makeCoExpressionPlot(one_month_org, "RBFOX1", "RBFOX3", title_prefix = "1 Month")  # No colocalization (1 RBFOX1 and 1 RBFOX3 specific clusters)
makeCoExpressionPlot(three_month_org, "RBFOX1", "RBFOX3", title_prefix = "3 Month")  # Colocalization in 2 cell groups
makeCoExpressionPlot(six_month_org, "RBFOX1", "RBFOX3", title_prefix = "6 Month")  # Colocalization in 1 cell type

# Expression distribution analysis -------------------------------------------
cat("\n=== EXPRESSION DISTRIBUTION ANALYSIS ===\n")

# Ridge plots - visualize single cell expression distributions in each cluster
# Hypothesis: bimodal distributions might result from differential isoform expression
# Note: NRXN2 expressed in (IN, CPN, oRG/Oligo progenitors cycling, Oligo progenitors, oRG cycling progenitors)
cat("Generating ridge plots for all timepoints...\n")
plot_ridge_features(one_month_org, ALL_FEATURES)
plot_ridge_features(three_month_org, ALL_FEATURES)
plot_ridge_features(six_month_org, ALL_FEATURES)

# Violin plots by chunks
cat("Generating violin plots by feature chunks...\n")
plot_feature_chunks(one_month_org, ALL_FEATURES, title_prefix = "One month")
plot_feature_chunks(three_month_org, ALL_FEATURES, title_prefix = "Three months")
plot_feature_chunks(six_month_org, ALL_FEATURES, title_prefix = "Six months")

# Pseudotime analysis -------------------------------------------------------
cat("\n=== PSEUDOTIME ANALYSIS ===\n")
plot_pseudotime_comparison(one_month_org, "1 month")
plot_pseudotime_comparison(three_month_org, "3 month") 
plot_pseudotime_comparison(six_month_org, "6 month")

# ============================================================================
# TODO: FUTURE ANALYSES
# ============================================================================

# TODO: Function that creates ridge plots for isoforms comprising a gene of interest
# (NRXN1, NRXN2, NRXN3)
# 
# Subfunction 1: get list of all transcripts expressed in seu.obj for a given gene list
# Create function that retains ridge plot for gene and then ridge plots of all its 
# isoforms that are expressed
#
# TODO: Source entropy functions for NRXN1, NRXN2, NRXN3, DLG4, PTBP1
#
# TODO: create function that allows me to rank and factor cell-types by their avg cell pseudotime values (smallest to largest) for a given organoid seurat object (with a pseudotime metadata column)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Next steps:\n")
cat("1. Implement isoform-specific ridge plots\n")
cat("2. Calculate expression entropy metrics\n")
cat("3. Perform differential isoform expression analysis\n")
cat("4. Generate publication-ready figures\n")