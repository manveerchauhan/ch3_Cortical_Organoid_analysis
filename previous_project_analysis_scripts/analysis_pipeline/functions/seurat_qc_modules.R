# =============================================================================
# Seurat QC Modules - Comprehensive Quality Control Workflow
# =============================================================================
# Author: Manveer Chauhan
# Description: Modular Seurat-based QC functions for dual-assay (RNA + iso)
#              single-cell analysis with dynamic thresholds and optimization
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(clustree)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

# =============================================================================
# FUNCTION 1: Initialize Seurat Objects with Dual Assays
# =============================================================================

#' Initialize Seurat object with both RNA (gene) and iso (isoform) assays
#'
#' @param sample_id Sample identifier
#' @param gene_counts Gene-level count matrix (ENSG IDs as rownames)
#' @param iso_counts Isoform-level count matrix (ENST IDs as rownames)
#' @param id_maps ID mapping object from create_id_mappings_from_flames()
#' @param min_cells Minimum cells required per feature (default: 5)
#' @param min_features Minimum features required per cell (default: 500)
#' @return Seurat object with RNA and iso assays
#' @export
initializeSeuratObjs <- function(sample_id, gene_counts, iso_counts, id_maps,
                                  min_cells = 5, min_features = 500) {

  log_message(sprintf("Initializing Seurat object for %s", sample_id), "DEBUG")

  # === Step 1: Convert gene counts (ENSG → Gene Symbol) ===
  log_message("Converting gene IDs to symbols", "DEBUG")
  gene_counts_converted <- apply_gene_symbol_conversion(
    count_matrix = gene_counts,
    gene_hash = id_maps$gene_hash,
    handle_duplicates = "sum"  # Sum counts for duplicate symbols
  )

  # === Step 2: Convert isoform counts (ENST → GeneName_TxID) ===
  log_message("Converting transcript IDs to GeneName_TxID format", "DEBUG")
  iso_counts_converted <- apply_transcript_label_conversion(
    count_matrix = iso_counts,
    tx_hash = id_maps$tx_hash
  )

  # === Step 3: Create Seurat object with gene assay ===
  log_message(sprintf("Creating Seurat object with min.cells=%d, min.features=%d",
                      min_cells, min_features), "DEBUG")

  seurat_obj <- CreateSeuratObject(
    counts = gene_counts_converted,
    project = sample_id,
    assay = "RNA",
    min.cells = min_cells,
    min.features = min_features
  )

  # Store initial cell count
  n_cells_initial <- ncol(seurat_obj)

  # === Step 4: Add isoform assay ===
  # Match cells between gene and isoform matrices
  shared_cells <- intersect(colnames(seurat_obj), colnames(iso_counts_converted))

  if (length(shared_cells) == 0) {
    stop("No shared cells between gene and isoform matrices!")
  }

  log_message(sprintf("Shared cells between assays: %d (gene: %d, iso: %d)",
                      length(shared_cells), ncol(seurat_obj), ncol(iso_counts_converted)), "DEBUG")

  # Subset isoform counts to matched cells and filter features
  iso_counts_matched <- iso_counts_converted[, shared_cells]

  # Filter isoforms: keep only those detected in min_cells
  iso_cell_counts <- Matrix::rowSums(iso_counts_matched > 0)
  iso_counts_filtered <- iso_counts_matched[iso_cell_counts >= min_cells, ]

  log_message(sprintf("Isoforms after filtering (min.cells=%d): %d",
                      min_cells, nrow(iso_counts_filtered)), "DEBUG")

  # Subset Seurat object to shared cells
  seurat_obj <- subset(seurat_obj, cells = shared_cells)

  # Add isoform assay
  seurat_obj[["iso"]] <- CreateAssayObject(counts = iso_counts_filtered)

  # === Step 5: Calculate QC metrics ===
  log_message("Calculating QC metrics", "DEBUG")

  # Mitochondrial genes (different patterns for gene vs isoform)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA")
  seurat_obj[["percent.mt.iso"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "iso")

  # Store sample metadata
  seurat_obj$sample_id <- sample_id
  seurat_obj$qc_status <- "unfiltered"

  log_message(sprintf("Seurat object created: %d cells, %d genes, %d isoforms",
                      ncol(seurat_obj),
                      nrow(seurat_obj[["RNA"]]),
                      nrow(seurat_obj[["iso"]])), "INFO")

  return(seurat_obj)
}

# =============================================================================
# FUNCTION 2: Preliminary QC Filtering with Dynamic Thresholds
# =============================================================================

#' Apply QC filtering with timepoint-specific dynamic thresholds
#'
#' @param seurat_obj Seurat object
#' @param timepoint Timepoint identifier (e.g., "1M_Org", "3M_Org")
#' @param qc_thresholds List of QC thresholds by timepoint
#' @return Filtered Seurat object with QC plots
#' @export
runSeuratPreliminaryFiltering <- function(seurat_obj, timepoint, qc_thresholds) {

  log_message(sprintf("Running preliminary QC filtering for %s", timepoint), "DEBUG")

  # Get timepoint-specific thresholds
  if (!timepoint %in% names(qc_thresholds)) {
    stop(sprintf("Timepoint %s not found in qc_thresholds", timepoint))
  }

  thresh <- qc_thresholds[[timepoint]]
  mt_cutoff <- thresh$mt_percent
  sd_mult <- thresh$sd_multiplier

  # === Step 1: Calculate dynamic thresholds ===
  log_message(sprintf("Calculating dynamic thresholds (mean ± %.1f SD, MT < %d%%)",
                      sd_mult, mt_cutoff), "DEBUG")

  thresholds <- calculate_dynamic_thresholds(seurat_obj, sd_multiplier = sd_mult, assay = "RNA")

  log_message(sprintf("Thresholds - nFeature: [%.0f, %.0f], nCount: [%.0f, %.0f]",
                      thresholds$nFeature_min, thresholds$nFeature_max,
                      thresholds$nCount_min, thresholds$nCount_max), "DEBUG")

  # === Step 2: Store pre-filter metrics ===
  n_cells_before <- ncol(seurat_obj)

  # === Step 3: Apply filters ===
  seurat_obj_filtered <- subset(
    seurat_obj,
    subset = nFeature_RNA >= thresholds$nFeature_min &
             nFeature_RNA <= thresholds$nFeature_max &
             nCount_RNA >= thresholds$nCount_min &
             nCount_RNA <= thresholds$nCount_max &
             percent.mt < mt_cutoff
  )

  n_cells_after <- ncol(seurat_obj_filtered)
  pct_retained <- round(100 * n_cells_after / n_cells_before, 1)

  log_message(sprintf("QC filtering complete: %d → %d cells (%.1f%% retained)",
                      n_cells_before, n_cells_after, pct_retained), "INFO")

  # === Step 4: Update QC status ===
  seurat_obj_filtered$qc_status <- "filtered"
  seurat_obj_filtered$qc_retained <- TRUE
  seurat_obj$qc_retained <- colnames(seurat_obj) %in% colnames(seurat_obj_filtered)

  # === Step 5: Generate QC comparison plots ===
  qc_plots <- plot_qc_comparison(
    seurat_obj_before = seurat_obj,
    seurat_obj_after = seurat_obj_filtered,
    sample_id = unique(seurat_obj$sample_id)[1],
    thresholds = thresholds,
    mt_cutoff = mt_cutoff
  )

  # Store QC plots in object metadata
  seurat_obj_filtered@misc$qc_plots <- qc_plots
  seurat_obj_filtered@misc$qc_thresholds <- thresholds
  seurat_obj_filtered@misc$qc_stats <- list(
    n_cells_before = n_cells_before,
    n_cells_after = n_cells_after,
    pct_retained = pct_retained
  )

  return(seurat_obj_filtered)
}

# =============================================================================
# FUNCTION 3: Normalization, Scaling, and Cell Cycle Scoring
# =============================================================================

#' Normalize, scale, and score cell cycle phases
#'
#' @param seurat_obj Filtered Seurat object
#' @param cell_cycle_genes Path to cell cycle gene list (cycle.rda)
#' @param n_pcs_test Number of PCs to test for elbow detection (default: 50)
#' @return Seurat object with normalization and optimal PC selection
#' @export
NormaliseScaleAndElbow <- function(seurat_obj, cell_cycle_genes, n_pcs_test = 50) {

  log_message("Running SCTransform normalization", "DEBUG")

  # === Step 1: SCTransform normalization ===
  seurat_obj <- SCTransform(
    seurat_obj,
    assay = "RNA",
    new.assay.name = "SCT",
    verbose = FALSE
  )

  # === Step 2: Load cell cycle genes ===
  log_message("Loading cell cycle genes", "DEBUG")

  if (!file.exists(cell_cycle_genes)) {
    warning(sprintf("Cell cycle gene file not found: %s. Skipping cell cycle scoring.", cell_cycle_genes))
    cc_genes <- NULL
  } else {
    load(cell_cycle_genes)  # Loads 'cycle' object
    cc_genes <- cycle
  }

  # === Step 3: Cell cycle scoring ===
  if (!is.null(cc_genes)) {
    log_message("Scoring cell cycle phases (will NOT regress)", "DEBUG")

    seurat_obj <- CellCycleScoring(
      seurat_obj,
      s.features = cc_genes$s.genes,
      g2m.features = cc_genes$g2m.genes,
      set.ident = FALSE
    )

    # Log cell cycle distribution
    cc_table <- table(seurat_obj$Phase)
    log_message(sprintf("Cell cycle distribution: G1=%d, S=%d, G2M=%d",
                        cc_table["G1"], cc_table["S"], cc_table["G2M"]), "DEBUG")
  } else {
    seurat_obj$Phase <- "Unknown"
    seurat_obj$S.Score <- NA
    seurat_obj$G2M.Score <- NA
  }

  # === Step 4: PCA and elbow plot ===
  log_message(sprintf("Running PCA (testing %d PCs)", n_pcs_test), "DEBUG")

  seurat_obj <- RunPCA(
    seurat_obj,
    assay = "SCT",
    npcs = n_pcs_test,
    verbose = FALSE
  )

  # === Step 5: Determine optimal PCs ===
  optimal_pcs <- determine_optimal_pcs(seurat_obj, method = "elbow")

  log_message(sprintf("Optimal PCs selected: %d", optimal_pcs), "INFO")

  # Store optimal PCs in metadata
  seurat_obj@misc$optimal_pcs <- optimal_pcs

  # Create elbow plot
  elbow_plot <- ElbowPlot(seurat_obj, ndims = n_pcs_test) +
    geom_vline(xintercept = optimal_pcs, linetype = "dashed", color = "red") +
    labs(title = sprintf("PCA Elbow Plot - Optimal PCs: %d", optimal_pcs)) +
    theme_publication()

  seurat_obj@misc$elbow_plot <- elbow_plot

  return(seurat_obj)
}

# =============================================================================
# FUNCTION 4: Optimal Clustering Resolution via Silhouette
# =============================================================================

#' Find optimal clustering resolution using silhouette analysis
#'
#' @param seurat_obj Seurat object with PCA
#' @param resolutions Vector of resolutions to test (default: c(0.1, 0.3, 0.5, 0.8, 1.0))
#' @return Optimal resolution value
#' @export
find.optimal.cluster.res <- function(seurat_obj, resolutions = c(0.1, 0.3, 0.5, 0.8, 1.0)) {

  log_message(sprintf("Testing %d clustering resolutions", length(resolutions)), "DEBUG")

  optimal_pcs <- seurat_obj@misc$optimal_pcs
  if (is.null(optimal_pcs)) {
    optimal_pcs <- 30
    warning("Optimal PCs not found in metadata, using default: 30")
  }

  # === Step 1: Build SNN graph ===
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "pca",
    dims = 1:optimal_pcs,
    verbose = FALSE
  )

  # === Step 2: Test all resolutions ===
  silhouette_scores <- data.frame(
    resolution = numeric(),
    silhouette_score = numeric(),
    n_clusters = numeric()
  )

  for (res in resolutions) {
    seurat_obj <- FindClusters(
      seurat_obj,
      resolution = res,
      verbose = FALSE
    )

    # Calculate silhouette score
    sil_score <- calculate_silhouette_score(seurat_obj, optimal_pcs = optimal_pcs)
    n_clusters <- length(unique(seurat_obj$seurat_clusters))

    silhouette_scores <- rbind(silhouette_scores, data.frame(
      resolution = res,
      silhouette_score = sil_score,
      n_clusters = n_clusters
    ))

    log_message(sprintf("Resolution %.1f: %d clusters, silhouette = %.3f",
                        res, n_clusters, sil_score), "DEBUG")
  }

  # === Step 3: Select optimal resolution ===
  optimal_idx <- which.max(silhouette_scores$silhouette_score)
  optimal_res <- silhouette_scores$resolution[optimal_idx]
  optimal_sil <- silhouette_scores$silhouette_score[optimal_idx]
  optimal_n_clusters <- silhouette_scores$n_clusters[optimal_idx]

  log_message(sprintf("Optimal resolution: %.1f (%d clusters, silhouette = %.3f)",
                      optimal_res, optimal_n_clusters, optimal_sil), "INFO")

  # === Step 4: Re-run clustering with optimal resolution ===
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = optimal_res,
    verbose = FALSE
  )

  # === Step 5: Generate clustree visualization ===
  log_message("Generating clustree visualization", "DEBUG")

  clustree_plot <- clustree(seurat_obj, prefix = "SCT_snn_res.") +
    theme_publication() +
    labs(title = sprintf("Clustering Resolution Tree (Optimal: %.1f)", optimal_res))

  # Store results in metadata
  seurat_obj@misc$optimal_resolution <- optimal_res
  seurat_obj@misc$silhouette_scores <- silhouette_scores
  seurat_obj@misc$clustree_plot <- clustree_plot

  return(optimal_res)
}

# =============================================================================
# FUNCTION 5: UMAP and Doublet Removal
# =============================================================================

#' Run UMAP and remove doublets with DoubletFinder
#'
#' @param seurat_obj Seurat object with clustering
#' @param optimal_pcs Number of PCs to use for UMAP
#' @param optimal_res Optimal clustering resolution
#' @param doublet_rate Expected doublet rate (default: 0.016 = 1.6%)
#' @return Seurat object with UMAP and doublets removed
#' @export
runUMAP.removeDoublets <- function(seurat_obj, optimal_pcs, optimal_res, doublet_rate = 0.016) {

  log_message("Running UMAP", "DEBUG")

  # === Step 1: Run UMAP ===
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca",
    dims = 1:optimal_pcs,
    verbose = FALSE
  )

  # === Step 2: DoubletFinder parameter estimation ===
  log_message("Running DoubletFinder", "DEBUG")

  # Estimate pK parameter
  sweep_res <- paramSweep_v3(seurat_obj, PCs = 1:optimal_pcs, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)

  # Select optimal pK
  optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  log_message(sprintf("DoubletFinder optimal pK: %.2f", optimal_pk), "DEBUG")

  # === Step 3: Run DoubletFinder ===
  # Calculate expected doublets
  n_cells <- ncol(seurat_obj)
  n_doublets_expected <- round(n_cells * doublet_rate)

  log_message(sprintf("Expected doublets: %d (%.1f%% of %d cells)",
                      n_doublets_expected, doublet_rate * 100, n_cells), "DEBUG")

  # Homotypic doublet proportion
  homotypic_prop <- modelHomotypic(seurat_obj$seurat_clusters)
  n_doublets_adj <- round(n_doublets_expected * (1 - homotypic_prop))

  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(
    seurat_obj,
    PCs = 1:optimal_pcs,
    pN = 0.25,
    pK = optimal_pk,
    nExp = n_doublets_adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )

  # === Step 4: Extract doublet classifications ===
  # DoubletFinder adds columns with dynamic names - find the classification column
  df_cols <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  doublet_col <- df_cols[length(df_cols)]  # Use most recent run

  # Rename for convenience
  seurat_obj$doublet_classification <- seurat_obj@meta.data[[doublet_col]]

  n_doublets <- sum(seurat_obj$doublet_classification == "Doublet")
  pct_doublets <- round(100 * n_doublets / n_cells, 1)

  log_message(sprintf("Doublets detected: %d (%.1f%%)", n_doublets, pct_doublets), "INFO")

  # === Step 5: Remove doublets ===
  seurat_obj_singlets <- subset(seurat_obj, subset = doublet_classification == "Singlet")

  n_cells_final <- ncol(seurat_obj_singlets)
  log_message(sprintf("Cells after doublet removal: %d", n_cells_final), "INFO")

  # === Step 6: Re-run UMAP on clean data ===
  log_message("Re-running UMAP on singlets", "DEBUG")

  seurat_obj_singlets <- RunUMAP(
    seurat_obj_singlets,
    reduction = "pca",
    dims = 1:optimal_pcs,
    verbose = FALSE
  )

  # Store doublet removal stats
  seurat_obj_singlets@misc$doublet_stats <- list(
    n_doublets = n_doublets,
    pct_doublets = pct_doublets,
    n_cells_before = n_cells,
    n_cells_after = n_cells_final
  )

  return(seurat_obj_singlets)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Calculate dynamic QC thresholds
#'
#' @param seurat_obj Seurat object
#' @param sd_multiplier SD multiplier for thresholds (default: 1.5)
#' @param assay Assay to use (default: "RNA")
#' @return List of thresholds
calculate_dynamic_thresholds <- function(seurat_obj, sd_multiplier = 1.5, assay = "RNA") {

  nFeature_mean <- mean(seurat_obj@meta.data[[paste0("nFeature_", assay)]])
  nFeature_sd <- sd(seurat_obj@meta.data[[paste0("nFeature_", assay)]])

  nCount_mean <- mean(seurat_obj@meta.data[[paste0("nCount_", assay)]])
  nCount_sd <- sd(seurat_obj@meta.data[[paste0("nCount_", assay)]])

  return(list(
    nFeature_min = max(200, nFeature_mean - sd_multiplier * nFeature_sd),
    nFeature_max = nFeature_mean + sd_multiplier * nFeature_sd,
    nCount_min = max(500, nCount_mean - sd_multiplier * nCount_sd),
    nCount_max = nCount_mean + sd_multiplier * nCount_sd
  ))
}

#' Determine optimal number of PCs using elbow method
#'
#' @param seurat_obj Seurat object with PCA
#' @param method Method for PC selection ("elbow" or "variance")
#' @return Optimal number of PCs
determine_optimal_pcs <- function(seurat_obj, method = "elbow") {

  pca_stdev <- seurat_obj@reductions$pca@stdev

  if (method == "elbow") {
    # Find elbow point using percentage change
    pct_change <- c(0, diff(pca_stdev)) / pca_stdev[-length(pca_stdev)]
    optimal_pc <- which(abs(pct_change) < 0.1)[1]

    if (is.na(optimal_pc)) optimal_pc <- 30  # Default fallback

  } else if (method == "variance") {
    # Use 90% variance explained
    cumsum_var <- cumsum(pca_stdev^2) / sum(pca_stdev^2)
    optimal_pc <- which(cumsum_var >= 0.90)[1]
  }

  return(optimal_pc)
}

#' Calculate silhouette score for clustering
#'
#' @param seurat_obj Seurat object with clusters
#' @param optimal_pcs Number of PCs to use
#' @return Average silhouette score
calculate_silhouette_score <- function(seurat_obj, optimal_pcs) {

  # Extract PCA coordinates
  pca_coords <- Embeddings(seurat_obj, reduction = "pca")[, 1:optimal_pcs]

  # Calculate distance matrix
  dist_mat <- dist(pca_coords)

  # Calculate silhouette
  sil <- cluster::silhouette(
    x = as.numeric(seurat_obj$seurat_clusters),
    dist = dist_mat
  )

  # Return average silhouette width
  return(mean(sil[, 3]))
}

#' Plot QC comparison before and after filtering
#'
#' @param seurat_obj_before Seurat object before filtering
#' @param seurat_obj_after Seurat object after filtering
#' @param sample_id Sample identifier
#' @param thresholds List of QC thresholds
#' @param mt_cutoff MT% cutoff
#' @return List of ggplot objects
plot_qc_comparison <- function(seurat_obj_before, seurat_obj_after, sample_id, thresholds, mt_cutoff) {

  # Prepare data
  df_before <- seurat_obj_before@meta.data %>%
    mutate(status = "Before Filtering")

  df_after <- seurat_obj_after@meta.data %>%
    mutate(status = "After Filtering")

  df_combined <- rbind(df_before, df_after)

  # Plot 1: nFeature_RNA
  p1 <- ggplot(df_combined, aes(x = status, y = nFeature_RNA, fill = status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white") +
    geom_hline(yintercept = thresholds$nFeature_min, linetype = "dashed", color = "red") +
    geom_hline(yintercept = thresholds$nFeature_max, linetype = "dashed", color = "red") +
    labs(title = paste(sample_id, "- Features Detected"), x = "", y = "nFeature_RNA") +
    theme_publication() +
    theme(legend.position = "none")

  # Plot 2: nCount_RNA
  p2 <- ggplot(df_combined, aes(x = status, y = nCount_RNA, fill = status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white") +
    geom_hline(yintercept = thresholds$nCount_min, linetype = "dashed", color = "red") +
    geom_hline(yintercept = thresholds$nCount_max, linetype = "dashed", color = "red") +
    labs(title = paste(sample_id, "- Total Counts"), x = "", y = "nCount_RNA") +
    theme_publication() +
    theme(legend.position = "none")

  # Plot 3: percent.mt
  p3 <- ggplot(df_combined, aes(x = status, y = percent.mt, fill = status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white") +
    geom_hline(yintercept = mt_cutoff, linetype = "dashed", color = "red") +
    labs(title = paste(sample_id, "- Mitochondrial %"), x = "", y = "Percent MT") +
    theme_publication() +
    theme(legend.position = "none")

  return(list(nFeature = p1, nCount = p2, percent_mt = p3))
}

log_message("Seurat QC modules loaded successfully", "INFO")
