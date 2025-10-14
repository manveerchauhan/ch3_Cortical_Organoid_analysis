# Generalized Script for Multiple Sample scRNA-seq QC
library(Seurat)
library(tidyverse)
library(reticulate)
library(DoubletFinder)
library(gridExtra)
library(grid)
library(clustree)
library(RColorBrewer)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis")
set.seed(5728)
theme_set(ggmin::theme_min())

## Load supporting resources------
load("/data/gpfs/projects/punim0646/manveer/cycle.rda")
source("/data/gpfs/projects/punim2251/scRNAseq_scripts/silhouette.R")
id.gene.dict <- read.csv("/data/gpfs/projects/punim2251/resoucres/v41_ENSG_ID_GENEsymbol.csv", row.names = 1)
id.gene.dict$geneID <- row.names(id.gene.dict)

### Define sample configuration------
# Create a configuration data frame for all your samples
sample_config <- data.frame(
  sample_id = c("org1A", "org1B", "org3A", "org3B", "org3C", "org6A", "org6B", "org6C"),
  timepoint = c("1month", "1month", "3month", "3month", "3month", "6month", "6month", "6month"),
  replicate = c("A", "B", "A", "B", "C", "A", "B", "C"),
  gene_matrix_file = c(
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org1A_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org1B_gene_counts.csv", 
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org3A_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org3B_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org3C_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org6A_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org6B_gene_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/gene_counts/org6C_gene_counts.csv"
  ),
  iso_matrix_file = c(
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org1A_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org1B_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org3A_transcript_counts.csv", 
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org3B_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org3C_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org6A_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org6B_transcript_counts.csv",
    "/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/Org_Timeseries_Analysis/previous_raw_data/isoform_counts/org6C_transcript_counts.csv"
  ),
  # Sample-specific parameters (can be customized per sample)
  max_counts = c(33000, 33000, 35000, 35000, 35000, 40000, 40000, 40000),
  mt_threshold = c(15, 15, 15, 15, 15, 20, 20, 20),
  optimal_pcs = c(15, 15, 20, 20, 20, 20, 20, 20),
  cluster_resolution = c(0.6, 0.6, 0.7, 0.7, 0.7, 0.8, 0.8, 0.8),
  stringsAsFactors = FALSE
)

### Keep all existing functions (unchanged)------
countMatrixConvertIDstoSymbols <- function(geneCountMatrix,
                                           idGeneDict = id.gene.dict) {
  geneCountMatrix$geneID <- row.names(geneCountMatrix)
  
  result <- left_join(geneCountMatrix, 
                      idGeneDict, 
                      by = "geneID")
  
  result <- result %>%
    select(genesymbol, everything())
  
  columns_to_sum <- setdiff(names(result), "genesymbol")
  result[columns_to_sum] <- lapply(result[columns_to_sum], as.numeric)

  result$totalCounts <- rowSums(result[columns_to_sum], na.rm = TRUE)
  
  result <- result %>%
    select(totalCounts, everything())
  
  non_unique_rows <- result[duplicated(result$genesymbol) | duplicated(result$genesymbol, fromLast = TRUE), ]
  duplicates <- non_unique_rows
  
  result <- result %>% 
    arrange(desc(totalCounts))
  
  unique_counts <- result %>% 
    distinct(genesymbol, .keep_all = TRUE)
  
  unique_counts$totalCounts <- NULL
  rownames(unique_counts) <- unique_counts$genesymbol
  unique_counts$genesymbol <- NULL
  
  message("Successfully converted ENSG IDs to Gene Symbols using the provided dictionary")
  return(unique_counts)
}

readInCountMatrices <- function(isoformMatrixFile,
                              geneMatrixFile,
                              IDstoSymbols = T){
  counts.isoformLvl <- read.csv(isoformMatrixFile, row.names = 1)

  counts.geneLvl.IDs <- read.csv(geneMatrixFile, row.names = 1)
  if(IDstoSymbols){
    counts.geneLvl <- countMatrixConvertIDstoSymbols(geneCountMatrix = counts.geneLvl.IDs)
    counts.geneLvl$geneID <- NULL
  } else{
    counts.geneLvl <- counts.geneLvl.IDs
  }
  
  col_names_gene <- names(counts.geneLvl)
  counts.isoformLvl <- counts.isoformLvl[col_names_gene]
  
  message("Successfully imported and formatted gene and isoform level count matrices.")
  
  function_output <- list(counts.geneLvl, counts.isoformLvl)
  return(function_output)
}

initializeSeuratObjs <- function(geneMatrixFilePath,
                                 isoMatrixFilePath,
                                 SAMPLE_ID = "Unlabelled_Sample",
                                 convertIDstoSymbols = T){
  
  formatted_matrices <- readInCountMatrices(isoformMatrixFile = isoMatrixFilePath,
                                            geneMatrixFile = geneMatrixFilePath,
                                            IDstoSymbols = convertIDstoSymbols)
  
  seurat.obj <- CreateSeuratObject(counts = formatted_matrices[1], project = SAMPLE_ID,
                                   min.cells = 3, min.features = 1)
  seurat.obj.isoLvl <- CreateSeuratObject(counts = formatted_matrices[2], project = SAMPLE_ID,
                                          min.cells = 3, min.features = 1)
  
  function.output <- list(seurat.obj, seurat.obj.isoLvl)
  message("Unfiltered Seurat objs initialised: returning gene-lvl seurat obj in index 1, and iso-lvl obj in index 2")
  
  return(function.output)
}

runSeuratPreliminaryFiltering <- function(geneLvl.seuratObj,
                                          isoLvl.seuratObj,
                                          SAMPLE_ID = "Unlabelled_Sample",
                                          MAX.GENES = 2500,
                                          MIN.GENES = 999999999,
                                          MIN.COUNTS = 500,
                                          MAX.COUNTS = 100000,
                                          MT_THRESHOLD = 15,
                                          CUSTOM.MAX.GENES = NULL,
                                          CUSTOM.MIN.GENES = NULL){
  table <- data.frame()
  cluster_resolution.figs <- list()
  
  seurat.obj <- geneLvl.seuratObj
  seurat.obj.isoLvl <- isoLvl.seuratObj
  
  table <- rbind(table, data.frame("Cells"=dim(seurat.obj)[2],
                                   "Median genes per cell"=median(seurat.obj$nFeature_RNA), 
                                   row.names = paste0('min genes > 0'),check.names = FALSE))
  initial.cell.num <- dim(seurat.obj)[2]
  
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  median.mito.content.before = median(seurat.obj@meta.data[["percent.mt"]])
  
  VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  association.plt.before <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + NoLegend() + labs(title = "Association between reads and \nunique genes per cell BEFORE filtering")
  association.plt.before
  
  vln1 <- VlnPlot(seurat.obj, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nBefore Filtering") + NoLegend()
  vln2 <- VlnPlot(seurat.obj, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nBefore Filtering") + NoLegend()
  vln3 <- VlnPlot(seurat.obj, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nBefore Filtering") + NoLegend()
  
  vlnplots.before.QC <- list(vln1, vln2, vln3)
  
  MAX.GENES <- round(mean(seurat.obj$nFeature_RNA) + (1.5 * sd(seurat.obj$nFeature_RNA)))
  MIN.GENES <- round(mean(seurat.obj$nFeature_RNA) - (1.5 * sd(seurat.obj$nFeature_RNA)))
  if (!is.null(CUSTOM.MIN.GENES)) {
    MIN.GENES <- CUSTOM.MIN.GENES
    message(paste("Custom lower bound gene cutoff being used:", MIN.GENES))
  }
  if (!is.null(CUSTOM.MAX.GENES)) {
    MAX.GENES <- CUSTOM.MAX.GENES
    message(paste("Custom upper bound gene cutoff being used:", MAX.GENES))
  }
  
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > MIN.GENES & nFeature_RNA < MAX.GENES 
                        & percent.mt < MT_THRESHOLD & nCount_RNA < MAX.COUNTS & nCount_RNA > MIN.COUNTS)
  
  table <- rbind(table, data.frame("Cells" = dim(seurat.obj)[2],
                                   "Median genes per cell" = median(seurat.obj$nFeature_RNA), 
                                   row.names = paste0(MAX.GENES, ' > genes > ', MIN.GENES), check.names = FALSE))
  
  association.plt.after <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm") + labs(title = "Association between reads \nand unique genes per cell AFTER filering") +
    NoLegend()
  
  VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  vln4 <- VlnPlot(seurat.obj, features = c("nFeature_RNA")) + labs(title = "Unique Genes per Cell\nAfter Filtering") + NoLegend()
  vln5 <- VlnPlot(seurat.obj, features = c("nCount_RNA")) + labs(title = "Reads per Cell\nAfter Filtering") + NoLegend()
  vln6 <- VlnPlot(seurat.obj, features = c("percent.mt")) + labs(title = "Mitochondrial Content (%)\nAfter Filtering") + NoLegend()
  vlnplots.after.QC <- list(vln4, vln5, vln6)
  grid.draw(association.plt.before | association.plt.after)
  
  seurat.obj.isoLvl <- subset(seurat.obj.isoLvl, cells = row.names(seurat.obj@meta.data))
  message("Poor quality cells also removed from isoform level seurat object")

  filename <- paste0(SAMPLE_ID, "-filter-figs-QC.pdf")
  summaryTitle <- paste("QC SCRIPT", SAMPLE_ID, "FILTERING INFORMATION")
  
  filter.figs.layout <- grid.arrange(vlnplots.before.QC[[1]], vlnplots.before.QC[[2]], vlnplots.before.QC[[3]],
                                     vlnplots.after.QC[[1]], vlnplots.after.QC[[2]], vlnplots.after.QC[[3]],
                                     association.plt.before, association.plt.after, tableGrob(table),
                                     nrow = 3, ncol = 3,
                                     top=textGrob(summaryTitle,
                                                  gp = gpar(fontsize = 12)))
  ggsave(filename, filter.figs.layout, width = 14, height = 12)
  
  function.output <- list(seurat.obj, seurat.obj.isoLvl)
  message(paste("Preliminary filtering performed, summary generated in :", filename))
  message("returning gene-lvl seurat obj in index [[1]] and iso-lvl obj in index [[2]]")
  return(function.output)
}

NormaliseScaleAndElbow <- function(SAMPLE_ID = "Unlabelled_Sample",
                                   seurat.obj,
                                   seurat.obj.iso,
                                   labelCellCycle = T,
                                   regressCellPhase = F){
  
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj.iso <- NormalizeData(seurat.obj.iso)
  
  seurat.obj <- FindVariableFeatures(seurat.obj, 
                                     selection.method = 'vst',
                                     nfeatures = 2000)
  seurat.obj.iso <- FindVariableFeatures(seurat.obj.iso, 
                                         selection.method = 'vst',
                                         nfeatures = 2000)
  
  top10 <- head(VariableFeatures(seurat.obj), 10)
  
  variable.feature.plot <- VariableFeaturePlot(seurat.obj) +
    labs(title = "Top 2000 Variable Genes across Cells") +
    NoLegend() +
    ggmin::theme_min()
  variable.feature.plot <- LabelPoints(plot = variable.feature.plot, points = top10, labels = top10,
                                       repel = TRUE, xnudge = 0, ynudge = 0)
  variable.feature.plot
  
  if(labelCellCycle){
    seurat.obj <- CellCycleScoring(seurat.obj,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)
  }

  all.features <- rownames(seurat.obj)
  all.features.iso <- rownames(seurat.obj.iso)
  
  if(regressCellPhase){
    seurat.obj$CC.Difference <- seurat.obj$S.Score - seurat.obj$G2M.Score
    seurat.obj <- ScaleData(seurat.obj, vars.to.regress = "CC.Difference", features = all.features)
    message("Regressed out cell cycle effect")
  } else{
    seurat.obj <- ScaleData(seurat.obj, features = all.features)
  }
  
  seurat.obj.iso <- ScaleData(seurat.obj.iso, features = all.features.iso)
  
  seurat.obj <- RunPCA(seurat.obj,
                       features = VariableFeatures(object = seurat.obj))
  seurat.obj.iso <- RunPCA(seurat.obj.iso,
                           features = VariableFeatures(object = seurat.obj.iso))
  
  elbow.plt <- ElbowPlot(seurat.obj) + 
    labs(title = 'Standard Deviation explained by each PC\nGene Lvl') +
    ggmin::theme_min()
  elbow.plt.iso <- ElbowPlot(seurat.obj.iso) + 
    labs(title = 'Standard Deviation explained by each PC\nIso Lvl') +
    ggmin::theme_min()
  
  if(labelCellCycle){
    cellCycle.pc1.pc2 <- DimPlot(seurat.obj, reduction = "pca",
                                 label = FALSE,  label.size = 3,
                                 repel = TRUE, group.by = 'Phase') + 
      labs(title = "Cell Cycle Effect on PC1 and PC2", color = "Cell Phase") +
      ggmin::theme_min()
    
    cellCycle.pc2.pc3 <- DimPlot(seurat.obj, reduction = "pca",
                                 dims = c(2,3),
                                 label = FALSE,  label.size = 3,
                                 repel = TRUE, group.by = 'Phase') + 
      labs(title = "Cell Cycle Effect on PC2 and PC3", color = "Cell Phase") +
      ggmin::theme_min()
  }
  
  if(regressCellPhase){
    filename <- paste0(SAMPLE_ID, "removed-cellCycle-PCs.pdf")
    summaryTitle <- paste("Cell Cycle (Regression) Effects + Elbow Plot", SAMPLE_ID, "Information\n")
  } else if(labelCellCycle){
    filename <- paste0(SAMPLE_ID, "-cellCycle-PCs.pdf")
    summaryTitle <- paste("Cell Cycle Effects + Elbow Plot", SAMPLE_ID, "Information\n")
  }else{
    filename <- paste0(SAMPLE_ID, "-elbowPlots.pdf")
    summaryTitle <- paste("Elbow Plots", SAMPLE_ID, "Information\n")
  }
  
  if(labelCellCycle){
    cellcycle.pc.figs.layout <- grid.arrange(cellCycle.pc1.pc2, cellCycle.pc2.pc3,
                                             elbow.plt, elbow.plt.iso,
                                             variable.feature.plot,
                                             nrow = 3, ncol = 2,
                                             top=textGrob(summaryTitle,
                                                          gp = gpar(fontsize = 12)))
    ggsave(filename, cellcycle.pc.figs.layout, width = 12, height = 12)
  } else{
    cellcycle.pc.figs.layout <- grid.arrange(elbow.plt, elbow.plt.iso,
                                             variable.feature.plot,
                                             nrow = 2, ncol = 2,
                                             top=textGrob(summaryTitle,
                                                          gp = gpar(fontsize = 12)))
    ggsave(filename, cellcycle.pc.figs.layout, width = 12, height = 12)
  }

  function.output <- list(seurat.obj, seurat.obj.iso)
  message(paste("Elbow and cell cycle plots generated in:", filename)) 
  message("Returning gene-lvl seurat obj in index 1, and iso-lvl obj in index 2")
  return(function.output)
}

find.optimal.cluster.res <- function(seurat.obj,
                                     SAMPLE_ID = "Unlabelled_Sample",
                                     npc = 20,
                                     desired.cluster.res = 1.2,
                                     commit.to.custom.res = F){
  cluster.resolution.figs <- list()
  
  initial.seurat.obj <- seurat.obj
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npc)
  
  avg.clusterSilVal.list <- as.data.frame(optimize_silhouette(sobject = seurat.obj,
                                                              test_res = seq(0.1, 1.2, by = 0.1),
                                                              summary_plot = TRUE,
                                                              reduction = "pca"))
  avg.clusterSilVal.list.filtered <- avg.clusterSilVal.list %>%
    dplyr::rename(avg_sil_vals = sil_vals) %>%
    group_by(num_clusters) %>%
    dplyr::slice(which.max(avg_sil_vals)) %>%
    ungroup() %>%
    arrange(desc(avg_sil_vals)) %>% 
    mutate(avg_sil_vals = round(avg_sil_vals, 4)) %>% 
    slice_head(n = 5)
  
  resolutions.to.plot <- c(avg.clusterSilVal.list.filtered$res_vals, desired.cluster.res)
  resolutions.to.test <- seq(0.1, 1.2, by = 0.1)
  
  seurat.obj <- FindClusters(seurat.obj, resolution = resolutions.to.test)
  
  for (i in 1:length(resolutions.to.plot)) {
    cluster.resolution.figs[[i]] <- DimPlot(seurat.obj,
                                            group.by = paste0("RNA_snn_res.", resolutions.to.plot[[i]]),
                                            label = T)
  }
  
  clusterTreePath <- paste0("./", SAMPLE_ID, "Clustree.pdf")
  cluster.tree <- clustree(seurat.obj) 
  ggsave(cluster.tree + scale_color_manual(values = brewer.pal(12, "Paired")) + 
           guides(color = guide_legend(ncol = 2)), 
         width = 10, height = 8, filename = clusterTreePath)
  
  if(commit.to.custom.res){
    message(paste("Committing to cluster res of:", desired.cluster.res))
    commited.res.seu.obj <- FindNeighbors(initial.seurat.obj, dims = 1:npc)
    commited.res.seu.obj <- FindClusters(commited.res.seu.obj, resolution = desired.cluster.res)
    
    return(c(commited.res.seu.obj, avg.clusterSilVal.list.filtered))
  } else{
    return(c(seurat.obj, avg.clusterSilVal.list.filtered))
  }
}

runUMAP.removeDoublets <- function(seurat.obj,
                                   npc = 20,
                                   SAMPLE_ID = "Unlabelled_sample",
                                   cluster.res = 0.6){
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc)
  
  sweep.res.list <- paramSweep(seurat.obj, PCs = 1:npc, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  BCmvn <- find.pK(sweep.stats)
  
  pK <- BCmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- seurat.obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  nExp_poi <- round(0.016*nrow(seurat.obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat.obj <- doubletFinder(seurat.obj, 
                              PCs = 1:20, 
                              pN = 0.25, 
                              pK = pK, 
                              nExp = nExp_poi.adj,
                              reuse.pANN = FALSE, sct = FALSE)
  
  colnames(seurat.obj@meta.data) <- sub("DF.classifications_.*$", "DF.classifications", colnames(seurat.obj@meta.data))
  
  statsDoublets <- seurat.obj@meta.data %>% 
    group_by(DF.classifications) %>%
    summarize(median_nCount_RNA = median(nCount_RNA),
              median_nFeature_RNA = median(nFeature_RNA),
              count = n())
  
  doublets.umap <- DimPlot(seurat.obj, reduction = 'umap', group.by = "DF.classifications") + 
    labs(title = "Doublets that were detected and removed") +
    ggmin::theme_min()
  doublets.umap
  
  seurat.obj <- subset(seurat.obj, subset = DF.classifications == 'Singlet')
  
  umap_title <- paste0(SAMPLE_ID, "\nResolution: ", cluster.res, ", Dimensions: ", npc)
  
  baseUMAPplot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, ) + 
    labs(color = "cluster \n(from PCA)", 
         title = umap_title) +
    ggmin::theme_min()
  baseUMAPplot
  
  return(c(seurat.obj, doublets.umap, baseUMAPplot))
}

### New wrapper function to process a single sample------
process_single_sample <- function(sample_row) {
  sample_id <- sample_row$sample_id
  gene_file <- sample_row$gene_matrix_file
  iso_file <- sample_row$iso_matrix_file
  max_counts <- sample_row$max_counts
  mt_thresh <- sample_row$mt_threshold
  n_pcs <- sample_row$optimal_pcs
  cluster_res <- sample_row$cluster_resolution
  
  message(paste("Processing sample:", sample_id))
  
  # Step 1: Initialize Seurat objects
  unfiltered.seu.objs <- initializeSeuratObjs(
    SAMPLE_ID = sample_id,
    geneMatrixFilePath = gene_file,
    isoMatrixFilePath = iso_file,
    convertIDstoSymbols = FALSE
  )
  
  # Step 2: Preliminary QC filtering  
  prelim.QCed.seu.objs <- runSeuratPreliminaryFiltering(
    SAMPLE_ID = sample_id,
    geneLvl.seuratObj = unfiltered.seu.objs[[1]],
    isoLvl.seuratObj = unfiltered.seu.objs[[2]],
    MAX.COUNTS = max_counts,
    MT_THRESHOLD = mt_thresh
  )
  
  # Step 3: Normalization and scaling
  normalised.scaled.seu.objs <- NormaliseScaleAndElbow(
    SAMPLE_ID = sample_id,
    seurat.obj = prelim.QCed.seu.objs[[1]],
    seurat.obj.iso = prelim.QCed.seu.objs[[2]],
    regressCellPhase = FALSE,
    labelCellCycle = FALSE
  )
  
  # Step 4: Find optimal clustering resolution
  opt.cluster.res <- find.optimal.cluster.res(
    seurat.obj = normalised.scaled.seu.objs[[1]],
    npc = n_pcs,
    desired.cluster.res = cluster_res,
    commit.to.custom.res = TRUE,
    SAMPLE_ID = sample_id
  )
  
  # Step 5: UMAP and doublet removal
  final.seuObj <- runUMAP.removeDoublets(
    seurat.obj = opt.cluster.res[[1]],
    npc = n_pcs,
    SAMPLE_ID = sample_id,
    cluster.res = cluster_res
  )
  
  # Step 6: Save final object
  output_file <- paste0("processed_", sample_id, ".rds")
  saveRDS(final.seuObj[[1]], file = output_file)
  
  message(paste("Sample", sample_id, "processing complete. Saved to:", output_file))
  
  # Return the final processed object
  return(list(
    sample_id = sample_id,
    seurat_obj = final.seuObj[[1]],
    doublet_plot = final.seuObj[[2]],
    umap_plot = final.seuObj[[3]]
  ))
}

### Main execution loop------
# Initialize storage for results
processed_samples <- list()
processing_summary <- data.frame()

# Create output directory for this batch
output_dir <- paste0("batch_processed_", Sys.Date())
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)

# Process each sample
for (i in 1:nrow(sample_config)) {
  tryCatch({
    message(paste("=== Processing sample", i, "of", nrow(sample_config), "==="))
    
    # Process single sample
    result <- process_single_sample(sample_config[i, ])
    
    # Store result
    processed_samples[[result$sample_id]] <- result
    
    # Track processing summary
    sample_summary <- data.frame(
      sample_id = result$sample_id,
      timepoint = sample_config[i, "timepoint"],
      replicate = sample_config[i, "replicate"],
      final_cell_count = ncol(result$seurat_obj),
      final_gene_count = nrow(result$seurat_obj),
      status = "Success",
      stringsAsFactors = FALSE
    )
    
    processing_summary <- rbind(processing_summary, sample_summary)
    
    # Optional: Clear memory between samples
    gc()
    
  }, error = function(e) {
    message(paste("ERROR processing sample", sample_config[i, "sample_id"], ":", e$message))
    
    # Log the error
    error_summary <- data.frame(
      sample_id = sample_config[i, "sample_id"],
      timepoint = sample_config[i, "timepoint"],
      replicate = sample_config[i, "replicate"],
      final_cell_count = NA,
      final_gene_count = NA,
      status = paste("Error:", e$message),
      stringsAsFactors = FALSE
    )
    processing_summary <- rbind(processing_summary, error_summary)
  })
}

# Save processing summary
write.csv(processing_summary, "processing_summary.csv", row.names = FALSE)
saveRDS(processed_samples, "all_processed_samples.rds")

message("=== Batch processing complete ===")
print(processing_summary)

### Optional: Generate batch summary plots------
# Create comparative plots across all samples
if (nrow(processing_summary[processing_summary$status == "Success", ]) > 1) {
  
  # Cell count comparison
  cell_count_plot <- ggplot(processing_summary[processing_summary$status == "Success", ], 
                           aes(x = sample_id, y = final_cell_count, fill = timepoint)) +
    geom_bar(stat = "identity") +
    labs(title = "Final Cell Counts per Sample", 
         x = "Sample", y = "Cell Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("batch_cell_counts.pdf", cell_count_plot, width = 10, height = 6)
}

### Extract individual processed objects for further analysis------
# Access individual samples like this:
# org1A_processed <- processed_samples[["org1A"]]$seurat_obj
# org1B_processed <- processed_samples[["org1B"]]$seurat_obj
# etc.

# Or loop through to assign them to individual variables
for (sample_name in names(processed_samples)) {
  if (processed_samples[[sample_name]]$sample_id %in% processing_summary$sample_id[processing_summary$status == "Success"]) {
    assign(paste0(sample_name, "_processed"), processed_samples[[sample_name]]$seurat_obj, envir = .GlobalEnv)
    message(paste("Assigned processed object to variable:", paste0(sample_name, "_processed")))
  }
}