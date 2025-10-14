library(Seurat)
library(tidyverse)
library(pheatmap)
library(presto)
library(gridExtra)
library(grid)

setwd("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis")

# Import seurat objects
one_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("/data/gpfs/projects/punim2251/Neurodevelopmental_Models_analysis/six_month_seurat.intergrated_harm.isofrom.rds")

Idents(one_month_org) <- "cluster_annotations"
Idents(three_month_org) <- "cluster_annotations"
Idents(six_month_org) <- "cluster_annotations"

# Define the trait-specific filtered database paths
trait_files <- list(
  #ASD = "Updated_Mitch_Filtered_Database/ASD_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  BP = "Updated_Mitch_Filtered_Database/BP_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  #MDD = "Updated_Mitch_Filtered_Database/MDD_filtered_5e-8_Sig_isoTWAS_developmental.csv",
  SCZ = "Updated_Mitch_Filtered_Database/SCZ_filtered_5e-8_Sig_isoTWAS_developmental.csv"
)

# Read in all trait-specific TWAS data and combine
cat("Loading trait-specific GWAS-filtered TWAS databases...\n")
isoform_TWAS_list <- list()
for(trait_name in names(trait_files)) {
  file_path <- trait_files[[trait_name]]
  if(file.exists(file_path)) {
    trait_data <- read.csv(file_path, stringsAsFactors = FALSE)
    cat("Loaded", nrow(trait_data), "records for", trait_name, "\n")
    isoform_TWAS_list[[trait_name]] <- trait_data
  } else {
    cat("Warning: File not found:", file_path, "\n")
  }
}

# Combine all trait-specific data into one dataframe
isoform_TWAS_df <- bind_rows(isoform_TWAS_list)
cat("Combined dataset contains", nrow(isoform_TWAS_df), "total records across", 
    length(unique(isoform_TWAS_df$Trait)), "traits\n")

# Create the by-condition lookup
isoform_TWAS_byCondition <- isoform_TWAS_df %>% 
  dplyr::select(c("Transcript", "Trait"))

# Get unique traits
traits <- unique(isoform_TWAS_df$Trait)
cat("Traits to analyze:", paste(traits, collapse = ", "), "\n")

# Function that returns a dataframe with DETs overlapping with TWAS dataframe
returnTWAS_DE_Features <- function(seu.obj){
  # calculate DEGs
  cat("Running FindAllMarkers for isoform-level analysis...\n")
  cellTypeMarkers_sc <- FindAllMarkers(seu.obj,
                                       assay = "iso",
                                       slot = "counts",
                                       test.use = "wilcox")
  
  cat("Found", nrow(cellTypeMarkers_sc), "total markers\n")
  
  # filter for p.adj <= 0.05 and effect size >= 0.5 log2FC, then only keep features in the reference TWAS df
  DE_iso_filtered <- cellTypeMarkers_sc %>%
    dplyr::filter(p_val_adj <= 0.05,
                  abs(avg_log2FC) >= 0.5,
                  !str_detect(gene, "Bambu")) %>%
    dplyr::mutate(trimmed_TxID = sub("\\..*", "", gene)) %>% 
    dplyr::filter(trimmed_TxID %in% isoform_TWAS_df$Transcript) %>% 
    left_join(isoform_TWAS_byCondition, by = c("trimmed_TxID" = "Transcript"))
  
  cat("After filtering:", nrow(DE_iso_filtered), "DETs overlap with TWAS transcripts\n")
  
  return(DE_iso_filtered)
}

# Create dfs of differentially expressed isoforms that overlap with TWAS transcripts----
cat("\nProcessing 6-month organoids...\n")
TWAS_DE_6M_Org <- returnTWAS_DE_Features(six_month_org)

cat("\nProcessing 3-month organoids...\n")
TWAS_DE_3M_Org <- returnTWAS_DE_Features(three_month_org)

cat("\nProcessing 1-month organoids...\n")
TWAS_DE_1M_Org <- returnTWAS_DE_Features(one_month_org)

# split by TWAS trait
splitDataframeByTrait <- function(TWAS_DE_df){
  TWAS_DE_df_splitByTrait <- TWAS_DE_df %>%
    group_by(Trait) %>%
    group_split() %>%
    setNames(unique(TWAS_DE_df$Trait))
  
  return(TWAS_DE_df_splitByTrait)
}

TWAS_DE_6M_Org.split <- splitDataframeByTrait(TWAS_DE_6M_Org)
TWAS_DE_3M_Org.split <- splitDataframeByTrait(TWAS_DE_3M_Org)
TWAS_DE_1M_Org.split <- splitDataframeByTrait(TWAS_DE_1M_Org)

cat("Data split by trait. Available traits in results:\n")
cat("6M:", names(TWAS_DE_6M_Org.split), "\n")
cat("3M:", names(TWAS_DE_3M_Org.split), "\n")
cat("1M:", names(TWAS_DE_1M_Org.split), "\n")

# Function that collapses the TWAS DE dataframe to gene level list
CollapseIsoDEdf2Gene <- function(TWAS_DE_df){
  IsoGene_TWAS_dict <- isoform_TWAS_df %>% 
    dplyr::select(c("Transcript", "HGNC"))

  TWAS_DE_df <- TWAS_DE_df %>% 
    left_join(IsoGene_TWAS_dict, by = c("trimmed_TxID" = "Transcript"))

  gene_lvl_DE_list <- unique(TWAS_DE_df$HGNC)
  
  return(gene_lvl_DE_list)
}

# Create heatmaps for each condition, celltype and timepoint----
generate_DE_heatmap <- function(seu.obj, 
                                DE_df, 
                                sampleID,
                                assay = "iso"){
    if(assay == "RNA") {
      featureList = CollapseIsoDEdf2Gene(DE_df)
    } else {
      featureList = DE_df$gene
    }
    
    cat("Generating", assay, "heatmap for", sampleID, "with", 
        length(featureList), "features\n")
    
    # Extract pseudobulked expression values for each cell type
    pseudobulk.Expression.mat <- AverageExpression(seu.obj, 
                                                   features = featureList,
                                                   assays = assay) %>% 
      as.data.frame()
    
    colnames(pseudobulk.Expression.mat) <- sub("^iso\\.", "", 
                                               colnames(pseudobulk.Expression.mat))
    colnames(pseudobulk.Expression.mat) <- sub("^RNA\\.", "", 
                                               colnames(pseudobulk.Expression.mat))
    

    # Create metadata containing the timepoint
    cluster_metadata <- data.frame(
      row.names = colnames(pseudobulk.Expression.mat)
    ) %>% 
      dplyr::mutate(Timepoint = sampleID)

  DEG_heatmap <- pheatmap::pheatmap(pseudobulk.Expression.mat,
                                    cluster_rows = TRUE,
                                    show_rownames = TRUE,
                                    annotation = cluster_metadata[, "Timepoint", 
                                                                  drop=FALSE], 
                                    border_color = NA, 
                                    fontsize = 10, 
                                    scale = "row", 
                                    fontsize_row = 10, 
                                    height = 20,
                                    annotation_names_col = FALSE,
                                    angle_col = 45)
  
  return(DEG_heatmap)
}

# Generate sample heatmaps for all data
cat("\nGenerating sample heatmaps for all TWAS-DE overlaps...\n")

DE_heatmapData_6M_Orgs <- generate_DE_heatmap(seu.obj = six_month_org,
                            DE_df = TWAS_DE_6M_Org,
                            sampleID = "6M_Orgs")

DE_heatmapData_3M_Orgs <- generate_DE_heatmap(seu.obj = three_month_org,
                                              DE_df = TWAS_DE_3M_Org,
                                              sampleID = "3M_Orgs")

DE_heatmapData_1M_Orgs <- generate_DE_heatmap(seu.obj = one_month_org,
                                              DE_df = TWAS_DE_1M_Org,
                                              sampleID = "1M_Orgs")

# Create pdf summary containing heatmaps by TWAS trait
createPDFSummary <- function(isoLvl = T, 
                             output_dir = "Updated_GWAS_Filtered_Heatmaps"){
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  for(trait_name in traits) {
    cat("Processing trait:", trait_name, "\n")
    
    # Check if trait exists in all timepoint data
    has_6M <- trait_name %in% names(TWAS_DE_6M_Org.split)
    has_3M <- trait_name %in% names(TWAS_DE_3M_Org.split) 
    has_1M <- trait_name %in% names(TWAS_DE_1M_Org.split)
    
    cat("  Data availability - 6M:", has_6M, "3M:", has_3M, "1M:", has_1M, "\n")
    
    if(!has_6M && !has_3M && !has_1M) {
      cat("  Skipping", trait_name, "- no data available\n")
      next
    }
    
    if(isoLvl) {
      assayType = "iso"
      fileName = file.path(output_dir, paste0(trait_name, "_GWAS_filtered_ISO.png"))
    } else {
      assayType = "RNA"
      fileName = file.path(output_dir, paste0(trait_name, "_GWAS_filtered_GENE.png"))
    }
    
    heatmaps_to_plot <- list()
    
    # Generate heatmaps for available timepoints
    if(has_6M && nrow(TWAS_DE_6M_Org.split[[trait_name]]) > 0) {
      DE_heatmapData_6M_Orgs <- generate_DE_heatmap(seu.obj = six_month_org,
                                                    DE_df = TWAS_DE_6M_Org.split[[trait_name]],
                                                    sampleID = "6M_Orgs",
                                                    assay = assayType)
      heatmaps_to_plot[["6M"]] <- DE_heatmapData_6M_Orgs[[4]]
    }
    
    if(has_3M && nrow(TWAS_DE_3M_Org.split[[trait_name]]) > 0) {
      DE_heatmapData_3M_Orgs <- generate_DE_heatmap(seu.obj = three_month_org,
                                                    DE_df = TWAS_DE_3M_Org.split[[trait_name]],
                                                    sampleID = "3M_Orgs",
                                                    assay = assayType)
      heatmaps_to_plot[["3M"]] <- DE_heatmapData_3M_Orgs[[4]]
    }
    
    if(has_1M && nrow(TWAS_DE_1M_Org.split[[trait_name]]) > 0) {
      DE_heatmapData_1M_Orgs <- generate_DE_heatmap(seu.obj = one_month_org,
                                                    DE_df = TWAS_DE_1M_Org.split[[trait_name]],
                                                    sampleID = "1M_Orgs",
                                                    assay = assayType)
      heatmaps_to_plot[["1M"]] <- DE_heatmapData_1M_Orgs[[4]]
    }
    
    # Only create plot if we have heatmaps
    if(length(heatmaps_to_plot) > 0) {
      cat("  Saving", length(heatmaps_to_plot), "heatmaps for", trait_name, "to", fileName, "\n")
      
      # Arrange available heatmaps
      heatmaps.grouped <- do.call(grid.arrange, c(heatmaps_to_plot, list(ncol = 1)))
      ggsave(fileName, heatmaps.grouped, width = 9, height = 8 * length(heatmaps_to_plot), dpi = 300)
    } else {
      cat("  No valid heatmaps to plot for", trait_name, "\n")
    }
  }
}

# Create trait-specific heatmaps using GWAS-filtered data
cat("\nCreating GWAS-filtered trait-specific heatmaps...\n")


cat("\nGenerating ISOFORM-level heatmaps...\n")
createPDFSummary(isoLvl = TRUE)

cat("\nAnalysis complete! GWAS-filtered heatmaps saved to Updated_GWAS_Filtered_Heatmaps/\n")

# Print summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Input data summary:\n")
for(trait_name in traits) {
  trait_count <- sum(isoform_TWAS_df$Trait == trait_name)
  cat("  ", trait_name, ":", trait_count, 
      "GWAS-filtered records\n")
}

cat("\nDifferential expression overlaps:\n")
cat("6M organoids:", nrow(TWAS_DE_6M_Org), 
    "TWAS-DE overlaps\n")
cat("3M organoids:", nrow(TWAS_DE_3M_Org), 
    "TWAS-DE overlaps\n") 
cat("1M organoids:", nrow(TWAS_DE_1M_Org), 
    "TWAS-DE overlaps\n")



#################

# ===== Write "isovis" CSVs with columns: gene_id, transcript_id, <CellType...> =====
# Values are pseudobulk expression (AverageExpression on 'iso') of isoforms that are DE & TWAS-overlapping.

# Output directory
isovis_outdir <- "Updated_Isovis"
if (!dir.exists(isovis_outdir)) dir.create(isovis_outdir, recursive = TRUE)

# --- helpers to parse feature names like "ENST00000335137.4-GENE1" ---
extract_txid_full <- function(x) sub("-.*$", "", x)                  # keep version; drop -GENE
extract_parent_symbol <- function(x) ifelse(grepl("-", x),
                                            sub(".*-", "", x),
                                            NA_character_)
# try to discover a gene_id column in your TWAS table (optional)
detect_gene_id_col <- function(df) {
  candidates <- c("gene_id", "EnsemblGene", "ENSG", "Ensembl_GeneID", "GeneID", "Ensembl", "Gene") # last resort
  intersect(candidates, names(df))[1]
}
twas_gene_id_col <- detect_gene_id_col(isoform_TWAS_df)  # may be NULL

# Build one timepoint isovis data.frame for a trait (expression, not logFC)
build_isovis_expr_one_timepoint <- function(seu.obj, twas_de_split_list, trait_name, timepoint_label) {
  if (!(trait_name %in% names(twas_de_split_list))) return(NULL)
  de_df <- twas_de_split_list[[trait_name]]
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
  
  # feature names present in the iso assay (e.g., "ENST...-GENE")
  feat_names <- unique(de_df$gene)
  
  # Average expression per cluster (pseudobulk)
  avg <- AverageExpression(seu.obj, features = feat_names, assays = "iso")$iso
  avg <- as.data.frame(avg)
  # clean column names (remove any "iso." prefix)
  colnames(avg) <- sub("^iso\\.", "", colnames(avg))
  
  # add transcript_id (full, with version) and gene_id (from TWAS if available, else parent symbol)
  df <- avg %>%
    tibble::rownames_to_column(var = "feature_name") %>%
    dplyr::mutate(
      transcript_id = extract_txid_full(feature_name),
      parent_symbol = extract_parent_symbol(feature_name)
    )
  
  # try to map to Ensembl gene_id using TWAS table (by transcript)
  # We normalize TWAS Transcript to transcript without version for a tolerant join;
  # if you have exact full IDs in TWAS, this still works because we coalesce.
  twas_map <- isoform_TWAS_df %>%
    dplyr::mutate(
      Transcript_full = Transcript,                              # as provided
      Transcript_nov  = sub("\\..*$", "", Transcript)            # trimmed version
    ) %>%
    {
      if (!is.null(twas_gene_id_col)) {
        dplyr::select(., !!twas_gene_id_col, Transcript_full, Transcript_nov, HGNC)
      } else {
        dplyr::select(., Transcript_full, Transcript_nov, HGNC)
      }
    } %>%
    dplyr::distinct()
  
  df <- df %>%
    dplyr::mutate(transcript_nov = sub("\\..*$", "", transcript_id)) %>%
    # left join twice to increase chances: full match first, then nov
    dplyr::left_join(twas_map, by = c("transcript_id" = "Transcript_full")) %>%
    dplyr::left_join(twas_map, by = c("transcript_nov" = "Transcript_nov"), suffix = c("", ".nov")) %>%
    dplyr::mutate(
      # pick gene_id: prefer Ensembl if available in TWAS; else fall back to parent symbol
      gene_id = dplyr::coalesce(!!rlang::sym(twas_gene_id_col %||% ""), !!rlang::sym(paste0(twas_gene_id_col, ".nov") %||% ""), parent_symbol)
    )
  
  # Final tidy table:
  # gene_id, transcript_id, then only cell-type columns (no timepoint prefix in headers)
  celltype_cols <- setdiff(names(df), c("feature_name","transcript_id","parent_symbol",
                                        "transcript_nov","HGNC","HGNC.nov","Transcript_full",
                                        "Transcript_nov","gene_id",
                                        if (!is.null(twas_gene_id_col)) c(twas_gene_id_col, paste0(twas_gene_id_col,".nov")) else character(0)))
  out <- df %>%
    dplyr::select(gene_id, transcript_id, dplyr::all_of(celltype_cols))
  
  # Optional: add a tiny epsilon to avoid exact zeros (if needed for log-scaling downstream)
  # out[celltype_cols] <- lapply(out[celltype_cols], function(v) ifelse(is.na(v), 0, v))
  
  attr(out, "timepoint_label") <- timepoint_label
  out
}

# Write per-trait, per-timepoint CSVs (columns: gene_id, transcript_id, <CellTypes...>)
write_isovis_for_trait <- function(trait_name) {
  cat("Creating expression isovis CSVs for", trait_name, "...\n")
  
  tab_6M <- build_isovis_expr_one_timepoint(six_month_org,   TWAS_DE_6M_Org.split, trait_name, "6M")
  tab_3M <- build_isovis_expr_one_timepoint(three_month_org, TWAS_DE_3M_Org.split, trait_name, "3M")
  tab_1M <- build_isovis_expr_one_timepoint(one_month_org,   TWAS_DE_1M_Org.split, trait_name, "1M")
  
  lst <- list(`6M`=tab_6M, `3M`=tab_3M, `1M`=tab_1M)
  for (nm in names(lst)) {
    tb <- lst[[nm]]
    if (is.null(tb) || nrow(tb) == 0) next
    # ensure numeric columns are numeric and NAs -> 0 (optional)
    ct_cols <- setdiff(names(tb), c("gene_id","transcript_id"))
    tb[ct_cols] <- lapply(tb[ct_cols], function(v) { v[is.na(v)] <- 0; as.numeric(v) })
    
    outfile <- file.path(isovis_outdir, paste0(trait_name, "_isovis_", nm, ".csv"))
    readr::write_csv(tb, outfile)
    cat("  Saved:", outfile, "with", nrow(tb), "isoforms and", length(ct_cols), "cell types.\n")
  }
  
  # Optional: also write an aggregated (across timepoints) file by taking the max per cell type
  if (!all(vapply(lst, is.null, logical(1)))) {
    # harmonize cell-type columns across timepoints by name, then take rowwise max
    # first, bind rows with a 'timepoint' tag, then aggregate
    bind_and_agg <- function(xlist) {
      xlist <- xlist[!vapply(xlist, is.null, logical(1))]
      if (length(xlist) == 0) return(NULL)
      # full join on gene_id + transcript_id so we don't lose isoforms unique to a TP
      out <- Reduce(function(a,b){
        dplyr::full_join(a, b, by = c("gene_id","transcript_id"), suffix = c(".a",".b"))
      }, xlist)
      
      # find cell-type columns (any non id columns)
      ct_cols <- setdiff(names(out), c("gene_id","transcript_id"))
      # If the same cell type appears multiple times (from different TPs), take rowwise max
      # Collapse by base name (strip any ".a"/".b" suffixes introduced by joins)
      base <- sub("\\.[ab]$", "", ct_cols)
      agg <- out[, c("gene_id","transcript_id"), drop = FALSE]
      for (nm in unique(base)) {
        cols <- ct_cols[base == nm]
        mat  <- as.matrix(out[, cols, drop = FALSE])
        agg[[nm]] <- apply(mat, 1, function(v) suppressWarnings(max(as.numeric(v), na.rm = TRUE)))
        agg[[nm]][!is.finite(agg[[nm]])] <- 0
      }
      agg
    }
    
    agg_tab <- bind_and_agg(lst)
    if (!is.null(agg_tab)) {
      outfile <- file.path(isovis_outdir, paste0(trait_name, "_isovis_aggregated.csv"))
      readr::write_csv(agg_tab, outfile)
      cat("  Saved aggregated:", outfile, "with", nrow(agg_tab), "isoforms and",
          ncol(agg_tab) - 2, "cell types.\n")
    }
  }
}

# Only run for traits present, and for your request focus on SCZ and BP
traits_of_interest <- intersect(c("SCZ", "BP"), traits)
for (tr in traits_of_interest) write_isovis_for_trait(tr)

cat("Done. Isovis CSVs are in:", isovis_outdir, "\n")

