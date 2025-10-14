library(Seurat)
library(tidyverse)
library(pheatmap)
library(presto)
library(readxl)
library(gridExtra)

setwd("/data/gpfs/projects/punim2251/Mitch")

# Import seurat objects
one_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/six_month_seurat.intergrated_harm.isofrom.rds")

Idents(one_month_org) <- "cluster_annotations"
Idents(three_month_org) <- "cluster_annotations"
Idents(six_month_org) <- "cluster_annotations"

View(six_month_org@meta.data)
# Read in TWAS data
isoform_TWAS_df <- read_excel("/data/gpfs/projects/punim2251/Mitch/FILTERED_Bhattacharya2023 isoTWAS Developmental.xlsx")
isoform_TWAS_byCondition <- isoform_TWAS_df %>% 
  dplyr::select(c("Transcript", "Trait"))
traits <- unique(isoform_TWAS_df$Trait)

# Function that returns a dataframe with DETs overlapping with TWAS dataframe
returnTWAS_DE_Features <- function(seu.obj){
  # calculate DEGs
  cellTypeMarkers_sc <- FindAllMarkers(seu.obj,
                                       assay = "iso",
                                       slot = "counts",
                                       test.use = "wilcox")
  
  # filter for p.adj <= 0.05, then only keep features in the reference TWAS df
  DE_iso_filtered <- cellTypeMarkers_sc %>%
    dplyr::filter(p_val_adj <= 0.05,
                  !str_detect(gene, "Bambu")) %>%
    dplyr::mutate(trimmed_TxID = sub("\\..*", "", gene)) %>% 
    dplyr::filter(trimmed_TxID %in% isoform_TWAS_df$Transcript) %>% 
    left_join(isoform_TWAS_byCondition, by = c("trimmed_TxID" = "Transcript"))
  
  return(DE_iso_filtered)
}

# Create dfs of differentially expressed isoforms that overlap with TWAS transcripts----
TWAS_DE_6M_Org <- returnTWAS_DE_Features(six_month_org)
TWAS_DE_3M_Org <- returnTWAS_DE_Features(three_month_org)
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
    
    # Extract pseudobulked expression values for each cell type
    pseudobulk.Expression.mat <- AverageExpression(seu.obj, 
                                                   features = featureList,
                                                   assays = assay) %>% 
      as.data.frame()
    
    colnames(pseudobulk.Expression.mat) <- sub("^iso\\.", "", colnames(pseudobulk.Expression.mat))
    colnames(pseudobulk.Expression.mat) <- sub("^RNA\\.", "", colnames(pseudobulk.Expression.mat))
    

    # Create metadata containing the timepoint
    cluster_metadata <- data.frame(
      row.names = colnames(pseudobulk.Expression.mat)
    ) %>% 
      dplyr::mutate(Timepoint = sampleID)

  DEG_heatmap <- pheatmap::pheatmap(pseudobulk.Expression.mat,
                                    cluster_rows = TRUE,
                                    show_rownames = FALSE,
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

DE_heatmapData_6M_Orgs <- generate_DE_heatmap(seu.obj = six_month_org,
                            DE_df = TWAS_DE_6M_Org,
                            sampleID = "6M_Orgs")

DE_heatmapData_3M_Orgs <- generate_DE_heatmap(seu.obj = three_month_org,
                                              DE_df = TWAS_DE_3M_Org,
                                              sampleID = "3M_Orgs")

DE_heatmapData_1M_Orgs <- generate_DE_heatmap(seu.obj = one_month_org,
                                              DE_df = TWAS_DE_1M_Org,
                                              sampleID = "1M_Orgs")
class(DE_heatmapData_1M_Orgs)


# Create pdf summary containing heatmaps by TWAS trait
createPDFSummary <- function(isoLvl = T){
  for(trait_name in traits) {
    if(isoLvl) {
      assayType = "iso"
      fileName = paste0(trait_name, "_ISO.png")
    } else {
      assayType = "RNA"
      fileName = paste0(trait_name, "_GENE.png")
    }
    
    DE_heatmapData_6M_Orgs <- generate_DE_heatmap(seu.obj = six_month_org,
                                                  DE_df =   TWAS_DE_6M_Org.split[[trait_name]],
                                                  sampleID = "6M_Orgs",
                                                  assay = assayType)
    
    DE_heatmapData_3M_Orgs <- generate_DE_heatmap(seu.obj = three_month_org,
                                                  DE_df =   TWAS_DE_3M_Org.split[[trait_name]],
                                                  sampleID = "3M_Orgs",
                                                  assay = assayType)
    
    DE_heatmapData_1M_Orgs <- generate_DE_heatmap(seu.obj = one_month_org,
                                                  DE_df =   TWAS_DE_1M_Org.split[[trait_name]],
                                                  sampleID = "1M_Orgs",
                                                  assay = assayType)
    
    # Then arrange all heatmaps
    heatmaps.grouped <- grid.arrange(DE_heatmapData_6M_Orgs[[4]],
                                     DE_heatmapData_3M_Orgs[[4]],
                                     DE_heatmapData_1M_Orgs[[4]], 
                                     ncol = 1)
    ggsave(fileName, heatmaps.grouped, width = 9, height = 24, dpi = 300)
  }
}

createPDFSummary(isoLvl = F)

createPDFSummary(isoLvl = T)

