# script to look at cell type proportions in seurat object
library(Seurat)
library(tidyverse)
library(gridExtra)
library(speckle)

theme_set(ggmin::theme_min())  # pre-set ggplot theme
##
#   Add cell-type proportion graphs using this: https://academic.oup.com/bioinformatics/article/38/20/4720/6675456
#   https://github.com/Oshlack/speckle
##

# import organoid seurat objects
one_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/one_month_seurat.intergrated_harm.isofrom.rds")
three_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/three_month_seurat.intergrated_harm.isofrom.rds")
six_month_org <- readRDS("/data/gpfs/projects/punim2251/Mitch/six_month_seurat.intergrated_harm.isofrom.rds")
View(one_month_org@meta.data)
# add additional metadata
one_month_org@meta.data$timepoint <- "1M_Org"
three_month_org@meta.data$timepoint <- "3M_Org"
six_month_org@meta.data$timepoint <- "6M_Org"

# function to create cell type proportion plot
MakeGlobalCellTypePropPlot <- function(seu.obj){
  plt <- plotCellTypeProps(clusters = seu.obj@meta.data$cluster_annotations, 
                    sample = seu.obj@meta.data$timepoint)
  
  return(plt)
}

# function to create cell type proportion plot with replicates
MakeReplicateCellTypePropPlot <- function(seu.obj){
  plt <- plotCellTypeProps(clusters = seu.obj@meta.data$cluster_annotations, 
                           sample = seu.obj@meta.data$orig.ident)
  
  return(plt)
}

Org_1M_global <- MakeGlobalCellTypePropPlot(one_month_org)
Org_3M_global <- MakeGlobalCellTypePropPlot(three_month_org)
Org_6M_global <- MakeGlobalCellTypePropPlot(six_month_org)

Org_1M_replicate <- MakeReplicateCellTypePropPlot(one_month_org)
Org_3M_replicate <- MakeReplicateCellTypePropPlot(three_month_org)
Org_6M_replicate <- MakeReplicateCellTypePropPlot(six_month_org)


org.global.prop.plts <- grid.arrange(Org_1M_global,
                                     Org_3M_global,
                                     Org_6M_global,
                                     ncol = 3)
org.replicate.prop.plts <- grid.arrange(Org_1M_replicate,
                                        Org_3M_replicate,
                                        Org_6M_replicate,
                                        ncol = 3)

ggsave("org_global_cellTypeProp.pdf", org.global.prop.plts, width = 30, 
       height = 9, dpi = 300)
ggsave("org_replicate_cellTypeProp.pdf", org.replicate.prop.plts, width = 30, 
       height = 9, dpi = 300)
