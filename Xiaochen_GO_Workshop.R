library(Seurat)
library(SeuratData)
library(cowplot)
library(Rsamtools)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(Signac)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

clusterProfiler::dotplot()
#SeuratData::InstallData("pbmcMultiome")

# load both modalities
pbmc.rna <- SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- SeuratData::LoadData("pbmcMultiome", "pbmc.atac")

# We use Seurat V5 object
pbmc.rna[["RNA"]] <- as(pbmc.rna[["RNA"]], Class = "Assay5")
# repeat QC steps performed in the WNN vignette
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")

# Process Seurat object
pbmc.rna <- pbmc.rna %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(., dims = 1:30) %>% 
  FindNeighbors(., dims = 1:30) %>% 
  FindClusters(., resolution = 0.5)

# Visualize the RNA data
DimPlot(pbmc.rna, group.by = "seurat_clusters", label = TRUE)

# Find marker genes for each cell subtype
Idents(pbmc.rna) <- "seurat_clusters"
marker.genes.pbmc.rna <- FindAllMarkers(pbmc.rna, only.pos = TRUE, 
                                        min.pct = 0.5, logfc.threshold = 0.5)

Cluster2.markers <- marker.genes.pbmc.rna[marker.genes.pbmc.rna$cluster==2 
                                          & marker.genes.pbmc.rna$p_val_adj < 0.05, ]
Cluster3.markers <- marker.genes.pbmc.rna[marker.genes.pbmc.rna$cluster==3 
                                          & marker.genes.pbmc.rna$p_val_adj < 0.05, ]
Cluster8.markers <- marker.genes.pbmc.rna[marker.genes.pbmc.rna$cluster==8
                                          & marker.genes.pbmc.rna$p_val_adj < 0.05, ]

Cluster2.gene_ids <- bitr(Cluster2.markers$gene, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Cluster3.gene_ids <- bitr(Cluster3.markers$gene, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Cluster8.gene_ids <- bitr(Cluster8.markers$gene, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform GO analysis for cluster 2
Cluster2.ego <- enrichGO(gene = Cluster2.gene_ids$ENTREZID, 
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP", # biological process
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

barplot(Cluster2.ego, showCategory=10, 
        title="GO Enrichment Analysis for Cluster 2")

p1 <- dotplot(Cluster2.ego, showCategory=20,
              title="GO Enrichment Analysis for Cluster 2")
p1 + theme(axis.text.y = element_text(size = 3))

Cluster2.ego <- pairwise_termsim(Cluster2.ego)
emapplot(Cluster2.ego, showCategory=20)

cnetplot(Cluster2.ego, categorySize="pvalue", 
         foldChange=Cluster2.gene_ids$ENTREZID, showCategory = 3)

# Perform GO analysis for all clusters----
cluster_gene_list <- split(marker.genes.pbmc.rna$gene, marker.genes.pbmc.rna$cluster)

cluster_gene_list <- lapply(cluster_gene_list, function(genes) {
  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(gene_ids$ENTREZID)
})

# GO enrichment analysis for BP
go_compare_bp <- compareCluster(geneCluster = cluster_gene_list, 
                                fun = "enrichGO", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "BP",  # "BP"、"MF" or "CC"
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)

go_compare_mf <- compareCluster(geneCluster = cluster_gene_list, 
                                fun = "enrichGO", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "MF",  # "BP"、"MF" or "CC"
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)

go_compare_cc <- compareCluster(geneCluster = cluster_gene_list, 
                                fun = "enrichGO", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "CC",  # "BP"、"MF" or "CC"
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)

# Dot plot for GO terms across all clusters
cluster.p1 <- dotplot(go_compare_bp, showCategory = 1,
                      title = "GO Enrichment (Biological Process)")
cluster.p1 + theme(axis.text.y = element_text(size = 5))

# Great clear plot -> good thing to print on A4 paper to show information across every cluster
trim_bp <- pairwise_termsim(go_compare_bp, showCategory = 50)
treeplot(trim_bp, showCategory = 5)

# Visualise GO enrichment results for other two ontologies-----
# (cellular component and molecular function)
cluster.p2 <- dotplot(go_compare_cc, showCategory = 1,
                      title = "GO Enrichment (Cellular Component)")
cluster.p2 + theme(axis.text.y = element_text(size = 5))
trim_cc <- pairwise_termsim(go_compare_cc, showCategory = 50)

cluster.p3 <- dotplot(go_compare_mf, showCategory = 1,
                      title = "GO Enrichment (Molecular Function)")
cluster.p3 + theme(axis.text.y = element_text(size = 5))

trim_mf <- pairwise_termsim(go_compare_mf, showCategory = 50)
treeplot(trim_mf, showCategory = 5)