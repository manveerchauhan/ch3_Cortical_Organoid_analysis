
###Radial Glia Cell#### 
FeaturePlot(topup_seurat.merged, features = c("HOPX","TNC", "PTPRZ1", "FAM107A","LIFR", "MOXD1", "PTPRZ1"), min.cutoff = 'q10') | DimPlot(topup_seurat.merged, label = TRUE, label.size = 3) # oRG Cluster  ### cluster 6
FeaturePlot(topup_seurat.merged, features = c("VIM", "NES", "PAX6", "HES1", "SLC1A3", "CDH2", "SOX2"), min.cutoff = 'q10') # vRG  # not conclusive 

#cortical SST interneuron 
FeaturePlot(topup_seurat.merged, features = c("SST", "VIP", "CALB2", "CALB1", "NPY")) | DimPlot(topup_seurat.merged, label = TRUE, label.size = 3) ### cluster 7

#Parvalbumin interneuron
FeaturePlot(topup_seurat.merged, features = c("PVALB", "SLIT2"))

#other interneuron Markers 
FeaturePlot(topup_seurat.merged, features = c("NOS1", "CCK", "SLC32A1"))

#Intermediate progenitor cells
FeaturePlot(topup_seurat.merged, features = c("NEUROD6","PPP1R17", "SOX4", "EOMES", "NEUROG2", "ASCL1", "DLL1", "DCX", "NES", "PAX6", "SOX2", "SOX9", "HES1"))

#Neuronal intermediate progenitor cel
FeaturePlot(topup_seurat.merged, features = c("EOMES","PPP1R17", "NEUROG1"))

#stem cells 
FeaturePlot(topup_seurat.merged, features = c("SOX2","NESTIN", "MYC", "PROM1", "POU5F1", "OLIG2", "PAX6", "SOX1")) # clsuter 1

#glutamtergic markers 
FeaturePlot(topup_seurat.merged, features = c("SLC17A6","SLC17A7", "SLC17A8"), min.cutoff = 'q10') # clsuter 8 some type of Glut precursor i think. TBR1 is expressed ehre to

#microGlial cells : NULL
FeaturePlot(topup_seurat.merged, features = c("CX3CR1","OBIF", "SP1999"), min.cutoff = 'q10') # clsuter 8 some type of Glut precursor i think. TBR1 is expressed ehre to

#Serotonergic Neurons: NULL
FeaturePlot(topup_seurat.merged, features = c("TPH2"), min.cutoff = 'q10') 
#Pyramidal Neurons:: NULL
FeaturePlot(topup_seurat.merged, features = c("CAMK2A", "GRIN2A"), min.cutoff = 'q10') #some expresstion of GRIN2A in cluster 6

#oligodendracytes
FeaturePlot(topup_seurat.merged, features = c("OLIG2", "OLIG1", "MBP", "SPG2", "CNP", "MOG"))

#neuron specific genes
FeaturePlot(topup_seurat.merged, features = c("MAPT", "GABBR2", "KLC1", "CACNA1C", "DOC2A", "RBFOX1", "RBFOX2", "CLCN3", "REST", "SRRM4")) # REST is the opposite 

#Astrocytes
FeaturePlot(topup_seurat.merged, features = c("GFAP", "S100B", "ALDH1L1", "GLAST", "SLC1A2", "SLC1A3"))

#explore sub stpes 
FeaturePlot(topup_seurat.merged, features = c("RELN", "SSTR1", "NPY", "CCK", "CR", "CRH", "LHX6", "SLC6A1"))
