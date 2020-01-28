library(scRNAseq)
library(Seurat)
library(scater)

setwd("Q:/BIML2020/data/3.Further/trajectory_data/kidney/")
# h17w_z1_counts/h17w_z2_counts

h17w_z1_counts <- readRDS("Q:/BIML2020/data/3.Further/trajectory_data/kidney/h17w_z1_counts.rds")
h17w_z2_counts <- readRDS("Q:/BIML2020/data/3.Further/trajectory_data/kidney/h17w_z2_counts.rds")
counts=cbind(h17w_z1_counts,h17w_z2_counts)

dim(counts[rowSums(counts)>0,])
counts=counts[rowSums(counts)>0,]
dim(counts)

# For comparison, we now apply sctransform normalization
kidney <- CreateSeuratObject(counts = counts)

# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
kidney <- SCTransform(object = kidney, verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
kidney <- RunPCA(object = kidney, verbose = FALSE)
kidney <- RunUMAP(object = kidney, dims = 1:30, verbose = FALSE)
kidney <- FindNeighbors(object = kidney, dims = 1:10, verbose = FALSE)
kidney <- FindClusters(object = kidney, verbose = FALSE, resolution = 1)

DimPlot(kidney, reduction = 'umap', label = T, label.size = 7.5)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

NPC_genes=c("TMEM100","CITED1","MEOX1")
podocyte_genes=c("NPHS2","MAFB","PTPRO","CLIC5","PODXL")
prolif_NPC_genes=c("TOP2A","SIX1","TMEM100")
prolif_committed_NP_genes=c("PCNA","LHX1","PAX2")
early_podocytes=c("PAX8","PCDH9","OLFM3")


FeaturePlot(kidney, features = c("CLIC5","MAFB","PODXL"), cols = c('#e0e0e0','#b2182b'), ncol = 3)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

# NPC to podocytes
cells=rownames(kidney@meta.data[kidney@active.ident %in% c(1,3,13,15,9,14,12,11,4,16,2,5),])
trajectory_counts=counts[,cells]
setwd("Q:/BIML2020/data/3.Further/trajectory_data/kidney/")
saveRDS(trajectory_counts, file='trajectory_counts.rds')

# For comparison, we now apply sctransform normalization
kidney_trajectory <- CreateSeuratObject(counts = trajectory_counts)

# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
kidney_trajectory <- SCTransform(object = kidney_trajectory, verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
kidney_trajectory <- RunPCA(object = kidney_trajectory, verbose = FALSE)
kidney_trajectory <- RunUMAP(object = kidney_trajectory, dims = 1:30, verbose = FALSE)
kidney_trajectory <- FindNeighbors(object = kidney_trajectory, dims = 1:10, verbose = FALSE)
kidney_trajectory <- FindClusters(object = kidney_trajectory, verbose = FALSE, resolution = 1)

NPC_genes=c("TMEM100","CITED1","MEOX1")
podocyte_genes=c("NPHS2","MAFB","PTPRO","CLIC5","PODXL")
prolif_NPC_genes=c("TOP2A","SIX1","TMEM100")
prolif_committed_NP_genes=c("PCNA","LHX1","PAX2")
early_podocytes=c("PAX8","PCDH9","OLFM3")

DimPlot(kidney_trajectory, reduction = 'umap', label = T, label.size = 7.5)+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

DotPlot(kidney_trajectory, 
        features = unique(c(NPC_genes, podocyte_genes, prolif_NPC_genes, prolif_committed_NP_genes, early_podocytes)))


# cluster-1,2,3: podocyte
# cluster-0,9: NPC
# cluster-7: prolif-NPC
# cluster-11: e-podocyte
# cluster-6,12,15: committed-NP

kidney_trajectory$cellstate="prolif_comm_NPC"
kidney_trajectory@meta.data[kidney_trajectory@active.ident %in% c(0,9),]$cellstate="NPC"
kidney_trajectory@meta.data[kidney_trajectory@active.ident %in% c(7),]$cellstate="prolif_NPC"
kidney_trajectory@meta.data[kidney_trajectory@active.ident %in% c(6,12,15),]$cellstate="comm_NPC"
kidney_trajectory@meta.data[kidney_trajectory@active.ident %in% c(10,11),]$cellstate="early_podocyte"
kidney_trajectory@meta.data[kidney_trajectory@active.ident %in% c(1:3),]$cellstate="podocyte"


DimPlot(kidney_trajectory, reduction = 'umap', label = T, label.size = 7.5, group.by = 'cellstate')+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

setwd("Q:/BIML2020/data/3.Further/trajectory_data/kidney/")
saveRDS(kidney_trajectory, file = 'kidney_trajectory.rds')