# BIML 2020 Tutorial
# Part 3. Batch correction
# 10x genomics or DropSeq-based data QC
hnscc_umiCounts <- readRDS("/BiO/sample/day1/advanced/BatchCorrection/hnscc_umiCounts.rds")

library(Seurat)
library(cowplot)

hnscc <- CreateSeuratObject(counts = hnscc_umiCounts)
hnscc <- NormalizeData(hnscc)
hnscc <- FindVariableFeatures(hnscc, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all.genes <- rownames(hnscc)
hnscc <- ScaleData(hnscc, features = all.genes)

# These are now standard steps in the Seurat workflow for visualization and clustering
hnscc <- RunPCA(object = hnscc, verbose = FALSE)
hnscc <- RunUMAP(object = hnscc, dims = 1:30, verbose = FALSE)
hnscc <- FindNeighbors(object = hnscc, dims = 1:20, verbose = FALSE)
hnscc <- FindClusters(object = hnscc, verbose = FALSE, resolution = 1)

# add patients information to metadata slot in seurat object
hnscc$sample_info <- unlist(lapply(strsplit(colnames(hnscc), "-", fixed=TRUE) ,FUN=function(x){paste(x[1])}))

DimPlot(hnscc, reduction = 'umap', label = F, label.size = 7.5, group.by = 'sample_info')
fig1=DimPlot(hnscc, reduction = 'umap', label = T, label.size = 7.5)
fig2=DimPlot(hnscc, reduction = 'umap', label = T, label.size = 7.5, split.by = 'sample_info')
plot_grid(fig1,fig2)

# splitting seurat object by each patients
# Normalization and feature selection (highly variable gene selection) - one more time by each patients
hnscc.list <- SplitObject(hnscc, split.by = "sample_info")
for (i in 1:length(hnscc.list)) {
  hnscc.list[[i]] <- NormalizeData(hnscc.list[[i]], verbose = FALSE)
  hnscc.list[[i]] <- FindVariableFeatures(hnscc.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

# Data integration step (Seurat V3. Integration)
hnscc.anchors <- FindIntegrationAnchors(object.list = hnscc.list, dims = 1:30)
hnscc.integrated <- IntegrateData(anchorset = hnscc.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(hnscc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
hnscc.integrated <- ScaleData(hnscc.integrated, verbose = FALSE)
hnscc.integrated <- RunPCA(hnscc.integrated, npcs = 30, verbose = FALSE)
hnscc.integrated <- RunUMAP(hnscc.integrated, reduction = "pca", dims = 1:30)

hnscc.integrated <- FindNeighbors(object = hnscc.integrated, dims = 1:30, verbose = FALSE)
hnscc.integrated <- FindClusters(object = hnscc.integrated, verbose = FALSE, resolution = 1)

fig1=DimPlot(hnscc.integrated, reduction = 'umap', label = T, label.size = 7.5)
fig2=DimPlot(hnscc.integrated, reduction = 'umap', label = T, label.size = 7.5, split.by = 'sample_info')
DimPlot(hnscc.integrated, reduction = 'umap', label = F, label.size = 7.5, group.by = 'sample_info')

p1 <- DimPlot(hnscc.integrated, reduction = "umap", group.by = "sample_info")
p2 <- DimPlot(hnscc.integrated, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2)


saveRDS(hnscc.integrated, file = 'hnscc.integrated.rds')
