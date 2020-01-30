# BIML 2020 Tutorial
# Part 3. Batch correction
# 10x genomics or DropSeq-based data QC
library(devtools)
library(harmony)
# install_github("immunogenomics/harmony") ### run this command, when your R-env don't have harmony library
library(dplyr)

hnscc_umiCounts <- readRDS("/BiO/sample/day1/advanced/BatchCorrection/hnscc_umiCounts.rds")

library(Seurat)
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
hnscc$sample_info <- unlist(lapply(strsplit(colnames(hnscc), "-", fixed=TRUE) ,FUN=function(x){paste(x[1])}))

DimPlot(hnscc, reduction = 'umap', label = F, label.size = 7.5, group.by = 'sample_info')
fig1=DimPlot(hnscc, reduction = 'umap', label = T, label.size = 7.5)
fig2=DimPlot(hnscc, reduction = 'umap', label = T, label.size = 7.5, split.by = 'sample_info')
plot_grid(fig1,fig2)

library(cowplot)
options(repr.plot.height = 2.5, repr.plot.width = 6)

set.seed(123456)
hnscc <- hnscc %>% 
  RunHarmony("sample_info", plot_convergence = F, max.iter.harmony = 100)

hnscc <- hnscc %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(hnscc, reduction = 'umap', label = T, label.size = 7.5, group.by = 'sample_info')
saveRDS(hnscc, file = 'hnscc_harmony.rds')
