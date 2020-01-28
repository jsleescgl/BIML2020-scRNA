library(DropletUtils)
library(Seurat)
library(sctransform)
library(dplyr)
library(scRNAseq)
library(ggplot2)

dirExp1='Q:\\BIML2020\\data\\2.Fundamental_pipeline\\exp1\\raw10x'
dirExp2='Q:\\BIML2020\\data\\2.Fundamental_pipeline\\exp2\\raw10x'

# experiment set1
exp1_counts <- counts(sce_filt)

library(Seurat)
rownames(exp1_counts)=gsub('_','-',rownames(exp1_counts))
lung_tme_exp1 <- CreateSeuratObject(counts = exp1_counts)

# Normalization
# Note that this single command replaces NormalizeData, ScaleData, and FindVariableFeatures.
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
lung_tme_exp1 <- SCTransform(object = lung_tme_exp1, verbose = FALSE)
# Perform dimensionality reduction by PCA and UMAP embedding

# These are now standard steps in the Seurat workflow for visualization and clustering
lung_tme_exp1 <- RunPCA(object = lung_tme_exp1, verbose = FALSE)
lung_tme_exp1 <- RunUMAP(object = lung_tme_exp1, dims = 1:30, verbose = FALSE)
lung_tme_exp1 <- FindNeighbors(object = lung_tme_exp1, dims = 1:20, verbose = FALSE)
lung_tme_exp1 <- FindClusters(object = lung_tme_exp1, verbose = FALSE, resolution = 1)

DimPlot(lung_tme_exp1, reduction = 'umap', label = T, label.size = 7.5)+
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


FeaturePlot(lung_tme_exp1, features = c('Ptprc','Cd3d','Cd3e','Cd3g','Cd19','Cd79a','Cd79b','Cd14','Fcgr3','Cd68','Irf8','Itgae','Itgax','Nkg7','S100a9','Cd200r3'), cols = c('#e0e0e0','#b2182b'), ncol = 4)+
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



# singler - using signature based cell type ID (correlation-based)
library(SingleR)
immgen.se <- ImmGenData()
singler_sce=SingleCellExperiment(list(counts=lung_tme_exp1@assays$SCT@counts))
logcounts(singler_sce)=lung_tme_exp1@assays$SCT@data

pred.lung_exp1 <- SingleR(test = singler_sce, ref = immgen.se, labels = immgen.se$label.main)
pred.lung_exp1$pruned.labels

identical(colnames(lung_tme_exp1), rownames(pred.lung_exp1))
lung_tme_exp1$singler_ft_label=pred.lung_exp1$labels
# pred.lung_exp1$singler_pr_label=pred.lung_exp1$pruned.labels

DimPlot(lung_tme_exp1, reduction = 'umap', label = F, label.size = 7.5, group.by = 'singler_ft_label')+
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

saveRDS(lung_tme_exp1, file='Q:\\BIML2020\\data\\mouse_lung_SCLC\\lung_tme_exp1_seuratset.rds')


lung_tme_exp1@meta.data
