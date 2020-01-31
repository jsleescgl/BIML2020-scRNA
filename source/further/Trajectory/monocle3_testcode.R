# devtools::install_github('cole-trapnell-lab/leidenbase', force = T)
# devtools::install_github('cole-trapnell-lab/monocle3', force = T)

library(Seurat)
library(scRNAseq)
library(monocle3)
library(scater)
library(leidenbase)
dir("/BiO/sample/day1/advanced/CelltypeID_Trajectory/")
setwd("/BiO/sample/day1/advanced/CelltypeID_Trajectory/")
trajectory_counts <- readRDS("/BiO/sample/day1/advanced/CelltypeID_Trajectory/trajectory_counts.rds")
dim(trajectory_counts)

kidney_trajectory <- readRDS("kidney_trajectory.rds")
kidney_trajectory$region <- unlist(lapply(strsplit(colnames(kidney_trajectory), "-", fixed=TRUE) ,FUN=function(x){paste(x[1])}))

identical(colnames(trajectory_counts), colnames(kidney_trajectory))
head(kidney_trajectory@meta.data)

cell_metadata=data.frame(celltype=kidney_trajectory$cellstate,
                         region=kidney_trajectory$region)
rownames(cell_metadata)=colnames(kidney_trajectory)

gene_annotation=data.frame(gene_short_name=rownames(trajectory_counts))
rownames(gene_annotation)=rownames(trajectory_counts)

library(scRNAseq)
library(SingleCellExperiment)
library(monocle3)

dim(trajectory_counts)
dim(cell_metadata)

cds <- new_cell_data_set(trajectory_counts,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)
plot_cells(cds)

cds <- reduce_dimension(cds, reduction_method="tSNE")

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
colData(cds)$region <- unlist(lapply(strsplit(rownames(colData(cds)), "-", fixed=TRUE) ,FUN=function(x){paste(x[1])}))

plot_cells(cds, color_cells_by="region")

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")

pseudotime_genes <- c("SIX1","SIX2","EYA1",
                      "LHX1","PAX8","FBLN2",
                      "PLA2R1","ARMH4","ANXA1")

pseudotime_genes <- c("SIX1","SIX2","EYA1",
                      "LHX1","PAX8","FBLN2")

plot_cells(cds,
           genes=pseudotime_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds)
colData(cds)
cds@clusters$UMAP$cluster_result$optim_res$membership

identical(rownames(colData(cds)), names(cds@clusters$UMAP$cluster_result$optim_res$membership))
colData(cds)$clusters=cds@clusters$UMAP$cluster_result$optim_res$membership
colData(cds)$clusters_type=paste("Cluster",colData(cds)$clusters,sep = '-')

plot_cells(cds,
           color_cells_by = "clusters_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
