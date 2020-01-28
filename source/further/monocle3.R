.libPaths( c( "C:\\Users\\jsleescgl\\tools" , .libPaths() ) )
.libPaths()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install('scRNAseq')

library(Seurat)
library(scRNAseq)
library(monocle3)

setwd("Q:/BIML2020/data/3.Further/trajectory_data/kidney/")
# h17w_z1_counts/h17w_z2_counts

counts=cbind(h17w_z1_counts,h17w_z2_counts)
dim(counts)
dim(counts[rowSums(counts)>0,])

library(scater)
counts=counts[rowSums(counts)>0,]
dim(counts)

cell_metadata=data.frame(celltype=rep('unknown',8519))
rownames(cell_metadata)=colnames(counts)

gene_annotation=data.frame(gene_short_name=rownames(counts))
rownames(gene_annotation)=rownames(counts)


cds <- new_cell_data_set(as(counts, "sparseMatrix"),
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
