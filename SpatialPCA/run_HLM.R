library(SpatialPCA)
library(ggplot2)
library(Matrix)
library(Seurat)

start_time <- Sys.time()

adata <- Load10X_Spatial('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/Human_Lymph_Node')

adata <- PercentageFeatureSet(adata, "^mt-", col.name = "percent_mito")

adata <- subset(adata, subset = nCount_Spatial > 6000)
gene_counts <- rowSums(GetAssayData(adata, slot = "counts") > 0)
keep_genes <- names(gene_counts[gene_counts >= 10])
adata <- subset(adata, features = keep_genes)

xy_coords <- adata@images$slice1@coordinates
xy_coords <- xy_coords[c('imagerow', 'imagecol')]
colnames(xy_coords) <- c('x_coord', 'y_coord')

count_sub <- adata@assays$Spatial@data
print(dim(count_sub))  # The count matrix
xy_coords <- as.matrix(xy_coords)
rownames(xy_coords) <- colnames(count_sub)