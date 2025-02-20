library(SpatialPCA)
library(ggplot2)
library(Matrix)
library(Seurat)

input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/Human_Ovarian_Cancer"
adata <- Load10X_Spatial(input_path)

adata <- PercentageFeatureSet(adata, "^mt-", col.name = "percent_mito")

adata <- subset(adata, subset = nCount_Spatial > 6000)
gene_counts <- rowSums(GetAssayData(adata, slot = "counts") > 0)
keep_genes <- names(gene_counts[gene_counts >= 10])
adata <- subset(adata, features = keep_genes)

# xy_coords <- adata@images$slice1@coordinates
xy_coords <- GetTissueCoordinates(adata)
xy_coords <- xy_coords[c('imagerow', 'imagecol')]
colnames(xy_coords) <- c('x_coord', 'y_coord')

# count_sub <- adata@assays$Spatial@data
count_sub <- GetAssayData(adata, assay = "Spatial", layer = "counts")
xy_coords <- as.matrix(xy_coords)
rownames(xy_coords) <- colnames(count_sub)

LIBD <- CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project="SpatialPCA", gene.type="spatial",
                               sparkversion="spark", numCores_spark=5, gene.number=3000, customGenelist=NULL,
                               min.loctions=20, min.features=20)
LIBD <- SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ", bandwidth.set.by.user=NULL)
LIBD <- SpatialPCA_EstimateLoading(LIBD,fast=TRUE,SpatialPCnum=20)
LIBD <- SpatialPCA_SpatialPCs(LIBD, fast=TRUE)

highres_ST <- SpatialPCA_highresolution(LIBD, platform='ST', newlocation = NULL)
expr_pred <- highres_ST@W %*% highres_ST@highPCs
# Save results
dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/Human_Ovarian_Cancer/SpatialPCA')
if (!dir.exists(dir.output)) {
  dir.create(dir.output, recursive = TRUE)
}
write.csv(expr_pred, file = file.path(dir.output, 'exp_mat.csv'), row.names = TRUE)
saveRDS(highres_ST, file = file.path(dir.output, 'resultObject.rds'))
