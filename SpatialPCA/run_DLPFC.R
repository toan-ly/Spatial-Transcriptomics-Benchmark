library(SpatialPCA)
library(Seurat)
library(ggplot2)
library(Matrix)
library(SPARK)
library(Seurat)
library(peakRAM)


batch_cluster_map <- list(
#   '151669' = 5, '151670' = 5, '151671' = 5, '151672' = 5,
#   '151673' = 7, '151674' = 7, '151675' = 7, '151676' = 7,
#   '151507' = 7, '151508' = 7, 
  '151509' = 7, 
  '151510' = 7
)

input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/DLPFC_spatialpca"

for (sample.name in names(batch_cluster_map)) {
    cat("Processing batch:", sample.name, "\n")
    cluster.number = batch_cluster_map[[sample.name]]
    dir.output <- file.path('/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/Results/DLPFC/SpatialPCA', sample.name)
    if (!dir.exists(dir.output)) {
    dir.create(dir.output, recursive = TRUE)
    }

    load(file.path(input_path, paste0('LIBD_sample', sample.name, ".RData")))
    xy_coords = as.matrix(xy_coords)
    rownames(xy_coords) = colnames(count_sub)

    LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)
    mem <- peakRAM({
        start_time <- Sys.time()
        LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL, sparseKernel = TRUE)
        LIBD = SpatialPCA_EstimateLoading(LIBD,fast=TRUE,SpatialPCnum=20)
        LIBD = SpatialPCA_SpatialPCs(LIBD, fast=TRUE)

        # clusterlabel <- walktrap_clustering(clusternum=as.numeric(cluster.number),latent_dat=LIBD@SpatialPCs,knearest=70 )

        highres_ST <- SpatialPCA_highresolution(LIBD, platform='Visium', newlocation = NULL)

        # cluster_SpatialPCA_high = walktrap_clustering(7, latent_dat=highres_ST@highPCs,200)
        # color_in=c(  "palegreen4", "chocolate1","plum1",  "#F0E442","mediumaquamarine","dodgerblue","lightblue2")
        # title_in="SpatialPCA High resolution"
        # plot_cluster(highres_ST@highPos, as.character(cluster_SpatialPCA_high), pointsize=2,text_size=20 ,title_in,color_in,legend="bottom")

        expr_pred <- highres_ST@W %*% highres_ST@highPCs
        end_time <- Sys.time()
        T = end_time - start_time
    })
    T
    # Save results
    write.csv(expr_pred, file = file.path(dir.output, 'exp_mat.csv'), row.names = TRUE)
    # saveRDS(highres_ST, file = file.path(dir.output, 'resultObject.rds'))
    gc()

}