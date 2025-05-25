library(PRECAST)
library(Seurat)
library(ggplot2)
library(bench)
library(clevr)  # For homogeneity, completeness, v-measure
library(cluster) # For silhouette score
library(aricode)
library(mclust)
library(dplyr)
library(BayesSpace)
library(SingleCellExperiment)

calculate_metrics <- function(ground_truth, clusters, data_matrix) {
  tryCatch({
    # Remove NA values
    valid_indices <- !is.na(clusters) & !is.na(ground_truth)
    clusters <- clusters[valid_indices]
    ground_truth <- ground_truth[valid_indices]
    data_matrix <- data_matrix[valid_indices, ]

    if (length(ground_truth) != length(clusters)) {
      stop("Mismatch in length between ground_truth and clusters.")
    }
    if (any(is.na(ground_truth)) || any(is.na(clusters))) {
      stop("Missing values found in ground_truth or clusters.")
    }

    # Convert factors to numeric for clustering metrics
    # if (is.factor(clusters)) {
    #   clusters <- as.numeric(as.character(clusters))
    # }
    # if (is.factor(ground_truth)) {
    #   ground_truth <- as.numeric(as.character(ground_truth))
    # }

    # # Convert to factors
    # ground_truth <- as.factor(ground_truth)
    # clusters <- as.factor(clusters)

    ARI <- adjustedRandIndex(ground_truth, clusters)
    AMI <- AMI(ground_truth, clusters)
    homogeneity <- clevr::homogeneity(ground_truth, clusters)
    completeness <- clevr::completeness(ground_truth, clusters)
    v_measure <- clevr::v_measure(ground_truth, clusters)

    # Calculate Silhouette Score (ASW)
    dist_matrix <- dist(data_matrix)
    sil <- silhouette(clusters, dist_matrix)
    ASW <- mean(sil[, 3])

    # Placeholder for CHAOS and PAS
    CHAOS <- NA  
    PAS <- NA    

    return(list(ARI = ARI, 
                AMI = AMI, 
                Homogeneity = homogeneity, 
                Completeness = completeness, 
                V_Measure = v_measure, 
                ASW = ASW, 
                CHAOS = CHAOS, 
                PAS = PAS))
  }, error = function(e) {
    warning("Error calculating metrics: ", e$message)
    return(list(ARI = NA, 
                AMI = NA, 
                Homogeneity = NA, 
                Completeness = NA, 
                V_Measure = NA, 
                ASW = NA, 
                CHAOS = NA, 
                PAS = NA))
  })
}

save_results <- function(sp_data, metrics_df, pca_data, dir.output) {
  # Save metrics
  write.csv(metrics_df, file = file.path(dir.output, 'metrics.csv'), row.names = FALSE)
          
  cluster_plot <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3) + 
    labs(title=paste("PRECAST (ARI = ", round(metrics_df$ARI, 2), ")", sep="")) +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_blank())
  
  ggsave(file.path(dir.output, 'clustering.pdf'), 
         plot = cluster_plot, width = 6, height = 6, dpi = 300)
  
  # Save reduced dimension data
  write.csv(pca_data,
            file = file.path(dir.output, 'low_dim_data.csv'),
            row.names = TRUE)
  
  # Save metadata
  write.csv(sp_data@meta.data,
            file = file.path(dir.output, 'cell_metadata.csv'),
            row.names = TRUE)
  
  # Save UMAP coordinates
  umap_coords <- as.data.frame(sp_data@reductions$UMAP@cell.embeddings)
  umap_coords$spot_id <- rownames(umap_coords)
  write.csv(umap_coords,
            file = file.path(dir.output, "spatial_umap_coords.csv"),
            row.names = FALSE)

#   umap_plot <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2, color = as.factor(sp_data@meta.data$cluster))) +
#     geom_point(size = 1, alpha = 0.8) +
#     scale_color_brewer(palette = "Set1") +
#     labs(title = "PRECAST", x = "UMAP 1", y = "UMAP 2") +
#     theme(plot.title = element_text(hjust = 0.5, size = 16),
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(color = "black"),
#         legend.title = element_blank())
#   ggsave(file.path(dir.output, 'umap.pdf'), plot = umap_plot,
#          width = 6, height = 6, dpi = 300, device = "pdf")
}


input_path <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/data/BRCA1"
dir.output <- "/Users/toanne/Desktop/Spatial-Transcriptomics-Benchmark/results3/BRCA1/PRECAST"

if (!dir.exists(dir.output)) {
  dir.create(dir.output, recursive = TRUE)
}

sample.name <- 'V1_Human_Breast_Cancer_Block_A_Section_1'
n_clusters <- 20

dir.input <- file.path(input_path, sample.name)



## ensure the row.names of metadata in metaList are the same as that of colnames count matrix
## in countList

githubURL <- "https://github.com/feiyoung/PRECAST/blob/main/vignettes_data/bc2.rda?raw=true"
download.file(githubURL,"bc2.rda",mode='wb')
load("bc2.rda")



bc2 <- list(bc2[[1]])

benchmark <- mark({
    ## Get the gene-by-spot read count matrices countList <- lapply(bc2, function(x)
    ## x[['RNA']]@counts)
    countList <- lapply(bc2, function(x) {
        assay <- DefaultAssay(x)
        GetAssayData(x, assay = assay, slot = "counts")

    })

    M <- length(countList)
    ## Get the meta data of each spot for each data batch
    metadataList <- lapply(bc2, function(x) x@meta.data)

    for (r in 1:M) {
        meta_data <- metadataList[[r]]
        all(c("row", "col") %in% colnames(meta_data))  ## the names are correct!
        head(meta_data[, c("row", "col")])
    }


    ## ensure the row.names of metadata in metaList are the same as that of colnames count matrix
    ## in countList

    for (r in 1:M) {
        row.names(metadataList[[r]]) <- colnames(countList[[r]])
    }
    bc2 <- Filter(function(x) x@project.name == sample.name, bc2)


    ## Create the Seurat list object

    seuList <- list()
    for (r in 1:M) {
        seuList[[r]] <- CreateSeuratObject(counts = countList[[r]], meta.data = metadataList[[r]], project = "BreastCancerPRECAST")
    }

    bc2 <- seuList
    rm(seuList)
    head(meta_data[, c("row", "col")])

    ## Create PRECASTObject.
    set.seed(2022)
    PRECASTObj <- CreatePRECASTObject(bc2, project = "BC2", gene.number = 2000, selectGenesMethod = "SPARK-X",
        premin.spots = 20, premin.features = 20, postmin.spots = 1, postmin.features = 10)

    PRECASTObj <- AddAdjList(PRECASTObj, platform = "Visium")
    PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, verbose = TRUE, maxIter = 30)
    PRECASTObj <- PRECAST(PRECASTObj, K = n_clusters)

    ## backup the fitting results in resList
    resList <- PRECASTObj@resList
    PRECASTObj <- SelectModel(PRECASTObj)

    seuInt <- IntegrateSpaData(PRECASTObj, species = "Human")

    seuInt <- AddUMAP(seuInt, n_comp=2)

    df_meta <- read.table(file.path(dir.input, 'metadata.tsv'), sep = '\t', header=TRUE)
    sp_data <- AddMetaData(seuInt, metadata = df_meta$fine_annot_type, col.name = 'fine_annot_type')

    gt <- sp_data@meta.data$fine_annot_type
    pred <- PRECASTObj@resList$cluster[[1]]
    pca_data <- PRECASTObj@resList$hZ[[1]]     
    metrics <- calculate_metrics(gt, pred, pca_data)
}, iterations=1L)

execution_time <- as.numeric(benchmark$time[[1]]) 
memory_usage <- as.numeric(benchmark$mem_alloc[[1]]) / (1024^2)

metrics_df <- data.frame(
    Sample = sample.name,
    ARI = metrics$ARI,
    AMI = metrics$AMI,
    Homogeneity = metrics$Homogeneity,
    Completeness = metrics$Completeness,
    V_Measure = metrics$V_Measure,
    ASW = metrics$ASW,
    Time = execution_time,
    Memory = memory_usage
)

sp_temp <- Load10X_Spatial(dir.input, filename = "filtered_feature_bc_matrix.h5")
sp_data[['Spatial']] <- sp_temp[['Spatial']]
sp_data@images <- sp_temp@images

save_results(sp_data, metrics_df, pca_data, dir.output)
