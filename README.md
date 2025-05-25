# Spatial Transcriptomics Benchmark

We benchmark 16 popular computational methods in spatial transcriptomics. This repository focuses on clustering and deconvolution tasks using public datasets. Each method is organized in its own folder with scripts to reproduce the results. 

## 📚 Table of Contents
- [Methods Included](#methods-included)
- [Datasets](#datasets)
- [Evaluation Metrics](#evaluation-metrics)
- [Installation](#installation)
- [Reference](#reference)

## Methods Included
| Category         | Method        | Folder |
|------------------|---------------|----------------|
| Clustering       | [BASS](https://github.com/zhengli09/BASS)    | [methods/BASS](./methods/BASS) |
| Clustering       | [BayesSpace](https://github.com/edward130603/BayesSpace)    | [methods/BayesSpace](./methods/BayesSpace) |
| Clustering       | [GraphST](https://github.com/JinmiaoChenLab/GraphST)       | [methods/GraphST](./methods/GraphST) |
| Clustering       | [STAGATE](https://github.com/zhanglabtools/STAGATE)       | [methods/STAGATE](./methods/STAGATE) |
| Clustering       | [PRECAST](https://github.com/feiyoung/PRECAST)       | [methods/PRECAST](./methods/PRECAST) |
| Clustering       | [SEDR](https://github.com/JinmiaoChenLab/SEDR)          | [methods/SEDR](./methods/SEDR) |
| Clustering       | [SpaGCN](https://github.com/jianhuupenn/SpaGCN)        | [methods/SpaGCN](./methods/SpaGCN) |
| Clustering       | [Seurat](https://github.com/satijalab/seurat)        | [methods/Seurat](./methods/Seurat) |
| Clustering      | [DeepST](https://github.com/JiangBioLab/DeepST)        | [methods/DeepST](./methods/DeepST) |
| Clustering      | [SpatialPCA](https://github.com/shangll123/SpatialPCA)    | [methods/SpatialPCA](./methods/SpatialPCA) |
| Clustering      | [SpaceFlow](https://github.com/hongleir/SpaceFlow)     | [methods/SpaceFlow](./methods/SpaceFlow) |
| Clustering      | [conST](https://github.com/ys-zong/conST)         | [methods/conST](./methods/conST) |
| Clustering      | [stLearn](https://github.com/BiomedicalMachineLearning/stLearn)       | [methods/stLearn](./methods/stLearn) |
| Deconvolution    | [CARD](https://github.com/YMa-lab/CARD)          | [methods/CARD](./methods/CARD) |
| Deconvolution    | [RCTD](https://github.com/dmcable/spacexr)          | [methods/RCTD](./methods/RCTD) |
| Deconvolution    | [SPOTlight](https://github.com/MarcElosua/SPOTlight)     | [methods/SPOTlight](./methods/SPOTlight) |

## Datasets
Please refer to [this github repo](https://github.com/OliiverHu/BenchmarkST_reproducibility/blob/main/docs/source/Data%20availability.rst) for the datasets

## Evaluation Metrics
We evaluate each method based on these metrics:

<table>
    <tr>
      <td colspan="4">Accuracy</td>   
      <td colspan="3">Continuity</td>
      <!-- <td colspan="2">Maker score</td> -->
      <td colspan="2">Scalability</td>
    </tr>
    <tr>
      <td>ARI</td>
      <td>AMI</td>
      <td>Homogeneity</td>
      <td>Completeness</td>
      <td>ASW</td>
      <td>CHAOS</td>
      <td>PAS</td>
      <!-- <td>Moran's I</td> -->
      <!-- <td>Geary's C</td> -->
      <td>Time (s)</td>
      <td>Memory (MB)</td>
    </tr>
</table>


## Installation
We highly recommend following the official Github repository of each method for installation instructions

## Reference
We also reference from the following public benchmarking repos:

- [BenchmarkST](https://github.com/maiziezhoulab/BenchmarkST
)
- [STdeconv_benchmark](https://github.com/leihouyeung/STdeconv_benchmark)
- [SDMBench](https://github.com/zhaofangyuan98/SDMBench)