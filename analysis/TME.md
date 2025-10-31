<!-- Written by Yun Yan -->

# Identifying cell states of the TME cell types

This vignette shows how we perform clustering and annotations on the TME cell types to identify the cell states. For discovery, we primarily use 2 approaches: 

- Leveraging the scRNA-seq data of the Human Normal Breast Atlas (HBCA, [PMID: 37380767](https://pubmed.ncbi.nlm.nih.gov/37380767/)) as the reference. 
- Unbiased clustering of our own data: using the [clustree](https://github.com/lazappi/clustree) package and performing DEGs to investigate if under-/over- clustering. If the resulting subclusters share too many DEGs or cannot detect DEGs, it indicates over-clustering. 

For validation, we use the Xenium as described in another vignette [here](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/xenium.cellstate.md). 

**Related figures**:

- Fig. 4

**Rscript file path**: 

- <kbd>analysis/scripts/sc_pp/prepare_sr.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/sc_pp/prepare_sr.R)). Making an analysis-ready Seurat object. 
- <kbd>analysis/scripts/sc_pp/integrate_seurat.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/sc_pp/integrate_seurat.R)). Integrate scRNA-seq data of multiple samples. 


**Synopsis**

```
prepare_sr.R FILE_PATH_TO_SEURAT DIR_RESULT CCC_PLAN N_PC N_PC_NN ASSAY
```

| Parameter           | Meaning                                                                                                                                                                            |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FILE_PATH_TO_SEURAT | file path to the Seurat object                                                                                                                                                     |
| DIR_RESULT          | directory to save all the results                                                                                                                                                  |
| CCC_PLAN            | cell cycling correction plan. 0: do not perform. 1: the standard correction which corrects all the cell cycling phases. 2: the alternative correction which keeps the G1/G0 phase. |
| N_PC                | Number of PCA components.                                                                                                                                                          |
| N_PC_NN             | Number of PCA components to run `RunUMAP` and `FindNeighbors` graph                                                                                                                    |
| ASSAY               | Assay to use. Default: RNA                                                                                                                                                         |

```
integrate_seurat.R FILE_PATH_TO_SEURAT N_PC_NN
```

| Parameter            | Meaning                                                                                                  |
| ------------------- | -------------------------------------------------------------------------------------------------------- |
| FILE_PATH_TO_SEURAT | file path to the merged Seurat object. The `patient` in the meta.data is used to specify the patients. |
| N_PC_NN             | Number of PCA components to run PCA and to run `FindNeighbors` and `RunUMAP`.                            |


**Output**

- An analysis-ready Seurat object, which contains subclusters in a various of clustering resolutions. 
- An integrated Seurat object that is also analysis-ready. 


**Visualization**

We utilize the above scripts repeatedly to identify cell types and also cell states in various tasks. For example, 1) we use `integrate_seurat.R` to integrate our data and the HBCA data for each cell type to propose cell states. 2) we also use `prepare_sr.R` to unbiased cluster our data to investigate the proper clusters with the help of the `clustree`. The `clustree` shows the relationships between the clusters of different clustering resolutions (**panel a**). The cluster-3 at resolution=0.05 splits into the cluster-8 and cluster-14 at the resolution=0.4, which correspond to the cDC1 and mDC cell states (**panel b and c**). 


<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/clustree_demo.png?raw=true" width="900">



$${\color{grey}\text{Written by Yun Yan}}$$
