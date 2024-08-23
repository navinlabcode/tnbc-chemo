# Identifying aneuploid cells

**Related figures**: Fig. S1. 


**Rscript file path**: <kbd>analysis/scripts/copykat_mix.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/copykat_mix.R)). Identifying aneuploid cells is standardized as this Rscript file which runs like a command line tool. 

**Synopsis**

``` terminal
Rscript analysis/scripts/copykat_mix.R PATH_TO_SEURAT_OBJECT OUTPUT_DIR SAMPLE_NAME
```

| Parameter             | Meaning                                                 |
| --------------------- | ------------------------------------------------------- |
| PATH_TO_SEURAT_OBJECT | file path to the Seurat object                          |
| OUTPUT_DIR            | output directory to save all the results                |
| SAMPLE_NAME           | sample name or any text of identifying samples/patients |


**Output**

| Tirosh-based esult                                                                                                                                                     | Leiden-based result                                                                                                                                        | Final-result                                                                                                                                        |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/identifying_aneuploid_cells/copykat_aneuploidy_prediction_tirosh.png?raw=true"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/identifying_aneuploid_cells/copykat_DNA_DR.category.png?raw=true" > | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/identifying_aneuploid_cells/copykat_heatmap4.png?raw=true" > |
| <kbd>Ref</kbd>=the mixed normal cells; <kbd>Obs-diploid</kbd>=a query cell inferred as diploid; <kbd>Obs-tumor</kbd>=a query cell inferred as aneuploid                | <kbd>ref</kbd>=the mixed normal cells; <kbd>obs</kbd>=query cells                                                                                          | Heatmap showing the result of 4 strategies (see details below)                                                                                      |

**Rationale**

Inferring copy number alternations (CNAs) from scRNA-seq data is straightforward by using CopyKAT or InferCNV. We used CopyKAT in our study. However, it is not straightforward to chose a proper threshold to robustly separate aneuploid cells and non-aneuploid cells. For example, over-expressions of nearby genes (e.g., HLA* genes) would lead to spurious focal CNA events. The gene expression baseline of the reference cells would also affect CNAs inference. To address this problem, we mixed normal cells (including normal epithelial, immune, and stromal cells) with each tumor sample to run CopyKAT. Then we refer to 4 strategies to collectively determine the aneuploid cells

1. The default copyKAT method. CopyKAT uses `hclust` to determine aneuploid and non-aneuploid cells. 
2. Tirosh-based method ([Tirosh, I et al 2016](https://doi.org/10.1126/science.aad0501)). It calculates whether a query cell has common CNA events with great magnitude. 
3. Leiden-based method. The mixed normal cells are 'negative control'. Any query cells that are clustered together with these normal cells are expected to be non-aneuploid cells. 
4. Livnat-based method ([Livnat Jerby-Arnon et al 2018](https://doi.org/10.1016/j.cell.2018.09.006)). By leveraging the public TCGA data, it calculates whether a query cell is similar to the breast cancer cell or the normal cell. 



