<!-- Written by Yun Yan -->

# Identifying cell types in Visium HD data

This vignette shows how we identified cell types (especially cancer cells) in Visium HD data. We primarily used the Robust Cell Type Decomposition (RCTD) method (https://github.com/dmcable/spacexr) to identify cell types. In addition, we used the tool 'Copykat' (https://github.com/navinlabcode/copykat) to infer copy number alterations which further helps identifying aneuploid cells. 

In detail, we performed the following 3 analysis to determine cell types, especially cancer cells. 

- **Analysis 1**. Running RCTD and using our own TNBC scRNA-seq data as the reference data to call cell types. 
- **Analysis 2**. Running RCTD and using the normal human breast cell atlas (HBCA) scRNA-seq data (PMID: 37380767) as the reference dataset to call cell types. 
- **Analysis 3**. Running CopyKat to identify the aneuploid cells. 

A cell is finalized as a 'cancer cell', if it is identified as a 'Tumor' cell in Analysis-1, a 'epithelial' cell in Analysis-2, and/or an 'aneuploid' cell in Analysis-3. 


**Related figures**:

- Extended Data Fig. 2i

**Rscript file path**: 

- <kbd>analysis/scripts/visiumHD/visium_HD.initiate_direct.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visium_HD.initiate_direct.R)). Initiate the Seurat object of reading the Visium HD data of each sample. 
- <kbd>analysis/scripts/visiumHD/visium_HD.prepare.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visium_HD.prepare.R)). Make the analysis-ready object. 
- <kbd>analysis/scripts/visiumHD/visium_HD.identify_celltype.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visium_HD.identify_celltype.R)). Identifying cell types using the scRNA-seq data of either the TNBC or the normal breast tissue as the reference data to run RCTD. 
- <kbd>analysis/scripts/visiumHD/visium_HD.copykat_mix.simple.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visium_HD.copykat_mix.simple.R)). Inferring the copy number alterations and identifying aneuploid cells. 
- <kbd>analysis/scripts/visiumHD/visium_HD.finalize_cancer_cells.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visium_HD.finalize_cancer_cells.R)). Finalizing the identities of cell types especially the cancer cells. 
- <kbd>analysis/scripts/visiumHD/visiumHD.addmodulescore_ForCancerCells.alt.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visiumHD/visiumHD.addmodulescore_ForCancerCells.alt.R)). Calculating module scores of any gene signatures on the cancer cells only. 

**Output**

- An analysis-ready Seurat object of Visium HD data, which contains the cell types identities, especially the cancer cells. 
- Module scores of any gene signatures on the cancer cells. 

**Visualization**

In this tissue example, the cancer cells harbor evident CNA events, and are predicted to be epithelial cells (i.e., luminal HR cells) based on the HBCA dataset and the tumor cells based on the TNBC dataset. 

| H&E image                                                                                                                                    | Cell types using HBCA                                                                                                                                | Cell types using TNBC                                                                                                                                | Aneuploid cells                                                                                                                                | Finalized cancer cells                                                                                                                             |
| -------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/spatial.image.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/spatial.celltype_HBCA.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/spatial.celltype_TNBC.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/spatial.copykat.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/spatial.Tumor_cells.pdf.png?raw=true" width="400"> |


| Heatmap                                                                                                                                                   |
| --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/visiumHD/copykat_heatmap3.celltypes_css.png?raw=true" width="800"> |


$${\color{grey}\text{Written by Yun Yan}}$$
