<!-- Written by Yun Yan -->

# Identifying cell types of Xenium data

This vignette shows how we identify cell types of Xenium data for each patient. The primary approach is Seurat's 'Label Transfer', which transfers cell labels (i.e., cell types in this vignette) of cells in the scRNA-seq data to the Xenium data. 


**Related figures**:

- Fig. 1g,h


## Rationale

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/xenium/css_celltyping.png?raw=true" width="800">


A common challenge when using the default Label Transfer workflow is that it is often difficult to justify the threshold of prediction probability that determines cell identity assignments. For example, if a cell is predicted to be an endothelial cell with a probability of 0.8, is this 0.8 high enough to be confident? What about 0.9 or 0.7? Another challenge is that it is not straightforward to evaluate the performance of Label Transfer on a given dataset, even though Seurat’s Label Transfer is, in itself, a well-designed and robust algorithm.

To address these issues, we implemented a cross-validation–like strategy for Label Transfer, as illustrated in panel a. 1) We first split the ~5,000 genes into 10 buckets, each containing genes with similar statistical properties of expression in the reference scRNA-seq dataset. 2) We then perform Label Transfer in a 9-fold cross-validation manner: Label Transfer is run 10 times in total, and in each run, 9 buckets of genes are used while leaving out one bucket. 3) After 10 iterations, each cell receives 10 votes for its predicted cell type. A cell must receive all 10 votes for the same cell type to be confidently classified; otherwise, it is designated as “low-confidence” and excluded from further analyses.

Because the Label Transfer algorithm transfers not only categorical labels but also gene expression predictions, the unused (withheld) gene bucket in each run can be used to directly compare the predicted versus actual gene expression levels. This enables a real-time evaluation of Label Transfer performance on each dataset. Our approach is inspired by the benchmarking study by [Li et al, Nature Methods, 2022 (PMID: 35577954)](https://www.nature.com/articles/s41592-022-01480-9), which systematically evaluated computational integration methods for single-cell and spatial transcriptomic data.


We decide to use 10 votes as the threshold, because we find that the cells having <10 votes have drastically smaller probability of prediction, smaller total molecule counts, and smaller total number of detected genes. This pattern is found in different cell types. This pattern is also found in multiple samples. 


## Workflow

**Rscript file path**: 

- <kbd>analysis/scripts/xenium/xenium.prepare_gene_portions.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.prepare_gene_portions.R)). Prepare the 10 gene buckets based on the scRNA-seq reference data. 
- <kbd>analysis/scripts/xenium/xenium.init.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.init.R)). Initiate the Seurat object of Xenium data for each patient sample. 
- <kbd>analysis/scripts/xenium/xenium.celltype.TransferLabel_EvalMode.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.celltype.TransferLabel_EvalMode.R)). Run the cross-validate-like workflow of Label Transfer. 
- <kbd>analysis/scripts/xenium/xenium.celltype.TransferLabel_EvalMode.consensus_decision.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.celltype.TransferLabel_EvalMode.consensus_decision.R)). Collect the 10 runs of the results of Label Transfer and determine cell types. 
- <kbd>analysis/scripts/xenium/xenium.clean_cells.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.clean_cells.R)). Clean cells that have lower votes and are located in outlier space based on the H&E image. 


**Output**

- An analysis-ready Seurat object for each sample. 
- Cell types for cells. 

**Visualization**

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/xenium/css_celltyping_example.png?raw=true" width="800">


$${\color{grey}\text{Written by Yun Yan}}$$
