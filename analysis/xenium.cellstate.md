<!-- Written by Yun Yan -->

# Identifying cell types of Xenium data

This vignette shows how we identify cell states of Xenium data for each patient. The primary approach is Seurat's 'Label Transfer', which transfers cell labels (i.e., cell states in this vignette) of cells in the scRNA-seq data to the Xenium data. 


**Related figures**:

- Extended Data Fig. 5j-m
- Extended Data Fig. 7j-n


## Rationale

Identifying cell states shares the same rationale and workflow as we identify cell types (See [here](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/xenium.celltype.md#rationale)). 

Moreover, there is a unique caveat: what if the marker genes of a novel cell state identified in the reference scRNA-seq data are missing from the default 5,000-gene panel used by Xenium? In such cases, the Label Transfer strategy — and likely any computational method — is inherently limited and will fail to correctly transfer the corresponding cell labels. Therefore, it is essential to first assess whether the default 5,000 genes in the Xenium panel are sufficient to capture all cell states present in the reference scRNA-seq dataset. If critical marker genes are missing, Xenium provides the option to include up to 100 customized genes, which can be strategically selected to recover the missing cell states. In our case, the default 5000 gene panel is sufficient to recover all the cell states, since the cell clusters are well separated, as illustrated in the following figure. 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/xenium/scRNA_sensitivity.png?raw=true" width="600">

Last, because some cell states may be rare, we take advantage of our large cohorts of Xenium samples and integrate samples for each cell type. Then we applied the Label Transfer to identify cell states for each cell type. 

**Rscript file path**: 

- <kbd>analysis/scripts/xenium/xenium.test_xenium_panel_power_on_sc.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.test_xenium_panel_power_on_sc.R)). Assess whether the default 5,000 genes in the Xenium panel are sufficient to capture all cell states present in the reference scRNA-seq dataset.
- <kbd>analysis/scripts/xenium/xenium.merge_objects.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.merge_objects.R)). Merge samples for each cell type. 
- <kbd>analysis/scripts/xenium/integrate_seurat.xenium.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/integrate_seurat.xenium.R)). Perform the CCA integration of the merged samples for each cell type. 
- <kbd>analysis/scripts/xenium/xenium.cellstate.MapQuery_worker.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.cellstate.MapQuery_worker.R)). Run Label Transfer for each cell type. 
- <kbd>analysis/scripts/xenium/xenium.cellstate.MapQuery.consensus_decision.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.cellstate.MapQuery.consensus_decision.R)). Collect the 10 runs of Label Transfer results to determine cell states for each cell type. 
- <kbd>analysis/scripts/xenium/xenium.cellstate.clean2.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xenium/xenium.cellstate.clean2.R)). QC cleaning the cells having low votes. 


**Output**

- Integrated Seurat object of all samples for each cell type. 
- Cell states. 
- Dotplot and other visualizations. 

**Visualization** 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/xenium/Endo.dotplot2.DEGs.cell_state_paper.pdf.png?raw=true" width="800">


$${\color{grey}\text{Written by Yun Yan}}$$
