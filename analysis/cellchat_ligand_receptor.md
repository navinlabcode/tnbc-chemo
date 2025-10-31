<!-- Written by Yun Yan -->

# Ligand-Receptor inference using CellChat

This vignette uses the ligand-receptor analysis to provide insights of cellular communications of cells within ecotypes. It performs the ligand-receptor inference to infer the cell-to-cell communications. Then it computes the ligand-receptor pairs of cells to compare signals within ecotypes and between ecotypes. 


**Related figures**:
- Extended Data Fig. 8

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cellchat/heatmap.markerLR_per_ecotrait.pdf.png?raw=true" width="400">


**Rscript file path**: 

- <kbd>analysis/scripts/cellchat/cellchat.s0a.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cellchat/cellchat.s0a.R)). Initiate CellChat object by dowsampling cells from the intact Seurat object. 
- <kbd>analysis/scripts/cellchat/cellchat.s1.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cellchat/cellchat.s1.R)). Run the CellChat workflow, including `computeCommunProb`. 
- <kbd>analysis/scripts/cellchat/cellchat.s2.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cellchat/cellchat.s2.R)). Compute the ligand-receptor signals of cells within ecotypes and between ecotypes. 
- <kbd>analysis/scripts/cellchat/cellchat.s2b.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cellchat/cellchat.s2b.R)). Detailed visualization of a specific ligand-receptor among cells. 
- <kbd>analysis/scripts/cellchat/cellchat_func.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cellchat/cellchat_func.R)). Helper functions. 


**Output**

- An analysis-ready CellChat object
- Ligand-receptor pairs that are specific to ecotypes. 



The ligand-receptor signals of cells that are between ecotypes and within ecotypes

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cellchat/barplot.markerLR_per_ecotrait.pdf.png?raw=true" width="800">


Specific examples of ligand-receptor within the cells. 

| Ecotype 3 (EMT-related)                                                                                                                                      | Ecotype 8 (IFN TME 'hot')                                                                                                                              |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cellchat/ecotrait3.LR_SPP1_ITGAV_ITGB1.pdf.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cellchat/ecotrait8.LR_HLA.F_CD8A.pdf.png?raw=true" width="300"> |

$${\color{grey}\text{Written by Yun Yan}}$$
