<!-- Written by Yun Yan -->

# Inferring cell cycling phases of cancer cells

**Related figures**: Fig. 3i

**Rscript file path**: 

- <kbd>analysis/scripts/default_seurat_pipe.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/default_seurat_pipe.R))
- <kbd>analysis/scripts/cell_cycle_phase_comparison.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/cell_cycle_phase_comparison.R))

**Synopsis**

1. Run <kbd>default_seurat_pipe.R</kbd>. Inferring cell-cycle phases is implemented and included in the default scRNA-seq processing pipeline <kbd>default_seurat_pipe.R</kbd>. In brief, it uses `CellCycleScoring` to computationally infer cell cycling phases. 

``` console
Rscript analysis/scripts/default_seurat_pipe.R PATH_TO_SEURAT_OBJECT OUTPUT_DIR \
	MODE_CELL_CYCLE_CORRECTION \
	N_PC N_PC_NN N_HVG \
	ASSAY
```

| Parameter                    | Meaning                                                                                                                                                                                                                    |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `PATH_TO_SEURAT_OBJECT`      | file path to the Seurat object                                                                                                                                                                                             |
| `OUTPUT_DIR`                 | output directory to save all the results                                                                                                                                                                                   |
| `MODE_CELL_CYCLE_CORRECTION` | strategy to perform cell cycling correction. 0='no'. 1='standard workflow'. 2='alternate workflow'. See details in Seurat's vignette [here](https://satijalab.org/seurat/articles/cell_cycle_vignette#alternate-workflow). |
| `N_PC`                       | number of PCA components to perform `RunPCA`                                                                                                                                                                               |
| `N_PC_NN`                    | number of PCA components to perform `RunUMAP` and `FindNeighbors`                                                                                                                                                          |
| `N_HVG`                      | number of highly variable genes (HVG) to perform `FindVariableFeatures`                                                                                                                                                    |
| `ASSAY`                      | name of assay to use.                                                                                                                                                                                                      |

2. Run <kbd>cell_cycle_phase_comparison.R</kbd> to determine frequencies of cell cycling phases and make comparisons of response groups. 

**Output**

| Frequencies of cell cycling phases                                                                                                                                 | Module score of S-phase genes                                                                                                                        | Module score of G2/M genes                                                                                                                             |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cell_cycle/barplot.cell_cycle_proportion.PCR.pdf.png?raw=true" width="500"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cell_cycle/barplot.S.Score.PCR.pdf.png?raw=true" width="250"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cell_cycle/barplot.G2M.Score.PCR.pdf.png?raw=true" width="250"> |


$${\color{grey}\text{Written by Yun Yan}}$$
