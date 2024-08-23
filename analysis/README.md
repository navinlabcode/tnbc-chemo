# Analysis


## Identifying aneuploid cells

| Analysis                             | Main method                                                                                                                                       | Related figures  |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| Identifying aneuploid cells [(link)]() | [CopyKAT (Gao, R. et al. 2021)](https://doi.org/10.1038/s41587-020-00795-2) and [Tirosh, I. et al. 2016](https://doi.org/10.1126/science.aad0501) | Fig. 1, Fig. S1. |

## Cancer-intrinsic archetypes

See detailed instruction: 

| Analysis                                                           | Main method                                                                       | Related figures |
| ------------------------------------------------------------------ | --------------------------------------------------------------------------------- | --------------- |
| Identifying archetypes [(link)]()                                  | non-negative matrix factorization at pseudo-bulk level                            | Fig. 2          |
| Identifying high cancer cell purity spots of Visium [(link)]()      | CopyKAT inferring copy number alterations for Visium spatial transcriptomics data | Fig. 2, Fig. S2 |
| Clinical association of archetypes and response groups [(link)]()  | Contingency table and chi-squared test                                            | Fig. 2          |
| Clinical association of archetypes and overall survival [(link)]() | Cox proportional-hazards model                                                    | Fig. 2          |



## Gene expression metaprograms of cancer cells

| Analysis                                                    | Main method                                            | Related figures |
| ----------------------------------------------------------- | ------------------------------------------------------ | --------------- |
| Identifying metaprograms [(link)]()                         | non-negative matrix factorization at single-cell level | Fig. 3          |
| Determining cellular frequencies of metaprograms [(link)]() | Module score                                           | Fig. 3          |
| Inferring cell phases of cancer cells [(link)]()            | Module score                                           | Fig. 3          |


## Cell states of immune and stromal cell types

| Analysis                                                            | Main method                             | Related figures |
| ------------------------------------------------------------------- | --------------------------------------- | --------------- |
| Identifying cell states of immune and stromal cell types [(link)]() | Differentially expressed genes analysis | Fig. 4          |
| Determining cellular frequencies of cell states [(link)]()          |                                         | Fig. 4          |


## Ecotypes

| Analysis                                            | Main method         | Related figures |
| --------------------------------------------------- | ------------------- | --------------- |
| Identifying ecotypes [(link)]()                     | Community detection | Fig. 5          |
| Comparing ecotypes across patient groups [(link)]() | Z-score             | Fig. 5          |


## Cell-based classifier of predicting response

| Analysis                                                  | Main method         | Related figures |
| --------------------------------------------------------- | ------------------- | --------------- |
| Training and testing the cell-based classifier [(link)]() | Logistic regression | Fig. 6a-c       |


## Gene-based classifier of predicting response

| Analysis                                      | Main method                    | Related figures |
| --------------------------------------------- | ------------------------------ | --------------- |
| Training the gene-based classifier [(link)]() | Logistic regression            | Fig. 6e-g       |
| Testing the gene-based classifier [(link)]()  | Cox proportional-hazards model | Fig. 6e-g       |



## Generic analysis

| Analysis                                                             | Main methods/packages |
| -------------------------------------------------------------------- | --------------------- |
| Gene enrichment analysis                                             | clusterProfiler       |
| Differentially expressed genes analysis for single-cell RNA-seq data | Seurat `wilcox.test`  |
| Differentially expressed genes analysis for bulk RNA-seq data        | DESeq2                |

