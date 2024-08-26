<!-- Written by Yun Yan -->

# Gene expression metaprograms of cancer cells

**Related figures**: Fig. 3, Fig. S3. 

This instruction will address the following two questions: 
1. How to identify meteprograms of cancer cells? 
2. How to determine the cell frequency of metaprograms? 


## Identifying metaprograms of cancer cells

In brief, we performed non-negative matrix factorization (NMF) to cancer cells for each sample. It identified NMF factors representing the heterogeneously expressed gene program in each sample. Then the similar gene programs recurrent in multiple samples were merged to result in metaprograms. 


### Step 1. Identifying heterogeneously expressed program in each sample

**Script files**: 

- <kbd>analysis/scripts/fastnmf.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/fastnmf.R)) 
- <kbd>analysis/scripts/metamodule.fnmf.wrapper.snk.py</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule.fnmf.wrapper.snk.py))

Performing NMF on cancer cells of each sample has been standardized as the Rscript <kbd>fastnmf.R</kbd>which runs as a command line tool. 

**Synopsis**

``` console
Rscript analysis/scripts/fastnmf.R PATH_TO_MATRIX NMF_RANK OUT_DIR 
```

| Parameter      | Meaning                                 |
| -------------- | --------------------------------------- |
| PATH_TO_MATRIX | file path to the expression matrix      |
| NMF_RANK       | a number specifying the rank to run NMF |
| OUT_DIR        | output directory to save results        |


In addition, to facilitate parallel computation, we provide a wrapper written in Python's [snakemake](https://snakemake.github.io/) to apply the above Rscript to samples. For example, the following command simultaneously processes 4 samples in parallel. 

``` console
snakemake -s analysis/scripts/metamodule.fnmf.wrapper.snk.py --cores 4
```

**Output**

The output include NMF results per rank per sample, in addition to several data visualization. The marker genes of each program are determined as previously described ([Moncada, R. et al. 2020](https://doi.org/10.1038/s41587-019-0392-8)), instead of choosing an arbitrary number of top genes. 

| Folder structure listing the results                                                                                                              | Marker genes of a NMF factor example                                                                                                                  |
| ------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/nmf_folders.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/program_example.png?raw=true" width="500"> |


### Step 2. Determining metaprograms

**Rscript files**: 

- <kbd>analysis/scripts/metamodule_fnmf.s1.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s1.R))
- <kbd>analysis/scripts/metamodule_fnmf.s2a.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s2a.R))
- <kbd>analysis/scripts/metamodule_fnmf.s2c.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s2c.R))
- <kbd>analysis/scripts/metamodule_fnmf.s2e.alt.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s2e.alt.R))
- <kbd>analysis/scripts/util.nmf.viz.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/util.nmf.viz.R)). This script contains the required helper functions related to metaprograms. 

**Synopsis**

1. Run <kbd>metamodule_fnmf.s1.R</kbd> to perform basic QC to remove NMF factors that are only found in 1 cell and contain less than 3 marker genes, because these NMF factors are probably noisy or computational artifact of using a high rank. 
2. Run <kbd>metamodule_fnmf.s2a.R</kbd> to create Jaccard similarity matrix of all the NMF factors (i.e., gene programs). 
3. Run <kbd>metamodule_fnmf.s2c.R</kbd> to propose merging similar NMF factors to form metaprograms. It performs hierarchical clustering and uses a series of clustering resolution to merge the similar NMF factors into tentative metaprograms.
4. Run <kbd>metamodule_fnmf.s2e.alt.R</kbd> to merge similar tentative metaprograms and delete noise-like metaprograms. 
5. Re-run <kbd>metamodule_fnmf.s2c.R</kbd> and <kbd>metamodule_fnmf.s2e.alt.R</kbd> to avoid under-clustering and over-clustering to finally determine the metaprograms. 


**Output**

To merge similar NMF factors and delete noisy tentative metaprograms, the R script <kbd>metamodule_fnmf.s2e.alt.R</kbd> outputs the intra-metaprogram and inter-metaprograms similarity. 

In this case of 30 tentative metaprograms, over-clustering is obvious. If the NMF factors of a tentative metaprogram have a low Jaccard similarity between each other, this metaprogram has a low intra-MP similarity, so it is noisy-like and will be removed. If two tentative metaprograms have a high similarity between each other (i.e., a high inter-MP similarity), these two metaprograms should be further merged to avoid over-clustering. 

| Tentative 30 metaprograms                                                                                                                                                                           | Intra-metaprograms similarity                                                                                                                                             | Inter-metaprograms similarity                                                                                                                                                                     |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/heatmap.NMFs_jaccard_to_modules.M30.by_hclust_ward.D2.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/lolipop.intra_MM_similarity.pdf.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/heatmap.meta_module_center_similarity_binarized.jaccard.pdf.png?raw=true" width="300"> |
| Intentionally start with an over-clustered result                                                                                                                                                   | Low intra-MP suggests noise-like MPs                                                                                                                                      | `No`: suggest to keep. `Yes`: suggest to merge                                                                                                                                                    |

After determining the metaprograms, the <kbd>metamodule_fnmf.s2c.R</kbd> and <kbd>metamodule_fnmf.s2e.alt.R</kbd> also generates the final data visualization of metaprograms. 

| Final metaprograms                                                                                                                                                                             | Patient number of metaprograms                                                                                                                                                          |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/heatmap.NMFs_jaccard_to_modules.M13.by_hclust_manual.pdf.png?raw=true" width="500"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/jaccard-manual.barplot.num_samples_per_module.pdf.png?raw=true" width="300"> |



### Step 3. Determining marker genes of metaprograms

In brief, for each metaprogram, we compute the average of the gene loading matrix of the included NMF factors. In this way, we create a gene loading matrix in a shape of gene by metaprograms, which applies the same strategy as we determine marker genes for each NMF factor in previous step. 


**Script files**: 

- <kbd>analysis/scripts/metamodule_fnmf.s3.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s3.R)). 

**Synopsis**

``` console
Rscript analysis/scripts/metamodule_fnmf.s3.R
```

**Output**

Marker genes of metaprograms. Top marker genes have higher contributions to the corresponding metaprogram. 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/metaprogram_marker_gene.png?raw=true" width="600">

## Determining cellular frequencies of metaprograms

In brief, for each metaprogram, we calculate module scores in cancer cells. Then we determine a cell is expressing the meteproram and termed the metaprogram-positive if the module score is great than 0.1. Last, cell frequency of a metaprogram in a sample is the ratio of the number of the metaprogram-positive cancer cells and the total cancer cell number in the sample. 

**Rscript files**: 

- <kbd>analysis/scripts/metamodule_cell_frequency.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_cell_frequency.R)). 

**Synopsis**

``` console
Rscript analysis/scripts/metamodule_cell_frequency.R
```

**Output**

- Module scores of metaprograms in single cancer cells
- Cellular frequencies of metaprograms in each sample
- Comparison of module scores and cell frequencies of metaprograms in pCR v.s. RD. 
- Comparison of cell frequencies of metaprograms across the four archetypes. 

| Module scores of metaprograms in an example                                                                                                                     | Cell frequency differences of metaprograms in pCR v.s. RD                                                                                                                | Cell frequency differences of metaprograms across archetypes                                                                                                         |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/demo_MM.ARTC100.score.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/lolipop.cell_fraction_diff.PCR.pdf.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/heatmap.cell_fraction.nmf4.pdf.png?raw=true" width="300"> |




| Boxplot comparing module scores of M5-IFN in pCR v.s. RD                                                                                                                                         | Boxplot comparing cell frequencies of M5-IFN in pCR v.s. RD                                                                                                                      |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/VlnPlot.bypcr.Signature_M05__Interferon_Seurat.patient.pdf.png?raw=true" height="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/VlnPlot.bypcr.M05__Interferon.cell_pct.pdf.png?raw=true" height="200"> |

Boxplots are also generated for the other metaprograms. M5-IFN is simply an example.


---

## Rationale

***Why not computationally integrate all cancer cells to identify metaprograms?***

Performing the classical 'integration-clustering' strategy on cancer cells could be an approach to identify metaprograms. Yes, it was once used in a paper ([Karaayvaz, M. et al. 2018](https://www.nature.com/articles/s41467-018-06052-0#Fig2)), which considered 'patients' as batches and computationally integrated datasets. However, patients could represent meaningful biological factors, for example, some subtypes of patients. It is not easy to justify whether computational integration over- or under-corrects the batch effects and the biological factors. Therefore, we do not recommend relying on computational integration. 

***Why not perform NMF on all single cells to identify metaprograms?***

Simply performing NMF on all single cells could be an approach to identify metaprograms. However, there are 3 problems. 1) It is not easy to justify the number of ranks of NMF. In theory, NMF is expected to reveal the possibly meaningful biological factors among the cells. If the rank is small (e.g., R=4), NMF would probably identify the most dominating factors separating cancer cells, which could be our archetypes or any other sub-typing categories. With the rank number is getting increased, NMF would tend to pick up other less dominating factors, which could be metaprograms or any other gene pathways. If the rank is very high (e.g., R=1000), many noisy gene signatures are expected to show up. Therefore, the 'sweat point' of ranks is not easy to identify. 2) It takes long time to perform NMF on all cancer cells even with a fixed rank. 3) It requires re-performing NMF on all cells if there are new patient data coming in. Token together, simply performing NMF on all cells has great challenges in scientific rational and practical operation. 


***Why not use the top 50 genes to define marker genes for each NMF factor?***

Based on the $$W$$ matrix after running NMF , choosing an arbitrary number of top genes (e.g., top 50) for each factor could be an approach to determine marker genes of each factor (i.e., program). But we did not use it because, for example, the top 37th gene in a factor $a$ could be the top 4th gene in another factor $b$, making the 37th gene failing the definition of marker gene for the factor $a$. Besides, it not easy to justify why top 50. What about top 30 or top 100 ? Therefore, we used the strategy as previously described ([Moncada, R. et al. 2020](https://doi.org/10.1038/s41587-019-0392-8)), which makes sure the marker genes of each NMF factor (i.e., program) have the highest contribution to the corresponding factor than the other factors. 

***Why not use consensus NMF ?***

Consensus NMF (cNMF) ([Kotliar, et al. 2019](https://doi.org/10.7554/eLife.43803)) appears similar to the metaprogram analysis as there are multiple runs of NMF, however, they are distinct. In general, NMF works like principal component analysis (PCA) to break down an input gene expression matrix (genes by cells) to two matrices: one is in a shape of genes by factors (or components) and the other is factors (or components) by cells. Unlike PCA, the result of NMF is not deterministic, meaning that different runs of NMF using the same input gene expression matrix end up with different results. In contrast, PCA is deterministic. For example, the first PCA components always explains the highest variation existing in the data. Re-running PCA with the same input always have the same results. Therefore, to address the problem of 'randomness' in NMF, consensus NMF is developed, which essentially performs NMF multiple times to find the consensus results of all NMF runs. 

***Why are not the cellular frequencies of metaprograms summed to 100%? ***

In contrast to immune cells that have specific cell state identities, cancer cells could be 'multi-tasking' in expressing multiple gene metaprograms. For example, it could simultaneously express hypoxia-related genes reflecting the local tumor environment and express interferon signaling-related genes to respond to the nearby immune cells. Therefore, in our study, cancer cells are not designated to a specific cell identity. The sum of cellular frequencies of cancer cells expressing metaprograms is not necessarily 100% for a tumor. 


>[!NOTE]
> ***Caveat of metaprogram analysis***
>
>Metaprogram analysis is good at resolving intra-tumoral heterogeneity. The **caveat** of metaprogram analysis is that it is less likely to identify any gene set uniformly expressed by all cells in a sample. For example, an AR-related gene set is uniformly expressed in a high level in all cancer cells in some patients (i.e., the luminal hormone responsive archetype). But the metaprogram analysis does not identify it, which is against a common expectation that metaprogram analysis should identify biologically meaningful gene sets. To clarify, these gene sets are actually diverse at the patient level, representing the **inter**-tumoral heterogeneity, whereas metaprogram analysis is designed to address **intra**-tumoral heerogeneity. This is also why we suggest the **archetype analysis** ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/archetype.md)) to address **inter**-tumor heterogeneity, making these two sets of analysis complementary to each other to better understand cancer cells. 
>
>


$${\color{grey}\text{Written by Yun Yan}}$$
