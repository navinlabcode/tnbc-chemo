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

### Step 3. Determining marker genes of metaprograms

In brief, for each metaprogram, we compute the average of the gene loading matrix of the included NMF factors. In this way, we create a gene loading matrix in a shape of gene by metaprograms, which applies the same strategy as we determine marker genes for each NMF factor in previous step. 


**Script files**: 

- <kbd>analysis/scripts/metamodule_fnmf.s3.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/metamodule_fnmf.s3.R)). 

**Synopsis**

``` console
Rscript analysis/scripts/metamodule_fnmf.s3.R
```

**Output**

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/metaprogram_marker_gene.png?raw=true" width="400">

## Determining cellular frequencies of metaprograms

**Rscript files**: 

- <kbd>analysis/scripts/xxx.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xxx.R)). 


**Synopsis**

``` console
Rscript 
```

| Parameter             | Meaning                                                 |
| --------------------- | ------------------------------------------------------- |



**Output**

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/yy.png?raw=true" width="400">

---

## Rationale

*Why not use the top 50 genes to define marker genes for each NMF factor?*

Based on the $$W$$ matrix after running NMF , choosing an arbitrary number of top genes (e.g., top 50) for each factor could be an approach to determine marker genes of each factor (i.e., program). But we did not use it because, for example, the top 37th gene in a factor $a$ could be the top 4th gene in another factor $b$, making the 37th gene failing the definition of marker gene for the factor $a$. Besides, it not easy to justify why top 50. What about top 30 or top 100 ? Therefore, we used the strategy as previously described ([Moncada, R. et al. 2020](https://doi.org/10.1038/s41587-019-0392-8)), which make sure the marker genes of each NMF factor (i.e., program) have the highest contribution to the corresponding factor. 


*Why not computationally integrate all cells?*

*Why not perform NMF on all single cells?*

*Why not use consensus NMF ([Kotliar, et al. 2019](https://doi.org/10.7554/eLife.43803))?*

>[!NOTE]
> *A possible limitation of metaprogram analysis*
>
>
>
>
>
>



$${\color{grey}\text{Written by Yun Yan}}$$
