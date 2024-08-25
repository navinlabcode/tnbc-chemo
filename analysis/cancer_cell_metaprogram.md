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


In addition, to facilitate parallel computation, we provide a wrapper written in Python's snakemake to call the above Rscript. 

``` console
snakemake -s analysis/scripts/metamodule.fnmf.wrapper.snk.py --cores 4
```

**Output**

The output include NMF results per rank per sample, in addition to several data visualization. The marker genes of each program are determined ([Moncada, R. et al. 2020](https://doi.org/10.1038/s41587-019-0392-8))). 

| Folder structure listing the results                                                                                                              | Marker genes of a NMF factor example                                                                                                                  |
| ------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/nmf_folders.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/cancer_metaprograms/program_example.png?raw=true" width="500"> |


### Step 2. Determining metaprograms

**Script files**: 

- <kbd>analysis/scripts/xxx.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xxx.R)). 
- 
**Synopsis**

``` console
Rscript 
```

| Parameter             | Meaning                                                 |
| --------------------- | ------------------------------------------------------- |

**Output**


### Step 3. Determining marker genes of metaprograms

**Script files**: 

- <kbd>analysis/scripts/xxx.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/xxx.R)). 
- 
**Synopsis**

``` console
Rscript 
```

| Parameter             | Meaning                                                 |
| --------------------- | ------------------------------------------------------- |


**Output**


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
