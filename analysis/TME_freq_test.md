<!-- Written by Yun Yan -->

# Comparing cell abundances in pCR and RD patients

We implements a standardized pipeline `std_cell_fraction_test.A_vs_B.R` to compare cell states percentages between the pCR and RD patients. 


**Related figures**:

- Fig. 4

**Rscript file path**: 

- <kbd>analysis/scripts/sc_pp/std_cell_fraction_test.A_vs_B.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/sc_pp/std_cell_fraction_test.A_vs_B.R)). 

**Synopsis**

```
std_cell_fraction_test.A_vs_B.R PATH_TO_DATA_FRAME N_CELLS_AT_LEAST ON_WHAT BY_WHAT PATIENT_IND A B
```

| Parameter          | Meaning                                                                                                   |
| ------------------ | --------------------------------------------------------------------------------------------------------- |
| PATH_TO_DATA_FRAME | data frame of single cells                                                                                |
| N_CELLS_AT_LEAST   | Exclude the sample if it has less than `N_CELLS_AT_LEAST` cells.                                          |
| ON_WHAT            | Cell label. It could be cell type, or cell state.                                                         |
| BY_WHAT            | Conditions of comparison. It could be 'therapy response' (pCR vs RD), or 'TIL level group' (low vs high). |
| PATIENT_IND        | Column name to specify patients. Default: 'patient'                                                       |
| A                  | Condition A. e.g., 'pCR' .                                                                                |
| B                  | Condition B. e.g., 'non-pCR'.                                                                             |



**Output**

- Test results with p-value corrected by the BH-method (i.e., FDR). 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/wilcoxon.png?raw=true" width="400">

- Boxplot

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/USE.boxplot.test.compact.pdf.Mye.100.pdf.png?raw=true" width="1000">

- Lolipop plot

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/USE.lolipop.cell_fraction_diff.PCR.colored.pdf.Mye.100.pdf.png?raw=true" width="400">

- Barplot

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/barplot.summaried.frac.pdf.png?raw=true" width="300">

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/tme_cellstates/barplot.frac.pdf.png?raw=true" width="600">

---

This pipeline performs many other statistical tests in additon to the default Wilcoxon test: 

- two-sided Kolmogorov-Smirnov test ([Jeong Seok Lee et al, Sci Immunol . 2020, PMID: 32651212](https://www.science.org/doi/10.1126/sciimmunol.abd1554))
- dunn test
- t test
- kruskal test

$${\color{grey}\text{Written by Yun Yan}}$$
