<!-- Written by Yun Yan -->

# Archetypes and Overall Survival

**Related figures**: Fig. 2i

**Rscript file path**: 

- <kbd>analysis/scripts/KM_METABRIC_ALT.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/KM_METABRIC_ALT.R)). 
- <kbd>analysis/scripts/util.ruok.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/util.ruok.R)). Helper functions. 


**Data**:

- We selected the patients data from the *public* METABRIC cohort that had TNBC and received chemotherapy. See the METABRIC patient IDs [here](https://github.com/navinlabcode/tnbc-chemo/tree/main/other/METABRIC).  

**Synopsis**

Run <kbd>KM_METABRIC_ALT.R</kbd> to calculate [ssGSEA score](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html#2_Introduction) ([Barbie, D. et al. 2009](https://www.nature.com/articles/nature08460)) of the archetypes genes on the pseudo-bulk data of METABRIC cohort and then investigate the association of archetype expressions and overall survival (OS). 

**Output**

| Multivariate Cox Proportional-Hazards Modeling result                                                                                                             |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/ggforest.ssgsea.OS_STATUS.zscored.pdf.png?raw=true" width="800"> |


$${\color{grey}\text{Written by Yun Yan}}$$
