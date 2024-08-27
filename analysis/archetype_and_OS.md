<!-- Written by Yun Yan -->

# Archetypes and Overall Survival

**Related figures**: Fig. 2i

**Rscript file path**: 

- <kbd>analysis/scripts/KM_METABRIC_ALT.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/KM_METABRIC_ALT.R)). 
- <kbd>analysis/scripts/util.ruok.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/util.ruok.R)). Helper functions. 
**Data**:

- We selected the patients data from the *public* METABRIC cohort that had TNBC and received chemotherapy. See the METABRIC patient IDs [here](https://github.com/navinlabcode/tnbc-chemo/tree/main/other/METABRIC).  

**Synopsis**

Run `KM_METABRIC_ALT.R` to calculate GSVA score of the archetypes genes on the pseudo-bulk data of METABRIC cohort and then investigate the association of archetype expressions and overall survival (OS). 

**Output**

| Multivariate Cox Proportional-Hazards Modeling result                                                                                                             |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/ggforest.ssgsea.OS_STATUS.zscored.pdf.png?raw=true" width="600"> |



**Rationale**


$${\color{grey}\text{Written by Yun Yan}}$$
