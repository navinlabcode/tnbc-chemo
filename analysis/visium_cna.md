<!-- Written by Yun Yan -->

# Identifying high tumor purity Visium spots

**Related figures**: Fig. 2, Fig. S2

**Rscript file path**: 
- <kbd>analysis/scripts/visium_st/1copykat.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/1copycat.R)). 
- <kbd>analysis/scripts/visium_st/2tumor_copykat_signal.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/2tumor_copykat_signal.R)). 
- <kbd>analysis/scripts/visium_st/3calculate_metatraits_score_byAUC.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/3calculate_metatraits_score_byAUC.R)). 
- <kbd>analysis/scripts/visium_st/4viz_features_across_samples.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/4viz_features_across_samples.R)). 
- <kbd>analysis/scripts/visium_st/5vln_test_features_across_conditions.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/5vln_test_features_across_conditions.R)). 
- <kbd>analysis/scripts/visium_st/scRNA_packages.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/visium_st/scRNA_packages.R)). Helper functions.

**Synopsis**

1. Run <kbd>1copykat.R</kbd> to use CopyKAT to infer copy number alterantions for the Visium ST data of each sample.
2. Run <kbd>2tumor_copykat_signal.R</kbd> to identify the high tumor purity spots. 
3. Run <kbd>3calculate_metatraits_score_byAUC.R</kbd> to quantify the AUCell score (similar to module scores) of archetype marker genes. 
4. Run <kbd>4viz_features_across_samples.R</kbd> and <kbd>5vln_test_features_across_conditions.R</kbd> to visualize and statistically compare archetype module scores. 

**Output**

| Justify the cutoff of CNA strength to propose the high tumor purity spots                                                                                    |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/msq_box_for_cutoff_selection.pdf.png?raw=true" width="600"> |


| CNA inference results in a sample                                                                                                                   | Spatial locations of the high tumor purity spots                                                                                                 |
| --------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/copykat_htmap_ART23.pdf.png?raw=true" width="500"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/ART23_tumor_area.pdf.png?raw=true" width="500"> |


| Comparing archetype module scores in a sample                                                                                                |
| -------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/ARC_ST_ART23.pdf.png?raw=true" width="800"> |


$${\color{grey}\text{Written by Yun Yan}}$$
