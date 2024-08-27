<!-- Written by Yun Yan -->

# Ecotype analysis

## Determining ecotypes

**Related figures**: Fig. 5, Fig. S8


**Rscript file path**: 

- <kbd>analysis/scripts/ecotype/ecotype_0_create_feature_corr.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_0_create_feature_corr.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_1_define.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_1_define.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_2_igraph_viz.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_2_igraph_viz.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_helper_functions.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_helper_functions.R)): Helper functions related to ecotype analysis. 

**Synopsis**

1. Run <kbd>ecotype_0_create_feature_corr.R</kbd> to create the correlation matrix of the cancer metaprogram and TME cell states. 
1. Run <kbd>ecotype_1_define.R</kbd> to define ecotypes. It also justifies the choices of clustering resolutions for community detection and number of ecotypes. 
1. Run <kbd>ecotype_2_igraph_viz.R</kbd> to create the graph laytout to visualize ecotypes. 


**Output**


| Justify the clustering resolution to define ecotypes                                                                                                          | Alleviate randomness of community detection                                                                                   |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/lineplot_jusitfy_number_ecohubs.pdf.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_correlations.diagnosis_consensus_clustering.pdf.png?raw=true" width="400"> |


| Co-occurrence (correlation matrix)                                                                                                                             | Graph of ecotypes                                                                                                                                                  |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/heatmap.feature_correlations.003.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_correlations.003.pval.pdf.png?raw=true" width="400"> |


## Comparing ecotypes across patient groups


**Related figures**: Fig. 5, Fig. S8


**Rscript file path**: 

- <kbd>analysis/scripts/ecotype/ecotype_3_contexts.archetype.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_3_contexts.archetype.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_3_contexts.response.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_3_contexts.response.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_4_compare_contexts.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_4_compare_contexts.R)). 
- <kbd>analysis/scripts/ecotype/ecotype_helper_functions.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ecotype/ecotype_helper_functions.R)). 


**Synopsis**

1. Run <kbd>ecotype_3_contexts.archetype.R</kbd> to create and visualize graphs of ecotypes of the four archetypes. 
1. Run <kbd>ecotype_3_contexts.response.R</kbd> to create and visualize graphs of ecotypes of the two response groups (pCR vs RD). 
1. Run <kbd>ecotype_4_compare_contexts.R</kbd> to compare the graphs of ecotypes 1) between archetypes and 2) between archetypes and response groups.


**Output**

| Graph of ecotypes in ARC1                                                                                                                                                                       | Graph of ecotypes in ARC2                                                                                                                                                                       | Graph of ecotypes in ARC3                                                                                                                                                                       | Graph of ecotypes in ARC4                                                                                                                                                                       |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges_top.by_archetype.001.pdf.png?raw=true" width="250"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges_top.by_archetype.002.pdf.png?raw=true" width="250"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges_top.by_archetype.003.pdf.png?raw=true" width="250"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges_top.by_archetype.004.pdf.png?raw=true" width="250"> |


| Graph of ecotypes in pCR                                                                                                                                                                     | Graph of ecotypes in RD                                                                                                                                                                      |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges.by_pCR_status.001.pdf.png?raw=true" width="450"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/igraph.feature_zscore_copresence.fixed_edges.by_pCR_status.002.pdf.png?raw=true" width="450"> |



| Graph similarity between archetypes                                                                                                                                             | Graph similarity between archetypes and response groups                                                                                                                            |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/heatmap_graph_similarity_within_archetypes_zscore.pdf.png?raw=true" width="300"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ecotype/heatmap_graph_similarity_response_x_archetype_zscore.pdf.png?raw=true" width="300"> |


$${\color{grey}\text{Written by Yun Yan}}$$
