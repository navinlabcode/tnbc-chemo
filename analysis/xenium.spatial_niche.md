<!-- Written by Yun Yan -->

# Identifying spatial niches using Xenium data


**Related figures**:

- Fig. 5d-k

**Rscript file path**: 

- <kbd>analysis/scripts/spatial_niche/xenium.spatialEcotype.step0.prepare_inputs.winner.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/xenium.spatialEcotype.step0.prepare_inputs.winner.R)). Make the data frame which saves the cell identities and spatial x-y locations.
- <kbd>analysis/scripts/spatial_niche/xenium.spatialEcotype.step1.make_ROIs.winner.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/xenium.spatialEcotype.step1.make_ROIs.winner.R)). For each cell, compute the cell states composition. Radius is 30 µm in our case. 
- <kbd>analysis/scripts/spatial_niche/xenium.spatialEcotype.step2a.initiate_niche_assay.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/xenium.spatialEcotype.step2a.initiate_niche_assay.R)). Create a Seurat object. 
- <kbd>analysis/scripts/spatial_niche/xenium.spatialEcotype.step2b.propose_niche_each_sample.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/xenium.spatialEcotype.step2b.propose_niche_each_sample.R)). Cluster the cell state compositions of single cells to identify niches in each sampple. 
- <kbd>analysis/scripts/spatial_niche/xenium.spatialEcotype.step2c.propose_MetaNiche_from_niches.across_samples.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/xenium.spatialEcotype.step2c.propose_MetaNiche_from_niches.across_samples.R)). Find the niches that are reoccuring in multiple patients. 
- <kbd>analysis/scripts/spatial_niche/util.ConsensusClusterPlus.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/spatial_niche/util.ConsensusClusterPlus.R)). Helper functions. 


**Output**

- Cell states composition of spatial niches. 
- Spatial visualizaiton of spatial niches. 

**Visualization**


Clustering on all cell communities found in all samples. 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/heatmap_niche_score_matrix_all_samples_css_clusters_niche_cluster.cut10.pdf.png?raw=true" width="800">

Heatmap showing the relative cell states composition of spatial niches. 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/heatmap_metaniche_scale_matrix_cut10.pdf.png?raw=true" width="800">

More examples of niches on samples

| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/spatial_dimplot_MetaNiche_10_ART10.pdf.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/spatial_dimplot_MetaNiche_10_ART23.pdf.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/spatial_dimplot_MetaNiche_10_ART236.pdf.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/spatial_dimplot_MetaNiche_10_ART272.pdf.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/spatial_niche/spatial_dimplot_MetaNiche_10_ART304.pdf.png?raw=true" width="200"> |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |




$${\color{grey}\text{Written by Yun Yan}}$$
