<!-- Written by Yun Yan -->

# Identifying archetype

**Related figures**: Fig. 2, Fig. S2. 


**Rscript file path**: 
- <kbd>analysis/scripts/psbulk.prepare.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/psbulk.prepare.R)). 
- <kbd>analysis/scripts/psbulk.fastNMF.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/psbulk.fastNMF.R)). 
- <kbd>analysis/scripts/util.nmf.viz.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/util.nmf.viz.R)). Helper functions related to NMF.
- <kbd>analysis/scripts/pkg_enricher.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/pkg_enricher.R)). Helper functions related to gene enrichment analysis. 

**Synopsis**

1. Run `psbulk.prepare.R` to create pseudo-bulk gene expression matrix. 
2. Run `psbulk.fastNMF.R` to run NMF using the pseudo-bulk gene expression matrix and perform various analysis including differentially expressed genes (DEGs) analysis, gene signature enrichment analysis (GSEA) and other analysis. 

**Output**

- Pseudo-bulk gene expression matrices that are log-normalized or normalized by using DESeq2. 
- Justifying the proper number of ranks, i.e., the number of archetypes. 
- Results of NMF with different ranks. 
	- Marker genes of archetypes.
	- Gene enrichment analysis for archetypes. 
	- Dictionary of designating patients to archetypes. 
-  DEGs across archetypes by using DESeq2. 
- Comparing various genes of interest (e.g., cancer cells, clinically druggable genes, epithelial markers, etc) across archetypes. 


| Computationally 'best' rank number                                                                                                                                     | Biologically 'best' rank number                                                                                                                                 |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/survey_ranks.silhouette.on_sample_corr.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/pheatmap.cut.sample_correlation.pdf.png?raw=true" width="400"> |


| Final archetypes                                                                                                                                  | GSEA                                                                                                                                                                         |
| ------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/expression.scaled.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/enricher.heatmap.anno_modules.MSigDBHallmark.pdf.png?raw=true" width="400"> |



| Clinically druggable genes of archetypes                                                                                                                                |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/expression.boxplot.clinical_genes_combo.pdf.png?raw=true" width="400"> |


$${\color{grey}\text{Written by Yun Yan}}$$
