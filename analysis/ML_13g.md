<!-- Written by Yun Yan -->

# The 13-gene classifier of predicting risk score

**Related figures**:

- Fig. 6d-g


## Building the 13-gene classifier

**Rscript file path**: 

- <kbd>analysis/scripts/ML13g/build_ML_model.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ML13g/build_ML_model.R)). 


**Output**

- The 13 genes that are finally used to build the model. 
- The major cell source of the 13 genes.  
- The result of the logistic regression model including the beta-efficient of the genes.


## Testing the 13-gene classifier


**Rscript file path**: 

- <kbd>analysis/scripts/ML13g/ML_13gene.apply_to_others.response.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ML13g/ML_13gene.apply_to_others.response.R)). Apply the model to the external cohorts (BrighTNess and I-SPY2) that have chemotherapy response data. 
- <kbd>analysis/scripts/ML13g/ML_13gene.apply_to_others.survival.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/ML13g/ML_13gene.apply_to_others.survival.R)). Apply the model to the external cohorts (METABRIC and SCAN-B) that have overall survival data. 


**Output**
- The predicted risk scores of non-responding to chemotherapy for each patient in the cohort. 
- Boxplot comparing the risk score between the responders and non-responders. 
- Survival analysis showing the association of the risk scores and overall survival. 

<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/ML13g/yyy.png?raw=true" width="400">


$${\color{grey}\text{Written by Yun Yan}}$$
