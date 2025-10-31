<!-- Written by Yun Yan -->

# Building the cell abundance-based classifier of predicting response


<img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/MLcell/scML_workflow.png?raw=true" width="700">


**Related figures**:

- Fig. 6a-c

**Rscript file path**: 

- <kbd>analysis/scripts/scML/scML_classifier.v3.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/scML/scML_classifier.v3.R)). 


**Output**

- The cell abundance-based models including the logistic regression model, the random forest model, and the Latent Dirichlet Allocation (LDA) model. The logistic regression model was used in the manuscript. 
- Feature importance of the cell states. 

| Training                                                                                                                                      | Testing                                                                                                                                         |
| --------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/MLcell/roc_zoo.training.pdf.png?raw=true" width="400"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/MLcell/roc_zoo.validation.pdf.png?raw=true" width="400"> |

$${\color{grey}\text{Written by Yun Yan}}$$
