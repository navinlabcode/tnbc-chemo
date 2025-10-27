<!-- Written by Yun Yan -->

# Predict archetypes for external TNBC patient cohorts

This script uses the patient data of our study as the internal control and mixes the external patient cohorts with them to run the same NMF procedure as described in the tutorial ['Identifying archetype'](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/archetype.md). 


**Related figures**:

- Fig. 2h
- Extended Data Fig. 2e


**Rscript file path**: 

- <kbd>scripts/archetype/predict_in_other_cohorts/predict_archetypes.other_cohorts.R</kbd> ([link](https://github.com/navinlabcode/tnbc-chemo/blob/main/analysis/scripts/archetype/predict_in_other_cohorts/predict_archetypes.other_cohorts.R)). 

**Synopsis**

``` Bash
Rscript predict_archetypes.other_cohorts.R STUDY_NAME TREATMENT_ARM
```

Parameters | Values
STUDY_NAME | Specify the cohort. It supports the curated TNBC cohorts of METABRIC, SCAN-B, ISPY-2, and BrighTNess. 
TREATMENT_ARM | Specifity the name of treatment arm in the cohort. 

**Output**

1. The predicted 4 NMF groups.
2. The contingency table between the predicted NMF groups and the existing archetype assignments of the internal patient data of our own cohort. 
3. The correctness of the prediction. 

For example, using the BrighTNess cohort, the unbiased analysis of NMF identified 4 groups F1-F4. 


Using the internal control patients, F1-F4 highly concordant with the archetypes ARC 1-4. F1 should be ARC1, F2 should be ARC2, F3 should be ARC3, and F4 should be ARC4. The correctness of this prediciton is (23+28+25+14)/97 = ~93%. 

| groups | ARC1 | ARC2 | ARC3 | ARC4 |
| ------ | ---- | ---- | ---- | ---- |
| F1     | 23   | 4    | 3    | 0    |
| F2     | 0    | 28   | 0    | 0    |
| F3     | 0    | 0    | 25   | 0    |
| F4     | 0    | 0    | 0    | 14   |


By this step, each patient in the external cohort has been assigned to archetypes. Therefore, it is straighforward to perform further comparisons (e.g., archetypes v.s. NAC response; archetypes v.s. archetypes expressions) as described in previous tutorials.

| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/arc_pred_brightness.response.png?raw=true" width="200"> | <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/rchetype/arc_pred_brightness.response.chitest.png?raw=true" width="200"> |
| -------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |

| <img src="https://github.com/navinlabcode/tnbc-chemo/blob/main/website_images/analysis/archetype/arc_pred_brightness.boxplot?raw=true" width="400"> |
| --------------------------------------------------------------------------------------------------------------------------------------------------- |



$${\color{grey}\text{Written by Yun Yan}}$$
