# along-vertex
Along-Vertex Analysis for DTI Data: TBI vs PTE Classification

This project implements statistical and machine learning analyses for along-vertex diffusion tensor imaging (DTI) data to classify TBI (Traumatic Brain Injury) and PTE (Post-Traumatic Epilepsy) groups. It includes preprocessing, linear modeling, Elastic Net Regression, and visualization of results.

Requirements
- R version 4.x or higher
- Libraries:
  - dplyr
  - ggplot2
  - glmnet
  - pROC
  - caret
 
Installation
install.packages(c("dplyr", "ggplot2", "glmnet", "pROC", "caret"))

Input data
- `patient_demo.csv`: Demographic data containing columns such as 'subject', 'age', 'sex', and 'outcome label (TBI vs PTE)'.
- `bundles.along.vertex.dti.map_FA_mean.csv`: DTI data for FA values along vertices.
- `outcome_data.csv`: Outcome data with classifications for TBI and PTE groups.

Output
- Significant results CSV files saved in the `models` directory.
- ROC curves saved as PNG files in the `roc` directory.
- Boxplots visualizing significant variables saved in the `plots` directory.

