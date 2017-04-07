# Predictive Modeling

## Background

Support Vector Machines (SVM) is a powerful supervised learning method used to recognize patterns in data. SVM has been used in many bioinformatics applications for classification, regression, and outlier detection because it avoids overfitting, it can account for nonlinear relationships, and is robust to noise [1]. 
‘Elastic net’ is a regularized regression method that combines the L1 and L2 penalties from LASSO and ridge regularization methods. In this way, Elastic net overcomes some of the limitations of LASSO and ridge. Elastic net regularized regression (GLMnet) has been used to create a multi-tissue age predictor based on DNA methylation before [2], indicating that this might be well-suited for our dataset.
We propose to compare SVM and Elastic net regression methods in generating a classifier that will predict ethnicity based on DNA methylation features (CpGs). 

## Step 1: Filtering

We chose to use an arbitrary threshold to prefilter CpGs to retain. CpGs with a standard deviation (SD) greater than 0.10, leaving 10775 CpGs for model building. We reason that only variable CpGs are likely to be able to be used to distinguish ethnicity, however we are uncertain with the viability of this strategy in setting a threshold.

![Before filtering:](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-1.png)

![After filtering:](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-2.png)

## Step 2: Cross Validation

In order to estimate the optimal parameters of these models (Elastic net: α, λ; SVM: penalty factor C) and to estimate prediction accuracy, data was randomly divided into 5 sets for cross validation (CV), using a testing and training set. The samples in the training set was used to train the model with a grid of tuning parameters. Five values of α, and λ (25 different combinations) were tested for elastic net. The parameters resulting in the model with the highest ROC on the testing sets were chosen. An outer layer of cross validation (nested) was used to estimate test error.

![](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/Cross%20validation.png)






