### This folder contains scripts for predictive modeling

We do two things within our predictive modeling analysis: 

1. [Compare](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/CV_TestPerformance.md) logistic regression with a elastic net regularization vs support vector machines
2. [Generate](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors.md) the final classifer with glmnet and analyze the results

Because of long computation time when training the models and running cross validation, we split the code into different scripts. Here is a list of them and a short description:

* [CV_TestPerformance.md](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/CV_TestPerformance.md): Here, we compare SVM and glmnet using a repeated nested cross validation strategy (folds = 5, repeats = 3). This script takes a very long time to run (~1-3hours).

* [BuildModel_AnalyzePredictors.Rmd](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors.md): This script is for training the final model and running the prediction on the unlabelled test set. We also do some small exploratory analyses on the predictors and compare the two datasets. Some of the analysis of the predictors is also in [Functional Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/Functional%20Analysis%20Write%20Up.md).