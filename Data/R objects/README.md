### The Methylation Badassays: R objects

This folder contains R objects that can be loaded into the analysis scripts to avoid long computation times associated with manipulating the large datasets. Please see individual scripts for the use of these objects.

The following R objects in this folder:
  * netFit5predictor.rds - glmnet fitted model resulting in 5 predictors
  * netFit_alpha55_lambda_077_51pred.rds - glmnet model with alpha = 0.55, and lambda = 0.077 (51 predictors)
  * netFit_alpha75_lambda077.rds - glmnet model with alpha = 0.75, and lambda = 0.077
  * netFit_alpha75_lambda25.rds - glmnet model with alpha = 0.75, and lambda = 0.25
  * netFitfinal.rds - glmnet model **load this one for 'BuildModel_AnalyzePredictors.Rmd' script*
  
  * y_predictNet.rds - contains predictions on dataset2 using netFitfinal.rds