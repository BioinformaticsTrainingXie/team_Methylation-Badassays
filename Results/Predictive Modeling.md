# Predictive Modeling

## Background

Support Vector Machines (SVM) is a powerful supervised learning method used to recognize patterns in data. SVM has been used in many bioinformatics applications for classification, regression, and outlier detection because it avoids overfitting, it can account for nonlinear relationships, and is robust to noise [1]. 
‘Elastic net’ is a regularized regression method that combines the L1 and L2 penalties from LASSO and ridge regularization methods. In this way, Elastic net overcomes some of the limitations of LASSO and ridge. Elastic net regularized regression (GLMnet) has been used to create a multi-tissue age predictor based on DNA methylation before [2], indicating that this might be well-suited for our dataset.
We propose to compare SVM and Elastic net regression methods in generating a classifier that will predict ethnicity based on DNA methylation features (CpGs). 

## Step 1: Filtering

We chose to use an arbitrary threshold to prefilter CpGs to retain. CpGs with a standard deviation (SD) greater than 0.10, leaving 10775 CpGs for model building. We reason that only variable CpGs are likely to be able to be used to distinguish ethnicity, however we are uncertain with the viability of this strategy in setting a threshold.

### Before filtering
![Before filtering:](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-1.png)

### After filtering
![After filtering:](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-2.png)

## Step 2: Cross Validation

In order to estimate the optimal parameters of these models (Elastic net: α, λ; SVM: penalty factor C) and to estimate prediction accuracy, data was randomly divided into 5 sets for cross validation (CV), using a testing and training set. The samples in the training set was used to train the model with a grid of tuning parameters. Five values of α, and λ (25 different combinations) were tested for elastic net. The parameters resulting in the model with the highest ROC on the testing sets were chosen. An outer layer of cross validation (nested) was used to estimate test error.

![](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/Cross%20validation.png)


### Tuning Alpha and Lambda in Elastic Net Regularization

Elastic net regression employs a penalty that is a combination of L1 and L2 norm controlled by α and λ tuning parameters. We wanted to use more L1-regularization to obtain a small panel of biomarkers (α = 0.75, λ = 0.25). We show the relationship between α, λ, and training error (AUC).

![tune alpha lambda](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/examine%20CV-1.png)
![heatmap](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/examine%20CV-2.png)

### Model Performance

The training performance was AUC = 0.981, 0.988 for glmnet and SVM, respectively. However, glmnet was more stable across repeats (we used repeatedcv, repeats = 3) during the estimation of test error (0.977 +- 0.024 vs 0.978 +- 0.05). Despite higher performance of SVM, we chose to build the final model with glmnet. The final tuning parameters used were α = 0.75, λ = 0.25.

![Test error](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/TestError.png)

The final model utilizes 11 CpGs as predictors.

## Step 3: Testing model on unlabeled test data

### Predicting ethnicity

Using the glmnet model, which utilizes 11 predictor CpGs to predict Ethnicity, we assessed the heterogeneity of our unlabeled test data. We found that all Samples were classified as Caucasian, with a probability ranging from 0.57 to 0.82.

### Clustering Analysis

![Train all CpG](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-1.png)
![Train just predictors](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-2.png)

Clustering combined Train and Test data using the predictor CpGs only results in separation of two main clusters, consisting of primarily Asians or Caucasians. The test training samples (unlabeled) primarily fall into the Caucasian cluster.

![Train and test](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/cluster%20both%20test%20and%20train-1.png)

## Step 4: Analyze predictors

### Plot all 11 predictors

![11 predictors](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/plot%20top%2035-1.png)

Next we plot the top and last predictor to see their difference in methylation. 

![1st CpG](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/plotting%20CpGs-1.png)
![11th CpG](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/plotting%20CpGs-2.png)

## Step 5: PCA on test set
Here we look at a PCA on merged test and training sets to show that the first PC is significantly different between the two datasets.

![PCA on merged data](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/PCA%20on%20testtrain.png)

We suspect that the classification is performing poorly on the test data because we doubt that it is truly entirely Caucasian. From this PCA plot on test and train, we can see that the first PC is much different in train vs test (bottom right panel). This indicates that our two datasets are very different, which might make the classifier unsuitable for the test.

## Follow-ups:

(All analysis done below were after the poster session and are explorations of our "future directions".)

### Homogenize Test and Training Data & Re-fit Elastic Net


As shown in the last step, the first PC for merged data differentiates training and test set, which means there are systematic differences between the two sets. We could discard the top PC first before predicting test set. This is done by:

1. Reconstruct the (centered and scaled) merged dataset by discarding the top PC;

2. Re-scale and re-center the correlation matrix for the reconstructed merged data.

After discarding the top PC for merged data, we separate training and test dataset again and use the homogenized training set to re-fit the logistic regression model with elastic net. Resulting probabilities for a test sample to be Asian is shown in the histogram below:

![Histogram on Predicted Asian Probability after Homogenization](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/predict%20redo-1.png)

### Weighted Classes and Up-sampling

We realized after the poster session that another reason our `glmnet` model predicts poorly is because we have imbalanced number of classes. We attempt some methods that remedies this problem and update our results:

* Many of the predictive models for classification have the ability to use case weights where each individual data point can be given more emphasis in the model training phase. One approach to rebalancing the training set would be to increase the weights for the samples in the minority classes. This can be interpreted as having identical duplicate data points with the exact same predictor values. Logistic regression, for example, can utilize "Asian" class weights in this way.

![Histogram on Predicted Asian Probability using weighted glmnet](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/fit%20weighted%20glmnet-1.png)


* Instead of having the model dealing with imbalanced ratio of classes, we can attempt to balance the class frequencies. There are post-hoc sampling approaches that can help attenuate the effects of the imbalance during model training. Two general post hoc approaches are down-sampling and up-sampling the data. Here we will try out up-sampling, which is a technique that simulates additional data points to improve balance across classes.

![Histogram on Predicted Asian Probability with up-sampling](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors_files/figure-markdown_github/fit%20up-sampling%20glmnet-1.png)

Good news is that both weighted glmnet and up-sampling boosted the number of predicted Asians to around 8 samples, which is slightly more believable. Still, more work can be done in the future.

#### References
[1] - Vapnik VN. Statistical Learning Theory. Wiley, New York (1998)

[2] - Horvath S. DNA methylation age of human tissues and cell types. Genome Biol. 2013;14:R115.

[3] - Kuhn, M. (2008). Caret package. Journal of Statistical Software, 28(5)


