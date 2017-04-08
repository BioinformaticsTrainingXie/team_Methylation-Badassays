ModelBuilding
================
Victor, Ming
April 2, 2017

-   [S0 Set up workspace](#s0-set-up-workspace)
    -   [Load packages](#load-packages)
    -   [Parallel processing](#parallel-processing)
    -   [Load data](#load-data)
-   [S1 Prefiltering Features](#s1-prefiltering-features)
    -   [1.1 Remove NAs](#remove-nas)
    -   [1.2 Reduce CpGs to match test and train](#reduce-cpgs-to-match-test-and-train)
    -   [1.3 Prefiltering cpgs (most variable)](#prefiltering-cpgs-most-variable)
    -   [2.2 SVM with linear kernel](#svm-with-linear-kernel)
-   [S3 Predict Ethnicity for external data Set](#s3-predict-ethnicity-for-external-data-set)
    -   [3.1 glmnet](#glmnet)
    -   [3.2 SVM](#svm)
-   [S4 Analysis of Predictors](#s4-analysis-of-predictors)
    -   [4.2 Plot CpG Predictors](#plot-cpg-predictors)
-   [S5 Tune alpha and lambda 10 x 10 grid](#s5-tune-alpha-and-lambda-10-x-10-grid)
-   [S6 Re-do Modeling and Prediction by Homogenizing Training and Test Sets](#s6-re-do-modeling-and-prediction-by-homogenizing-training-and-test-sets)
    -   [6.1 Merge Test Set with Training Set](#merge-test-set-with-training-set)
    -   [6.2 Unsupervised Clustering PCA](#unsupervised-clustering-pca)
    -   [6.3 Re-do Elastic Net Logistic Regression](#re-do-elastic-net-logistic-regression)
-   [step 6: Trying weighted cases and up-sampling](#step-6-trying-weighted-cases-and-up-sampling)
    -   [Unequal Class Weights](#unequal-class-weights)
    -   [up-sampling](#up-sampling)

S0 Set up workspace
===================

Load packages
-------------

``` r
#source("https://bioconductor.org/biocLite.R")
#biocLite('e1071')                                    # required for glmnet in caret
#biocLite('pROC')
library(pROC)
```

    ## Warning: package 'pROC' was built under R version 3.3.3

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 3.3.3

``` r
library(limma)
```

    ## Warning: package 'limma' was built under R version 3.3.3

``` r
library(caret)
```

    ## Warning: package 'caret' was built under R version 3.3.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 3.3.3

``` r
library(parallel)
library(doParallel)
```

    ## Warning: package 'doParallel' was built under R version 3.3.3

    ## Warning: package 'foreach' was built under R version 3.3.3

    ## Warning: package 'iterators' was built under R version 3.3.3

``` r
library(readxl)
```

    ## Warning: package 'readxl' was built under R version 3.3.3

Parallel processing
-------------------

Set up parallel processing to speed up trian() Make sure to specify in trainControl()

``` r
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
```

Load data
---------

Read in pre-processed data: \*Make sure the pre-processed data (data.txt, which is in data.zip) is present in the ../processed\_data/ directory.

``` r
# load data (pre-processed training set)
train.data <- read.table('../../data/Processed Data/data.txt')
str(train.data)                                   ## row names are CpG sites, column names are sample names
```

    ## 'data.frame':    464923 obs. of  45 variables:
    ##  $ PM104: num  0.9219 0.569 0.0435 0.0635 0.0368 ...
    ##  $ PM112: num  0.8546 0.6886 0.0724 0.0856 0.0364 ...
    ##  $ PM114: num  0.8319 0.6578 0.0906 0.0909 0.0393 ...
    ##  $ PM115: num  0.9021 0.6146 0.0665 0.0839 0.0389 ...
    ##  $ PM119: num  0.896 0.6565 0.0591 0.0987 0.0336 ...
    ##  $ PM120: num  0.8658 0.6338 0.0982 0.1254 0.0421 ...
    ##  $ PM123: num  0.8769 0.6214 0.0571 0.0582 0.0477 ...
    ##  $ PM124: num  0.8481 0.6759 0.0711 0.1022 0.0437 ...
    ##  $ PM130: num  0.821 0.6326 0.1116 0.1038 0.0449 ...
    ##  $ PM136: num  0.9024 0.569 0.0674 0.0671 0.0493 ...
    ##  $ PM139: num  0.8259 0.6098 0.0556 0.0628 0.0371 ...
    ##  $ PM142: num  0.8099 0.6414 0.1094 0.1451 0.0435 ...
    ##  $ PM153: num  0.8911 0.6525 0.0939 0.0888 0.044 ...
    ##  $ PM155: num  0.8384 0.5506 0.0802 0.0984 0.0458 ...
    ##  $ PM158: num  0.835 0.6604 0.1136 0.1122 0.0393 ...
    ##  $ PM167: num  0.8618 0.5825 0.0807 0.0774 0.0424 ...
    ##  $ PM181: num  0.9024 0.6273 0.0885 0.1101 0.0416 ...
    ##  $ PM20 : num  0.8925 0.677 0.1326 0.1507 0.0538 ...
    ##  $ PM205: num  0.851 0.641 0.139 0.142 0.125 ...
    ##  $ PM226: num  0.8822 0.693 0.0972 0.1074 0.0367 ...
    ##  $ PM227: num  0.9049 0.6118 0.0462 0.0719 0.0343 ...
    ##  $ PM233: num  0.8022 0.6165 0.0744 0.1352 0.0448 ...
    ##  $ PM243: num  0.8551 0.6229 0.0896 0.1001 0.0376 ...
    ##  $ PM249: num  0.9153 0.6164 0.0654 0.0946 0.0345 ...
    ##  $ PM29 : num  0.9114 0.6393 0.0399 0.0364 0.0296 ...
    ##  $ PM30 : num  0.8835 0.5735 0.0885 0.1435 0.0407 ...
    ##  $ PM4  : num  0.8463 0.5852 0.0733 0.0869 0.0408 ...
    ##  $ PM40 : num  0.8385 0.6786 0.0586 0.0718 0.0364 ...
    ##  $ PM41 : num  0.8149 0.6576 0.0663 0.1045 0.0485 ...
    ##  $ PM44 : num  0.8089 0.5596 0.1227 0.1471 0.0378 ...
    ##  $ PM46 : num  0.9026 0.6467 0.0751 0.1048 0.0372 ...
    ##  $ PM47 : num  0.8491 0.6345 0.0653 0.0967 0.0468 ...
    ##  $ PM52 : num  0.8891 0.5681 0.0782 0.1018 0.0362 ...
    ##  $ PM53 : num  0.8566 0.6803 0.0866 0.1076 0.0466 ...
    ##  $ PM54 : num  0.8513 0.7054 0.0998 0.1685 0.0427 ...
    ##  $ PM55 : num  0.8317 0.5979 0.0954 0.1153 0.0354 ...
    ##  $ PM58 : num  0.8322 0.6706 0.0754 0.1237 0.0492 ...
    ##  $ PM66 : num  0.9138 0.6296 0.0859 0.1162 0.041 ...
    ##  $ PM71 : num  0.7972 0.5837 0.1344 0.1746 0.0477 ...
    ##  $ PM72 : num  0.8706 0.6164 0.1114 0.1194 0.0396 ...
    ##  $ PM74 : num  0.8611 0.5987 0.0988 0.0997 0.037 ...
    ##  $ PM76 : num  0.8401 0.5985 0.0728 0.1062 0.0409 ...
    ##  $ PM84 : num  0.8955 0.6762 0.2148 0.1548 0.0433 ...
    ##  $ PM9  : num  0.8784 0.6027 0.1118 0.1305 0.0441 ...
    ##  $ PM98 : num  0.8279 0.629 0.0827 0.0837 0.0418 ...

``` r
# transpose our data to have rows as observations, which is more convenient later on for building models
train.data <- as.data.frame(t(train.data))

# load metadata
train.design <- read.csv("../../data/Processed Data/des.txt", sep="\t", header=TRUE)
str(train.design)
```

    ## 'data.frame':    45 obs. of  5 variables:
    ##  $ Samplename  : Factor w/ 45 levels "PM104","PM112",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Sample_Group: Factor w/ 3 levels "CONTROL","IUGR",..: 1 1 1 3 3 1 2 1 2 1 ...
    ##  $ ga          : num  40.7 38.9 38.6 41.1 37.1 38 35.7 40 36.9 38.6 ...
    ##  $ sex         : Factor w/ 2 levels "F","M": 2 1 2 2 2 1 1 2 2 1 ...
    ##  $ Ethnicity   : Factor w/ 2 levels "Asian","Caucasian": 2 1 2 2 1 2 2 2 2 2 ...

``` r
row.names(train.data) == train.design$Samplename               # check that the samples are in same order
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [29] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [43] TRUE TRUE TRUE

Read in test data:

``` r
# read pre-processed test data
test.data <- read.table("../../Data/Processed Data/Test data/Matrix.processed.betas.placenta.txt", row.names = 1, header = T)
test.data <- as.data.frame(t(test.data))   #transpose data

# meta data for test data
test.design <-  read_excel("../../data/Processed Data/Test data/metadata.GA_illumina_methylation.xls", 
    sheet = "Metadata", skip = 28)

# subset only columns we need and rename them
test.design <- test.design[test.design$`Sample name` %in% rownames(test.data),]
test.design <- test.design[,c(1,7,8,10)]
colnames(test.design)[1] <- "Samplename"
colnames(test.design)[3] <- "sex"
colnames(test.design)[4] <- "ga"

str(test.design)
```

    ## Classes 'tbl_df', 'tbl' and 'data.frame':    52 obs. of  4 variables:
    ##  $ Samplename                 : chr  "FT21_v" "FT36_v" "FT26_v" "NTD8_v" ...
    ##  $ characteristics: NTD status: chr  "spina bifida" "control" "anencephaly" "spina bifida" ...
    ##  $ sex                        : chr  "M" "F" "M" "M" ...
    ##  $ ga                         : num  21.1 15 19.3 20 16.7 23.7 17 21.8 22 17 ...

S1 Prefiltering Features
========================

Reducing the number of features can significantly reduce computational time, which is desirable when the dataset is large. However, we must be careful not remove potentially 'interesting' features that have a high chance of being useful in building a classifier.

1.1 Remove NAs
--------------

We should remove any sites with NAs in them or else predictions cannot be generated for these samples if the CpG site is chosen as a predictor.

``` r
#remove sites with NAs
sum(is.na(test.data)) # 52000 total entries that are NA
```

    ## [1] 52000

``` r
test.rmna <- test.data[, colSums(is.na(test.data)) == 0]  # remove columns with NAs present
```

1.2 Reduce CpGs to match test and train
---------------------------------------

Some CpGs were removed in the Test dataset from preprocessing and QC. These need to be removed, or errors might occur when trying to predict.

``` r
# this isn't necessary if the test data didn't have CpGs removed (as a result of QC/preprocessing)
train.data <- train.data[,colnames(train.data) %in% colnames(test.rmna)]
```

1.3 Prefiltering cpgs (most variable)
-------------------------------------

The goal of this prefiltering section is to reduce computational time without compromising detecting interesting features.

``` r
train.sd <- apply(as.matrix(train.data), MARGIN = 2,FUN = sd) #caculate SD for each feature
hist(train.sd)                    # histogram of the s.d.'s
abline(v = mean(train.sd)) 
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-1.png)

``` r
# filter CpG sites with low s.d: only keep those with s.d higher than the average s.d across all CpG sites
train.gsd <- subset(train.sd, train.sd > 0.10)
hist(train.gsd)
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/prefiltering%20based%20on%20SD-2.png)

``` r
# subset training data to only highly variable features
train.data.gsd <- train.data[,colnames(train.data) %in% names(train.gsd)]
```

We reduced the \# of features to 10775 to reduce computation time. `train.data.gsd` is the working dataset. \# S2 Supervised classification: We decided to try two different models for building our classifer: elastic net logistic regression (`glmnet`) and support vector machines (SVM). Both of these models have been used in the literature to build predictive models based on 450k DNA methylation data (Horvath 2013, De Carli et al 2017), indicating that they may be well-suited for our dataset. \#\# 2.1 logistic regression with elastic net regularization

``` r
#renamed training data for standard coding
x.train <- train.data.gsd
y.train <- train.design$Ethnicity 
```

``` r
k = 5
M = 3

fitControl <- trainControl(method = "repeatedcv", 
                                                     number = k,                 # Number of folds
                                                     repeats = M,
                                                     ## Estimate class probabilities
                                                     classProbs = TRUE,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = twoClassSummary,
                                                     allowParallel = TRUE # allow parallel processing
                                                     )  

netGrid <- expand.grid(alpha = c(0.75),
                           lambda = c(0.077, 0.25))
```

We specify the model to be built using repeated cross validation with a fold = 5, and repeats = 3. We tune the model holding alpha constant (alpha = 0.75), keeping alpha high to favour L1 norm to achieve a small panel of biomarkers. Lambda, the magnitude of the penalty, is tested at 0.077, and 0.25.

``` r
set.seed(2017)                                         # training models requires the use of random #s. Setting (set.seed()) the randomness ensures reproducibility

system.time(netFit <- train(x = x.train,   # samples need to be in rows, features need to be columns
                                y = y.train,                  
                                method = "glmnet",                     # glmnet model
                                trControl = fitControl,                # use fitControl to specify cross validation
                                tuneGrid = netGrid,
                                preProcess = c( "center", "scale"),    # Center and Scale the data
                                metric = 'Accuracy')                        # ROC because distribution is slightly skewed
)
```

    ## Loading required package: glmnet

    ## Warning: package 'glmnet' was built under R version 3.3.3

    ## Loading required package: Matrix

    ## Loaded glmnet 2.0-5

    ## 
    ## Attaching package: 'glmnet'

    ## The following object is masked from 'package:pROC':
    ## 
    ##     auc

    ## Warning in train.default(x = x.train, y = y.train, method = "glmnet",
    ## trControl = fitControl, : The metric "Accuracy" was not in the result set.
    ## ROC will be used instead.

    ##    user  system elapsed 
    ##   12.95    0.10   94.92

``` r
netFit
```

    ## glmnet 
    ## 
    ##    45 samples
    ## 10775 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10775), scaled (10775) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 37, 36, 35, 36, 36, 36, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   lambda  ROC        Sens       Spec
    ##   0.077   0.9809524  0.6333333  1   
    ##   0.250   0.9809524  0.4333333  1   
    ## 
    ## Tuning parameter 'alpha' was held constant at a value of 0.75
    ## ROC was used to select the optimal model using  the largest value.
    ## The final values used for the model were alpha = 0.75 and lambda = 0.25.

``` r
netFit$results
```

    ##   alpha lambda       ROC      Sens Spec     ROCSD    SensSD SpecSD
    ## 1  0.75  0.077 0.9809524 0.6333333    1 0.0424012 0.2831232      0
    ## 2  0.75  0.250 0.9809524 0.4333333    1 0.0424012 0.3073181      0

``` r
#saveRDS(netFit, './Data/Processed Data/netFitfinal.rds')
```

``` r
predictorsNet <- predictors(netFit)
length(predictorsNet)
```

    ## [1] 11

``` r
#write.table(predictorsNet, './Data/Processed Data/predictorsGlmnet.txt')
```

Our glmnet-built model has chosen 11 CpGs that can be used to predict ethnicity.

2.2 SVM with linear kernel
--------------------------

This section is for building the model using SVM with a linear kernel (i.e. penalty parameter C = 1). However, because computational time is long, this section is excluded when ran, since we have chosen the glmnet model to be our final model.

``` r
svmControl <- trainControl(method="repeatedcv",   
                           number = 5,
                           repeats=3,           
                           summaryFunction=twoClassSummary, # Use AUC to pick the best model
                           classProbs=TRUE,
                           allowParallel = TRUE)

system.time(svmFit <- train(x=x.train,
                            y= y.train,
                            method = "svmLinear",
                            preProc = c("center","scale"),
                            metric="ROC",
                            trControl= svmControl)  )
```

    ## Loading required package: kernlab

    ## 
    ## Attaching package: 'kernlab'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha

    ##    user  system elapsed 
    ##   14.83    0.08   99.48

``` r
svmFit
```

    ## Support Vector Machines with Linear Kernel 
    ## 
    ##    45 samples
    ## 10775 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10775), scaled (10775) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 36, 36, 36, 36, 36, 35, ... 
    ## Resampling results:
    ## 
    ##   ROC  Sens       Spec     
    ##   1    0.9111111  0.9904762
    ## 
    ## Tuning parameter 'C' was held constant at a value of 1
    ## 

S3 Predict Ethnicity for external data Set
==========================================

Next, we use the models we built and run it on an external data set, where there is no ethnicity information.

3.1 glmnet
----------

Using the `predict()` function we can obtain both binary class prediction results as well as probabilities of being "Asian" under logistic regression model.

``` r
#subset x.test down to the sites used for training (after prefilter)
x.test <- test.data[,colnames(test.data) %in% names(x.train)]

#class predictions
y.predictNet <- predict(netFit,  x.test)
y.predictNet
```

    ##  [1] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ##  [8] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [15] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [22] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [29] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [36] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [43] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [50] Caucasian Caucasian Caucasian
    ## Levels: Asian Caucasian

``` r
#class probabilities
y.predictNetProb <- predict(netFit, x.test, type = 'prob')
y.predictNetProb
```

    ##        Asian Caucasian
    ## 1  0.3251351 0.6748649
    ## 2  0.3578469 0.6421531
    ## 3  0.2633016 0.7366984
    ## 4  0.1780938 0.8219062
    ## 5  0.2850498 0.7149502
    ## 6  0.3327671 0.6672329
    ## 7  0.3038513 0.6961487
    ## 8  0.2379285 0.7620715
    ## 9  0.3429558 0.6570442
    ## 10 0.2889463 0.7110537
    ## 11 0.2718087 0.7281913
    ## 12 0.3105621 0.6894379
    ## 13 0.2860582 0.7139418
    ## 14 0.2853955 0.7146045
    ## 15 0.2793128 0.7206872
    ## 16 0.1953898 0.8046102
    ## 17 0.1998339 0.8001661
    ## 18 0.2355024 0.7644976
    ## 19 0.1949862 0.8050138
    ## 20 0.3643206 0.6356794
    ## 21 0.3421961 0.6578039
    ## 22 0.3202268 0.6797732
    ## 23 0.3054885 0.6945115
    ## 24 0.3763022 0.6236978
    ## 25 0.3917575 0.6082425
    ## 26 0.3064885 0.6935115
    ## 27 0.3628590 0.6371410
    ## 28 0.3245240 0.6754760
    ## 29 0.3675905 0.6324095
    ## 30 0.2889059 0.7110941
    ## 31 0.2168600 0.7831400
    ## 32 0.4541071 0.5458929
    ## 33 0.1972951 0.8027049
    ## 34 0.2246739 0.7753261
    ## 35 0.2807834 0.7192166
    ## 36 0.3326410 0.6673590
    ## 37 0.2260130 0.7739870
    ## 38 0.3023590 0.6976410
    ## 39 0.2964268 0.7035732
    ## 40 0.2992843 0.7007157
    ## 41 0.4258332 0.5741668
    ## 42 0.2806238 0.7193762
    ## 43 0.3687338 0.6312662
    ## 44 0.2317289 0.7682711
    ## 45 0.2065376 0.7934624
    ## 46 0.3363971 0.6636029
    ## 47 0.1820875 0.8179125
    ## 48 0.2686947 0.7313053
    ## 49 0.3029359 0.6970641
    ## 50 0.2364834 0.7635166
    ## 51 0.3123047 0.6876953
    ## 52 0.3035512 0.6964488

``` r
#saveRDS(y.predictNet, './data/Processed Data/y_predictNet.rds')
```

It looks like our model classifies the entire external dataset is Caucasian. This is suspicious, as we believe the samples to come from a relatively heterogenous population. However, due to time constraints, we decided to move ahead and perform downstream analysis. If there was more time, we might think about where we can change our modelling process to produce more sensible results.

#### Some explanations for this result:

-   It's possible that the data set is truly all Caucasian.

-   There are just too many predictors for elastic net to select, some features can be significantly different by ethnicity in training data by chance but were still selected as predictors, which introduces overfitting;

-   The test dataset is too 'different' to have the classifier ran on. (mostly differentially methylated between test and training set across all sites)

-   The self-reported ethnicities in the training data is unreliable.

3.2 SVM
-------

``` r
y.predictSVM <- predict(svmFit,  x.test)

y.predictSVM
```

    ##  [1] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ##  [8] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [15] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [22] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [29] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [36] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [43] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [50] Caucasian Caucasian Caucasian
    ## Levels: Asian Caucasian

S4 Analysis of Predictors
=========================

Here we pull out the CpG sites and look at them more closely. First we will see if clustering with only the predictors separates asians and caucasians \#\# 4.1 Clustering

``` r
library(ggdendro)
```

    ## Warning: package 'ggdendro' was built under R version 3.3.3

``` r
library(sparcl) # ColorDendrogram
library(dendextend)
```

    ## Warning: package 'dendextend' was built under R version 3.3.3

``` r
#without all CpGs used to train
hclust <- hclust(dist(x.train, method = 'euclidean'))

#swap labels with ethnicity
swaplabels <- function(hclust, des){     # des is a design matrix containing 'Samplename' and 'Ethnicity' col
  labels <- data.frame(labels(hclust))   # pulls out current labels (samplename)
  colnames(labels) <- 'Samplename'
  labels <- labels %>% left_join(select(des, Samplename, Ethnicity), by = 'Samplename')
  labels(hclust) <- as.character(labels$Ethnicity)
  return(hclust)
}

hclust <- swaplabels(hclust, train.design)
y1 = cutree(hclust, 3)
ColorDendrogram(hclust, y = y1, labels = names(y1), branchlength = 1.0, main = 'Clustering train on all CpGs')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-1.png)

``` r
#with predictors only
x.train.predictors <- x.train[,colnames(x.train) %in% predictorsNet]
hclust2 <- hclust(dist(x.train.predictors, method = 'euclidean'))
hclust2 <- swaplabels(hclust2, train.design)          #swap labels with ethnicity
y2 = cutree(hclust2, 2)
ColorDendrogram(hclust2, y = y2, labels = names(y2), branchlength = 0.3, main = 'Clustering train with predictors only')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-2.png) We see that clustering with the predictors extracted from our classifier, our training data clusters into two homogenous groups consisting of Asians and Caucasians. This might indicate overfitting, as there were 0 missclassifications.

``` r
# Hierarchical clustering of predicted data, distance measured by Euclidean distance, average linkage
hclust3 <- hclust(dist(x.test, method = 'euclidean'))
y3 = cutree(hclust3, 2)
ColorDendrogram(hclust3, y=y3, labels = names(y3), branchlength = 2, main = 'Clustering Test with all CpGs')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/Clustering%20test%20data-1.png)

``` r
# clustering only with the predictors
x.test.predictors <- x.test[,colnames(x.test) %in% predictorsNet]
hclust4 <- hclust(dist(x.test.predictors, method = 'euclidean'))
y4 = cutree(hclust4, 2)
ColorDendrogram(hclust4, y=y4, labels = names(y4), branchlength = 0.25, main = 'Clustering Test with predictors only')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/Clustering%20test%20data-2.png) We can see that clustering the external data set (test) does not improve the separation of the test data into two main clusters, indicating that the classifier is not producing a heterogenous set of predictions.

``` r
# add samplename column to match on
x.test.predictors <- x.test.predictors %>% 
                        tibble::rownames_to_column('Samplename')
x.train.predictors <- x.train.predictors %>%
                        tibble::rownames_to_column('Samplename') 

# replace sample name with true ethnicity info
#x.train.predictors <- x.train.predictors %>% left_join(select(train.design, Samplename, Ethnicity), by = 'Samplename') 
#x.train.predictors$Samplename <- x.train.predictors$Ethnicity 

# combine train and test
x.test.train.predictors <- full_join(x.train.predictors, x.test.predictors) %>%
                            tibble::column_to_rownames('Samplename')
```

    ## Joining, by = c("Samplename", "cg24673385", "cg15486123", "cg22853943", "cg13921903", "cg08704934", "cg14581129", "cg05393297", "cg25025879", "cg16329197", "cg12011926", "cg12602405")

``` r
# clustering
hclust5 <- hclust(dist(x.test.train.predictors, method = 'euclidean'))
labels5 <- data.frame(labels(hclust5))   # pulls out current labels (samplename)
colnames(labels5) <- 'Samplename'
labels5 <- labels5 %>% left_join(select(train.design, Samplename, Ethnicity), by = 'Samplename')
#replace train samples with ethnicity labels
#labels5$Samplename[!is.na(labels5$Ethnicity)] <- as.character(labels5$Ethnicity[!is.na(labels5$Ethnicity)])

labels(hclust5) <- labels5$Samplename


hclust5 <- swaplabels(hclust5, train.design)
labels(hclust5)
```

    ##  [1] "Caucasian" "Caucasian" "Caucasian" "Caucasian" "Caucasian"
    ##  [6] "Caucasian" "Caucasian" "Caucasian" "Caucasian" "Caucasian"
    ## [11] "Caucasian" "Caucasian" "Caucasian" "Caucasian" "Caucasian"
    ## [16] "Caucasian" "Caucasian" "Caucasian" "Caucasian" "Caucasian"
    ## [21] "Caucasian" "Caucasian" "Caucasian" NA          NA         
    ## [26] "Caucasian" NA          NA          NA          NA         
    ## [31] NA          NA          NA          "Caucasian" NA         
    ## [36] "Caucasian" "Caucasian" "Caucasian" "Caucasian" "Caucasian"
    ## [41] NA          NA          "Caucasian" "Caucasian" NA         
    ## [46] "Asian"     NA          "Caucasian" NA          NA         
    ## [51] NA          NA          NA          NA          NA         
    ## [56] NA          NA          NA          NA          NA         
    ## [61] NA          NA          NA          NA          NA         
    ## [66] NA          NA          NA          NA          NA         
    ## [71] NA          "Asian"     "Asian"     NA          NA         
    ## [76] NA          NA          NA          NA          NA         
    ## [81] NA          NA          NA          NA          "Asian"    
    ## [86] "Asian"     "Asian"     "Asian"     "Asian"     "Asian"    
    ## [91] "Asian"     "Asian"     "Asian"     NA          NA         
    ## [96] NA          NA

``` r
y5 = cutree(hclust5, 5)
ColorDendrogram(hclust5, y = y5, labels = names(y5), branchlength = 0.3, main = 'Clustering train with predictors only')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/cluster%20both%20test%20and%20train-1.png) When we perform hierarchical clustering with the entire train and test set, we can see that Caucasians and Asians mainly separate into the two largest clusters, with the majority of the test set (unlabeled branches) clustering closer to the Caucasians samples.

4.2 Plot CpG Predictors
-----------------------

``` r
glmImp <- varImp(netFit, scale = F) # gives the t-statistic for all CpGs in the dataset
plot(glmImp, top = length(predictors(netFit)))
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/plot%20top%2035-1.png) Here we plot the 11 predictor CpGs against 'importance' which is calculated based on their relative t-statistic score.

``` r
# For training data set
cpg1 <- x.train.predictors %>% select(Samplename, cg16329197) %>% 
                                left_join(train.design, 'Samplename')
ggplot(cpg1, aes(x=Ethnicity, y=cg16329197))+
  geom_boxplot()+
  ggtitle('Top CpG predictor methylation in Training data is differentially
          methylated')+
  ylab('cg16329197 methylation')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/plotting%20CpGs-1.png)

``` r
# Pick 11th ranked CpG
cpg2 <- x.train.predictors %>% select(Samplename, cg22853943) %>% 
                                left_join(train.design, 'Samplename')
ggplot(cpg2, aes(x=Ethnicity, y=cg22853943))+
  geom_boxplot()+
  ggtitle('11th ranked CpG predictor methylation in Training data is
          differentially methylated')+
  ylab('cg22853943 methylation')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/plotting%20CpGs-2.png) We can see that the 1st and 11th ranked predictor CpG are both obviously differentially methylated in the training dataset between the asians and caucasians. This is a good sign that the model has chosen 'useful' CpG sites. However, perhaps these CpGs fit our training data too well.

S5 Tune alpha and lambda 10 x 10 grid
=====================================

This section is to tune parameters across 100 different combinations of alpha and lambda

``` r
netGrid100 <-  expand.grid(alpha = c(0.10, 0.20, 0.30, 0.40, 0.50, 
                                     0.60, 0.70, 0.80, 0.90, 1.00),
                           lambda = c(0.025, 0.050, 0.075, 0.10, 0.15, 
                                      0.20, 0.25, 0.30, 0.40, 0.50))

set.seed(2017)                              

system.time(netFit100 <- train(x = x.train, 
                                y = y.train,                  
                                method = "glmnet",                     # glmnet model
                                trControl = fitControl,                # use fitControl to specify cross validation
                                tuneGrid = netGrid100,
                                preProcess = c( "center", "scale"),    # Center and Scale the data
                                metric = 'ROC')                        # ROC because distribution is slightly skewed
)
```

    ##    user  system elapsed 
    ##   13.14    0.06  820.76

``` r
netFit100
```

    ## glmnet 
    ## 
    ##    45 samples
    ## 10775 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10775), scaled (10775) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 37, 36, 35, 36, 36, 36, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   alpha  lambda  ROC        Sens        Spec     
    ##   0.1    0.025   0.9920635  0.90000000  1.0000000
    ##   0.1    0.050   0.9920635  0.87777778  1.0000000
    ##   0.1    0.075   0.9920635  0.87777778  1.0000000
    ##   0.1    0.100   0.9920635  0.85555556  1.0000000
    ##   0.1    0.150   0.9920635  0.82222222  1.0000000
    ##   0.1    0.200   0.9920635  0.82222222  1.0000000
    ##   0.1    0.250   0.9920635  0.82222222  1.0000000
    ##   0.1    0.300   0.9968254  0.82222222  1.0000000
    ##   0.1    0.400   0.9968254  0.82222222  1.0000000
    ##   0.1    0.500   0.9968254  0.82222222  1.0000000
    ##   0.2    0.025   0.9920635  0.87777778  1.0000000
    ##   0.2    0.050   0.9920635  0.85555556  1.0000000
    ##   0.2    0.075   0.9920635  0.82222222  1.0000000
    ##   0.2    0.100   0.9920635  0.82222222  1.0000000
    ##   0.2    0.150   0.9920635  0.82222222  1.0000000
    ##   0.2    0.200   0.9920635  0.82222222  1.0000000
    ##   0.2    0.250   0.9920635  0.80000000  1.0000000
    ##   0.2    0.300   0.9920635  0.80000000  1.0000000
    ##   0.2    0.400   0.9920635  0.77777778  1.0000000
    ##   0.2    0.500   0.9920635  0.65555556  1.0000000
    ##   0.3    0.025   0.9920635  0.82222222  1.0000000
    ##   0.3    0.050   0.9920635  0.82222222  1.0000000
    ##   0.3    0.075   0.9920635  0.82222222  1.0000000
    ##   0.3    0.100   0.9920635  0.80000000  1.0000000
    ##   0.3    0.150   0.9888889  0.80000000  1.0000000
    ##   0.3    0.200   0.9888889  0.77777778  1.0000000
    ##   0.3    0.250   0.9888889  0.74444444  1.0000000
    ##   0.3    0.300   0.9888889  0.65555556  1.0000000
    ##   0.3    0.400   0.9936508  0.60000000  1.0000000
    ##   0.3    0.500   0.9936508  0.51111111  1.0000000
    ##   0.4    0.025   0.9888889  0.82222222  1.0000000
    ##   0.4    0.050   0.9888889  0.80000000  1.0000000
    ##   0.4    0.075   0.9888889  0.80000000  1.0000000
    ##   0.4    0.100   0.9888889  0.80000000  1.0000000
    ##   0.4    0.150   0.9888889  0.74444444  1.0000000
    ##   0.4    0.200   0.9888889  0.68888889  1.0000000
    ##   0.4    0.250   0.9936508  0.65555556  1.0000000
    ##   0.4    0.300   0.9936508  0.63333333  1.0000000
    ##   0.4    0.400   0.9904762  0.56666667  1.0000000
    ##   0.4    0.500   0.9904762  0.32222222  1.0000000
    ##   0.5    0.025   0.9888889  0.80000000  0.9904762
    ##   0.5    0.050   0.9888889  0.80000000  1.0000000
    ##   0.5    0.075   0.9888889  0.74444444  1.0000000
    ##   0.5    0.100   0.9888889  0.74444444  1.0000000
    ##   0.5    0.150   0.9888889  0.68888889  1.0000000
    ##   0.5    0.200   0.9904762  0.65555556  1.0000000
    ##   0.5    0.250   0.9904762  0.60000000  1.0000000
    ##   0.5    0.300   0.9857143  0.56666667  1.0000000
    ##   0.5    0.400   0.9857143  0.32222222  1.0000000
    ##   0.5    0.500   0.9809524  0.04444444  1.0000000
    ##   0.6    0.025   0.9888889  0.80000000  0.9904762
    ##   0.6    0.050   0.9888889  0.74444444  1.0000000
    ##   0.6    0.075   0.9888889  0.71111111  1.0000000
    ##   0.6    0.100   0.9857143  0.71111111  1.0000000
    ##   0.6    0.150   0.9857143  0.63333333  1.0000000
    ##   0.6    0.200   0.9857143  0.60000000  1.0000000
    ##   0.6    0.250   0.9857143  0.56666667  1.0000000
    ##   0.6    0.300   0.9857143  0.45555556  1.0000000
    ##   0.6    0.400   0.9809524  0.13333333  1.0000000
    ##   0.6    0.500   0.9761905  0.00000000  1.0000000
    ##   0.7    0.025   0.9772487  0.71111111  0.9904762
    ##   0.7    0.050   0.9777778  0.71111111  1.0000000
    ##   0.7    0.075   0.9809524  0.68888889  1.0000000
    ##   0.7    0.100   0.9809524  0.63333333  1.0000000
    ##   0.7    0.150   0.9809524  0.60000000  1.0000000
    ##   0.7    0.200   0.9809524  0.56666667  1.0000000
    ##   0.7    0.250   0.9809524  0.47777778  1.0000000
    ##   0.7    0.300   0.9761905  0.41111111  1.0000000
    ##   0.7    0.400   0.9761905  0.00000000  1.0000000
    ##   0.7    0.500   0.9650794  0.00000000  1.0000000
    ##   0.8    0.025   0.9777778  0.68888889  0.9904762
    ##   0.8    0.050   0.9809524  0.68888889  1.0000000
    ##   0.8    0.075   0.9809524  0.63333333  1.0000000
    ##   0.8    0.100   0.9772487  0.60000000  1.0000000
    ##   0.8    0.150   0.9809524  0.56666667  1.0000000
    ##   0.8    0.200   0.9761905  0.51111111  1.0000000
    ##   0.8    0.250   0.9698413  0.43333333  1.0000000
    ##   0.8    0.300   0.9746032  0.22222222  1.0000000
    ##   0.8    0.400   0.9650794  0.00000000  1.0000000
    ##   0.8    0.500   0.5000000  0.00000000  1.0000000
    ##   0.9    0.025   0.9809524  0.63333333  0.9904762
    ##   0.9    0.050   0.9724868  0.63333333  1.0000000
    ##   0.9    0.075   0.9693122  0.60000000  1.0000000
    ##   0.9    0.100   0.9693122  0.60000000  1.0000000
    ##   0.9    0.150   0.9637566  0.51111111  1.0000000
    ##   0.9    0.200   0.9682540  0.48888889  1.0000000
    ##   0.9    0.250   0.9626984  0.41111111  1.0000000
    ##   0.9    0.300   0.9682540  0.07777778  1.0000000
    ##   0.9    0.400   0.9449735  0.00000000  1.0000000
    ##   0.9    0.500   0.5000000  0.00000000  1.0000000
    ##   1.0    0.025   0.9624339  0.63333333  1.0000000
    ##   1.0    0.050   0.9656085  0.60000000  1.0000000
    ##   1.0    0.075   0.9619048  0.62222222  1.0000000
    ##   1.0    0.100   0.9544974  0.58888889  1.0000000
    ##   1.0    0.150   0.9558201  0.56666667  1.0000000
    ##   1.0    0.200   0.9650794  0.51111111  1.0000000
    ##   1.0    0.250   0.9756614  0.41111111  1.0000000
    ##   1.0    0.300   0.9656085  0.00000000  1.0000000
    ##   1.0    0.400   0.5000000  0.00000000  1.0000000
    ##   1.0    0.500   0.5000000  0.00000000  1.0000000
    ## 
    ## ROC was used to select the optimal model using  the largest value.
    ## The final values used for the model were alpha = 0.1 and lambda = 0.5.

``` r
trellis.par.set(caretTheme())
ggplot(netFit100)
```

    ## Warning: Ignoring unknown aesthetics: shape

    ## Warning: The shape palette can deal with a maximum of 6 discrete values
    ## because more than 6 becomes difficult to discriminate; you have
    ## 10. Consider specifying shapes manually if you must have them.

    ## Warning: Removed 40 rows containing missing values (geom_point).

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/examine%20CV-1.png)

``` r
#heatmap of results
plot(netFit100, metric = "ROC", plotType = "level",
     scales = list(x = list(rot = 90)))
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/examine%20CV-2.png)

``` r
glmImp100 <- varImp(netFit100, scale = F) # gives the t-statistic for all CpGs in the dataset
plot(glmImp100, top = length(predictors(netFit100)))
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/examine%20CV-3.png)

``` r
length(predictors(netFit100))
```

    ## [1] 200

S6 Re-do Modeling and Prediction by Homogenizing Training and Test Sets
=======================================================================

Following the discussions in step 3, we aim to re-do our model fitting and prediction by first assessing how different is the test set from our training data using PCA. Then we attempt to homogenize the two datasets by removing the top PC. After homogenization, we re-fit elastic net logistic regression with the transformed training set and re-predict the transformed test set.

6.1 Merge Test Set with Training Set
------------------------------------

Merge the two datasets, only CpGs present in both sets are kept.

``` r
merged.all <- rbind(train.data, test.data[, colnames(test.data) %in% colnames(train.data)])
merged.design <- rbind(train.design[,c("Samplename","ga","sex")], test.design[,c("Samplename","ga","sex")])
merged.design$Group = relevel(
    factor(c(rep("Train",nrow(train.data)),rep("Test", nrow(test.data)))), 
    ref = "Train")
```

6.2 Unsupervised Clustering PCA
-------------------------------

First we perform PCA on the merged dataset to spot potential systematic differences between test vs. training set. PCA is performed for the correlation matrix.

``` r
pc.merged <- prcomp(merged.all, center=T, scale = T) # perform PCA on correlation matrix
PC1to5 <- data.frame(pc.merged$x[,1:5])              # Take out first 5 PCs
PC1to5 <- PC1to5 %>% tibble::rownames_to_column('Samplename') %>%             # Put sample names into a column 
                    left_join(merged.design, 'Samplename')                         # Join the metadata info 
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factor and character vector, coercing into character vector

``` r
summary(pc.merged)
```

    ## Importance of components:
    ##                             PC1      PC2       PC3       PC4       PC5
    ## Standard deviation     354.6554 163.0598 137.09167 117.92856 106.27552
    ## Proportion of Variance   0.3151   0.0666   0.04708   0.03484   0.02829
    ## Cumulative Proportion    0.3151   0.3817   0.42876   0.46360   0.49189
    ##                             PC6      PC7     PC8      PC9     PC10
    ## Standard deviation     87.71157 82.54326 70.9159 69.07545 66.53617
    ## Proportion of Variance  0.01927  0.01707  0.0126  0.01195  0.01109
    ## Cumulative Proportion   0.51116  0.52823  0.5408  0.55278  0.56387
    ##                            PC11     PC12     PC13     PC14     PC15
    ## Standard deviation     64.35782 61.32775 59.96083 59.46560 58.37218
    ## Proportion of Variance  0.01038  0.00942  0.00901  0.00886  0.00854
    ## Cumulative Proportion   0.57424  0.58366  0.59267  0.60153  0.61006
    ##                            PC16    PC17     PC18     PC19     PC20
    ## Standard deviation     57.86636 56.1464 55.14601 54.18909 53.26386
    ## Proportion of Variance  0.00839  0.0079  0.00762  0.00736  0.00711
    ## Cumulative Proportion   0.61845  0.6263  0.63397  0.64132  0.64843
    ##                            PC21     PC22     PC23     PC24     PC25
    ## Standard deviation     52.98262 52.15731 51.88641 51.42077 51.00811
    ## Proportion of Variance  0.00703  0.00681  0.00674  0.00662  0.00652
    ## Cumulative Proportion   0.65546  0.66227  0.66902  0.67564  0.68216
    ##                            PC26     PC27     PC28     PC29     PC30
    ## Standard deviation     50.51731 50.25923 49.64607 49.25210 49.12714
    ## Proportion of Variance  0.00639  0.00633  0.00617  0.00608  0.00605
    ## Cumulative Proportion   0.68855  0.69488  0.70105  0.70713  0.71318
    ##                            PC31     PC32     PC33     PC34     PC35
    ## Standard deviation     48.85641 48.48633 48.28366 47.97271 47.73080
    ## Proportion of Variance  0.00598  0.00589  0.00584  0.00576  0.00571
    ## Cumulative Proportion   0.71915  0.72504  0.73088  0.73665  0.74236
    ##                            PC36     PC37    PC38     PC39     PC40
    ## Standard deviation     47.58246 47.02875 46.8765 46.28328 46.11304
    ## Proportion of Variance  0.00567  0.00554  0.0055  0.00537  0.00533
    ## Cumulative Proportion   0.74803  0.75357  0.7591  0.76444  0.76976
    ##                            PC41     PC42     PC43     PC44     PC45
    ## Standard deviation     45.93975 45.80212 45.36270 45.27876 44.99385
    ## Proportion of Variance  0.00529  0.00526  0.00515  0.00514  0.00507
    ## Cumulative Proportion   0.77505  0.78031  0.78546  0.79060  0.79567
    ##                            PC46     PC47     PC48     PC49    PC50
    ## Standard deviation     44.70837 44.32468 44.14109 44.00117 43.7699
    ## Proportion of Variance  0.00501  0.00492  0.00488  0.00485  0.0048
    ## Cumulative Proportion   0.80067  0.80560  0.81048  0.81533  0.8201
    ##                            PC51     PC52     PC53     PC54    PC55
    ## Standard deviation     43.64809 43.40574 43.18623 43.01068 42.8638
    ## Proportion of Variance  0.00477  0.00472  0.00467  0.00463  0.0046
    ## Cumulative Proportion   0.82490  0.82962  0.83429  0.83892  0.8435
    ##                            PC56     PC57     PC58     PC59     PC60
    ## Standard deviation     42.65923 42.31732 42.20239 42.03615 42.00739
    ## Proportion of Variance  0.00456  0.00449  0.00446  0.00443  0.00442
    ## Cumulative Proportion   0.84808  0.85257  0.85703  0.86146  0.86588
    ##                            PC61     PC62     PC63     PC64     PC65
    ## Standard deviation     41.64098 41.60680 41.40311 41.31415 41.05270
    ## Proportion of Variance  0.00434  0.00434  0.00429  0.00428  0.00422
    ## Cumulative Proportion   0.87022  0.87456  0.87885  0.88313  0.88735
    ##                            PC66     PC67     PC68     PC69     PC70
    ## Standard deviation     40.80753 40.71364 40.53636 40.34016 40.08043
    ## Proportion of Variance  0.00417  0.00415  0.00412  0.00408  0.00402
    ## Cumulative Proportion   0.89152  0.89567  0.89979  0.90387  0.90789
    ##                           PC71     PC72     PC73     PC74     PC75
    ## Standard deviation     39.9449 39.87343 39.67109 39.55220 39.51245
    ## Proportion of Variance  0.0040  0.00398  0.00394  0.00392  0.00391
    ## Cumulative Proportion   0.9119  0.91587  0.91981  0.92373  0.92764
    ##                            PC76     PC77     PC78    PC79     PC80
    ## Standard deviation     39.33753 39.11535 38.72265 38.4267 38.33379
    ## Proportion of Variance  0.00388  0.00383  0.00376  0.0037  0.00368
    ## Cumulative Proportion   0.93152  0.93535  0.93911  0.9428  0.94649
    ##                            PC81     PC82     PC83     PC84     PC85
    ## Standard deviation     38.22757 38.08962 37.87122 37.63844 37.60782
    ## Proportion of Variance  0.00366  0.00363  0.00359  0.00355  0.00354
    ## Cumulative Proportion   0.95015  0.95378  0.95737  0.96092  0.96447
    ##                           PC86     PC87     PC88     PC89     PC90
    ## Standard deviation     37.3855 37.26789 37.05874 36.89020 36.63627
    ## Proportion of Variance  0.0035  0.00348  0.00344  0.00341  0.00336
    ## Cumulative Proportion   0.9680  0.97145  0.97489  0.97829  0.98166
    ##                            PC91     PC92     PC93     PC94     PC95
    ## Standard deviation     36.17745 35.95192 35.38269 34.90980 34.56520
    ## Proportion of Variance  0.00328  0.00324  0.00314  0.00305  0.00299
    ## Cumulative Proportion   0.98494  0.98817  0.99131  0.99436  0.99735
    ##                            PC96      PC97
    ## Standard deviation     32.49466 1.395e-12
    ## Proportion of Variance  0.00265 0.000e+00
    ## Cumulative Proportion   1.00000 1.000e+00

``` r
# first PC can explain 31.5% variance

# scree plot
plot(pc.merged, type = "l", main = "PCA Scree Plot for Merged Data")
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/merged%20PCA-1.png)

``` r
# scatter plot matrix for the first 5 PCs
splom(PC1to5[,c(2:6,9)], raster = TRUE)
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/merged%20PCA-2.png)

``` r
# first PC is significantly different between training vs. test
```

The first PC differentiates training and test set, which means there are systematic differences between the two sets. We must discard the top PC first before predicting test set. This is done by:

1.  Reconstruct the (centered and scaled) merged dataset by discarding the top PC;

2.  Re-scale and re-center the correlation matrix for the reconstructed merged data;

``` r
# discard the first PC
Xhat<- pc.merged$x[,-1] %*% t(pc.merged$rotation[,-1]) 

# back-scale features to original center and scale
merged.trunc <- scale(Xhat, center = F, scale = 1/pc.merged$scale)
merged.trunc <- scale(merged.trunc, center = -1 * pc.merged$center, scale = FALSE)
str(merged.trunc)
```

    ##  num [1:97, 1:399206] 0.897 0.836 0.821 0.886 0.878 ...
    ##  - attr(*, "scaled:scale")= Named num [1:399206] 29 41.3 35.4 16.7 12.8 ...
    ##   ..- attr(*, "names")= chr [1:399206] "cg13869341" "cg12045430" "cg20826792" "cg20253340" ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:97] "PM104" "PM112" "PM114" "PM115" ...
    ##   ..$ : chr [1:399206] "cg13869341" "cg12045430" "cg20826792" "cg20253340" ...
    ##  - attr(*, "scaled:center")= Named num [1:399206] -0.848 -0.093 -0.119 -0.359 -0.312 ...
    ##   ..- attr(*, "names")= chr [1:399206] "cg13869341" "cg12045430" "cg20826792" "cg20253340" ...

``` r
# verify that PC truncation was sucessful
pc.trunc <- prcomp(merged.trunc, center=T, scale = T)

PC1to5.trunc <- data.frame(pc.trunc$x[,1:5])              # Take out first 5 PCs
PC1to5.trunc <- PC1to5.trunc %>% tibble::rownames_to_column('Samplename') %>%             # Put sample names into a column
                    left_join(merged.design, 'Samplename')                         # Join the metadata info
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factor and character vector, coercing into character vector

``` r
# scatter plot matrix for the first 5 PCs
splom(PC1to5.trunc[,c(2:6,9)], raster = TRUE)
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/discard%20the%20first%20pc-1.png)

``` r
# no more separation between test and training
```

6.3 Re-do Elastic Net Logistic Regression
-----------------------------------------

After discarding the top PC for merged data, we separate training and test dataset again and use the homogenized training set to re-fit the logistic regression model with elastic net.

``` r
# subset only training samples
x.train.redo <- as.data.frame(merged.trunc[1:45,])

# filter low variances features
train.sd.redo <- apply(as.matrix(x.train.redo), MARGIN = 2,FUN = sd)
train.gsd.redo <- subset(train.sd.redo, train.sd.redo > 0.10)
train.data.gsd.redo <- x.train.redo[,colnames(x.train.redo) %in% names(train.gsd.redo)]

x.train.redo <- train.data.gsd.redo
```

``` r
netGrid <- expand.grid(alpha = (1:9)/10,
                          lambda = seq(0.05,0.5,length.out = 9))
netGrid <- expand.grid(alpha = c(0.75),
                           lambda = c(0.077, 0.25))
```

Re-fit the `glmnet` model:

``` r
set.seed(2017)                                         # training models requires the use of random #s. Setting (set.seed()) the randomness ensures reproducibility


system.time(netFit.redo <- train(x = x.train.redo,   # samples need to be in rows, features need to be columns
                                y = y.train,                  
                                method = "glmnet",                     # glmnet model
                                trControl = fitControl,                # use fitControl to specify cross validation
                                tuneGrid = netGrid,
                                preProcess = c( "center", "scale"),    # Center and Scale the data
                                metric = 'ROC')                        # ROC because distribution is slightly skewed
)
```

    ##    user  system elapsed 
    ##   12.30    0.22   89.53

``` r
netFit.redo
```

    ## glmnet 
    ## 
    ##    45 samples
    ## 10551 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10551), scaled (10551) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 37, 36, 35, 36, 36, 36, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   lambda  ROC        Sens       Spec
    ##   0.077   0.9809524  0.7111111  1   
    ##   0.250   0.9841270  0.4666667  1   
    ## 
    ## Tuning parameter 'alpha' was held constant at a value of 0.75
    ## ROC was used to select the optimal model using  the largest value.
    ## The final values used for the model were alpha = 0.75 and lambda = 0.25.

``` r
length(predictors(netFit.redo))
```

    ## [1] 10

``` r
# Histogram of probability(Asian) for training set
probTest.redo <- predict(netFit.redo, x.train.redo, type = 'prob')
ethProb.redo <- probTest.redo[,'Asian']
hist(ethProb.redo)
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/tune%20glmnet%20parameters%20redo-1.png)

Predict the homogenized test set:

``` r
# subset the test set
x.test.redo <- merged.trunc[46:97,]

x.test.redo <- x.test.redo[, colnames(x.test.redo) %in% colnames(x.train.redo)]

# classification

y.predict.redo <- predict(netFit.redo,  x.test.redo)

# predicted probability to be Asian
y.predict.redo <- predict(netFit.redo,  x.test.redo, type = "prob")
y.predict.redo[,"Asian"]
```

    ##  [1] 0.1493926 0.2254065 0.2156075 0.1963437 0.1767613 0.3002691 0.2228806
    ##  [8] 0.1719623 0.1965631 0.1739325 0.1980987 0.1475837 0.2171292 0.2061519
    ## [15] 0.4719445 0.1769546 0.3045401 0.4689060 0.2103207 0.1860418 0.1329820
    ## [22] 0.1852908 0.1902585 0.5018671 0.2469115 0.1894678 0.3748567 0.1702613
    ## [29] 0.1675490 0.1438982 0.2534604 0.3014062 0.3243579 0.2091487 0.2097832
    ## [36] 0.3609782 0.2081737 0.2138131 0.1930854 0.2625372 0.1621245 0.4572762
    ## [43] 0.3080319 0.1920813 0.3763593 0.2627948 0.2399595 0.1859450 0.2413379
    ## [50] 0.3578074 0.1388270 0.4826315

``` r
# histogramfor the prob. to be Asian
hist(y.predict.redo[,"Asian"], main = "Predicted Probability to be Asian")
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/predict%20redo-1.png)

step 6: Trying weighted cases and up-sampling
=============================================

We realized after the poster session that another reason our `glmnet` model predicts poorly is because we have imbalanced number of classes. We attempt some methods that remedies this problem and update our results:

Unequal Class Weights
---------------------

Many of the predictive models for classification have the ability to use case weights where each individual data point can be given more emphasis in the model training phase. One approach to rebalancing the training set would be to increase the weights for the samples in the minority classes. This can be interpreted as having identical duplicate data points with the exact same predictor values. Logistic regression, for example, can utilize "Asian" class weights in this way.

-   From: [How do I handle an unbalanced training set?](https://www.researchgate.net/post/In_classification_how_do_i_handle_an_unbalanced_training_set)

Here The weights ratio between an Asian and a Caucasian is n\_Caucasian/n\_Asian. Also, we would use the homogenized data.

``` r
# tuning parameter grid
netGrid <- expand.grid(alpha = (5:9)/10,
                          lambda = seq(0.1,0.5,length.out = 5))

# Create model weights (they sum to one)

model_weights <- ifelse(y.train == "Asian",
                        (1/table(y.train)["Asian"]) * 0.5,
                        (1/table(y.train)["Caucasian"]) * 0.5)

# Build weighted model

weighted_fit <- train(x = x.train.redo,
                      y = y.train,
                      method = "glmnet",
                      weights = model_weights,
                      metric = "ROC",
                      trControl = fitControl,
                                    tuneGrid = netGrid,
                                    preProcess = c( "center", "scale"))

weighted_fit
```

    ## glmnet 
    ## 
    ##    45 samples
    ## 10552 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10551), scaled (10551) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 36, 36, 36, 36, 36, 35, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   alpha  lambda  ROC        Sens       Spec     
    ##   0.5    0.1     0.9952381  0.8666667  1.0000000
    ##   0.5    0.2     0.9952381  0.8666667  1.0000000
    ##   0.5    0.3     0.9952381  0.8333333  1.0000000
    ##   0.5    0.4     0.9952381  0.8111111  1.0000000
    ##   0.5    0.5     0.9952381  0.8000000  0.9904762
    ##   0.6    0.1     0.9952381  0.8111111  1.0000000
    ##   0.6    0.2     0.9952381  0.8111111  1.0000000
    ##   0.6    0.3     0.9904762  0.7888889  0.9904762
    ##   0.6    0.4     0.9952381  0.7666667  0.9904762
    ##   0.6    0.5     0.9915344  0.7666667  0.9904762
    ##   0.7    0.1     0.9952381  0.8111111  1.0000000
    ##   0.7    0.2     0.9904762  0.7555556  1.0000000
    ##   0.7    0.3     0.9867725  0.7666667  0.9904762
    ##   0.7    0.4     0.9915344  0.7333333  0.9904762
    ##   0.7    0.5     0.9915344  0.7666667  0.9619048
    ##   0.8    0.1     0.9952381  0.7777778  0.9888889
    ##   0.8    0.2     0.9904762  0.7555556  0.9793651
    ##   0.8    0.3     0.9878307  0.7333333  0.9793651
    ##   0.8    0.4     0.9915344  0.7333333  0.9809524
    ##   0.8    0.5     0.9804233  0.8444444  0.8761905
    ##   0.9    0.1     0.9904762  0.7555556  0.9793651
    ##   0.9    0.2     0.9962963  0.7333333  0.9793651
    ##   0.9    0.3     0.9962963  0.7555556  0.9904762
    ##   0.9    0.4     0.9878307  0.8000000  0.9523810
    ##   0.9    0.5     0.5915344  0.6000000  0.4000000
    ## 
    ## ROC was used to select the optimal model using  the largest value.
    ## The final values used for the model were alpha = 0.9 and lambda = 0.3.

``` r
# prediction, classification results

y.predict.weight <- predict(weighted_fit,  x.test.redo)
y.predict.weight
```

    ##  [1] Caucasian Caucasian Caucasian Caucasian Caucasian Asian     Caucasian
    ##  [8] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [15] Asian     Caucasian Caucasian Asian     Caucasian Caucasian Caucasian
    ## [22] Caucasian Caucasian Asian     Caucasian Caucasian Asian     Caucasian
    ## [29] Caucasian Caucasian Caucasian Caucasian Asian     Caucasian Caucasian
    ## [36] Asian     Caucasian Caucasian Caucasian Caucasian Caucasian Asian    
    ## [43] Asian     Caucasian Asian     Caucasian Caucasian Caucasian Caucasian
    ## [50] Asian     Caucasian Asian    
    ## Levels: Asian Caucasian

``` r
y.predict.weight.des <- data.frame(Samplename = rownames(x.test.redo), Ethnicity = paste("Test",y.predict.weight,sep = "_")) # to be used for dendrogram

# predicted probability to be Asian
y.predict.weight <- predict(weighted_fit,  x.test.redo, type = "prob")
y.predict.weight[,"Asian"]
```

    ##  [1] 0.3393053 0.4034113 0.4285244 0.3969305 0.3728601 0.5068674 0.4035799
    ##  [8] 0.3784845 0.4093938 0.3766687 0.3645585 0.3315011 0.3714280 0.3838460
    ## [15] 0.6155602 0.3667729 0.4818798 0.6081624 0.4180080 0.3405599 0.2964691
    ## [22] 0.3768089 0.3980019 0.5903607 0.4607889 0.3655339 0.5535791 0.3462776
    ## [29] 0.3668982 0.3267600 0.4225762 0.4892697 0.5100308 0.3826325 0.4274563
    ## [36] 0.5013843 0.4179416 0.4179803 0.3883662 0.4244360 0.3540207 0.5804047
    ## [43] 0.5067055 0.3717680 0.5481965 0.4581657 0.4335487 0.3776910 0.4292777
    ## [50] 0.5213675 0.3275139 0.6189497

``` r
# histogram for the prob. to be Asian
hist(y.predict.weight[,"Asian"], main = "Predicted Probability to be Asian")
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/fit%20weighted%20glmnet-1.png)

Now map predicted results to dendrogram from hierarchical clustering of the merged data:

``` r
# make a design matrix containing info on whether it's training or test data, and ethnicity info

weight.des <- rbind(train.design[,c("Samplename","Ethnicity")],y.predict.weight.des)

# clustering for (centered and scaled) merged data
hclust.weight <- hclust(dist(scale(merged.trunc,center = T,scale = T), method = 'euclidean'), method = "average")

labels.weight <- data.frame(labels(hclust.weight))   # pulls out current labels (samplename)
colnames(labels.weight) <- 'Samplename'
labels.weight <- labels.weight %>% left_join(weight.des, by = 'Samplename')
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
labels(hclust.weight) <- labels.weight$Samplename


hclust.weight <- swaplabels(hclust.weight, weight.des)
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
labels(hclust.weight)
```

    ##  [1] "Caucasian"      "Caucasian"      "Caucasian"      "Test_Asian"    
    ##  [5] "Test_Asian"     "Asian"          "Test_Caucasian" "Test_Asian"    
    ##  [9] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Asian"    
    ## [13] "Test_Asian"     "Test_Caucasian" "Test_Caucasian" "Test_Caucasian"
    ## [17] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Caucasian"
    ## [21] "Test_Caucasian" "Test_Asian"     "Test_Caucasian" "Caucasian"     
    ## [25] "Test_Caucasian" "Caucasian"      "Caucasian"      "Test_Caucasian"
    ## [29] "Caucasian"      "Caucasian"      "Test_Asian"     "Test_Asian"    
    ## [33] "Test_Caucasian" "Caucasian"      "Test_Caucasian" "Test_Caucasian"
    ## [37] "Test_Caucasian" "Test_Caucasian" "Test_Asian"     "Caucasian"     
    ## [41] "Test_Asian"     "Asian"          "Asian"          "Test_Caucasian"
    ## [45] "Test_Caucasian" "Caucasian"      "Test_Caucasian" "Test_Caucasian"
    ## [49] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Caucasian"     
    ## [53] "Asian"          "Caucasian"      "Caucasian"      "Caucasian"     
    ## [57] "Caucasian"      "Asian"          "Test_Caucasian" "Caucasian"     
    ## [61] "Caucasian"      "Test_Caucasian" "Test_Asian"     "Caucasian"     
    ## [65] "Asian"          "Caucasian"      "Asian"          "Caucasian"     
    ## [69] "Test_Caucasian" "Caucasian"      "Caucasian"      "Caucasian"     
    ## [73] "Caucasian"      "Caucasian"      "Caucasian"      "Caucasian"     
    ## [77] "Asian"          "Asian"          "Test_Caucasian" "Test_Caucasian"
    ## [81] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Test_Caucasian"
    ## [85] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Asian"         
    ## [89] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Asian"    
    ## [93] "Test_Caucasian" "Caucasian"      "Asian"          "Caucasian"     
    ## [97] "Asian"

``` r
dendro.weight = cutree(hclust.weight, 2)
ColorDendrogram(hclust.weight, y = dendro.weight, labels = names(dendro.weight), branchlength = 0.3, main = 'Clustering train and test with labels from weighted glmnet')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/dendrogram%20weighted%20glmnet-1.png)

up-sampling
-----------

instead of having the model dealing with imbalanced ratio of classes, we can attempt to balance the class frequencies. There are post-hoc sampling approaches that can help attenuate the effects of the imbalance during model training. Two general post hoc approaches are down-sampling and up-sampling the data. Here we will try out up-sampling, which is a technique that simulates additional data points to improve balance across classes.

-   From: [How do I handle an unbalanced training set?](https://www.researchgate.net/post/In_classification_how_do_i_handle_an_unbalanced_training_set)

``` r
# Build up-sampled model

fitControl$sampling <- "up"

up_fit <- train(x = x.train.redo,
                      y = y.train,
                      method = "glmnet",
                      metric = "ROC",
                      trControl = fitControl,
                                    tuneGrid = netGrid,
                                    preProcess = c( "center", "scale"))

up_fit
```

    ## glmnet 
    ## 
    ##    45 samples
    ## 10551 predictors
    ##     2 classes: 'Asian', 'Caucasian' 
    ## 
    ## Pre-processing: centered (10551), scaled (10551) 
    ## Resampling: Cross-Validated (5 fold, repeated 3 times) 
    ## Summary of sample sizes: 35, 35, 37, 36, 37, 36, ... 
    ## Addtional sampling using up-sampling prior to pre-processing
    ## 
    ## Resampling results across tuning parameters:
    ## 
    ##   alpha  lambda  ROC        Sens       Spec     
    ##   0.5    0.1     1.0000000  0.8888889  1.0000000
    ##   0.5    0.2     1.0000000  0.8555556  1.0000000
    ##   0.5    0.3     1.0000000  0.8333333  1.0000000
    ##   0.5    0.4     1.0000000  0.8333333  0.9904762
    ##   0.5    0.5     1.0000000  0.8111111  0.9904762
    ##   0.6    0.1     1.0000000  0.8777778  1.0000000
    ##   0.6    0.2     1.0000000  0.8000000  1.0000000
    ##   0.6    0.3     0.9968254  0.7777778  1.0000000
    ##   0.6    0.4     0.9920635  0.7222222  1.0000000
    ##   0.6    0.5     0.9645503  0.7222222  0.9698413
    ##   0.7    0.1     1.0000000  0.7666667  1.0000000
    ##   0.7    0.2     0.9936508  0.7333333  1.0000000
    ##   0.7    0.3     0.9936508  0.7111111  0.9904762
    ##   0.7    0.4     0.9698413  0.7222222  0.9603175
    ##   0.7    0.5     0.9386243  0.6444444  0.9269841
    ##   0.8    0.1     0.9904762  0.7444444  1.0000000
    ##   0.8    0.2     0.9772487  0.7111111  0.9904762
    ##   0.8    0.3     0.9576720  0.7555556  0.9603175
    ##   0.8    0.4     0.9370370  0.6666667  0.9492063
    ##   0.8    0.5     0.9256614  0.6444444  0.9571429
    ##   0.9    0.1     0.9931217  0.7444444  1.0000000
    ##   0.9    0.2     0.9629630  0.7000000  0.9888889
    ##   0.9    0.3     0.9457672  0.7000000  0.9666667
    ##   0.9    0.4     0.9161376  0.6666667  0.9666667
    ##   0.9    0.5     0.6497354  0.4888889  0.7111111
    ## 
    ## ROC was used to select the optimal model using  the largest value.
    ## The final values used for the model were alpha = 0.5 and lambda = 0.5.

``` r
# prediction, classification results

y.predict.up <- predict(up_fit,  x.test.redo)
y.predict.up
```

    ##  [1] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ##  [8] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [15] Asian     Caucasian Caucasian Asian     Caucasian Caucasian Caucasian
    ## [22] Caucasian Caucasian Asian     Caucasian Caucasian Asian     Caucasian
    ## [29] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [36] Asian     Caucasian Caucasian Caucasian Caucasian Caucasian Asian    
    ## [43] Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian Caucasian
    ## [50] Asian     Caucasian Asian    
    ## Levels: Asian Caucasian

``` r
y.predict.up.des <- data.frame(Samplename = rownames(x.test.redo), Ethnicity = paste("Test",y.predict.up,sep = "_")) # to be used for dendrogram

# predicted probability to be Asian
y.predict.up <- predict(up_fit,  x.test.redo, type = "prob")
y.predict.up[,"Asian"]
```

    ##  [1] 0.3626702 0.4096970 0.4085382 0.3862106 0.3924155 0.4289514 0.4343024
    ##  [8] 0.3774437 0.3928205 0.3971509 0.3971190 0.3523198 0.4183526 0.3991514
    ## [15] 0.5361824 0.3596419 0.4899164 0.5383241 0.4204269 0.3870706 0.3313562
    ## [22] 0.3827367 0.3920565 0.5820861 0.4249275 0.4104848 0.5012311 0.3761821
    ## [29] 0.3471587 0.3244149 0.4026396 0.4978440 0.4641937 0.4133080 0.4068565
    ## [36] 0.5372756 0.4184331 0.4209291 0.4091264 0.4582510 0.3701803 0.5185645
    ## [43] 0.4762741 0.3767989 0.4645193 0.4610973 0.4447341 0.3824280 0.4215725
    ## [50] 0.5050892 0.3320359 0.5475099

``` r
# histogramfor the prob. to be Asian
hist(y.predict.up[,"Asian"], main = "Predicted Probability to be Asian")
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/fit%20up-sampling%20glmnet-1.png)

Now map predicted results to dendrogram from hierarchical clustering of the merged data:

``` r
# make a design matrix containing info on whether it's training or test data, and ethnicity info

up.des <- rbind(train.design[,c("Samplename","Ethnicity")],y.predict.up.des)

# clustering for (centered and scaled) merged data
hclust.up <- hclust(dist(scale(merged.trunc,center = T,scale = T), method = 'euclidean'), method = "average")

labels.up <- data.frame(labels(hclust.up))   # pulls out current labels (samplename)
colnames(labels.up) <- 'Samplename'
labels.up <- labels.up %>% left_join(up.des, by = 'Samplename')
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
labels(hclust.up) <- labels.up$Samplename


hclust.up <- swaplabels(hclust.up, up.des)
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
labels(hclust.up)
```

    ##  [1] "Caucasian"      "Caucasian"      "Caucasian"      "Test_Caucasian"
    ##  [5] "Test_Asian"     "Asian"          "Test_Caucasian" "Test_Caucasian"
    ##  [9] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Asian"    
    ## [13] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Caucasian"
    ## [17] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Caucasian"
    ## [21] "Test_Caucasian" "Test_Asian"     "Test_Caucasian" "Caucasian"     
    ## [25] "Test_Caucasian" "Caucasian"      "Caucasian"      "Test_Caucasian"
    ## [29] "Caucasian"      "Caucasian"      "Test_Asian"     "Test_Asian"    
    ## [33] "Test_Caucasian" "Caucasian"      "Test_Caucasian" "Test_Caucasian"
    ## [37] "Test_Caucasian" "Test_Caucasian" "Test_Asian"     "Caucasian"     
    ## [41] "Test_Asian"     "Asian"          "Asian"          "Test_Caucasian"
    ## [45] "Test_Caucasian" "Caucasian"      "Test_Caucasian" "Test_Caucasian"
    ## [49] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Caucasian"     
    ## [53] "Asian"          "Caucasian"      "Caucasian"      "Caucasian"     
    ## [57] "Caucasian"      "Asian"          "Test_Caucasian" "Caucasian"     
    ## [61] "Caucasian"      "Test_Caucasian" "Test_Caucasian" "Caucasian"     
    ## [65] "Asian"          "Caucasian"      "Asian"          "Caucasian"     
    ## [69] "Test_Caucasian" "Caucasian"      "Caucasian"      "Caucasian"     
    ## [73] "Caucasian"      "Caucasian"      "Caucasian"      "Caucasian"     
    ## [77] "Asian"          "Asian"          "Test_Caucasian" "Test_Caucasian"
    ## [81] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Test_Caucasian"
    ## [85] "Test_Caucasian" "Test_Caucasian" "Caucasian"      "Asian"         
    ## [89] "Test_Caucasian" "Test_Caucasian" "Test_Caucasian" "Test_Asian"    
    ## [93] "Test_Caucasian" "Caucasian"      "Asian"          "Caucasian"     
    ## [97] "Asian"

``` r
dendro.up = cutree(hclust.up, 2)
ColorDendrogram(hclust.up, y = dendro.up, labels = names(dendro.up), branchlength = 0.3, main = 'Clustering train and test with labels from up-sampling glmnet')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/dendrogram%20up-sampling-1.png)

It is disappointing that the dendrograms cannot separate Asians from Caucasians that well, as Euclidean distance is not well-suited for calculating the distance for high-dimensional data. Good news is that both weighted glmnet and up-sampling boosted the number of predicted Asians to around 8 samples, which is slightly more believable. Still, more work can be done in the future.
