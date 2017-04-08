CV\_TestPerformance
================
Ming, Victor
April 4, 2017

(Repeated CV takes a long time to run for `glmnet`. Hence to ensure knitting finishes fast, we excluded CV part for comparing performance between glmnet vs. SVM in "BuildModel\_AnalyzePredictors.Rmd". CV is done in this file and CV results were saved as objects in data/R objects folder.

Step 1&2 are the same as in "BuildModel\_AnalyzePredictors.Rmd", only step 3 is new.)

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
setwd('../')                                           # note: all of these relative file path calls work only for knitting

# load data (pre-processed training set)
train.data <- read.table('../data/Processed Data/data.txt')
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
train.design <- read.csv("../data/Processed Data/des.txt", sep="\t", header=TRUE)
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

![](CV_TestPerformance_files/figure-markdown_github/prefiltering%20based%20on%20SD-1.png)

``` r
# filter CpG sites with low s.d: only keep those with s.d higher than the average s.d across all CpG sites
train.gsd <- subset(train.sd, train.sd > 0.10)
hist(train.gsd)
```

![](CV_TestPerformance_files/figure-markdown_github/prefiltering%20based%20on%20SD-2.png)

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
    ##   13.08    0.05   95.03

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
    ##   14.92    0.05   99.46

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

S3 Compare Model Performance Using Cross Validation
===================================================

``` r
# this is adapted from Amrit's code from lec 19
#use repeated CV to estimate test performance
set.seed(2018)

# list of lists containing fold ids
folds <- lapply(1:M, function(i) createFolds(y.train, k = k))

system.time(netTesterror <- lapply(folds, function(i){
  lapply(i, function(j){
    # tune parameters with CV
    set.seed(2019)
    fitControl <- trainControl(method = "repeatedcv", 
                                                     number = k,                 
                                                     repeats = M,
                                                     classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                                                     savePredictions = TRUE,
                                                     allowParallel = T)
    
    
    # build elastic net classifier
    netFit <- train(x =  x.train[-j,],   
                                y = y.train[-j],                  
                                method = "glmnet",                     
                                trControl = fitControl,
                                preProcess = c( 'center', 'scale'),
                                metric = 'ROC')   
    
    # Estimate probabilities of test predictions
    probTest <- predict(netFit, x.train[j,], type = 'prob')
    ethProb <- probTest[,'Asian']
    ethProb
  })

})
)
```

    ##    user  system elapsed 
    ##  457.75    0.60 4036.49

``` r
# netTesterror
#saveRDS(netTesterror,"netTesterror.rds")
```

``` r
# Computer classification performance measures
# enet
Performance <- mapply(function(x, y){
  auc <- pROC::roc(y.train[unlist(x)], unlist(y),
                   direction ='<',
                   levels = c('Caucasian', 'Asian'),
                   percent = TRUE)
  list(tpr = auc$sensitivities,
       fpr = 100 - auc$specificities,
       auc = round(auc$auc, 2))
}, x = folds, y = netTesterror)
Performance
```

    ##     [,1]       [,2]       [,3]      
    ## tpr Numeric,46 Numeric,46 Numeric,46
    ## fpr Numeric,46 Numeric,46 Numeric,46
    ## auc 100        100        98.23

``` r
# plot ROC curve

plot(Performance['tpr',][[1]] ~ Performance['fpr',][[1]],
     type = 'l', col = 1, xlab = '100 - sensitivity',
     ylab = 'Sensitivity', main = 'Enet')
for(i in length(folds)){
  points(Performance['tpr',][[i]] ~ Performance['fpr',][[i]],
         type = 'l', col = 2)
}
text(x = 60, y = 40, labels =
       paste0('mean AUC = ', round(mean(unlist(Performance['auc',])), 1),
              '+/-', round(sd(unlist(Performance['auc',])), 1), '%'))
```

![](CV_TestPerformance_files/figure-markdown_github/ROC%20curve-1.png)

Support Vector Machine
----------------------

We do not need repeated CV for linear kernel SVM as it does not have tuning parameters. Regular CV is conducted for SVM instead.

``` r
# Linear kernel SVM with CV
    fitControl <- trainControl(method = "cv", 
                                                     number = k,                 
                                                     classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                                                     savePredictions = TRUE,
                                                     allowParallel = T)
```

``` r
# this is adapted from Amrit's code from lec 19
#use repeated CV to estimate test performance
set.seed(2018)

# list of lists containing fold ids
folds <- createFolds(y.train, k = k)

system.time(svmTesterror <- lapply(folds, function(j){
    # tune parameters with CV
    fitControl <- trainControl(method = "cv", 
                                                     number = k,                 
                                                     classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                                                     savePredictions = TRUE, allowParallel = T)
    
    
    # build SVM classifier
    svmFit <- train(x=x.train[-j,],
                        y= y.train[-j],
                        method = "svmLinear",
                        preProc = c("center","scale"),
                        metric="ROC",
                        trControl=fitControl)
    
    # Estimate probabilities of test predictions
    probTest <- predict(svmFit, x.train[j,], type = 'prob')
    ethProb <- probTest[,'Asian']
    ethProb

})
)
```

    ##    user  system elapsed 
    ##  100.96    0.11  271.05

``` r
#svmTesterror
```

``` r
# Computer classification performance measures
# enet
Performance <- mapply(function(x, y){
  auc <- pROC::roc(y.train[unlist(x)], unlist(y),
                   direction ='<',
                   levels = c('Caucasian', 'Asian'),
                   percent = TRUE)
  list(tpr = auc$sensitivities,
       fpr = 100 - auc$specificities,
       auc = round(auc$auc, 2))
}, x = folds, y = svmTesterror)
Performance
```

    ##     Fold1      Fold2      Fold3      Fold4     Fold5     
    ## tpr Numeric,10 Numeric,10 Numeric,10 Numeric,9 Numeric,11
    ## fpr Numeric,10 Numeric,10 Numeric,10 Numeric,9 Numeric,11
    ## auc 100        88.89      100        100       100

``` r
# plot ROC curve

plot(Performance['tpr',][[1]] ~ Performance['fpr',][[1]],
     type = 'l', col = 1, xlab = '100 - sensitivity',
     ylab = 'Sensitivity', main = 'Linear Kernel SVM')
for(i in length(folds)){
  points(Performance['tpr',][[i]] ~ Performance['fpr',][[i]],
         type = 'l', col = 2)
}
text(x = 60, y = 40, labels =
       paste0('mean AUC = ', round(mean(unlist(Performance['auc',])), 1),
              '+/-', round(sd(unlist(Performance['auc',])), 1), '%'))
```

![](CV_TestPerformance_files/figure-markdown_github/svm%20ROC%20curve-1.png)
