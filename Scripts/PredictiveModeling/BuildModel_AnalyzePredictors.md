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
    -   [2.2 SVM linear](#svm-linear)
-   [S3 Predict Ethnicity for external data Set](#s3-predict-ethnicity-for-external-data-set)
    -   [3.1 glmnet](#glmnet)
    -   [3.2 SVM](#svm)
-   [S4 Analysis of Predictors](#s4-analysis-of-predictors)
    -   [4.2 Plot CpG Predictors](#plot-cpg-predictors)
-   [S5 Tune alpha and lambda 10 x 10 grid](#s5-tune-alpha-and-lambda-10-x-10-grid)

S0 Set up workspace
===================

Load packages
-------------

``` r
#source("https://bioconductor.org/biocLite.R")
#biocLite('e1071')                                    # required for glmnet in caret
#biocLite('pROC')
library(pROC)
library(ggplot2)
library(limma)
library(caret)
library(dplyr)
library(parallel)
library(doParallel)
```

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
design <- read.csv("../data/Processed Data/des.txt", sep="\t", header=TRUE)
str(design)
```

    ## 'data.frame':    45 obs. of  5 variables:
    ##  $ Samplename  : Factor w/ 45 levels "PM104","PM112",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ Sample_Group: Factor w/ 3 levels "CONTROL","IUGR",..: 1 1 1 3 3 1 2 1 2 1 ...
    ##  $ ga          : num  40.7 38.9 38.6 41.1 37.1 38 35.7 40 36.9 38.6 ...
    ##  $ sex         : Factor w/ 2 levels "F","M": 2 1 2 2 2 1 1 2 2 1 ...
    ##  $ Ethnicity   : Factor w/ 2 levels "Asian","Caucasian": 2 1 2 2 1 2 2 2 2 2 ...

``` r
row.names(train.data) == design$Samplename               # check that the samples are in same order
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [29] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [43] TRUE TRUE TRUE

Read in test data:

``` r
setwd('../') # for knitting
# read pre-processed test data
test.data <- read.table("../Data/Processed Data/Test data/Matrix.processed.betas.placenta.txt", row.names = 1, header = T)
test.data <- as.data.frame(t(test.data))   #transpose data
```

S1 Prefiltering Features
========================

Reducing the number of features can significantly reduce computational time, which is desirable when the dataset is large. However, we must be careful not remove potentially 'interesting' features that have a high chance of being useful in building a classifier.

1.1 Remove NAs
--------------

We should remove any sites with NAs in them or else predictions cannot be generated for these samples if the CpG site is chosen as a predictor

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
sd(train.sd)
```

    ## [1] 0.02742435

``` r
hist(train.sd)                    # histogram
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
train.data.gsd <- train.data[,colnames(train.data) %in% names(train.gsd)]
```

We reduced the \# of features to 'r ncol(train.data.gsd) to reduce computation time. train.data.gsd is the working dataset \# S2 Supervised classification: We decided to try two different models for building our classifer: elastic net logistic regression (glmnet) and support vector machines (SVM). Both of these models have been used in the literature to build predictive models based on 450k DNA methylation data (Horvath 2013, De Carli et al 2017), indicating that they may be well-suited for our dataset. \#\# 2.1 logistic regression with elastic net regularization

``` r
#renamed just so that I can copy Amrit's code
x.train <- train.data.gsd
y.train <- design$Ethnicity 
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
                                                     allowParallel = TRUE
                                                     )  

netGrid <- expand.grid(alpha = c(0.75),
                           lambda = c(0.077, 0.25))
```

We specify the model to be built using repeated cross validation with a fold = 5, and repeats = 3. We tune the model holding alpha constant (alpha = 0.75), keeping alpha high to favour L1 norm to achieve a small panel of biomarkers. Lambda, the magnitude of the penalty, is tested at 0.077, and 0.25.

``` r
set.seed(2017)                                         # training models requires the use of random #s. Setting (set.seed()) the randomness ensures reproducibility

#netGrid <- expand.grid(.alpha = seq(.05, 1, length = 15),
 #                                                   .lambda = c((1:5)/10)) # grid of tuning parameters to try out

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
    ##   24.95    1.99  247.61

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

Looks like our glmnet-built model has chosen 'r length(predictorsNet)' CpGs that can be used to predict ethnicity.

2.2 SVM linear
--------------

This section is for building the model using SVM. However, because computational time is long, this section is can be excluded when ran, since we have chosen the glmnet model to be our final model.

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
    ##   25.75    0.16  256.42

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

Next, we use the model we built and run it on an external data set, where there is no ethnicity information.

3.1 glmnet
----------

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

It looks like our model classifies the entire external dataset is Caucasian. This is suspicious, as we believe the samples to come from a relatively heterogenous population. However, due to time constraints, we decided to move ahead and perform downstream analysis. If there was more time, we might think about where we can change our model tuning process to produce more sensible results.

#### Some explanations for this result:

-   It's possible that the data set is truly all Caucasian.
-   The dataset is too 'different' to have the classifier ran on. (too much noise)
-   The self-reported ethnicities in the training data is too unreliable

3.2 SVM
-------

``` r
y.predictSVM <- predict(svmFit,  x.test)
#throws a warning
y.predictSVM
```

S4 Analysis of Predictors
=========================

Here we pull out the CpG sites and look at them more closely. First we will see if clustering with only the predictors separates asians and caucasians \#\# 4.1 Clustering

``` r
library(ggdendro)
library(sparcl) # ColorDendrogram
library(dendextend)
```

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

hclust <- swaplabels(hclust, design)
y1 = cutree(hclust, 3)
ColorDendrogram(hclust, y = y1, labels = names(y1), branchlength = 1.0, main = 'Clustering train on all CpGs')
```

![](BuildModel_AnalyzePredictors_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-1.png)

``` r
#with predictors only
x.train.predictors <- x.train[,colnames(x.train) %in% predictorsNet]
hclust2 <- hclust(dist(x.train.predictors, method = 'euclidean'))
hclust2 <- swaplabels(hclust2, design)          #swap labels with ethnicity
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
#x.train.predictors <- x.train.predictors %>% left_join(select(design, Samplename, Ethnicity), by = 'Samplename') 
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
labels5 <- labels5 %>% left_join(select(design, Samplename, Ethnicity), by = 'Samplename')
#replace train samples with ethnicity labels
#labels5$Samplename[!is.na(labels5$Ethnicity)] <- as.character(labels5$Ethnicity[!is.na(labels5$Ethnicity)])

labels(hclust5) <- labels5$Samplename


hclust5 <- swaplabels(hclust5, design)
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
                                left_join(design, 'Samplename')
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
                                left_join(design, 'Samplename')
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
    ##   22.61    6.06 2014.65

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
