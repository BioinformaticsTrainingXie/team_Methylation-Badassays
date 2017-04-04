ModelBuilding
================
Victor
April 2, 2017

-   [Parallel processing](#parallel-processing)
-   [Load data](#load-data)
-   [Section 1.0 Prefiltering Features](#section-1.0-prefiltering-features)
    -   [1.1 Remove NAs](#remove-nas)
    -   [1.2 Reduce CpGs to match test and train](#reduce-cpgs-to-match-test-and-train)
    -   [1.3 Prefiltering cpgs (most variable)](#prefiltering-cpgs-most-variable)
    -   [2.1 logistic regression with elastic net regularization](#logistic-regression-with-elastic-net-regularization)
    -   [2.2 SVM linear](#svm-linear)
-   [Section 3.0 Predict Ethnicity for Test Set](#section-3.0-predict-ethnicity-for-test-set)
    -   [3.1 glmnet](#glmnet)
    -   [3.2 SVM](#svm)
-   [Section 4.0 Analysis of Predictors](#section-4.0-analysis-of-predictors)
    -   [4.2 Plot CpG Predictors](#plot-cpg-predictors)
-   [Section 5.0](#section-5.0)

``` r
#source("https://bioconductor.org/biocLite.R")
#biocLite('e1071')                                    # required for glmnet in caret
#biocLite('pROC')
library(pROC)
```

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
library(ggplot2)
library(limma)
library(caret)
```

    ## Loading required package: lattice

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(parallel)
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

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

Section 1.0 Prefiltering Features
=================================

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

![](PredictiveModelingVY_files/figure-markdown_github/prefiltering%20based%20on%20SD-1.png)

``` r
# filter CpG sites with low s.d: only keep those with s.d higher than the average s.d across all CpG sites
train.gsd <- subset(train.sd, train.sd > 0.10)
hist(train.gsd)
```

![](PredictiveModelingVY_files/figure-markdown_github/prefiltering%20based%20on%20SD-2.png)

``` r
train.data.gsd <- train.data[,colnames(train.data) %in% names(train.gsd)]
```

We reduced the \# of features to 'r ncol(train.data.gsd)' to reduce computation time. train.data.gsd is the working dataset \# Section 2.0 Supervised classification:

2.1 logistic regression with elastic net regularization
-------------------------------------------------------

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
    ##   22.26    0.19  231.31

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

Looks like our glmnet-built model has chosen 'r length(predictors)' CpGs that can be used to predict ethnicity.

2.2 SVM linear
--------------

This section is for building the model using SVM. However, because computational time is long, this section is can be excluded if chosen (specify eval = false)

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
    ##   25.23    0.22  240.09

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

Section 3.0 Predict Ethnicity for Test Set
==========================================

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

3.2 SVM
-------

``` r
y.predictSVM <- predict(svmFit,  x.test)
#throws a warning
y.predictSVM
```

Section 4.0 Analysis of Predictors
==================================

Here we pull out the CpG sites and look at them more closely. First we will see if clustering with only the predictors separates asians and caucasians \#\# 4.1 Clustering

``` r
library(ggdendro)
library(sparcl) # ColorDendrogram
library(dendextend)
```

    ## 
    ## ---------------------
    ## Welcome to dendextend version 1.5.2
    ## Type citation('dendextend') for how to cite the package.
    ## 
    ## Type browseVignettes(package = 'dendextend') for the package vignette.
    ## The github page is: https://github.com/talgalili/dendextend/
    ## 
    ## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
    ## Or contact: <tal.galili@gmail.com>
    ## 
    ##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
    ## ---------------------

    ## 
    ## Attaching package: 'dendextend'

    ## The following object is masked from 'package:ggdendro':
    ## 
    ##     theme_dendro

    ## The following object is masked from 'package:stats':
    ## 
    ##     cutree

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

![](PredictiveModelingVY_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-1.png)

``` r
#with predictors only
x.train.predictors <- x.train[,colnames(x.train) %in% predictorsNet]
hclust2 <- hclust(dist(x.train.predictors, method = 'euclidean'))
hclust2 <- swaplabels(hclust2, design)          #swap labels with ethnicity
y2 = cutree(hclust2, 2)
ColorDendrogram(hclust2, y = y2, labels = names(y2), branchlength = 0.3, main = 'Clustering train with predictors only')
```

![](PredictiveModelingVY_files/figure-markdown_github/clustering%20train%20based%20on%20predictors-2.png)

``` r
# Hierarchical clustering of predicted data, distance measured by Euclidean distance, average linkage
hclust3 <- hclust(dist(x.test, method = 'euclidean'))
y3 = cutree(hclust3, 2)
ColorDendrogram(hclust3, y=y3, labels = names(y3), branchlength = 2, main = 'Clustering Test with all CpGs')
```

![](PredictiveModelingVY_files/figure-markdown_github/Clustering%20test%20data-1.png)

``` r
# clustering only with the predictors
x.test.predictors <- x.test[,colnames(x.test) %in% predictorsNet]
hclust4 <- hclust(dist(x.test.predictors, method = 'euclidean'))
labels(hclust4) <- ''
```

    ## Warning in `labels<-.hclust`(`*tmp*`, value = ""): The lengths of the new
    ## labels is shorter than the number of leaves in the hclust - labels are
    ## recycled.

``` r
y4 = cutree(hclust4, 1)
ColorDendrogram(hclust4, y=y4, labels = names(y4), branchlength = 0.25, main = 'Clustering Test with predictors only')
```

![](PredictiveModelingVY_files/figure-markdown_github/Clustering%20test%20data-2.png)

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
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
#replace train samples with ethnicity labels
#labels5$Samplename[!is.na(labels5$Ethnicity)] <- as.character(labels5$Ethnicity[!is.na(labels5$Ethnicity)])

labels(hclust5) <- labels5$Samplename
labels5
```

    ##    Samplename Ethnicity
    ## 1       PM104 Caucasian
    ## 2       PM115 Caucasian
    ## 3        PM20 Caucasian
    ## 4        PM52 Caucasian
    ## 5        PM55 Caucasian
    ## 6       PM123 Caucasian
    ## 7        PM58 Caucasian
    ## 8       PM181 Caucasian
    ## 9        PM74 Caucasian
    ## 10      PM124 Caucasian
    ## 11       PM66 Caucasian
    ## 12        PM4 Caucasian
    ## 13      PM167 Caucasian
    ## 14      PM205 Caucasian
    ## 15      PM142 Caucasian
    ## 16       PM44 Caucasian
    ## 17       PM72 Caucasian
    ## 18      PM227 Caucasian
    ## 19        PM9 Caucasian
    ## 20      PM155 Caucasian
    ## 21      PM249 Caucasian
    ## 22      PM136 Caucasian
    ## 23      PM243 Caucasian
    ## 24     FT35_v      <NA>
    ## 25     FT67_v      <NA>
    ## 26       PM30 Caucasian
    ## 27     FT28_v      <NA>
    ## 28     FT29_v      <NA>
    ## 29    PL148_v      <NA>
    ## 30     FT23_v      <NA>
    ## 31     FT65_v      <NA>
    ## 32     NTD1_v      <NA>
    ## 33     FT52_v      <NA>
    ## 34      PM120 Caucasian
    ## 35     FT13_v      <NA>
    ## 36      PM153 Caucasian
    ## 37      PM130 Caucasian
    ## 38       PM46 Caucasian
    ## 39      PM158 Caucasian
    ## 40       PM71 Caucasian
    ## 41     NTD6_v      <NA>
    ## 42     NTD9_v      <NA>
    ## 43      PM114 Caucasian
    ## 44       PM84 Caucasian
    ## 45     FT59_v      <NA>
    ## 46       PM29     Asian
    ## 47     NTD3_v      <NA>
    ## 48       PM54 Caucasian
    ## 49     FT85_v      <NA>
    ## 50     FT82_v      <NA>
    ## 51    NTD17_v      <NA>
    ## 52     FT60_v      <NA>
    ## 53     FT54_v      <NA>
    ## 54      FT3_v      <NA>
    ## 55     FT57_v      <NA>
    ## 56     FT64_v      <NA>
    ## 57     FT40_v      <NA>
    ## 58     NTD8_v      <NA>
    ## 59     FT22_v      <NA>
    ## 60     FT84_v      <NA>
    ## 61     FT16_v      <NA>
    ## 62     FT21_v      <NA>
    ## 63     FT47_v      <NA>
    ## 64     FT38_v      <NA>
    ## 65     FT78_v      <NA>
    ## 66     FT27_v      <NA>
    ## 67     FT72_v      <NA>
    ## 68      FT5_v      <NA>
    ## 69     FT75_v      <NA>
    ## 70     FT34_v      <NA>
    ## 71     FT36_v      <NA>
    ## 72      PM119     Asian
    ## 73       PM98     Asian
    ## 74    mt4.5_v      <NA>
    ## 75     FT41_v      <NA>
    ## 76    NTD16_v      <NA>
    ## 77     FT79_v      <NA>
    ## 78    PL118_v      <NA>
    ## 79     FT73_v      <NA>
    ## 80     FT74_v      <NA>
    ## 81     FT26_v      <NA>
    ## 82     NTD2_v      <NA>
    ## 83     FT62_v      <NA>
    ## 84    NTD14_v      <NA>
    ## 85      PM226     Asian
    ## 86       PM47     Asian
    ## 87       PM76     Asian
    ## 88       PM41     Asian
    ## 89       PM53     Asian
    ## 90      PM112     Asian
    ## 91      PM139     Asian
    ## 92      PM233     Asian
    ## 93       PM40     Asian
    ## 94     FT58_v      <NA>
    ## 95     FT39_v      <NA>
    ## 96     FT42_v      <NA>
    ## 97    PL149_v      <NA>

``` r
hclust5 <- swaplabels(hclust5, design)
```

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
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

![](PredictiveModelingVY_files/figure-markdown_github/cluster%20both%20test%20and%20train-1.png)

4.2 Plot CpG Predictors
-----------------------

``` r
glmImp <- varImp(netFit, scale = F) # gives the t-statistic for all CpGs in the dataset
plot(glmImp, top = 11)
```

![](PredictiveModelingVY_files/figure-markdown_github/plot%20top%2035-1.png)

Section 5.0
===========

This section is to tune parameters across 100 different combinations of alpha and lambda

``` r
netGrid100 <-  expand.grid(alpha = c(0.20, 0.40, 0.60, 0.80, 1.00),
                           lambda = c(0.05, 0.10, 0.15, 0.20, 0.25))

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
    ##   22.31    0.16 1031.61

``` r
trellis.par.set(caretTheme())
ggplot(netFit100)
#heatmap of results
plot(netFit100, metric = "ROC", plotType = "level",
     scales = list(x = list(rot = 90)))
length(predictors(netFit100))
glmImp100 <- varImp(netFit100, scale = F) # gives the t-statistic for all CpGs in the dataset
plot(glmImp100, top = 50)
```
