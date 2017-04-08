## The Methylation Badassays

### Introduction:

DNA methylation (DNAm) is linked to many diseases like cancer and autism. However DNAm marks can change due to environmental stimuli, cell types, gender, etc. In recent years, several DNAm studies have suggested that a large portion of DNAm variability is associated with genetic ancestry and is heritable, making DNAm a potential confounding factor which is not given enough consideration in the context of DNA methylation analysis. Differentially methylated CpG sites associated with pathology can be confounded by CpGs associated with genetic ancestry causing suprious results. Therefore, to invetigate how DNA methylation affects prenatal health, it is important for us to identify genetic ancestry-associated CpGs to figure out true positives. This DNAm variability in the placenta due to genetic ancestry needs to be accounted for in large scale DNAm studies, or else no meaningful interpretation of results can be done to assess prenatal health. In this project, we are going to investigate if DNA in placental tissue is differentially methylated across populations of different ancestry. 

**Hypothesis**: DNA in placental tissue is differentially methylated across populations of different ancestries.

We will first find methylation profiles in subjects from our [dataset 1](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) and the genetic ancestry (Asians or Caucasians) of our subjects is known. These profiles will then serve as a basis to cluster methylation data in our [dataset 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69502) in which genetic ancestry is not known. For both dataset 1 and dataset 2, DNAm was measured by 450K microarray from Illumina.

For the details of the project ideas, dataset and methods we used for this project, please check our [project proposal](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/project_proposal.md). 

------

### Team Member:

|  Name  | Department/Program  |GitHub ID |
|-------|---------------------|----------|
| Victor Yuan| Genome Science and Technology| @wvictor14 |
| Michael Yuen|Medical Genetics|@myuen89|
|Nivretta Thatra|Bioinformatics|@nivretta|
|Ming Wan|Statistics|@MingWan10|
|Anni Zhang|Genome Science and Technology|@annizubc|


------
### Workflow

This figure summarizes the workflow of our project:

![workflow](https://cloud.githubusercontent.com/assets/24922214/24690299/deb7dc8c-1980-11e7-9554-ec0ca92f4038.png)

------

**For all our processing steps below, please see the [Results](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Results) folder for a more detailed write up on our findings from our analyses. If you're interested in the code, see the markdown files in the [Scripts](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts) folder.**

### Preprocessing and Normalization

We first used this [script](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Preprocessing/PreprocessQC.md) to process (via quality control, filtering, and normalization) the [raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) of dataset 1 into to our [processed data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Processed%20Data). For detailed information of dataset 1, please see [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Data/Raw%20Data/samplesheet.csv).


### Exploratory Analysis

We [explored](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory.md) our data by generating sample-sample correlation heatmaps, plotting a few random CpGs and plotting the first few principal components.

### Differential Methylation Analysis

We used the R package [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) to identify differentially methylated probes between Asian and Caucasian samples. Please see our differential DNA methylation analysis [script](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma.md) in the [scripts folder](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts) for the code and details. Limma prioritized 13 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry using a cutoff off p value = 0.01.

### Building an Ancestry Classifer

To build the DNA methylation ancestry classifer, we compare [SVM](http://ca.wiley.com/WileyCDA/WileyTitle/productCd-0471030031.html) and [elastic net logistic regression (glmnet)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115) models. We ended up choosing glmnet for building the final model, and used a nested cross validation strategy to tune the penalization parameters, and for estimating the test error. After generating the final model, we analyzed the predictors, and examined the results of the predictions on the secondary unlabelled dataset. Please see the subdirectory [predictive modeling](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts/PredictiveModeling) for the markdown files and details. 

### Brief Functional Analysis

We looked the 13 CpG sites prioritized by limma and the 11 CpG sites prioritized by glmnet, in this [script](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/FunctionalAnalysis/FunctionalAnalysis.md). Using the COHCAP (City of Hope CpG Island Analysis Pipeline) package, the CpGs we mapped to chromosome, location, gene name and CpG island information. Each gene was annotated with its GO term using the package mygene.

### Summary

* Please see our [poster](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/poster.pdf)! :smile:

* SVM performed slightly better than glmnet (for both training and testing error)

* Final model used 11 CpG predictors and was built with glmnet with a AUC of 0.981 and 0.977+-0.024 for training and testing error respectively (α = 0.75, λ = 0.25).

* The classifier predicted all of the unlabeled test set to Caucasian, which we doubt is the true case.

* We suspect the test set is too ‘different’ from the training data set for the classifier to perform accurately on the test set

### Future Direction

* Normalizing and QCing the test and training datasets together may be necessary for DNA methylation classifiers to perform well.

* Using MDS ancestry coordinates from population stratification meta-analyses may provide ‘labels’ to assess classifier performance or iprove model building. (self-reported ancestry can be unreliable)


------
### Table of contents:

1. [Project proposal](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/project_proposal.md): includes the introduction to the ideas, dataset and methods we used in this project.

2. [Progress report](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/progress_report.md): contains the progress about our project.

3. [Data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data) folder contains metadata, raw data, and processed data.

    * [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Data/Raw%20Data/samplesheet.csv)
      + human placental tissue from 45 subjects with self reported ancestry
      + columns correspond to subject ancestry, name, sex, gestational age and what complications they had in pregnancy (none, intrauterine growth (IUGR) restriction, or late onset preeclampsia (LOPET), neither of which affect DNAm)
      + columns for Sentrix ID and position correspond to the sample’s batch ID and position on the Illumina microarray 
      + each row is one subject.
      
    * [Raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) for dataset 1.
    
    * [Processed data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Processed%20Data) this folder contains the processed data processed from raw data.

4. [Scripts](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts) folder contains the script for:

    * [Preprocessing](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Preprocessing/PreprocessQC.md): processing the raw data
    
    * [Exploratory Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory.md): explores our processed training data to see if there are any obvious underlying structure. 
    
    * [Differential Methylation Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma.md): done using limma on the processed data.
    
    * [Building the classifer](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/BuildModel_AnalyzePredictors.md): This script is for building the ancestry classifier, as well as for the analysis of the resulting predictor CpGs. This folder also contains the script to run the classifier  and analyze those results on the [second dataset](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-016-0054-8), whose genetic ancestry is unknown.
    
    * [Comparing SVM vs glmnet](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/PredictiveModeling.md): This script was used to compare glmnet and SVM. 

    * [Functional Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/FunctionalAnalysis/FunctionalAnalysis.md): for the functional analysis of the CpG sites prioritized by glmnet and limma.

5. [Results](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Results) contain a summary of our main findings. 

    * [Exploratory Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/Exploratory%20Analysis%20Write%20Up.md)
    * [Differential Methylation](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/Differential%20Methylation%20Analysis%20Write%20Up.md)
    * [Predictive Modeling](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/Predictive%20Modeling.md)
    * [Functional Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/Functional%20Analysis%20Write%20Up.md)
    
6. [Poster](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/poster.pdf)