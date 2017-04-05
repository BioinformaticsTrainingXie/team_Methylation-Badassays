## The Methylation Badassays

### Introduction:

DNA methylation (DNAm) is linked to many diseases like cancer and autism. However DNAm marks can change due to environmental stimuli, cell types and etc. In recent years, several DNAm studies have suggested that a large portion of DNAm variability is associated with genetic ancestry and is heritable, making DNAm a potential confounding factor which is not given enough consideration in the context of DNA methylation analysis. Differentially methylated CpG sites associated with pathology can ben confounded by CpGs associated with genetic ancestry causing suprious reuslts. Therefore, to invetigate how DNA methylation affects prenatal health, it is important for us to identify genetic ancestry-associated CpGs to figure out true positives. This DNAm variability in the placenta due to genetic ancestry needs to be accounted for in large scale DNAm studies, or else no meaningful interpretation of results can be done to assess prenatal health. In this project, we are going to investigate if DNA in placental tissue is differentially methylated across populations of different ancestry. 

**Hypothesis**: DNA in placental tissue is differentially methylated across populations of different ancestries.

We will first find methylation profiles in subjects from our [dataset 1](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) and the genetic ancestry (Asians or Caucasians) of our subjects is known. These profiles will then serve as a basis to cluster methylation data in our [dataset 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69502) which genetic ancestry is not known. For both dataset 1 and dataset 2, DNAm was measured by 450K microarray from Illumina.

For the details of the project ideas, dataset and methods we used for this project, please check the [project proposal](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/project_proposal.md). 

### Workflow

This figure summaries the workflow of our project:

![workflow](https://cloud.githubusercontent.com/assets/24922214/24690299/deb7dc8c-1980-11e7-9554-ec0ca92f4038.png)


We first used those [codes](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Preprocessing/PreprocessQC.md) to process [raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) of dataset 1 to [processed data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Processed%20Data) and we also performed [quality control](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Preprocessing/PreprocessQC.md). For the details information of dataset 1, please see [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Data/Raw%20Data/samplesheet.csv).

Then, we used a linear model to identify differentially methyalted probes with limma. Please see [Limma](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma.md) in scripts [folder](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts) for details. In [functional analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/FunctionalAnalysis/FunctionalAnalysis.md), we used the COHCAP (City of Hope CpG Island Analysis Pipeline) package, CpGs identified as predictors by glmnet and prioritized by limma were annotated with chromosome mapping, location, gene name and CpG island information.

We also used the processed data to perform [exploratory analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory.md).

We compare [SVM](http://ca.wiley.com/WileyCDA/WileyTitle/productCd-0471030031.html) and [elastic net logistic regression (glmnet)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115) to generate a claaifier that will predict ancestry based on DNA methylation fatures. Please see [predictive modeling](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/PredictiveModeling.Rmd) for details.


### Table of contents:

1. [Project proposal](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/project_proposal.md): project proposal includes the introduction to the ideas, dataset and methods we used in this project.

2. [Progress report](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/progress_report.md): Progress report contains the progress about our project.

3. [Data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data) folder contains metadata, raw data, and processed data.

* [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Data/Raw%20Data/samplesheet.csv)
  + human placental tissue from 45 subjects with self reported ethnicity
  + columns correspond to subject ethnicity, name, sex, gestational age and what complications they had in pregnancy (none, intrauterine growth (IUGR) restriction, or late onset preeclampsia (LOPET), neither of which affect DNAm)
  + columns for Sentrix ID and position correspond to the sample’s batch ID and position on the Illumina microarray 
  + each row is one subject.
  
* [Raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) for dataset 1.

* [Processed data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Processed%20Data) this folder contains the processed data processed from raw data.

4. [Scripts](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts): scripts folder contains the script for:
    * [Preprosessing](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Preprocessing/PreprocessQC.md):This [folder](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts/Preprocessing) contains scripts process raw data to processed data.
    
    * [Exploratory Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory.md): This [folder](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts/ExploratoryAnalysis) contains scripts to identify differentially methylated CpG sites between Caucasians and Asians in placental tissue. 
    
    * [Predictive Modeling](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/PredictiveModeling/PredictiveModeling.md): This [folder](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Scripts/PredictiveModeling) contains scripts to use the identified CpG sites from [Exploratory Analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory.md) to determine the genetic ancestry of [a second dataset](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-016-0054-8) whose genetic ancestry is unknown.

5. [Results](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/results) folder contains the results. 

6. [docs](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/docs)

7. [Poster](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/poster.pdf)