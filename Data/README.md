### The Methylation Badsassays: Data

This folder contains raw data, processed data and metadata for the training dataset. The data for the unlabelled (w/ ethnicity) placental DNA methylation samples can be downloaded, since it is too large to keep on github.
  
* [Raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Raw%20Data) for dataset 1, training data.
  + .idat formats

* [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Data/Raw%20Data/samplesheet.csv)
  + human placental tissue from 45 subjects with self reported ethnicity
  + columns correspond to subject ethnicity, name, sex, gestational age and what complications they had in pregnancy (none, intrauterine growth (IUGR) restriction, or late onset preeclampsia (LOPET), neither of which affect DNAm)
  + columns for Sentrix ID and position correspond to the sample’s batch ID and position on the Illumina microarray 
  + each row is one subject.
  + use this sample sheet to read in .idats
  
* [Processed data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/Data/Processed%20Data) this folder contains the resulting data after processing from raw data.

* Information about dataset 2 can be found in this [link](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-016-0054-8), and can be downloaded [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69502). There are 52 placental tissue DNA methylation (450k) samples that we use for our analysis.
