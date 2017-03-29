## The Methylation Badsassays: Data

Table of Contents:

1. [Metadata](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/data/Raw%20Data/supplementary%20clinical%20info) 

  + human placental tissue from 45 subjects with self reported ethnicity
  + columns correspond to subject ethnicity, name, sex, gestational age and what complications they had in pregnancy (none, intrauterine growth (IUGR) restriction, or late onset preeclampsia (LOPET), neither of which affect DNAm)
  + columns for Sentrix ID and position correspond to the sample’s batch ID and position on the Illumina microarray 
  + each row is one subject.

2. [Raw data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/data/Raw%20Data/IDATS)

  + each subject has two .idat files 
  + one .idat file contains the methylated intensity profiles for all 450k CpG sites
  + one .idat file contains the unmethylated intensity profiles for all 450K CpG sites
  + Methylation is determined by taking the ratio of the two intensities at each site
  
3. [Processed Data](https://github.com/STAT540-UBC/team_Methylation-Badassays/tree/master/data/processed_data)
