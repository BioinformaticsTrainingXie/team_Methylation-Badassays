# Differential methylation analysis (of Processed Data) with limma

Here we used a linear model to identify differentially methylated probes with limma, as shown in this older 540 [seminar](http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar08_methylation.html).

We first fit a linear model with ethnicity as the only covariate to obtain top differentially methylated CpG sites (testing the null hypothesis that samples in the two ethnic groups are drawn from populations with the same mean CpG site methylation). Using a cutoff of p value = 0.01, we identified 106 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry. 

## Plot of top 1000 Limma hit with ancestry as only covariate

![plot](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits-1.png)

However, DNA methylation is known to associate with gender, so the interaction effect ethnicity and gender was accounted for in another linear model. Using a cutoff of FDR = 0.01, we identified just [13 CpG sites](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/limma_pval0.01_ancestry_accountingforGender.txt) that are differentially methylated between Caucasian and Asian genetic ancestry, when the interaction effect of ethnicity and gender was accounted for. 

## Plot of [top 100 Limma](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/limma_top100_ethnicity_accountingforGender.txt) hits when accounting for Ethnicity and gender interaction effect

![plot](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits-1.png)

## Check overlap with GLMnet predictions

There is an overlap of 5 CpG sites (cg05393297, cg12011926, cg14581129, cg16329197
and cg25025879) between the those detected by [GLMnet](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/predictorsGlmnet.txt) and the [ones detected in linear regression analysis](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Results/limma_pval0.01_ancestry_accountingforGender.txt). 

### Plots of 3 of these CpG sites prioritized by both limma and glmnet:

![x](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite1.png)

![y](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite2.png)

![z](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite3.png)
  
  
### Plots of random CpG site: 
  
![i](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite1.png)
![j](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite2.png)
![k](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite3.png)
  
  


