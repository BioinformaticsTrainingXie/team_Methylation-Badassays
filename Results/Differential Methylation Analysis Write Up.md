# Differential methylation analysis (of Processed Data) with limma

Here we used a linear model to identify differentially methylated probes with limma, as shown in this older 540 [seminar](http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar08_methylation.html).

We first fit a linear model with ethnicity as the only covariate to obtain top differentially methylated CpG sites (testing the null hypothesis that samples in the two ethnic groups are drawn from populations with the same mean CpG site methylation). Using a cutoff of p value = 0.01, we identified 106 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry. 

##Plot of top 1000 Limma hit with ancestry as only covariate

![plot](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits-1.png)

However, DNA methylation is known to associate with gender, so the interaction effect ethnicity and gender was accounted for in another linear model. Using a cutoff of FDR = 0.01, we identified just 13 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry, when the interaction effect of ethnicity and gender was accounted for. Here are 6 of those sites:
  
  logFC	AveExpr	t	P.Value	adj.P.Val	B
cg16329197	0.5368513	0.4996451	9.546477	0	0.0000020	17.319139
cg25025879	0.4343004	0.4535298	9.205678	0	0.0000028	16.272412
cg05393297	0.4229273	0.6325893	8.211144	0	0.0000427	13.142162
cg14581129	0.2265901	0.5049442	6.992750	0	0.0016940	9.185289
cg26513180	-0.0294624	0.0339633	-6.732052	0	0.0025085	8.327817
cg19041462	0.1018915	0.8764931	6.689685	0	0.0025085	8.188292

## Plot of top 100 Limma hits when accounting for Ethnicity and gender interaction effect

![plot](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits-1.png)

## Check overlap with GLMnet predictions

There is an overlap of 5 CpG sites (cg05393297, cg12011926, cg14581129, cg16329197
and cg25025879) between the those detected by GLMnet and the ones detected in linear regression analysis. 

###Plots of 3 of these CpG sites:

![x](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite1.png)

![y](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite2.png)

![z](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/importantsite3.png)
  
  
###Plots of random CpG site: 
  
![i](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite1.png)
![j](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite2.png)
![k](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/Limma/Limma_files/figure-markdown_github/randosite3.png)
  
  


