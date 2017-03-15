# Progress report 

### What has changed based on the final proposal 

- **Did your dataset change? If so, why?**

  Went from 485512 to 464923 CpG sites.

  We removed probes with a bad detection p value > 0.01, and a beadcount of < 3.  ‘Bad detection p value > 0.01’ probes are probes that failed to be significantly different than the technical negative control probes on the Array provided by Illummina. Beadcount < 3 probes are those probes that were measured by less than 3 beads and are therefore unreliable measurements of the true methylation quantity. We also normalized Type I and Type II probes. See below for discussion and graphing. No ‘bad’ samples (major outliers) were identified. 


- **Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?**

- **Are there any changes in task assignments of group members?**

  There are no major changes to team members task. Some tasks will be assigned that are more specific than what was initially proposed. Such tasks are outlined below (and subject to change as the project progresses):
  
* Victor - Predictive Modeling Analysis (research on model selection, workflow, and coding)
* Ming - Predictive Modeling Analysis (research on model selection, workflow, and coding)
* Nivi - Exploratory Analysis - Github organization 
* Anni - Exploratory Analysis - Github organization / r markdown annotation
* Michael - Accessory (will help out wherever is needed)

 Everybody will help out on the poster.




### What is the progress of the analyses 

- **Since your initial proposal, you should have decided more concretely on what methods to use for each step of your analyses, and employed some of those methods.**

Because we already discussed the changes to our plan for the analysis in the above section, in brief, here is what our current (still deciding) plan for the methodology:

- Decide on a model-(logistic regression, glmnet, knn, SVM)
     - Leaning towards glmnet (See Horvath et al. 2013 - built an age-predictor on 27k) 

- Will check out caret and glmnet R packages to implement this in R

Please see above sections for our thoughts on the analysis plan.


- **Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?**


- **Can you write down your rationale for the normalization? Refer to the density plots. Can you put a picture of the plots here from the knitted html file? **

  There are a couple different normalization methods used in DNA methylation analysis; however, there isn't a consensus on which method is the best and metrics to evaluate how good normalization methods perform are vague and unclear. So we tried different available normalization methods to see which one works better. The normalization methods we tried are noob, functional normalization and quantile normalization. 

  We first tried noob background subtraction method with dye-bias normalization (the function preprocessNoob in Minfi)  which estimates background noise from the out-of-band probes and remove it for each sample separately, while the dye-bias normalization utilizes a subset of the control probes to estimate the dye bias (1). We then tried functional normalization (the preprocessFunnorm function in Minfi) which uses the internal control probes present on the array to infer between-array technical variation (2). In default first step of preprocessFunnorm function is background subtraction using preprocessNoob and then uses the first two principal components of the control probes to infer the unwanted variation. This functional normalization has been shown particularly useful when global changes are expected. Finally, we used quantile normalization (3) which is applied to the methylation and demethylation intensities separately. In this normalization method, the function preprocessQuantile first normalizes the type II probes across samples and then interpolates a reference distribution to which the type I probes were normalized; therefore, the distribution of type I and type II signals is forced to be the same. The quantile normalization is suggested for datasets where small differences rather than global changes are expected.

  The plots we obtained by applying the three normalization methods are:
  
  
![Noob](https://cloud.githubusercontent.com/assets/24922214/23965730/2b3f9984-0976-11e7-8e82-5268c1f0173c.png)


![funNorm_noob](https://cloud.githubusercontent.com/assets/24922214/23965740/3467b6e0-0976-11e7-8806-c8f9deea0b51.png)


![Quantile](https://cloud.githubusercontent.com/assets/24922214/23965746/38498f36-0976-11e7-8dc7-c840d19c0b9b.png)
  
  A good preprocessing method should make the peaks of type 1 & 2 probe distributions close together, so functional normalization and quantile normalization appears to be better at this task with our dataset. We chose functional normalization, but are unsure whether Quantile would be better, or if the difference is negligible.  

References:

1. Timothy J Triche, Daniel J Weisenberger, David Van Den Berg, Peter W Laird, and Kimberly D Siegmund. 2013. “Low-level processing of Illumina Infinium DNA Methylation BeadArrays.” Nucleic Acids Research 41 (7): e90. doi:10.1093/nar/gkt090.

2. Fortin, Jean-Philippe, and Kasper D. Hansen. 2015. “Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data.” BioRxiv. doi:10.1101/019000.

3. Nizar Touleimat, and Jörg Tost. 2012. “Complete pipeline for Infinium Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation.” Epigenomics 4 (3): 325–41. doi:10.2217/epi.12.21.









- **What R packages or other tools are you using for your analyses?**

  We used mainly the ‘Minfi’ package for preprocessing. Dplyr for data organizing / wrangling, and ‘ggplot2’ for graphing.


- **Provide the links to any markdown reports within your repo to refer to the relevant analysis.**

The link to scripts for processing the raw data to data/ processed data.

[Scripts](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC.Rmd)





**References**

**Minfi**:

Martin J Aryee, Andrew E Jaffe, Hector Corrada Bravo, Christine Ladd-Acosta, Andrew P Feinberg, Kasper D Hansen, and Rafael A Irizarry. 2014. “Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays.” Bioinformatics 30 (10): 1363–69. doi:10.1093/bioinformatics/btu049.

**wateRmelon**:

Pidsley, Ruth, Wong Y, C C, Volta, Manuela, Lunnon, Katie, Mill, Jonathan, Schalkwyk and C L (2013). “A data-driven approach to preprocessing Illumina 450K methylation array data.” BMC Genomics, 14, pp. 293. doi: 10.1186/1471-2164-14-293.

We plan to use **Caret** for the predictive modeling aspects of our analyses:

C. K. Williams, A. Engelhardt, T. Cooper, Z. Mayer, A. Ziem, L. Scrucca, Y. Tang, C. Candan, M. M. Kuhn, "Package ‘caret", 2015.




