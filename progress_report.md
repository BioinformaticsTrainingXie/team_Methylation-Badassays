Progress Report
================
Team Methylation Badassays
March 15, 2017

Progress report (7 pts)
-----------------------

### Section 1: What has changed based on the final proposal (1 pt.)

##### S.1.1: Did your dataset change? If so, why?

Went from 485512 to 464923 CpG sites.

We removed probes with a bad detection p value &gt; 0.01, and a beadcount of &lt; 3. ‘Bad detection p value &gt; 0.01’ probes are probes that failed to be significantly different than the technical negative control probes on the Array provided by Illummina. Beadcount &lt; 3 probes are those probes that were measured by less than 3 beads and are therefore unreliable measurements of the true methylation quantity. We also normalized Type I and Type II probes. See below for discussion on how we chose which normalization method to use. No ‘bad’ samples (major outliers) were identified.

##### S.1.2: Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?

One way we can think of, in order to identify CpG sites that are helpful in predicting ethnicity, is using classification methods like logistic regression:

1.  Training step:

    1.  Fit logistic regression with Prob(ethnicity ="Asian") as the response and all CpG sites as predictors, use regularization techniques (LASSO/elastic net) to choose which CpG sites should be kept. We plan to try glmnet model from "caret" package (there are other functions we can try as well)

    2.  Cross-validate our model (5/10 fold CV);

    3.  Alternatively, if time permits, we can try other classification methods like KNN/SVM, repeat training-CV procedure and compare performance across different models

2.  Apply our trained classifier to the test data set, although we cannot fully "test" our results; However, we propose to validate these results with the following methods:

3.  To double-check our findings, we can experiment with unsupervised clustering methods like PCA:

    1.  Training set only: use PCA to visualize which of the PCs are good classifiers of ethnicity given the labels we have, by examining correlations between PC vs. ethnicity;

    2.  Extract these PC loadings (“eigen-CpG sites”) and apply to the testing data. As a verification to results from step \#2, we expect our labels from step 2 to form distinguishable clusters in the PC plots;

4.  If time permits, we can look for differentially methylated sites between healthy placenta versus those with neural tube defects in our test data, adjusted for ethnicity.

##### S.1.3: Are there any changes in task assignments of group members?

There are no major changes to team members task. Some tasks will be assigned that are more specific than what was initially proposed. Such tasks are outlined below (and subject to change as the project progresses):

Victor - Predictive Modeling Analysis (research on model selection, workflow, and coding)

Ming - Predictive Modeling Analysis (research on model selection, workflow, and coding)

Nivi - Exploratory Analysis - Github organization

Anni - Exploratory Analysis - Github organization / r markdown annotation

Michael - Accessory (will help out wherever is needed)

Everybody will help out on the poster.

### Section 2: [What is the progress of the analyses](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC.md) (4 pts.)

##### S.2.1: Since your initial proposal, you should have decided more concretely on what methods to use for each step of your analyses, and employed some of those methods.

Because we already discussed the changes to our plan for the analysis in the above section, in brief, here is what our current (still deciding) plan for the methodology:

-   Decide on a model - (logistic regression, glmnet, knn, SVM)

    -   Leaning towards glmnet (See Horvath et al. 2013 - built an age-predictor on 27k
-   Will check out caret and glmnet R packages to implement this in R

Please see above sections for our thoughts on the analysis plan.

##### S.2.2: Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?

In summary, we have preprocessed our data (took longer than we anticipated!) into a format usable for analysis. Now, what remains is building the classifier and predicting on our testing set. We have progressed in these aspects mainly on the theoretical level - we have discussed it and are narrowing down the workflow and important decisions to be made (model selection). We have also identified a few R packages that we can use (caret, glmnet).

We also are doing some preliminary exploratory analysis to see the global heterogeneity and see if our samples are relatively represented on the additional covariates (gestational age, sex).

##### S.2.3: Can you write down your rationale for the normalization? Refer to the density plots. Can you put a picture of the plots here from the knitted html file?

It is important to normalize 450k data because of the technical biases introduced by the two probe types that populate these arrays \[4\]. The two types of probes found on the 450k array are classified as ‘Type I’ and ‘Type II’. Due to the inherent differences of the technology behind these designs, each class has a particular bias towards measurement of methylation. If we plot the beta value (measurement of DNAm) of all type I (infinium I) and type II (Infinium II) probes (figure 1), we can see that their distributions are dissimilar.

Figure 1)

![Raw](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC_files/figure-html/unnamed-chunk-12-1.png)

There are a couple different normalization methods used in DNA methylation analysis; however, there isn't a consensus on which method is the best and metrics to evaluate how good normalization methods perform are vague and unclear. So we tried different available normalization methods to see which one works better. The normalization methods we tried are noob, functional normalization and quantile normalization.

We first tried noob background subtraction method with dye-bias normalization (the function preprocessNoob in Minfi (Fig. 2a) which estimates background noise from the out-of-band probes and remove it for each sample separately, while the dye-bias normalization utilizes a subset of the control probes to estimate the dye bias (1). We then tried functional normalization (the preprocessFunnorm function in Minfi) (Fig. 2b) which uses the internal control probes present on the array to infer between-array technical variation (2). In default first step of preprocessFunnorm function is background subtraction using preprocessNoob and then uses the first two principal components of the control probes to infer the unwanted variation. This functional normalization has been shown particularly useful when global changes are expected. Finally, we used quantile normalization (3) which is applied to the methylation and demethylation intensities separately. In this normalization method, the function preprocessQuantile (Fig. 2c) first normalizes the type II probes across samples and then interpolates a reference distribution to which the type I probes were normalized; therefore, the distribution of type I and type II signals is forced to be the same. The quantile normalization is suggested for datasets where small differences rather than global changes are expected.

The plots we obtained by applying the three normalization methods are:

Figure 2a)

![Noob](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC_files/figure-html/unnamed-chunk-13-1.png)

Figure 2b)

![funNorm\_noob](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC_files/figure-html/unnamed-chunk-13-2.png)

Figure 2c)

![Quantile](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC_files/figure-html/unnamed-chunk-13-3.png)

A good preprocessing method should make the peaks of type 1 & 2 probe distributions close together, so functional normalization and quantile normalization appears to be better at this task with our dataset. We chose functional normalization, but are unsure whether Quantile would be better, or if the difference is negligible.

##### S.2.4: What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.

We used mainly the ‘Minfi’ package for preprocessing. Dplyr for data organizing / wrangling, and ‘ggplot2’ for graphing.

##### S.2.5: Provide the links to any markdown reports within your repo to refer to the relevant analysis.

See [preprocessing work](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC.md).

See [exploratory analysis attempt](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/Exploratory.md).

### References:

1.  Timothy J Triche, Daniel J Weisenberger, David Van Den Berg, Peter W Laird, and Kimberly D Siegmund. 2013. “Low-level processing of Illumina Infinium DNA Methylation BeadArrays.” Nucleic Acids Research 41 (7): e90. <doi:10.1093/nar/gkt090>.

2.  Fortin, Jean-Philippe, and Kasper D. Hansen. 2015. “Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data.” BioRxiv. <doi:10.1101/019000>.

3.  Nizar Touleimat, and Jörg Tost. 2012. “Complete pipeline for Infinium Human Methylation 450K BeadChip data processing using subset quantile normalization for accurate DNA methylation estimation.” Epigenomics 4 (3): 325–41. <doi:10.2217/epi.12.21>.

4.  Wu MC, Joubert BR, Kuan P, et al. A systematic assessment of normalization approaches for the Infinium 450K methylation platform. Epigenetics. 2014;9(2):318-329. <doi:10.4161/epi.27119>.

##### References: R packages

1.  Martin J Aryee, Andrew E Jaffe, Hector Corrada Bravo, Christine Ladd-Acosta, Andrew P Feinberg, Kasper D Hansen, and Rafael A Irizarry. 2014. “Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays.” Bioinformatics 30 (10): 1363–69. <doi:10.1093/bioinformatics/btu049>.

2.  Pidsley, Ruth, Wong Y, C C, Volta, Manuela, Lunnon, Katie, Mill, Jonathan, Schalkwyk and C L (2013). “A data-driven approach to preprocessing Illumina 450K methylation array data.” BMC Genomics, 14, pp. 293. doi: 10.1186/1471-2164-14-293.

3.  C. K. Williams, A. Engelhardt, T. Cooper, Z. Mayer, A. Ziem, L. Scrucca, Y. Tang, C. Candan, M. M. Kuhn, "Package ‘caret", 2015.

### Section 3: Results (2 pts.)

Unfortunately we were unable to produce significant results at this point in time. The preprocessing took longer than what was expected. Please see our [preprocessing .md file](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/PreprocessQC.md) to see what work we've done. We realize we probably shouldn't have spent so much time on the preprocessing but as Meg mentioned in the lecture that the preprocessing was very important, we focused a lot of our attention here. We have, however, thought a lot about how we are going to do the analysis, and it is becoming a lot clearer where to go from here. See section 1 for a brief outline of the plan (thanks @rbalshaw @farnushfarhadi @singha53 for the guidance!). It's still subject to change but the skeleton is there. 

<<<<<<< HEAD
##### Were you able to answer your hypothesis?

##### Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results. - List some challenges that you encountered? How will you address them?
=======
##### What are your primary results? We tried to do some preliminary exploratory analysis: [Here is an md describing snafus](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/Exploratory.md) and the start of an [Rmd](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/data/Scripts/Exploratory2.Rmd) working off of the pre processing script.
