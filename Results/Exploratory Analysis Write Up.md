# Exploratory Analysis of Processed Data

## Step 1: Exploratory Plots

### 1.1 Explore random CpG Sites with plotting

![Random CpG site plots](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Random%20CpG%20Site%20Plots-1.png)

![Random CpG site plots2](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Random%20CpG%20Site%20Plots-2.png)

### 1.2 Sample to Sample Correlations

![plot1](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Sample%20to%20Sample%20Correlations-1.png)

Sample group (Control, LOPET, IUGR) does not seem to affect the clustering of the samples. This fits with what it is known in the literature (control, LOPET, IUGR do not affect DNAm). Most of the samples seem to correlate evenly with each other.

Now let's look at sample to sample correlation, ordered by ethnicity then gender.
 
![plot2](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Sample%20to%20Sample%20Correlations2-1.png)

It's difficult to tell if gender affects affects the clustering of the samples. We need to look at the data further.


## Step 2: Unsupervised clustering


### 2.1: PCA on training data:
 
The plot() function returns a plot of the variances (y-axis) associated with the PCs (x-axis), which is useful to decide how many PCs to retain for further analysis.
summary(pc.train)

![plot3](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/pca-1.png)
 
Next we can see from plotting the first three principal components that our groups (Asian, Caucasian) do not seem to separate. This indicates that the main drivers of the variance in the data is something else.

![plot4](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/plot%20PCs-1.png)
![plot5](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/plot%20PCs-2.png)
![plot6](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/plot%20PCs-3.png)
 
With further plotting, it's not clear that our other variables are driving the variance in the data (sex, gestational age, and sample group).
 
Plotting scatter plots of the top 5 PCs against ethnicity, none of the PCs can clearly separate samples by ethnicity, disappointing. But this suggests we need to move on to other analysis, like differential methylation and glmnet.

![plot7](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Plot%20other%20metadata-1.png)
![plot8](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Plot%20other%20metadata-2.png)
![plot9](https://github.com/STAT540-UBC/team_Methylation-Badassays/blob/master/Scripts/ExploratoryAnalysis/Exploratory_files/figure-markdown_github/Plot%20other%20metadata-3.png)
