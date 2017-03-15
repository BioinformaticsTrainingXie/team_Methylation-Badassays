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
     - Leaning towards glmnet (See Horvath et al. 2013 - built an age-predictor on 27k 

- Will check out caret and glmnet R packages to implement this in R

Please see above sections for our thoughts on the analysis plan.



- **Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?**

- **What R packages or other tools are you using for your analyses?**

- **Provide the links to any markdown reports within your repo to refer to the relevant analysis.**










![Noob](https://cloud.githubusercontent.com/assets/24922214/23965730/2b3f9984-0976-11e7-8e82-5268c1f0173c.png)

![funNorm_noob](https://cloud.githubusercontent.com/assets/24922214/23965740/3467b6e0-0976-11e7-8806-c8f9deea0b51.png)

![Quantile](https://cloud.githubusercontent.com/assets/24922214/23965746/38498f36-0976-11e7-8dc7-c840d19c0b9b.png)


### Results 



**References**
