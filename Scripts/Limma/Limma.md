Differential methylation analysis (of Processed Data) with limma
================
Nivretta Thatra

Step 0: Load Packages and Data
==============================

### Installed the necessary dependencies by:

``` r
#install.packages("gplots")

library(gplots)
library(limma)
```

### Load in pre-processed data from dropbox.

``` r
#setwd("../Limma")


#download preprocessed data from this dropbox link directly into Limma folder
#data.txt is in gitignore, be sure not to try to push to github
#download.file("https://www.dropbox.com/s/s4xv0k0vsl0ka0r/data.txt?dl=1", "data.txt")

#using dropbox because google drive links do not work

data <- read.table("data.txt")

design <- read.csv("des.txt", sep="\t", header=TRUE)
```

Step 1: ID differentially methylated sites by Ethnicity only
============================================================

We can use a linear model to identify differentially methylated probes with `limma`, as shown in this older 540 [seminar](http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar08_methylation.html).

Create a design matrix looking at just the effect of ethnicity.

``` r
des.mat <- model.matrix(~Ethnicity, design)
```

Fit model to obtain top differentially methylated CpG sites.

``` r
DMRfit <- lmFit(data, des.mat)
DMRfitEb <- eBayes(DMRfit)

DMR <- topTable(DMRfitEb, coef = "EthnicityCaucasian", number = Inf, p.value = 0.01, adjust.method = "BH")

knitr::kable(head(DMR))  # top hits 
```

|            |       logFC|    AveExpr|          t|  P.Value|  adj.P.Val|         B|
|------------|-----------:|----------:|----------:|--------:|----------:|---------:|
| cg16329197 |   0.4916578|  0.4996451|  11.656642|        0|    0.0e+00|  23.94872|
| cg25025879 |   0.3917075|  0.4535298|  10.947130|        0|    0.0e+00|  21.86189|
| cg05393297 |   0.4004352|  0.6325893|  10.444695|        0|    0.0e+00|  20.34470|
| cg16808927 |  -0.2639234|  0.1481974|  -9.190328|        0|    1.0e-06|  16.41929|
| cg06903451 |   0.1702557|  0.5767780|   9.038860|        0|    1.2e-06|  15.93273|
| cg14581129 |   0.2138548|  0.5049442|   8.955025|        0|    1.3e-06|  15.66234|

So using a cutoff of FDR = 0.01, we identified 106 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry. Now we can make some plots to check these hits.

Step 3: Plotting
================

``` r
#get the top 100 differentially methylated CpG sites & values for all samples
DMR100 <- topTable(DMRfitEb, coef = "EthnicityCaucasian", number = 100)
DMR.CpG <- t(as.matrix(subset(data, rownames(data) %in% rownames(DMR100))))
str(DMR.CpG, max.level = 0)
```

    ##  num [1:45, 1:100] 0.625 0.445 0.669 0.678 0.443 ...
    ##  - attr(*, "dimnames")=List of 2

``` r
limma_top100_ethinicity <- list(rownames(DMR100))

#write_rds(list(rownames(DMR100)), path = "limma_top100_ethinicity.rds")

#rename the samples before plotting
meta <- cbind(design, DMR.CpG)
rownames(DMR.CpG) <- paste(rownames(DMR.CpG), sub(c(rownames(DMR.CpG)), "", meta$Ethnicity), sep = "_")
```

    ## Warning in sub(c(rownames(DMR.CpG)), "", meta$Ethnicity): argument
    ## 'pattern' has length > 1 and only the first element will be used

``` r
#plot heatmap
col <- c(rep("darkgoldenrod1", times = nrow(DMR.CpG)))
col[grepl("Asian", rownames(DMR.CpG))] <- "forestgreen"
op <- par(mai = rep(0.5, 4))
heatmap.2(DMR.CpG, col = redblue(256), RowSideColors = col, density.info = "none", 
    trace = "none", Rowv = TRUE, Colv = TRUE, labCol = FALSE, labRow = FALSE, 
    dendrogram = "row", margins = c(1, 5), main = "Differentially expressed genes 
    by Ethnicity")
legend("topright", c("Caucasian", "Asian"), col = c("darkgoldenrod1", "forestgreen"), 
    pch = 15)
```

![](Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits-1.png)

``` r
par(op)
```

Step 4: Check for interaction effects
-------------------------------------

Create a design matrix looking at the effect of ethnicity when taking into account ethnicity and sex. The Sample\_groups (control, lopet, IUGR) are shown in literature to not affect CpG site methylation, so let's just check for gender.

``` r
des.mat.gender <- model.matrix(~Ethnicity*sex, design)
```

Fit model to obtain top differentially methylated CpG sites.

``` r
DMRfit.gender <- lmFit(data, des.mat.gender)
DMRfitEb.gender <- eBayes(DMRfit.gender)

DMR.gender <- topTable(DMRfitEb.gender, coef = "EthnicityCaucasian", number = Inf, p.value = 0.01, adjust.method = "BH")

knitr::kable(head(DMR.gender))  # top hits 
```

|            |       logFC|    AveExpr|          t|  P.Value|  adj.P.Val|          B|
|------------|-----------:|----------:|----------:|--------:|----------:|----------:|
| cg16329197 |   0.5368513|  0.4996451|   9.546477|        0|  0.0000020|  17.319139|
| cg25025879 |   0.4343004|  0.4535298|   9.205678|        0|  0.0000028|  16.272412|
| cg05393297 |   0.4229273|  0.6325893|   8.211144|        0|  0.0000427|  13.142162|
| cg14581129 |   0.2265901|  0.5049442|   6.992750|        0|  0.0016940|   9.185289|
| cg26513180 |  -0.0294624|  0.0339633|  -6.732052|        0|  0.0025085|   8.327817|
| cg19041462 |   0.1018915|  0.8764931|   6.689685|        0|  0.0025085|   8.188292|

So using a cutoff of FDR = 0.01, we identified just 13 CpG sites that are differentially methylated between Caucasian and Asian genetic ancestry.

Step 5: Plot top 100 Limma hits when accounting for Ethnicity and gender interaction effect
-------------------------------------------------------------------------------------------

``` r
#get the top 100 differentially methylated CpG sites & values for all samples
DMR100.gender <- topTable(DMRfitEb.gender, coef = "EthnicityCaucasian", number = 100)
DMR.CpG.gender <- t(as.matrix(subset(data, rownames(data) %in% rownames(DMR100.gender))))


limma_top100_ethinicity.gender <- list(rownames(DMR100.gender))

#write_rds(list(rownames(DMR100)), path = "limma_top100_ethinicity.rds")

#rename the samples before plotting
meta2 <- cbind(design, DMR.CpG.gender)
rownames(DMR.CpG.gender) <- paste(rownames(DMR.CpG.gender), sub(c(rownames(DMR.CpG.gender)), "", meta2$Ethnicity), sep = "_")
```

    ## Warning in sub(c(rownames(DMR.CpG.gender)), "", meta2$Ethnicity): argument
    ## 'pattern' has length > 1 and only the first element will be used

``` r
#plot heatmap
col <- c(rep("darkgoldenrod1", times = nrow(DMR.CpG.gender)))
col[grepl("Asian", rownames(DMR.CpG.gender))] <- "forestgreen"
op <- par(mai = rep(0.5, 4))
heatmap.2(DMR.CpG.gender, col = redblue(256), RowSideColors = col, density.info = "none", 
    trace = "none", Rowv = TRUE, Colv = TRUE, labCol = FALSE, labRow = FALSE, 
    dendrogram = "row", margins = c(1, 5), main = "Differentially expressed genes by Ethnicity, 
    accounting for Gender")
legend("topright", c("Caucasian", "Asian"), col = c("darkgoldenrod1", "forestgreen"), 
    pch = 15)
```

![](Limma_files/figure-markdown_github/heatmap%20of%20beta%20values%20of%20top%20100%20Limma%20hits%20with%20interaction-1.png)

``` r
par(op)
```

Step 6: Check with GLMnet predictions
-------------------------------------

Hopefully there is overlap between the CpG sites detected by GLMnet and the ones detected in linear regression analysis. Let's check.

``` r
#load in the CpG sites predicted with GLMnet
interestingSites <- read.table("../../Data/Processed Data/predictorsGlmnet.txt")
colnames(interestingSites) <- c("SiteID")

#list of 13 differentially methylated sites
DMR.gender.sites <- data.frame(rownames(DMR.gender))
colnames(DMR.gender.sites) <- c("SiteID")

#find overlaps
overlaps <- merge(interestingSites, DMR.gender.sites)
overlaps
```

    ##       SiteID
    ## 1 cg05393297
    ## 2 cg12011926
    ## 3 cg14581129
    ## 4 cg16329197
    ## 5 cg25025879

There are 5 CpG sites that overlap - which is promising!
