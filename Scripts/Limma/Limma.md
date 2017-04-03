Differential methylation analysis (of Processed Data) with limma
================
Nivretta

Step 0: Load Packages and Data
==============================

### Installed the necessary dependencies by:

``` r
#install.packages("gplots")

library(gplots)
library(limma)
```

### Load in preprocessed data from dropbox.

``` r
#setwd("../Limma")

#download preprocessed data from this dropbox link directly into Limma folder
#data.txt is in gitignore, be sure not to try to push to github
#download.file("https://www.dropbox.com/s/s4xv0k0vsl0ka0r/data.txt?dl=1", "data.txt")


data <- read.table("data.txt")

design <- read.csv("des.txt", sep="\t", header=TRUE)
```

Step 1: ID differentially methylated sites
==========================================

We can use a linear model to identify differentially methylated probes with `limma`, as shown in this older 540 [seminar](http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar08_methylation.html).

Create a design matrix.

``` r
des.mat <- model.matrix(~Ethnicity, design)
```

Fit model to obtain top differentially methylated CpG sites.

``` r
DMRfit <- lmFit(data, des.mat)
DMRfitEb <- eBayes(DMRfit)
cutoff <- 0.01
DMR <- topTable(DMRfitEb, coef = "EthnicityCaucasian", number = Inf, p.value = cutoff)
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
#plot heatmap
###LABELS WRONG AT THE MOMENT...
#col <- c(rep("darkgoldenrod1", times = nrow(DMR.CpG)))
#col[grepl("Caucasian", rownames(DMR.CpG))] <- "forestgreen"
#op <- par(mai = rep(0.5, 4))
#heatmap.2(DMR.CpG, col = redblue(256), RowSideColors = col, density.info = "none", 
    #trace = "none", Rowv = TRUE, Colv = TRUE, labCol = FALSE, labRow = FALSE, 
    #dendrogram = "row", margins = c(1, 5))
#legend("topright", c("Caucasian", "Asian"), col = c("darkgoldenrod1", "forestgreen"), 
    #pch = 15)
```
