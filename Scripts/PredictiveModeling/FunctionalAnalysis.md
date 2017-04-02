FunctionalAnalysis
================
Nivretta Thatra

Step 0: Load Packages and Data
==============================

Load required packages:

``` r
#getting cpg site functional info
#source("https://bioconductor.org/biocLite.R")
#biocLite("COHCAP")

library("COHCAP")
```

    ## Loading required package: WriteXLS

    ## Loading required package: COHCAPanno

    ## Loading required package: RColorBrewer

    ## Loading required package: gplots

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(plyr)
```

    ## -------------------------------------------------------------------------

    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)

    ## -------------------------------------------------------------------------

    ## 
    ## Attaching package: 'plyr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize

Load a set of CpG sites of interest

``` r
#I've saved a set of interesting "test" sites in testCpGsites.txt

interestingSites <- read.table("testCpGsites.txt")
```

Step 1: Annotation to chromosome, gene, and CpG Islands
=======================================================

The [COHCAP](https://www.bioconductor.org/packages/devel/bioc/manuals/COHCAP/man/COHCAP.pdf) package has annotations available for 450k-UCSC, 450k-HMM and 27k array probes.

``` r
colnames(interestingSites) <- c("SiteID")

data(COHCAP.450k.HMM)
    
annotated <- join(interestingSites, COHCAP.450k.HMM)
```

    ## Joining by: SiteID

``` r
knitr::kable(head(annotated)) 
```

| SiteID     | Chr |     Loc| Gene   | Island          |
|:-----------|:----|-------:|:-------|:----------------|
| cg13869341 | 1   |   15865| WASH5P | NA              |
| cg14008030 | 1   |   18827| WASH5P | NA              |
| cg03130891 | 1   |   91550| NA     | NA              |
| cg24335620 | 1   |  135252| NA     | NA              |
| cg16162899 | 1   |  449076| NA     | NA              |
| cg24669183 | 1   |  534242| NA     | 1:523025-524193 |
