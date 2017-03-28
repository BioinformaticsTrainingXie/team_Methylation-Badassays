Exploratory
================
Nivretta Thatra
March 15, 2017

Installed the necessary dependencies by:

``` r
#source("https://bioconductor.org/biocLite.R")
#biocLite("wateRmelon")
#biocLite('limma')
#biocLite('minfi')
library(wateRmelon)
```

    ## Warning: package 'limma' was built under R version 3.3.2

    ## Warning: package 'matrixStats' was built under R version 3.3.2

    ## Warning: package 'scales' was built under R version 3.3.2

    ## Warning: package 'reshape2' was built under R version 3.3.2

    ## Warning: package 'ggplot2' was built under R version 3.3.2

    ## Warning: package 'GenomicFeatures' was built under R version 3.3.2

    ## Warning: package 'S4Vectors' was built under R version 3.3.2

    ## Warning: package 'GenomeInfoDb' was built under R version 3.3.2

    ## Warning: package 'GenomicRanges' was built under R version 3.3.2

    ## Warning: package 'AnnotationDbi' was built under R version 3.3.2

    ## Warning: package 'minfi' was built under R version 3.3.2

    ## Warning: package 'Biostrings' was built under R version 3.3.2

    ## Warning: package 'foreach' was built under R version 3.3.2

    ## Warning: package 'iterators' was built under R version 3.3.2

    ## Warning: package 'locfit' was built under R version 3.3.2

    ## Warning: package 'lumi' was built under R version 3.3.2

``` r
library(minfi)
library(dplyr)
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.3.2

    ## Warning: package 'tidyr' was built under R version 3.3.2

``` r
library(readr)
```

Load data
---------

Tried to load preprocessed data in attempt to work from it directly rather than having to regenerate it from preprocessing script...didn't work.

``` r
setwd("../processed_data/")

#data import is not working correctly, read_delim is the only function that gets close
#read.delim, read.csv don't work
#read_delim problem: row names are assigned to as first sample values, fine sample not imported
#data <- read_delim("data.txt"," ", escape_double = FALSE, trim_ws = TRUE)

design <- read.csv("des.txt", sep="\t", header=TRUE)

#trying to get metadata and data into one dataframe to do exploratory visualizations, not possible because data is too big
design$Samplename
```

    ##  [1] PM104 PM112 PM114 PM115 PM119 PM120 PM123 PM124 PM130 PM136 PM139
    ## [12] PM142 PM153 PM155 PM158 PM167 PM181 PM20  PM205 PM226 PM227 PM233
    ## [23] PM243 PM249 PM29  PM30  PM4   PM40  PM41  PM44  PM46  PM47  PM52 
    ## [34] PM53  PM54  PM55  PM58  PM66  PM71  PM72  PM74  PM76  PM84  PM9  
    ## [45] PM98 
    ## 45 Levels: PM104 PM112 PM114 PM115 PM119 PM120 PM123 PM124 PM130 ... PM98

``` r
#colnames(data) <- c(as.character(design$Samplename))

#full <- cbind(design, t(data))
```

Plotting
--------

Ideally we'd plot methylation at a random CpG sight and see how it varies with covariates. Below is sample code for how we'd that, if data loading was working.

``` r
#Sample code

#random cpg site
probe_row <- 2000

#get the site name 
probe_name <- colnames(full)[probe_row]

ggplot(full, aes(x = as.factor(time_num), y = full[probe_row])) + 
  geom_boxplot(aes(fill=Ethnicity), show.legend = TRUE) + 
  geom_jitter(width = 0.3) + 
  facet_wrap(~Sample_Group) + 
  xlab("Gestational Age") + 
  ylab("Beta values") +
  ggtitle(paste("Beta values for CpG site", probe_name)) +
  stat_summary(fun.y = mean, geom="point", colour="darkred", size= 3)
```
