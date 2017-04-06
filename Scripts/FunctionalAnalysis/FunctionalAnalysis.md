Functional analysis of sites prioritized by limma and glmnet
================
Nivretta Thatra

Step 0: Load Packages and Data
==============================

Load required packages:

``` r
#getting cpg site functional info
#source("https://bioconductor.org/biocLite.R")
#biocLite("COHCAP")
#biocLite("GO.db")
#biocLite("mygene")

library(COHCAP)
```

    ## Loading required package: WriteXLS

    ## Warning: package 'WriteXLS' was built under R version 3.3.2

    ## Loading required package: COHCAPanno

    ## Loading required package: RColorBrewer

    ## Loading required package: gplots

    ## Warning: package 'gplots' was built under R version 3.3.2

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

    ## Warning: package 'tidyverse' was built under R version 3.3.2

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr

    ## Warning: package 'ggplot2' was built under R version 3.3.2

    ## Warning: package 'tibble' was built under R version 3.3.3

    ## Warning: package 'tidyr' was built under R version 3.3.2

    ## Warning: package 'readr' was built under R version 3.3.3

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(plyr)
```

    ## Warning: package 'plyr' was built under R version 3.3.2

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

``` r
library(mygene)
```

    ## Loading required package: GenomicFeatures

    ## Warning: package 'GenomicFeatures' was built under R version 3.3.3

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colnames,
    ##     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
    ##     sort, table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Loading required package: S4Vectors

    ## Warning: package 'S4Vectors' was built under R version 3.3.3

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:gplots':
    ## 
    ##     space

    ## The following objects are masked from 'package:base':
    ## 
    ##     colMeans, colSums, expand.grid, rowMeans, rowSums

    ## Loading required package: IRanges

    ## Warning: package 'IRanges' was built under R version 3.3.3

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     desc

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     reduce, simplify

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, regroup, slice

    ## Loading required package: GenomeInfoDb

    ## Warning: package 'GenomeInfoDb' was built under R version 3.3.2

    ## Loading required package: GenomicRanges

    ## Warning: package 'GenomicRanges' was built under R version 3.3.3

    ## Loading required package: AnnotationDbi

    ## Warning: package 'AnnotationDbi' was built under R version 3.3.2

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 
    ## Attaching package: 'mygene'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     query

Load a set of CpG sites of interest

``` r
interestingSites_net <- read.table("CopyOfpredictorsGlmnet.txt")

interestingSites_limma <- read.table("CopyOflimma_pval0.01_ancestry_accountingforGender.txt")
```

Step 1: Glmnet sites' annotation to chromosome, location, gene and CpG Islands
==============================================================================

The [COHCAP](https://www.bioconductor.org/packages/devel/bioc/manuals/COHCAP/man/COHCAP.pdf) (City of Hope CpG Island Analysis Pipeline) package has annotations available for 450k-UCSC, 450k-HMM and 27k array probes. Specifically, these annotations contain which chromosome, location, gene and CpG island each CpG site maps to.

``` r
colnames(interestingSites_net) <- c("SiteID")

data(COHCAP.450k.HMM)
data(COHCAP.450k.UCSC)
annotated_net <- join(interestingSites_net, COHCAP.450k.HMM)
```

    ## Joining by: SiteID

``` r
annotatedUCSC <- join(interestingSites_net, COHCAP.450k.UCSC)
```

    ## Joining by: SiteID

``` r
knitr::kable(annotated_net) 
```

| SiteID     | Chr |        Loc| Gene    | Island                |
|:-----------|:----|----------:|:--------|:----------------------|
| cg24673385 | 1   |   21051448| SH2D5   | NA                    |
| cg15486123 | 1   |  152881815| IVL     | NA                    |
| cg22853943 | 2   |    4703822| NA      | NA                    |
| cg13921903 | 2   |  241290446| NA      | 2:240939049-240939317 |
| cg08704934 | 3   |  194826585| C3orf21 | 3:196307679-196307875 |
| cg14581129 | 12  |   53358946| NA      | 12:51645202-51645774  |
| cg05393297 | 12  |   53359155| NA      | 12:51645202-51645774  |
| cg25025879 | 12  |   53359317| NA      | 12:51645202-51645774  |
| cg16329197 | 12  |   53359506| NA      | 12:51645202-51645774  |
| cg12011926 | 15  |   29037099| NA      | NA                    |
| cg12602405 | 17  |   79456824| NA      | NA                    |

``` r
knitr::kable(annotatedUCSC)
```

| SiteID     | Chr |        Loc| Gene    | Island                   |
|:-----------|:----|----------:|:--------|:-------------------------|
| cg24673385 | 1   |   21051448| SH2D5   | NA                       |
| cg15486123 | 1   |  152881815| IVL     | NA                       |
| cg22853943 | 2   |    4703822| NA      | NA                       |
| cg13921903 | 2   |  241290446| NA      | chr2:241293841-241294408 |
| cg08704934 | 3   |  194826585| C3orf21 | NA                       |
| cg14581129 | 12  |   53358946| NA      | chr12:53359192-53359507  |
| cg05393297 | 12  |   53359155| NA      | chr12:53359192-53359507  |
| cg25025879 | 12  |   53359317| NA      | chr12:53359192-53359507  |
| cg16329197 | 12  |   53359506| NA      | chr12:53359192-53359507  |
| cg12011926 | 15  |   29037099| NA      | chr15:29033878-29034710  |
| cg12602405 | 17  |   79456824| NA      | chr17:79454734-79455823  |

Summary:
--------

The 11 CpG sites are found on chromosomes 1, 2, 3, 12, 15 or 17. The only difference between the two annotations is the CpG islands to which they map.

The CpG sites map to three genes: SH2D5, IVL and C3orf21. Now let's look at which GO terms the three identified genes map to.

Step 2: Glmnet CpG site-gene GO terms
=====================================

``` r
querySH2D5 <- mygene::query("SH2D5", fields='go', species='human')$hits
(SH2D5 <- lapply(querySH2D5, as.list))
```

    ## $`_id`
    ## $`_id`[[1]]
    ## [1] "400745"
    ## 
    ## 
    ## $`_score`
    ## $`_score`[[1]]
    ## [1] 434.2402
    ## 
    ## 
    ## $go
    ## $go$CC
    ## $go$CC[[1]]
    ##   evidence         id                  term
    ## 1      ISS GO:0014069  postsynaptic density
    ## 2      IEA GO:0030054         cell junction
    ## 3      IEA GO:0045211 postsynaptic membrane

``` r
queryIVL <- mygene::query("IVL", fields='go', species='human')$hits
(IVL <- lapply(queryIVL, as.list))
```

    ## $`_id`
    ## $`_id`[[1]]
    ## [1] "3713"
    ## 
    ## 
    ## $`_score`
    ## $`_score`[[1]]
    ## [1] 399.0674
    ## 
    ## 
    ## $go
    ## $go$BP
    ## $go$BP[[1]]
    ##   evidence         id   pubmed
    ## 1      IDA GO:0010224 16639001
    ## 2      IDA GO:0018149 10908733
    ## 3      TAS GO:0018153  1601889
    ## 4      IDA GO:0030216 10908733
    ## 5      TAS GO:0070268       NA
    ##                                                       term
    ## 1                                         response to UV-B
    ## 2                                    peptide cross-linking
    ## 3 isopeptide cross-linking via N6-(L-isoglutamyl)-L-lysine
    ## 4                             keratinocyte differentiation
    ## 5                                            cornification
    ## 
    ## 
    ## $go$CC
    ## $go$CC[[1]]
    ##   evidence         id          pubmed                  term
    ## 1      IDA GO:0001533        10908733    cornified envelope
    ## 2      TAS GO:0001533            NULL    cornified envelope
    ## 3      IDA GO:0005737 42494, 16639001             cytoplasm
    ## 4      IDA GO:0005813            NULL            centrosome
    ## 5      IDA GO:0005829            NULL               cytosol
    ## 6      TAS GO:0005829            NULL               cytosol
    ## 7      IDA GO:0016604            NULL          nuclear body
    ## 8      IDA GO:0070062        19199708 extracellular exosome
    ## 
    ## 
    ## $go$MF
    ## $go$MF[[1]]
    ##   evidence         id          pubmed                         term
    ## 1      IDA GO:0005198 42494, 10908733 structural molecule activity
    ## 2      IPI GO:0005515        21044950              protein binding
    ## 3      IDA GO:0030674        10908733    protein binding, bridging

``` r
queryC3orf21 <- mygene::query("C3orf21", fields='go', species='human')$hits
(C3orf21 <- lapply(queryC3orf21, as.list))
```

    ## $`_id`
    ## $`_id`[[1]]
    ## [1] "152002"
    ## 
    ## 
    ## $`_score`
    ## $`_score`[[1]]
    ## [1] 1.243485
    ## 
    ## 
    ## $go
    ## $go$BP
    ##   evidence         id            pubmed                term
    ## 1      IDA GO:0016266 8982869, 22117070 O-glycan processing
    ## 
    ## $go$CC
    ##   evidence         id   pubmed
    ## 1      IDA GO:0030176 22117070
    ##                                                   term
    ## 1 integral component of endoplasmic reticulum membrane
    ## 
    ## $go$MF
    ## $go$MF[[1]]
    ##   evidence         id            pubmed                            term
    ## 1      IDA GO:0000287           8982869           magnesium ion binding
    ## 2      IDA GO:0030145           8982869           manganese ion binding
    ## 3      IDA GO:0035252 8982869, 22117070 UDP-xylosyltransferase activity

``` r
#queryMany, a function from the mygene package which allows many genes to be queried at once, does not seem to be working
```

Summary:
--------

Gene SH2D5 is involved in postsynaptic density (presumably in neurons), cell junction, and the postsynaptic membrane. Gene IVL is involved in cornfication (hard layer of skin formation), peptide cross-linking, and other protein binding. Gene C3orf21 is a key component of ER membrane, is involved in enzyme transport, and in ion binding.

Step 3: Limma sites' annotation to chromosome, location, gene and CpG Islands
=============================================================================

Repeat step 1 but for sites prioritized by limma.

``` r
colnames(interestingSites_limma) <- c("SiteID")

annotated_limma <- join(interestingSites_limma, COHCAP.450k.HMM)
```

    ## Joining by: SiteID

``` r
knitr::kable(annotated_limma) 
```

| SiteID     | Chr |       Loc| Gene      | Island               |
|:-----------|:----|---------:|:----------|:---------------------|
| cg16329197 | 12  |  53359506| NA        | 12:51645202-51645774 |
| cg25025879 | 12  |  53359317| NA        | 12:51645202-51645774 |
| cg05393297 | 12  |  53359155| NA        | 12:51645202-51645774 |
| cg14581129 | 12  |  53358946| NA        | 12:51645202-51645774 |
| cg26513180 | 16  |  89883248| FANCA     | 16:88409682-88410944 |
| cg19041462 | 8   |  17146201| VPS37A    | NA                   |
| cg12011926 | 15  |  29037099| NA        | NA                   |
| cg05795554 | 16  |  89883212| FANCA     | 16:88409682-88410944 |
| cg16808927 | 16  |    710921| WDR90     | 16:650649-651454     |
| cg06903451 | 1   |   1323192| CCNL2     | 1:1312116-1313099    |
| cg12419862 | 22  |  24373484| LOC391322 | 22:22702803-22703712 |
| cg01172526 | X   |   2826829| ARSD      | NA                   |
| cg10265016 | 2   |  18061319| KCNS3     | NA                   |

Summary:
--------

The 13 CpG sites are found on chromosomes 1, 2, 8, 12, 15, 16, 22 or the X chromosome.

The CpG sites map to seven genes: FANCA, VPS37A, WDR90, CCNL2, LOC391322, ARSD and KCNS3. Now let's look at which GO terms the seven identified genes map to.

Step 4: Limma CpG site-gene GO terms
====================================

Hiding these results because they are quite lengthy. Please see summary below.

``` r
queryFANCA <- mygene::query("FANCA", fields='go', species='human')$hits
SH2D5 <- lapply(queryFANCA, as.list)

queryVPS37A <- mygene::query("VPS37A", fields='go', species='human')$hits
VPS37A <- lapply(queryVPS37A, as.list)

queryWDR90 <- mygene::query("WDR90", fields='go', species='human')$hits
WDR90 <- lapply(queryWDR90, as.list)

queryCCNL2 <- mygene::query("CCNL2", fields='go', species='human')$hits
CCNL2 <- lapply(queryCCNL2, as.list)

queryLOC391322 <- mygene::query("LOC391322", fields='go', species='human')$hits
LOC391322 <- lapply(queryLOC391322, as.list)

queryARSD <- mygene::query("ARSD", fields='go', species='human')$hits
ARSD <- lapply(queryARSD, as.list)

queryKCNS3 <- mygene::query("KCNS3", fields='go', species='human')$hits
KCNS3 <- lapply(queryKCNS3, as.list)
```

Summary:
--------

FANCA is involved in many processes, most notably in DNA repair, gonad development, and inflammation VPs37A is involved in protein transportation and viral life cycles. WDR90 is involved in protein binding. CCNL2 is involved in transcription and regulation of RNA polymerases. LOC391322 did not have any GO terms. ARSD is involved in many process, most notably lipid metabolic processin and ion binding. KCN3 is involved in ion transport and voltage-gated channel regulation.

There does not seem to be any clear patterns in the functions of genes identified by the two methods. Interestingly, limma found one important gene associated with gonad development, suggesting the importance of the interaction effect between gender and ancestry.

Notes:
======

The easiest way to get all annotation information (gene, chromosome, GO ID, GO term, etc) would be to use the package `IlluminaHumanMethylation450k.db` - but I can't get it to download correctly. I get this error: "ERROR: loading failed \* removing ‘/home/nivretta/R/x86\_64-redhat-linux-gnu-library/3.3/IlluminaHumanMethylation450k.db’ The downloaded source packages are in ‘/tmp/RtmpK8whuD/downloaded\_packages’ installation path not writeable, unable to update packages: cluster, lattice, Matrix, mgcv, nlme, survival Warning message: In install.packages(pkgs = doing, lib = lib, ...) : installation of package ‘IlluminaHumanMethylation450k.db’ had non-zero exit status".

And this seems to be [a documented problem](https://support.bioconductor.org/p/62068/), so I gave up and moved on.
