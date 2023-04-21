IP_MS_ANALYSIS
================
Finn
, Last edited on 21 April 2023

- <a href="#myd88" id="toc-myd88">MYD88</a>
- <a href="#irak4" id="toc-irak4">IRAK4</a>
- <a href="#irak1" id="toc-irak1">IRAK1</a>

### MYD88

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.1     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
    ##     yday, year
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
library(janitor)
```

    ## 
    ## Attaching package: 'janitor'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
library(ggrepel)
library(knitr)
library(svglite)
library(readxl)
library(UniprotR)
library(plotly)
```

    ## 
    ## Attaching package: 'plotly'
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

``` r
################################################################################
knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo = FALSE, fig.height=14, fig.width=21)
options(stringsAsFactors = FALSE)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

################################################################################

### DEFINE PATHS
MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/"), 
                 yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/",
                 no  =  "~/Documents/Github/transferGit/")
```

![](IP_MS_ANALYSIS_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### IRAK4

### IRAK1
