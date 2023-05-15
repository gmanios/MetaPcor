# MetaPcor

MetaPcor : A package to conduct meta-analysis of co-expression networks using partial correlation as effect size

Web version of MetaPcor available at : http://rs.dib.uth.gr:3839/metapcor/metapcor2023/




## Installation Guide

First of all, make sure the required packages that are metnioned below are installed. 

Install the required packages

```R
library(devtools)

install.packages(c("readxl", "igraph", "visNetwork",
                   "plotly", "GGally", "ggrepel",
                   "gprofiler2", "matrixcalc", 
                   "ff", "bit", "DescTools",
                   "forcats", "stringr",
                   "purrr", "readr",
                   "tidyr", "tibble","ggplot2", 
                   "tidyverse", "dplyr",
                   "plyr", "metafor", "metadat", 
                   "Matrix", "meta", "GeneNet",
                   "fdrtool", "longitudinal", 
                   "data.table", "corpcor"), 
                 dependencies = TRUE)
                 
                 
# Install geomnet from GitHub
devtools::install_github("sctyner/geomnet")

# Install space from GitHub
devtools::install_github("cran/space")     

# Install Bioconductor packages
BiocManager::install("GEOquery")
BiocManager::install("DExMA")

```
Install MetaPcor from GitHub
```R
devtools::install_github("gmanios/MetaPcor") 
```
## Examples 

```R 
rm(list=ls())  

source("metapcor_functions.R")

#Run MetaPcor 
pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANHT_12_V4'), option=2, method="sparse", meta_method= "random", pvalue_thres = 0.01,l1  = 0.6 ,l2 = 0)

#Run MetaPcor with option 4 (DE Analysis)
#pcor <- meta_pcor(folder_path = 'studies/' , option=4, method="sparse", meta_method= "random",l1  = 0.6 ,l2 = 0)


# Encrichment Analysis with gProfiler
ea_results <- enrichment_analysis(as.data.frame(pcor))

# Manhattan plot
ea_results

# Network plot with GraphViz
network_plot(pcor)

# Volcano Plot
volc_plot<- volc_plot_plotly(as.data.frame(pcor),pval_thres = 0.5,coeff_thres = 0.1)
volc_plot


```

## Requirements

* readxl 1.4.1
* igraph 1.3.5
* geomnet 0.3.1
* visNetwork 2.1.2
* plotly 4.10.1
* GGally 2.1.2
* ggrepel 0.9.2
* gprofiler2 0.2.1
* matrixcalc 1.0-6
* ff 4.0.7
* bit 4.0.4
* DescTools 0.99.47
* forcats 0.5.2
* stringr 1.4.1
* purrr 0.3.5
* readr 2.1.3
* tidyr 1.2.1
* tibble 3.1.8
* ggplot2 3.4.0
* tidyverse 1.3.2
* dplyr 1.0.10
* plyr 1.8.8
* metafor 3.8-1
* metadat 1.2-0
* Matrix 1.5-1
* meta 6.0-0
* GeneNet 1.2.16
* fdrtool 1.2.17
* longitudinal 1.1.13
* data.table 1.14.4
* space 0.1-1.1
* corpcor 1.6.10
* GEOquery 2.66.0
* Biobase 2.58.0
* BiocGenerics 0.44.0

