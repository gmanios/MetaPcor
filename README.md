[![GitHub issues](https://img.shields.io/github/issues/gmanios/MetaPcor?color=green)](https://github.com/gmanios/MetaPcor/issues/new)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fgmanios%2FMetaPcor&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
# MetaPcor

MetaPcor : A package to conduct meta-analysis of co-expression networks using partial correlation as effect size

Web version of MetaPcor available at : http://rs.dib.uth.gr:3839/metapcor/metapcor2023/
or http://195.251.108.211:3839/metapcor/metapcor2023/

 

## Installation Guide

 
Install MetaPcor from GitHub
```R
library(devtools)

devtools::install_github("gmanios/MetaPcor") 
```


## Example 1 (Partial correlation meta-analysis with thresholds)

```R 
library(MetaPcor)

#Run MetaPcor with option 2 (Partial correlation meta-analysis with thresholds)
pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=2, 
method="sparse", meta_method= "random", 
pvalue_thres = 0.01,l1  = 0.7 ,l2 = 0)

 
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
## Example 2 (Partial correlation meta-analysis without thresholds)

```R 
library(MetaPcor)

pcor <-  meta_pcor(folder_path = 'demo_files/GSE_DEMO/' , option=3, method="sparse", meta_method= "random",l1  = 0.8 ,l2 = 0)


#Uncoment to run
#pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANREF_8_V3'), option=3,method="sparse", meta_method= "random",l1  = 0.7 ,l2 = 0)

 
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
## Example 3 (Differential Expression Analysis and Partial correlation meta-analysis)

```R
library(MetaPcor)

#Run MetaPcor with option 4 (DE Analysis and Partial correlation meta-analysis)
 
pcor <- meta_pcor(folder_path = 'demo_files/DE/' , option=4, method="sparse", meta_method= "random",l1  = 0.6 ,l2 = 0)
 
# Encrichment Analysis with gProfiler
ea_results <- enrichment_analysis(as.data.frame(pcor))

# Manhattan plot
ea_results

# Network plot with GraphViz
network_plot(pcor)
 
# Volcano Plot
volc_plot<- volc_plot_plotly(as.data.frame(pcor),pval_thres = 0.01,coeff_thres = 0.2)
volc_plot
```
