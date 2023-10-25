# MetaPcor

MetaPcor : A package to conduct meta-analysis of co-expression networks using partial correlation as effect size

Web version of MetaPcor available at : http://rs.dib.uth.gr:3839/metapcor/metapcor2023/




## Web Version
The web version of MetaPcor is available [here](http://rs.dib.uth.gr:3839/metapcor/metapcor2023/).

 

## Installation Guide

 
Install MetaPcor from GitHub
```R
library(devtools)

devtools::install_github("gmanios/MetaPcor") 
```


## Example 1 (Partial correlation meta-analysis with thresholds)

```R 

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

pcor <-  meta_pcor(folder_path = 'studies/' , option=3, method="sparse", meta_method= "random",l1  = 0.8 ,l2 = 0)


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

#Run MetaPcor with option 4 (DE Analysis and Partial correlation meta-analysis)
 
pcor <- meta_pcor(folder_path = 'studies/' , option=4, method="sparse", meta_method= "random",l1  = 0.6 ,l2 = 0)


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
