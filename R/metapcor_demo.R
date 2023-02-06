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




