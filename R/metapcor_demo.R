rm(list=ls())  

source("metapcor_functions.R")



#Run MetaPcor 
pcor <- meta_pcor(GEO_names=c("GSE76427") ,target_namespace = c('ILLUMINA_HUMANHT_12_V4'), option=2, method="sparse", meta_method= "random", pvalue_thres = 0.01,l1  = 0.6 ,l2 = 0)

# Encrichment Analysis with gProfiler
ea_results <- enrichment_analysis(as.data.frame(pcor))

# Manhattan plot
ea_results

# Network plot with GraphViz
network_plot(pcor)

# Volcano Plot
volc_plot<- volc_plot_plotly(as.data.frame(pcor),pval_thres = 0.5,coeff_thres = 0.1)
volc_plot


#Run Meta-Analysis and partial correlation meta-analysis
pcor<- meta_pcor(folder_path = 'studies_breast_cancer_excel/', option=4, method="sparse", meta_method= "random", pvalue_thres = NULL, fdr_thres = NULL, coef_thres = NULL,l1  = 0.8 ,l2 = 0)

write.table(pcor,'breast_cancer_meta_analysis.tsv',sep = '\t')

pcor[pcor$pval<0.05]
# Encrichment Analysis with gProfiler
ea_results <- enrichment_analysis(as.data.frame(pcor))

# Manhattan plot
ea_results

# Network plot with GraphViz
network_plot(pcor)



expr_mat = list()
annotations_df = list()
col_names = c()
for (i in 1:length(list_of_studies)){
  
  # Read the file
  study <- as.data.frame(t(read.table(list_of_studies[i])))
  genes<- names(study)<- study[1,]
  study <- study[-c(1),]

  study<- as.data.frame(lapply(study,as.numeric))
  names(study)<- genes
  
  expr_mat<- append(expr_mat,list(as.data.table(study)))
}


# Import stat sig genes from meta-analysis
 
sig_genes <- read.table('studies2/stat_significant_genes.txt')
sig_genes<- sig_genes$X0


for (i in 1:length(expr_mat)){

  
  expr_mat[[i]]<-  expr_mat[[i]][ ,  sig_genes, with = FALSE] 
  }



#

gc()

folder_path = 'studies3/'
#Give the folder of the path to read the studies
list_of_studies = list.files(path=folder_path, pattern='.txt', all.files=FALSE, full.names=TRUE)
gc()

expr_mat = list()
annotations_df = list()
col_names = c()
for (i in 1:length(list_of_studies)){
  
  # Read the file
  study <- as.data.frame(t(read.table(list_of_studies[i])))
  genes<- names(study)<- study[1,]
  study <- study[-c(1),]
  
  study<- as.data.frame(lapply(study,as.numeric))
  names(study)<- genes
  
  expr_mat<- append(expr_mat,list(as.data.table(study[,1:10000])))
}


gc()

pcor_list5 <- pcor_neighborhood(list_of_files = expr_mat, significant= FALSE , pvalue_thres = 0.001,coef_thres= 0.01,l1= 0.6 , l2 = 0)      
gc()
pcor_test2 <- pcor_neighborhood(list_of_files = expr_mat[2], significant= FALSE , pvalue_thres = 0.001,coef_thres= 0.01,l1= 0.6 , l2 = 0)      
gc()

library(plyr)
# Merge with full join

pcor_final <- join_all(list(as.data.frame(pcor_test1),as.data.frame(pcor_test2)),type = 'full')

pcor_test <- pcor_list5

pcor_test$from <- as.character(pcor_test$from)
pcor_test$to <- as.character(pcor_test$from)

correlations = group_by(pcor_test,from)
meta_cor5  <- my_meta(correlations = as.data.table(pcor_final), method ="random")

pcor_test %>% 
  group_by(from)


write.csv(pcor_test2,"pcor2.csv",row.names = FALSE)


read.table("pcor2.csv",sep = ',',header = T)



