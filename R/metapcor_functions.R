library(GEOquery)
library(corpcor)
library(space)
library(data.table)
library(GeneNet)
library(meta)
library(metafor)
library(dplyr)
library(tidyverse)
library(purrr)
library(DescTools)
library(ggplot2)
library(ff)
library(matrixcalc)
library(gprofiler2)
library(visNetwork)
library(geomnet)
library(igraph)
library(plotly)
library(readxl)
library(DExMA)
memory.limit(99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999)


# User choose folder path of expression data or geo names for the meta- analysis
readfiles <- function(path){
  x<- read.table(path,sep ='\t',header = TRUE)
  
  study_t <-data.frame(x)
  # study_t <- study_t[1:100,]
  study_t <- study_t[-1,]
  
  study_t <- t(study_t)
  colnames(study_t) <- study_t[1,]
  study_t<- study_t[-1,]
  study_t<- as.data.table(study_t)[, lapply(.SD, as.numeric)]
  
  
  
  # # Convert each column from character to numeric
  # xt <- xt[, lapply(.SD, as.numeric)]
  
  
  return (study_t)
  
}

# User choose folder path of expression data or geo names for the meta- analysis
load_GSE_data <- function(folder_path=NULL, GEO_names=NULL, GPL_list= NULL, target_namespace = NULL){   
  if (!is.null(folder_path) & is.null(GEO_names)){
    
    file_names<-dir(path = folder_path)
    setwd(folder_path)
    list_of_files<-lapply(file_names, readfiles)
    
  }
  else if(!is.null(GEO_names) & is.null(folder_path)) {
    list_of_files<-lapply(GEO_names, getGEO)
    for(i in 1:length(list_of_files)){
      
      if(is.null(list_of_files[[i]][[paste(GEO_names[i],"_series_matrix.txt.gz", sep = "")]]@featureData@data[["Gene Symbol"]])==FALSE){
        namecol<- list_of_files[[i]][[paste(GEO_names[i],"_series_matrix.txt.gz", sep = "")]]@featureData@data[["Gene Symbol"]]
        list_of_files[[i]]<-data.table(t(list_of_files[[i]][[paste(GEO_names[i],"_series_matrix.txt.gz", sep = "")]]@assayData[["exprs"]]))
        colnames(list_of_files[[i]])<-namecol
        list_of_files[[i]]<-list_of_files[[i]][,which(colnames(list_of_files[[i]])!=""), with = FALSE]
        list_of_files[[i]]<-list_of_files[[i]][,which(duplicated(colnames(list_of_files[[i]]))==FALSE), with = FALSE]
      }
      else{
        print("Converting Probes to gene IDs with gConvert")
        list_of_files[[i]]<-data.table(t(list_of_files[[i]][[paste(GEO_names[i],"_series_matrix.txt.gz", sep = "")]]@assayData[["exprs"]]))
        
        # list_of_files[[i]] <- probe2geneID(list_of_files = list_of_files[[i]],GPL_list = GPL_list[[i]])
        # list_of_files[[i]] <- Filter(NROW, gene_annot_list)
        probe_list <- names(list_of_files[[i]])
        
        name_space = target_namespace[i]
        new_cols <- gconvert(probe_list, organism = 'hsapiens', target = name_space)
        
        
        #returns probes that g convert found
        probes_found <- intersect(probe_list,new_cols$input) 
        geo_df = data.frame(list_of_files[[i]])
        
        # Keep the found probes in the initial dataset
        geo_df_cropped = geo_df[,probes_found]
        
        # Keep the common probes in the gconvert output
        new_cols <- new_cols[new_cols$input %in% probes_found,]
        
        
        
        new_cols <- new_cols %>% distinct(input, .keep_all=TRUE)
        
        names(geo_df_cropped) <- new_cols$name
        list_of_files[[i]] <- as.data.table(geo_df_cropped) 
        list_of_files[[i]]<-list_of_files[[i]][,which(duplicated(colnames(list_of_files[[i]]))==FALSE), with = FALSE]
      }
      #list_of_files[[i]] <- list_of_files[[i]][,1:10000]
    }
  }
  else{
    stop("Error: choose files")
  }
  return(list_of_files)
}

# Calculate partial correlation with shrinkage method for each study
pcor_shrinkage <-function (list_of_files, significant=FALSE, pvalue_thres=NULL, fdr_thres=NULL, coef_thres=NULL){ 
  p_cor_list<-list()
  for (i in 1:length(list_of_files)){
    pcor1<-pcor.shrink(list_of_files[[i]])
    pcor1<-pcor1[1:nrow(pcor1),1:nrow(pcor1)]
    
    pcor1<-adjmatrix_to_edgelist(pcor1)
    
    #pcor1<-pcor1[which(pcor1$value!=1),]
    #names(pcor1) <- c("from","to", "ri")
    pcor1$ni <-nrow(list_of_files[[i]])
    
    #keep only the significant correlations for the meta-analysis
    
    if(significant==TRUE){
      
      degrees<- pcor1$ni-3
      test<- pcor1$ri*sqrt((degrees)/(1-(pcor1$ri^2)))
      pcor1$pvalue<-2*pt(q=abs(test), df=degrees, lower.tail=FALSE)
      
      #cutoff with pvalue threshold
      if(is.null(pvalue_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(fdr_thres)==TRUE){
        
        pcor1<-pcor1[which(pcor1$pvalue<pvalue_thres),]
      }
      #cutoff with pvalue threshold and partial correlation coefficient threshold
      else if(is.null(pvalue_thres)!=TRUE & is.null(coef_thres)!=TRUE & is.null(fdr_thres)==TRUE){
        
        pcor1<-pcor1[which(pcor1$pvalue<pvalue_thres & abs(pcor1$ri) > coef_thres ), ]
      }
      
      #cutoff with fdr threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(pvalue_thres)==TRUE){
        
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres),]
      }
      #cutoff with fdr threshold and partial correlation coefficient threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)!=TRUE & is.null(pvalue_thres)==TRUE){
        
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & abs(pcor1$ri) > coef_thres),]
      }
      #cutoff with partial correlation coefficient threshold
      else if(is.null(fdr_thres)==TRUE & is.null(coef_thres)!=TRUE & is.null(pvalue_thres)==TRUE){
        pcor1<-pcor1[which(abs(pcor1$ri) > coef_thres),]
      }
      #cutoff with fdr threshold and pvalue threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(pvalue_thres)!=TRUE){
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & pcor1$pvalue<pvalue_thres),]
      }
      #cutoff with fdr threshold, pvalue threshold and partial correlation coefficient threshold
      else{
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & pcor1$pvalue<pvalue_thres & abs(pcor1$ri) > coef_thres),]
      }
      pcor1<-pcor1[,c("from","to","ri","ni")]
    }
    
    p_cor_list <- append(p_cor_list,list(pcor1))
    
  }
  # Merge the studies
  correlations = p_cor_list %>% purrr::reduce(full_join, by = c("from","to"))
  return(correlations)
}  

# Calculate sparse partial correlation estimation with neighborhood selection approach method for each study
pcor_neighborhood<-function(list_of_files, l1, l2=0, significant=TRUE, pvalue_thres=NULL, fdr_thres=NULL, coef_thres=NULL){
  pcor_list<-list()
  
  for(i in 1:length(list_of_files)){
    
    pcor1<-space.neighbor(list_of_files[[i]], lam1=l1, lam2=l2)$ParCor
    
    colnames(pcor1)<-rownames(pcor1)<- colnames(list_of_files[[i]])
    #pcor1<-melt(pcor1)
    pcor1<-adjmatrix_to_edgelist(pcor1)
    
    #pcor1<-pcor1[which(pcor1$value!=1),]
    #pcor1<-pcor1[which(pcor1$value!=0),]
    #names(pcor1)<-c("from","to", "ri")
    pcor1$ni <-nrow(list_of_files[[i]])
    
    
    #keep only the significant correlations for the meta-analysis
    
    if(significant==TRUE){
      degrees<- pcor1$ni-3
      test<- pcor1$ri*sqrt((degrees)/(1-(pcor1$ri^2)))
      pcor1$pvalue<-2*pt(q=abs(test), df=degrees, lower.tail=FALSE)
      
      #cutoff with pvalue threshold
      if(is.null(pvalue_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(fdr_thres)==TRUE){
        
        pcor1<-pcor1[which(pcor1$pvalue<pvalue_thres),]
      }
      #cutoff with pvalue threshold and partial correlation coefficient threshold
      else if(is.null(pvalue_thres)!=TRUE & is.null(coef_thres)!=TRUE & is.null(fdr_thres)==TRUE){
        
        pcor1<-pcor1[which(pcor1$pvalue<pvalue_thres & abs(pcor1$ri) > coef_thres ), ]
      }
      
      #cutoff with fdr threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(pvalue_thres)==TRUE){
        
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres),]
      }
      #cutoff with fdr threshold and partial correlation coefficient threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)!=TRUE & is.null(pvalue_thres)==TRUE){
        
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & abs(pcor1$ri) > coef_thres),]
      }
      #cutoff with partial correlation coefficient threshold
      else if(is.null(fdr_thres)==TRUE & is.null(coef_thres)!=TRUE & is.null(pvalue_thres)==TRUE){
        pcor1<-pcor1[which(abs(pcor1$ri) > coef_thres),]
      }
      #cutoff with fdr threshold and pvalue threshold
      else if(is.null(fdr_thres)!=TRUE & is.null(coef_thres)==TRUE & is.null(pvalue_thres)!=TRUE){
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & pcor1$pvalue<pvalue_thres),]
      }
      #cutoff with fdr threshold, pvalue threshold and partial correlation coefficient threshold
      else{
        pcor1$fdrs<-p.adjust(p = pcor1$pvalue,method = 'fdr')
        pcor1<-pcor1[which(pcor1$fdrs<fdr_thres & pcor1$pvalue<pvalue_thres & abs(pcor1$ri) > coef_thres),]
      }
      pcor1<-pcor1[,c("from","to","ri","ni")]
    }
    
    pcor_list = append(pcor_list,list(pcor1))
  }
  # Merge the studies
  correlations = pcor_list %>% purrr::reduce(full_join, by = c("from","to"))
  return(correlations)
}

my_meta<-function (correlations, method= "random"){
  fisherZ<-correlations
  studies<-ncol(fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ])
  
  no_meta<-fisherZ[which(rowSums(!is.na(fisherZ))==4),]
  fisherZ<-fisherZ[which(rowSums(!is.na(fisherZ))>4),]
  
  for (i in 1:studies){
    no_meta[[paste0("zi",i)]]<-FisherZ(no_meta[, grepl( "ri" , names( no_meta) ) , with = FALSE][[i]])
    fisherZ[[paste0("zi",i)]]<-FisherZ(fisherZ[, grepl( "ri" , names( fisherZ) ) , with = FALSE][[i]])
  }
  
  no_meta<-no_meta[, !grepl( "ri" , names( no_meta)), with = FALSE]
  fisherZ<-fisherZ[, !grepl( "ri" , names( fisherZ)), with = FALSE]
  
  for (i in 1:studies){
    up<-(((no_meta[, grepl( "ni" , names( no_meta) ), with = FALSE ][[i]])-3) * no_meta[, grepl( "zi" , names( no_meta) ) , with = FALSE][[i]]) #/ ((no_meta[, grepl( "n" , names( no_meta) ) ][i])-3)
    down<- ((no_meta[, grepl( "ni" , names( no_meta) ) , with = FALSE][[i]])-3)
    st<-which(is.na(data.table(down))==FALSE)
    summ<-up/down
    no_meta[st,"Z_tot"]<-summ[st]
    no_meta[st,"SE_fe"]<-sqrt(1/ down)[st]
  }
  
  no_meta$z_score<-no_meta$Z_tot/no_meta$SE
  no_meta$pval<-2*pnorm(q=abs(no_meta$z_score),lower.tail=FALSE)
  no_meta$cor<-FisherZInv(no_meta$Z_tot)
  
  no_meta<-no_meta[,c("from","to","cor","pval")]
  
  
  up<-0
  down<-0
  for (i in 1:studies){
    up<-up+replace((((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3) * fisherZ[, grepl( "zi" , names( fisherZ) ) , with = FALSE][[i]]), is.na((((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3) * fisherZ[, grepl( "zi" , names( fisherZ) ) , with = FALSE][[i]])), 0)
    down<- down+ replace(((fisherZ[, grepl( "ni" , names( fisherZ) ) , with = FALSE][[i]])-3), is.na(((fisherZ[, grepl( "ni" , names( fisherZ) ) , with = FALSE][[i]])-3)),0)
  }
  
  fisherZ$Z_tot<-up/down
  fisherZ$SE_fe<-sqrt(1/ down)
  fisherZ$z_score_fe<-fisherZ$Z_tot/fisherZ$SE_fe
  
  if(method=="fixed"){
    fisherZ$pval<-2*pnorm(q=abs(fisherZ$z_score_fe),lower.tail=FALSE)
    fisherZ$cor<-FisherZInv(fisherZ$Z_tot)
    fisherZ<-fisherZ[,c("from","to","cor","pval")]
    fisherZ<-full_join(fisherZ, no_meta)
    fisherZ$fdrs<-p.adjust(p = fisherZ$pval,method = 'fdr')
    return(fisherZ)
  }
  
  else if(method=="random"){
    
    q<-0
    sw<-0
    sw2<-0
    for (i in 1:studies){
      q_stat<-((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3)* ((fisherZ[, grepl( "zi" , names( fisherZ) ) , with = FALSE][[i]] - fisherZ$Z_tot)^2)
      q_stat <- replace(q_stat, is.na(q_stat), 0)
      q<-q+q_stat
      sw<-sw+replace((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3, is.na((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3),0)
      sw2<-sw2+ replace((((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3)^2), is.na((((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3)^2)),0)
    }
    
    c<- sw - (sw2/sw)
    t2<- (q-(studies-1))/c
    t2<-replace(t2, t2<0 , 0)
    up<-0
    down<-0
    for (i in 1:studies){
      down<- down+ replace(((1/ ((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3))+t2)^(-1), is.na(((1/ ((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3))+t2)^(-1)),0)
      up<-up+replace((((1/ ((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3))+t2)^(-1)) *fisherZ[, grepl( "zi" , names( fisherZ) ) , with = FALSE][[i]], is.na((((1/ ((fisherZ[, grepl( "ni" , names( fisherZ) ), with = FALSE ][[i]])-3))+t2)^(-1)) *fisherZ[, grepl( "zi" , names( fisherZ) ) , with = FALSE][[i]]),0)
      
    }
    
    negative<-which(data.table(down)<0)
    positive<-which(data.table(down)>0)
    if(is_empty(negative)==FALSE){
      fisherZ[negative,"Z_tot_re"]<-fisherZ[negative,"Z_tot"]
      fisherZ[negative,"SE_re"]<-fisherZ[negative,"SE_fe"]
      fisherZ[negative,"z_score_re"]<-fisherZ[negative,"z_score_fe"]
      
      if(is_empty(positive)==FALSE){
        fisherZ[positive,"Z_tot_re"]<-up[positive]/down[positive]
        fisherZ[positive,"SE_re"]<-sqrt(1/ down[positive])
        fisherZ[positive,"z_score_re"]<-fisherZ[positive,"Z_tot_re"]/fisherZ[positive,"SE_re"]
      }
    }
    else{
      fisherZ$Z_tot_re<-up/down
      fisherZ$SE_re<-sqrt(1/ down)
      fisherZ$z_score_re<-fisherZ$Z_tot_re/fisherZ$SE_re
      
    }
    fisherZ$pval<-2*pnorm(q=abs(fisherZ$z_score_re),lower.tail=FALSE)
    fisherZ$cor<-FisherZInv(fisherZ$Z_tot_re)
    fisherZ<-fisherZ[,c("from","to","cor","pval")]
    fisherZ<-full_join(fisherZ, no_meta)
    fisherZ$fdrs<-p.adjust(p = fisherZ$pval,method = 'fdr')
    return(fisherZ)
  } 
  
}

# Calculate Pearson's correlation coefficients and
# perform  meta-analysis
calc_pearson_metacor <- function(list_of_files, meta_method= "random", pvalue_thres=NULL, fdr_thres=NULL, coef_thres=NULL){
  
  if (is.null(fdr_thres) == TRUE & is.null(pvalue_thres)==TRUE){
    stop("Error: You should use pvalue or fdr threshold")
  }
  
  cor_list <-  list()
  for (i in 1:length(list_of_files)){
    
    # Calculating Pearson's correlations coefficients using cor() and melt
    pearson <- cor(list_of_files[[i]])
    pearson<-adjmatrix_to_edgelist(pearson)
    
    # Remove self correlations
    #pearson<-pearson[which(pearson$value!=1),]
    
    #Rename the columns
    #names(pearson) <- c("from","to", "ri")
    if(is.null(coef_thres)==FALSE){
      pearson<-pearson[which(abs(pearson$ri) >= coef_thres),]
    }
    pearson$ni <-nrow(list_of_files[[i]])
    #print(pearson)
    cor_list <- append(cor_list ,list(pearson) )
    
  }
  # Merge the studies
  corr_list_data <- cor_list %>% purrr::reduce(full_join, by = c("from","to"))
  
  # Keep correlations above a threshold
  #meta_cor  <- run_metacor(correlations = corr_list_data)
  meta_cor  <- my_meta(correlations = corr_list_data , method= meta_method)
  
  # Run metacor with p-value threshold  or FDR
  if (is.null(fdr_thres) == FALSE & is.null(pvalue_thres)==TRUE ){
    cropped_meta_cor<-meta_cor[which(meta_cor$fdrs<fdr_thres),]
    
  }
  else if(is.null(fdr_thres) == TRUE & is.null(pvalue_thres)==FALSE){
    cropped_meta_cor<-meta_cor[which(meta_cor$pval<pvalue_thres),]
  }
  else if(is.null(fdr_thres) == FALSE & is.null(pvalue_thres)==FALSE){
    cropped_meta_cor<-meta_cor[which(meta_cor$pval<pvalue_thres & meta_cor$fdrs<fdr_thres),]
  }
  
  return (cropped_meta_cor)
  
  
}


# Make a function which will reduce the columns of the initial studies 
# The second argument needs to be a a correlation data.table 

keep_sig_pairs <-function (list_of_studies,sig_data){   
  
  cropped_studies <- list()
  sig_genes<-union(unique(sig_data$from), unique(sig_data$to))
  for (i in 1:length(list_of_studies)){
    
    cropped_studies[[i]]<-list_of_studies[[i]][, ..sig_genes]
    
  }
  return(cropped_studies)
  
}

probe2geneID <- function(list_of_files, GPL_list){
  
  #Load GPL platform 
  ids = idmap(GPL_list,destdir=tempdir())
  x = t(filterEM(t(list_of_files), ids))
  x <- data.table(x)
  
  return (x)
  
}


#create edge list from matrix
# if symmetric duplicates are removed

adjmatrix_to_edgelist<-function(mat){
  mat <-upper.triangle(mat)
  
  obj<-as.data.frame.table(mat)
  
  names(obj)<-c("from","to","ri")
  obj<-obj[which(obj$ri!=1),]
  obj<-obj[which(obj$ri!=0),]
  return(data.table(obj))
}


run_metacor <- function (correlations){
  
  meta_cor<-correlations[,1:2]
  
  
  ri<- correlations[,.SD,.SDcols=grep("ri",names(correlations),value =TRUE)]
  ni<- correlations[,.SD,.SDcols=grep("ni",names(correlations),value =TRUE)]
  
  for(i in 1:nrow(correlations)){
    m1 <- metacor(as.numeric(ri[i]),as.numeric(ni[i]),sm="ZCOR",method.tau = "HE")
    meta_cor$cor.fe[i]=transf.ztor(m1[["TE.fixed"]])
    meta_cor$cor.re[i]=transf.ztor(m1[["TE.random"]])
    meta_cor$pval[i] = m1[["pval.random"]]
  }
  meta_cor$fdrs<-p.adjust(p = meta_cor$pval,method = 'fdr')
  
  return (meta_cor)
  
}

volc_plot<- function (meta_data,coef_thres,pval_thres) {
  
  # Make plots 
  # The basic scatter plot: x is "partial correlation", y is "-log10(p-value)"
  
  # if partial correlation coefficient > coef_thres and pvalue < pval_thres, set as "UP Significant" 
  # else if partial correlation coefficient < coef_thres and pvalue < pval_thres, set as "DOWN Significant" 
  # else set as "Non Significant"
  
  meta_data$expression = ifelse(meta_data$pval < pval_thres & abs(meta_data$cor) >= coef_thres , ifelse(meta_data$cor> coef_thres ,'Up Significant','Down Significant'), 'Non Significant')
  
  ggplot(data = meta_data,
         aes(x = cor,
             y = -log10(pval), 
             colour=expression)) +
    geom_point() +
    # geom_text_repel(aes(label = from), 
    #                 min.segment.length = 0, 
    #                 seed = 42, 
    #                 box.padding = 0.1) +
    scale_color_manual(values=c("blue", "grey","red"))+
    # xlim(c(-pval_thres, pval_thres)) +
    # geom_vline(xintercept=c(-coef_thres,coef_thres),lty=4,col="black",lwd=0.8) +
    # geom_hline(yintercept = -log10(pval_thres),lty=4,col="black",lwd=0.8)+
    labs(x="partial correlation coefficient",
         y="-log10(p-value)",
         title="Volcano plot")
  
}


# Graphviz Simple Network Visualization 
network_plot <- function (file) {
  file<- file[order(-file$cor),]
  file <- head(file,5000)
  
  nodes = data.frame(id =unique(c(file$from,file$to)) , label =unique(c(file$from,file$to)) )
  # colnames(nodes) <- c("id", "label")
  
  edges = as.data.frame(file[,c(1,2,3)])
  
  colnames(edges) <- c("from", "to", "width")
  
  
  # in circle ?
  
  net_plot<- visNetwork(nodes, edges) %>%
    visIgraphLayout(layout = "layout_nicely") %>%
    visNodes(size = 5) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T), 
               nodesIdSelection = T)
  
  net_plot
  
}

# Enrichment Analysis and Manhattan Plot
enrichment_analysis <-function (file){
  
  
  # using subset function
  newdata <- subset(file, pval<0.05,
                    select=c(from,to))
  
  DEGs_list_from = newdata$from
  DEGs_list_to = newdata$to
  
  print(DEGs_list_to)
  DEGs = unique(c(DEGs_list_to,DEGs_list_from))
  print(DEGs)
  # enrichment analysis
  # gp_from = gost(DEGs_list_from, organism = "hsapiens")
  # gp_to = gost(DEGs_list_to, organism = "hsapiens") 
  gp_DEGs = gost(as.vector(DEGs), organism = "hsapiens") 
  
  #Manhattan Plot (Interactive)
  p<- gostplot(
    gp_DEGs,
    capped = TRUE,
    interactive = TRUE,
    pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
              "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99",
            HPA ="#6633cc", CORUM = "#66aa00", HP = "#990099")
  )
  p
  
  return (list(gp_DEGs[['result']],p))
  
}


meta_pcor <- function(folder_path=NULL, GEO_names=NULL, GPL_list= NULL, target_namespace = NULL, my_list=NULL, option, method, meta_method= "random", pvalue_thres = NULL, fdr_thres = NULL, coef_thres = NULL,l1  = NULL ,l2 = NULL ){
  
  value1 = pvalue_thres
  value2 = fdr_thres
  value3 = coef_thres
  value4 = l1
  value5 = l2
  value6 = folder_path
  value7 = target_namespace
  
  # Meta-analysis and partial correlation meta-analysis (folderpath needed)
  if (option == 4 ){
    
    # Meta-analysis with the DEGs of the studies
    
    DEG_simple <- DE_analysis(folder_path = value6 ,case = 'CASE', control = 'CONTROL', l1 = 0.6 , fold_threshold = 0.0, p_value_threshold = 0.05 )
    

    
    return(DEG_simple)
    
  }
  
  list_of_files<-load_GSE_data(folder_path, GEO_names, GPL_list, value7)
  
  if (option ==1){
    ##################################################
    #                   Option 1                     #
    ##################################################
    
    
    # Calculate Pearson's correlation coefficients and
    # perform  meta-analysis
    pearson_meta_cor <- calc_pearson_metacor(list_of_files = list_of_files,pvalue_thres = value1 ,fdr_thres =  value2 ,coef_thres = value3 )
    
    
    # Return a cropped list of studies (with significant correlations pairs which occurred either from meta-analysis or cut-offs)
    # to calculate partial correlations or conduct meta-analysis with the metacor package
    result <- keep_sig_pairs (list_of_studies = list_of_files, sig_data = pearson_meta_cor)
    
    # Calculate sparse partial correlation estimation with neighborhood selection approach method for each study
    
    if (method == 'sparse') {
      # We do not want any cut offs when significant == FALSE
      pcor_list <- pcor_neighborhood(list_of_files = result, l1 = value4, l2 = value5, significant=FALSE, pvalue_thres=NULL, fdr_thres=NULL, coef_thres=NULL)
    }
    else if (method == 'shrinkage'){
      pcor_list <- pcor_shrinkage(list_of_files = result, significant= FALSE , pvalue_thres = NULL, fdr_thres = NULL , coef_thres = NULL )
      
    }
    
    # Perform meta-analysis of correlations coefficients
    meta_cor  <- my_meta(correlations = pcor_list, method=meta_method)
    
    return(meta_cor)
    
  }
  
  
  if (option ==2) {
    
    # Calculate partial correlation with thresholds 
    if (method == 'sparse'){
      
      pcor_list <- pcor_neighborhood(list_of_files = list_of_files, l1 = value4, l2 = value5, significant=TRUE, pvalue_thres=value1, fdr_thres=value2, coef_thres=value3)      
      
    }
    else if (method == 'shrinkage') {
      pcor_list <- pcor_shrinkage(list_of_files = list_of_files, significant=TRUE, pvalue_thres=value1, fdr_thres=value2, coef_thres=value3)      
      
    }
    
    #print(pcor_list)
    
    # Perform meta-analysis of correlations coefficients
    meta_cor  <- my_meta(correlations = pcor_list, method=meta_method)
    
    return(meta_cor)
    
  }
  
  
  if (option ==3) {
    
    # Calculate partial correlation with thresholds 
    if (method == 'sparse'){
      
      pcor_list <- pcor_neighborhood(list_of_files = list_of_files, l1 = value4, l2 = value5, significant=FALSE, pvalue_thres=value1, fdr_thres=value2, coef_thres=value3)      
      
    }
    else if (method == 'shrinkage') {
      pcor_list <- pcor_shrinkage(list_of_files = list_of_files, significant=FALSE, pvalue_thres=value1, fdr_thres=value2, coef_thres=value3)      
      
    }
    
    #Perform meta-analysis of correlations coefficients
    
    meta_cor  <- my_meta(correlations = pcor_list, method=meta_method)
    
    return(meta_cor)
    
  }

  
}


run_metacor <- function (correlations){
  
  meta_cor<-correlations[,1:2]
  
  
  ri<- correlations[,.SD,.SDcols=grep("ri",names(correlations),value =TRUE)]
  ni<- correlations[,.SD,.SDcols=grep("ni",names(correlations),value =TRUE)]
  
  for(i in 1:nrow(correlations)){
    m1 <- metacor(as.numeric(ri[i]),as.numeric(ni[i]),sm="ZCOR",method.tau = "HE")
    meta_cor$cor.fe[i]=transf.ztor(m1[["TE.fixed"]])
    meta_cor$cor.re[i]=transf.ztor(m1[["TE.random"]])
    meta_cor$pval[i] = m1[["pval.random"]]
  }
  meta_cor$fdrs<-p.adjust(p = meta_cor$pval,method = 'fdr')
  
  return (meta_cor)
  
}


volc_plot_plotly<- function(pcor,pval_thres,coeff_thres){
  pcor$group ='Non Significant'
  pcor$group[pcor$pval<pval_thres& abs(pcor$cor)>coeff_thres] ='Up Significant'
  pcor$group[pcor$pval<pval_thres& abs(pcor$cor)<coeff_thres] ='Down Significant'
  
  
  
  # make the Plot.ly plot
  p <- plot_ly(type = 'scatter', x = pcor$cor, y = -log10(pcor$pval), text = paste(pcor$from,pcor$to), mode = "markers", color = pcor$group, colors=c('red','gray','blue')) %>% 
    layout(title ="Volcano Plot") 
  
  return (p)
  
}


DE_analysis <- function(folder_path,case,control, fold_threshold,p_value_threshold,l1 ){
  
  expr_mat = list()
  
  list_of_studies = list.files(path = folder_path, pattern ='.txt', all.files=FALSE, full.names=TRUE)
  
  l1_value = l1
  for (i in 1:length(list_of_studies)){
    # Read the file
    # inFile <- list_of_studies[i]
    # old_name <- inFile$datapath
    # dirstr <- dirname(inFile$datapath)
    # new_name <- paste(dirstr, inFile$name,sep="/")
    # file.rename(old_name, new_name)
    print("Reading files")
    print(as.vector(list_of_studies[[i]]))
    study <- as.data.frame(read.table(as.vector(list_of_studies[[i]]),sep ='\t',header = TRUE))
    
    study_t <- study[-1,]
    study_t <- t(study_t)
    
    colnames(study_t) <- study_t[1,]
    study_t<- study_t[-1,]
    rownames(study_t) <- 1:nrow(study_t)
    study_t <- as.data.table(study_t)
    
    study_t <- lapply(study_t,as.numeric)
    study_t <- as.data.table(study_t)
    
    study <- as.data.frame(study)
    
    # Split cases and controls
    cases <- study [,(study[1,]) == 'CASE' ]
    controls <- study [,(study[1,]) == 'CONTROL']
    
    
    
    
    
    # Delete first row (we don't need the annotation row now )
    cases <- cases[-1,]
    rownames(cases) <- colnames(study_t)
    cases<- t(cases)
    # colnames(cases) <- cases[1,]
    # cases<- cases[-1,]
    # rownames(cases) <- 1:nrow(cases)
    # cases <- as.data.frame(cases)
    #
    
    
    controls <- controls[-1,]
    rownames(controls) <- colnames(study_t)
    
    controls<- t(controls)
    
    
    cases <- matrix(as.numeric(cases),    # Convert to numeric matrix
                    ncol = ncol(cases))
    colnames(cases) <- colnames(study_t)
    
    
    controls <- matrix(as.numeric(controls),    # Convert to numeric matrix
                       ncol = ncol(controls))
    colnames(controls) <- colnames(study_t)
    
    
    
    
    # Apply log2 scale in cases and controls
    controls <- log2(controls)
    cases <- log2(cases)
    
    print(cases)
    print(controls)
    # Calculate the means of each group
    
    group1 <- apply(t(cases), 1, mean)
    group2 <- apply(t(controls), 1, mean)
    
    
    

    
    # Compute fold-change for biological significance
    # Difference between the means of the two conditions
    fold = group1 - group2
    #fold
    
    
    
    
    
    #pvalue
    pvalue = NULL # Empty list for the p-values
    tstat = NULL # Empty list of the t test statistics
    
    
    
    
    for (i in 1: length(colnames(cases))){
      t = t.test(cases[,i], controls[,i])
      # Put the current p-value in the pvalues list
      pvalue[i] = t$p.value
      # Put the current t-statistic in the tstats list
      tstat[i] = t$statistic
    }
    
    #give your fold cutoff
    fold_cutoff = fold_threshold
    #give your pvalue cutoff
    pvalue_cutoff = p_value_threshold
    
    
    # Fold-change filter for biological significance
    filter1 = abs(fold) >= fold_threshold
    
    
    # P-value filter for statistical significance
    filter2 = pvalue <= pvalue_cutoff
    
    names(filter2) <- names(filter1)
    
    cbind(filter1,filter2)
    
    # Combined filter
    filter3 = filter1 & filter2
    
    
    #filtered = study_t[filter3,]
    sig = as.data.frame(filter3)
    filtered = row.names(sig)[which(sig$filter3==TRUE)]
    
    # Finally we filtrer the dataset with the DE genes
    final = t(t(study_t)[rownames(t(study_t)) %in% filtered,])
    
    final <- as.data.table(final)
    expr_mat <- append(expr_mat,list(final))
    
  }
  
  #Sparse Correlation
  pcor_list <- pcor_neighborhood(list_of_files = expr_mat, l1 = l1_value, l2 = 0, significant=FALSE)
  
  
  
  #Perform meta-analysis of correlations coefficients
  
  meta_cor  <- my_meta(correlations = pcor_list, method='random')
  return(meta_cor)
}



