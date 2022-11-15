library(shiny)
library(AnnoProbe)
library(igraph)
library(networkD3)
library(matrixcalc)

memory.limit(size = 100000)
#Limit of each file for server : 80MB for each file
options(shiny.maxRequestSize=80*1024^2) 

source('fun_test.R')
# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  
  output$oid1 <- renderPrint({
    cat("MetaPcor: A package for meta-analysis of Gene Expression studies with partial correlation as effect size\n\n")
    x <- unlist(strsplit(input$GEO_names,","))
    cat("GEO Series:\n")
    print(x)

    cat("Analysis Option:\n")
    print(input$option)
    cat("Method:\n")
    print(input$method)
    cat("Significant:\n")
    print(input$significant)
    cat("P-value threshold:\n")
    print(input$pvalue_thres)
    cat("FDR threshold:\n")
    print(input$fdr_thres)
    cat("Coefficient threshold:\n")
    print(input$coef_thres)
    
    if (is.na(input$l1) == FALSE & is.na(input$l2) == FALSE){
      
      cat("l1:\n")
      print(input$l1)
      cat("l2:\n")
      print(input$l2)
      
      cat("Input Files\n")
      print(input$upload)
      
      
    }
    
  }
  
  )
  
  output$files <- renderTable(input$upload)
  
  
  # When the button is pressed
  observeEvent(input$execute,{
    inFile <- input$upload
    file_list <- list()
    # Read TAB delimited files
    if (is.null(inFile)!= TRUE) {
      for (i in 1:length(inFile$datapath)) {
        x<- inFile$datapath[[i]]
        file_list<- append(file_list, list(x))
      }
      
      if (input$option == "DE and partial correlation meta-analysis")
      {
        
        # Insert DE_analysis
        final <- DE_analysis(file_list,case ='CASE',control= 'CONTROL', fold_threshold = 0.05, p_value_threshold = 0.05 )
        
  
        
        print(final)
        output$final_head <- DT::renderDataTable({final})
        
        # Encrichment Analysis with gProfiler
        
        if (input$ea_option =="YES"){
          ea_results <- enrichment_analysis(as.data.frame(final))
          output$ea_results <- DT::renderDataTable({ea_results[[1]]})
          
          
          
          output$downloadEAData <- downloadHandler(
            filename = function() {
              paste("ea_results", ".csv",sep = '')
            },
            content = function(file) {
              write.csv(ea_results, file)
            }
          )
          
        }
        
        
        
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("meta_pcor_results", ".csv",sep = '')
          },
          content = function(file) {
            write.csv(final, file)
          }
        )
        # final_net<- reactive({data.frame(final)})
        
        
        output$network<- renderVisNetwork({
          network_plot(as.data.frame(final))
        })
        
        
        
        
        net_file <- network_plot(as.data.frame(final))
        
        # Volcano plot 
        output$volc_plot <- renderPlotly({volc_plot_plotly(pcor  = as.data.frame(final), pval_thres = input$pval_volc_plot,coeff_thres = input$coeff_volc_plot)})
        
        
        output$downloadPlot <- downloadHandler(
          filename = function() {
            paste("meta_pcor_network", ".html",sep = '')
          },
          content = function(file) {
            #write.csv(output$network, file)
            saveNetwork(net_file, file, selfcontained = TRUE)
          }
        )
        
        
        # Manhattan plot 
        output$manhattan_plot <- renderPlotly(  #Manhattan Plot (Interactive)
          p<- gostplot(
            ea_results[[2]],
            capped = TRUE,
            interactive = TRUE,
            pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
                      "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99",
                    HPA ="#6633cc", CORUM = "#66aa00", HP = "#990099")
          ))
        
        
        
        output$downloadVolcPlot <- downloadHandler(
          filename = function() {
            paste("meta_pcor_volc_plot", ".png",sep = '')
          },
          content = function(file) {
            ggsave(file,volc_plot(meta_data = final,coef_thres = input$coeff_volc_plot, pval_thres = pval_volc_plot) )
          }
        )

        
        
      }
      
      
    }
    
    else{
      x <- unlist(strsplit(input$GEO_names,","))
      
      # list_of_studies<- load_GSE_data(GEO_names = x)
      # Run the main function
      # GPL_list  <- unlist(strsplit(input$GPL_names,","))
      # gene_annot_list <- probe2geneID(list_of_files = list_of_studies,GPL_list = GPL_list)
      # new_gene_annot_list <- Filter(NROW, gene_annot_list)
      if(input$significant == 'FALSE'){
        bool = FALSE
      }
      else{
        bool = TRUE
      }
      # final <- meta_pcor(option = input$option , method = input$method, list_of_files = list_of_studies, pvalue_thres = input$pvalue_thres ,fdr_thres =  input$fdr_thres ,coef_thres = input$coef_thres, l1 = 0.8,l2 =0 )

      print(input$option)
      
      
  
      final<- meta_pcor(GEO_names=x, GPL_list= NULL, option=as.character(input$option) , method=input$method, meta_method= "random", pvalue_thres = input$pvalue_thres, fdr_thres = input$fdr_thres, coef_thres = input$coef_thres,l1  = input$l1 ,l2 = input$l2)
      
      print(final)
      output$final_head <- DT::renderDataTable({final})
      
      # Encrichment Analysis with gProfiler
      
      if (input$ea_option =="YES"){
        ea_results <- enrichment_analysis(as.data.frame(final))
        output$ea_results <- DT::renderDataTable({ea_results[[1]]})
        
        
        
        output$downloadEAData <- downloadHandler(
          filename = function() {
            paste("ea_results", ".csv",sep = '')
          },
          content = function(file) {
            write.csv(as.data.table(ea_results), file)
          }
        )
        
      }
    
      
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("meta_pcor_results", ".csv",sep = '')
        },
        content = function(file) {
          write.csv(final, file)
        }
      )
      # final_net<- reactive({data.frame(final)})
      
      
      output$network<- renderVisNetwork({
        network_plot(as.data.frame(final))
      })
      
      
     
      
      
      net_file <- network_plot(as.data.frame(final))
      
      # Volcano plot 
      output$volc_plot <- renderPlotly({volc_plot_plotly(pcor  = as.data.frame(final), pval_thres = input$pval_volc_plot,coeff_thres = input$coeff_volc_plot)})
      
      
      output$downloadPlot <- downloadHandler(
        filename = function() {
          paste("meta_pcor_network", ".html",sep = '')
        },
        content = function(file) {
          #write.csv(output$network, file)
          saveNetwork(net_file, file, selfcontained = TRUE)
        }
      )
      
      
      # Manhattan plot 
      output$manhattan_plot <- renderPlotly(  #Manhattan Plot (Interactive)
        p<- gostplot(
          ea_results[[2]],
          capped = TRUE,
          interactive = TRUE,
          pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", KEGG =
                    "#dd4477", REAC = "#3366cc", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99",
                  HPA ="#6633cc", CORUM = "#66aa00", HP = "#990099")
        ))
      
      
      
      output$downloadVolcPlot <- downloadHandler(
        filename = function() {
          paste("meta_pcor_volc_plot", ".png",sep = '')
        },
        content = function(file) {
          ggsave(file,volc_plot(meta_data = final,coef_thres = input$coeff_volc_plot, pval_thres = pval_volc_plot) )
        }
      )
   
      
    }# end else
    
    
    
  })
  

 
  
  
  
}
