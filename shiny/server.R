library(GEOquery)
library(corpcor)
library(space)
library(data.table)
library(GeneNet)
library(meta)
library(metafor)
library(dplyr)
library(tidyverse)
library(DescTools)
library(ggplot2)
library(ff)
library(matrixcalc)
library(gprofiler2)
library(ggrepel)
library(GGally)
library(plotly)
library(visNetwork)
library(geomnet)
library(igraph)
library(purrr)
library(readxl)
library(matrixcalc)
library(Cairo)
library(networkD3)
options(bitmapType='cairo')
memory.limit(size = 100000)
#Limit of each file for server : 800MB for each file
options(shiny.maxRequestSize=800*1024^2)

source('fun_test.R')
# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {

  output$download_example <- downloadHandler(
    filename = "GSE_demo.txt",
    content = function(file) {
      # Path to the example data file
      example_file <- "GSE_demo.txt"

      # Copy the example data file to the download destination
      file.copy(example_file, file)
    }
  )


  output$download_example_de <- downloadHandler(
    filename = "GSE_demo_DE.txt",
    content = function(file) {
      # Path to the example data file
      example_file <- "GSE_demo_DE.txt"

      # Copy the example data file to the download destination
      file.copy(example_file, file)
    }
  )

  # # Define a function to create a zip file
  # createZip <- function() {
  #   zip_file <- "data_examples.zip"
  #   dir_path <- "data_examples"
  #   zip(zip_file, files = list.files(dir_path, full.names = TRUE))
  #   return(zip_file)
  # }
  # output$download_examples <- downloadHandler(
  #   filename = function() {
  #     "example_files.zip"
  #   },
  #   content = function(file) {
  #     file.copy(createZip(), file)
  #   }
  # )

  output$oid1 <- renderPrint({
    cat("MetaPcor: A package for meta-analysis of Gene Expression studies with partial correlation as effect size\n\n")

    cat("Analysis Option:\n")
    print(input$option)

    cat("Significant:\n")
    print(input$significant)
    cat("P-value threshold:\n")
    print(input$pvalue_thres)
    cat("FDR threshold:\n")
    print(input$fdr_thres)
    cat("Coefficient threshold:\n")
    print(input$coef_thres)


    cat("l1:\n")
    print(input$l1)


    cat("Input Files\n")
    print(input$upload)




  }

  )
  new_dropdown_value <- reactive({
    input$norm_method
  })
  output$files <- renderTable(input$upload[1:2])
  observeEvent(input$norm_data, {
    if(input$norm_data == "YES") {
      output$new_dropdown <- renderUI({
        selectInput(inputId = "norm_method",
                    label = "Normalization Method:",
                    choices = c('RPKM', 'QUANT_NORM', 'Z_SCOR_STAND','LOG2'))
      })
    } else {
      output$new_dropdown <- renderUI({})
    }
  })


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
      print (input$option)
      if (input$option == "DE and partial correlation meta-analysis")
      {

        # Insert DE_analysis
        final <- DE_analysis(file_list,case ='CASE',control= 'CONTROL', fold_threshold = 0.05, p_value_threshold = 0.05, l1 = input$l1)

      }
      else{

        if (input$model == "Random Effects model")
        {
          model = 'random'
        }
        else{
          model = 'fixed'
        }
        final<- meta_pcor(file_names = file_list, option=as.character(input$option), method = "sparse", meta_method= model, pvalue_thres = input$pvalue_thres, fdr_thres = input$fdr_thres, coef_thres = input$coef_thres,l1  = input$l1 ,l2 = 0, norm_data = input$norm_data, norm_method = input$norm_method )
        print(final)

      }

      print(final)
      # Hide the loading message
      Sys.sleep(5)
      tags$script(HTML("$('#loading-message').hide();"))
      # tags$script(HTML("$('#loading-message').hide();"))
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
            fwrite(ea_results[[1]], file)

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





  )






}
