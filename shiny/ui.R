library(shiny)
source('fun_test.R')
library(igraph)
library(shinycssloaders)
library(shinyjs)
#library(networkD3)
# Define UI for dataset viewer app ----
ui <- fluidPage(

  # App title ----
  titlePanel("MetaPcor"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(



      #Sidebar controls
      # Input: Text for providing a caption ----
      # Note: Changes made to the caption in the textInput control
      # are updated in the output area immediately as you type

      # Change sidebar color
      tags$style(".well {background-color:#EBECF0;}"),

      helpText("Upload your file(s) here to proceed"),
      # helpText("For partial correlation meta-analysis, this is the proposed input"),
      #
      #
      # helpText("For DE and partial correlation meta-analysis, this is the proposed input"),
      #

      # textInput(inputId = "caption",
      #           label = "Caption:",
      #           value = "Data Summary"),


      fileInput("upload", NULL, buttonLabel = "Upload GSE Matrix file(s)", multiple = TRUE),
      tableOutput("files"),
      # Add a download button for the example file
      helpText("Download a GSE Matrix file"),

      downloadButton("download_example", "GSE_matrix example"),

      helpText("Download a GSE Matrix for DE and partial correlation meta-analysis"),

      # Add a download button for the example file
      downloadButton("download_example_de", "GSE_matrix DE example"),
      # downloadButton("download_examples", "Download Example Files (.zip)"),
      # Input: Selector for choosing dataset ----
      selectInput(inputId = "norm_data",
                  label = "Normalize Data:",
                  choices = c("NO", "YES")),
      uiOutput("new_dropdown"),




      # Input: Selector for choosing dataset ----
      selectInput(inputId = "option",
                  label = "Choose an option:",
                  choices = c("Partial correlation meta-analysis","Partial correlation meta-analysis with thresholds" ,"DE and partial correlation meta-analysis")),



      # selectInput(inputId = "significant",
      #             label = "Significant:",
      #             choices = c("FALSE","TRUE")),

      conditionalPanel(condition = "input.option == 'Partial correlation meta-analysis with thresholds'",
                       numericInput(inputId = "pvalue_thres",
                                    label = "p-value threshold:",
                                    value = 0.01),
                       numericInput(inputId = "fdr_thres",
                                    label = "FDR threshold:",
                                    value = NULL),
                       numericInput(inputId = "coef_thres",
                                    label = "Coefficient threshold:",
                                    value = NULL)),


      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "l1",
                   label = "l1 (Norm Penalty):",
                   value = 0.8),



      # Input: Selector for choosing dataset ----
      selectInput(inputId = "ea_option",
                  label = "Enrichment Analysis:",
                  choices = c("NO","YES")),
      helpText("Volcano plot options. When the Volcano plot is produced, you can change these values dynamically"),

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "coeff_volc_plot",
                   label = "Coefficient threshold for Volcano plot",
                   value = 0),

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "pval_volc_plot",
                   label = "p-value threshold for Volcano plot",
                   value = 0.01),

      actionButton("execute", "Execute",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),

      #submitButton("Insert GEO accession numbers", icon("play"))





    ),

    # Main panel for displaying outputs ----

    mainPanel(

      h4('You entered: '),
      verbatimTextOutput("oid1"),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div(id = "loading-message",
                                style = "position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);
                                                    z-index: 10000; background-color: rgba(255, 255, 255, 0.7);
                                                    padding: 20px; border-radius: 10px;",
                                tags$p("Loading...", style = "font-size: 18px; font-weight: bold; text-align: center;"))),

      # downloadButton("downloadData", "Download meta-analysis results"),
      conditionalPanel(condition = "input.execute > 0",
                       h4('Meta-analysis results')),
      DT::dataTableOutput("final_head"),


      # downloadButton("downloadData", "Download meta-analysis results"),
      conditionalPanel(condition = "input.execute > 0",
                       downloadButton("downloadData", "Download meta-analysis results")),

      # downloadButton("downloadData", "Download meta-analysis results"),
      conditionalPanel(condition = "input.ea_option =='YES' && input.execute > 0",
                       h4('Enrichment analysis results')),

      DT::dataTableOutput("ea_results"),

      # downloadButton("downloadEAData", "Download Enrichment analysis results"),
      # downloadButton("downloadEAData", "Download Enrichment analysis results",
      #                conditionalPanel(condition = "input.execute > 0")),

      conditionalPanel(condition = "input.ea_option =='YES' && input.execute > 0",

                       downloadButton("downloadEAData", "Download Enrichment analysis results")),
      conditionalPanel(condition = "input.execute > 0",


                      h4('Plots'),
                      tabsetPanel(
                        tabPanel("Network",
                                 downloadButton("downloadPlot", "Download Network"),
                                 visNetworkOutput('network')


                        ),
                        tabPanel("Volcano Plot",
                                 plotlyOutput("volc_plot")
                        ),

                        tabPanel("Manhattan Plot",
                                 plotlyOutput("manhattan_plot")

                        ),

                      ),
                      # conditionalPanel(
                      # condition = "input.execute > 0",
                      # downloadButton("downloadPlot", "Download Network"))
                      ),

                      ),


      # downloadButton("downloadPlot", "Download Network"),

      # downloadButton("downloadPlot", "Download Network",
      #                conditionalPanel(condition = "input.plotTabs == 'Network'")),
      # conditionalPanel(
      # condition = "input.tabset == 'Network'",
      # shinyjs::hidden(downloadButton("downloadPlot", "Download Network"))),
      # downloadButton("downloadVolcPlot", "Download Volcano plot"),
      # Output: Verbatim text for data summary ----
      # Output: HTML table with requested number of observations ----
    )
  )
