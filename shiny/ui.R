library(shiny)
source('fun_test.R')
library(igraph)
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
      helpText("Choose your analysis settings to proceed"),

      # textInput(inputId = "caption",
      #           label = "Caption:",
      #           value = "Data Summary"),


      fileInput("upload", NULL, buttonLabel = "Upload .xlsx file(s)", multiple = TRUE),
      tableOutput("files"),




      # Input: Selector for choosing dataset ----
      selectInput(inputId = "option",
                  label = "Choose an option:",
                  choices = c("Partial correlation meta-analysis with thresholds", "Partial correlation meta-analysis","DE and partial correlation meta-analysis")),



      # Input: Selector for choosing dataset ----
      selectInput(inputId = "significant",
                  label = "Significant:",
                  choices = c("TRUE", "FALSE")),



      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "pvalue_thres",
                   label = "p-value threshold:",
                   value = 0.01),


      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "fdr_thres",
                   label = "FDR threshold:",
                   value = NULL),
      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "coef_thres",
                   label = "Coefficient threshold:",
                   value = NULL),

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "l1",
                   label = "l1 (Norm Penalty):",
                   value = 0.8),



      # Input: Selector for choosing dataset ----
      selectInput(inputId = "ea_option",
                  label = "Enrichment Analysis:",
                  choices = c("YES", "NO")),

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "coeff_volc_plot",
                   label = "Coefficient threshold for Volcano plot",
                   value = 0),

      # Input: Numeric entry for number of obs to view ----
      numericInput(inputId = "pval_volc_plot",
                   label = "p-value threshold for Volcano plot",
                   value = 0.01),

      actionButton("execute", "Execute the Main Function",style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),

      #submitButton("Insert GEO accession numbers", icon("play"))





    ),

    # Main panel for displaying outputs ----

    mainPanel(

      h4('You entered: '),
      verbatimTextOutput("oid1"),
      #verbatimTextOutput("oid2")

      DT::dataTableOutput("final_head"),

      downloadButton("downloadData", "Download meta-analysis results"),

      DT::dataTableOutput("ea_results"),

      downloadButton("downloadEAData", "Download Enrichment analysis results"),


      tabsetPanel(
        tabPanel("Network",
                 visNetworkOutput('network')
        ),

        tabPanel("Volcano Plot",
                 plotlyOutput("volc_plot")
        ),

        tabPanel("Manhattan Plot",
                 plotlyOutput("manhattan_plot")

        ),

      ),

      downloadButton("downloadPlot", "Download Network"),
      # downloadButton("downloadVolcPlot", "Download Volcano plot"),
      # Output: Verbatim text for data summary ----
      # Output: HTML table with requested number of observations ----

    )
  )
)
