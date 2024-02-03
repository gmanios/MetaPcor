library(shiny)
library(shinythemes)  # Load the shinythemes package
source('fun_test.R')
library(igraph)
library(shinycssloaders)
library(shinyjs)

setDTthreads(150)
ui <- navbarPage(

  # App title ----
  title = "MetaPcor",

  # Apply the Yeti theme to the entire app
  theme = shinytheme("yeti"),  # Use the Yeti theme
  # Add custom CSS to change the navbar color
  tags$head(
    tags$style(
      HTML(
        "
      /* Navbar background color */
      .navbar-default {
        background-color: #2D396A;
      }

      /* Hover or click background color for links */
      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:hover,
      .navbar-default .navbar-nav > .active > a:focus,
      .navbar-default .navbar-nav > li > a:hover,
      .navbar-default .navbar-nav > li > a:focus,
      .navbar-default .navbar-nav > li > a:active,
      .navbar-default .navbar-nav > li > a.active,
      .navbar-default .navbar-nav > .open > a,
      .navbar-default .navbar-nav > .open > a:hover,
      .navbar-default .navbar-nav > .open > a:focus {
        background-color: #1A223E;
      }

      /* Loading bar background color */
      .shiny-progress-bar {
        background-color: lightgreen;
      }

      /* Button background color */
      .btn-primary {
        background-color: blue;
        border-color: blue;
      }
      "
      )
    )
  ),

  # Add Navbar with tabs ----
  tabPanel(
    "Welcome",
    fluidPage(
      # Your Welcome tab content goes here
      h2("Welcome to MetaPcor! "),
      br(),
      br(),
      "MetaPcor : A package to conduct meta-analysis of co-expression networks using partial correlation as effect size",
      br(),
      br(),
      "Correspondence : pbagos@compgen.org",
      br(),
      br(),
      "Download the source code of MetaPcor from",
      tags$a(href="https://github.com/gmanios/MetaPcor", "GitHub"),
      br(),
      br(),
      tags$a(href="http://www.compgen.org","http://www.compgen.org"),
      br(),
      br(),
      "Main developers",
      tags$li("Georgios A. Manios - University of Thessaly"),
      tags$li("Ioanna V. Sasilioglou - University of Thessaly")



    )
  ),
  tabPanel(
    "Analysis",
    fluidPage(
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
          
          conditionalPanel(condition = "input.option == 'DE and partial correlation meta-analysis'",
                           selectInput(inputId = "log2_norm",
                                       label = "Apply log tranformation:",
                                       choices = c( "TRUE","FALSE")),
                           numericInput(inputId = "log2FC_threshold",
                                        label = "log2FC threshold:",
                                        value = 1.0)),
          
          # Input: Selector for choosing dataset ----
          selectInput(inputId = "model",
                      label = "Choose meta-analysis model:",
                      choices = c("Random Effects model","Fixed Effects model")),

          # Input: Numeric entry for number of obs to view ----
          numericInput(inputId = "l1",
                       label = "l1 (Norm Penalty):",
                       value = 0.7),



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

          #verbatimTextOutput("oid1"),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                           tags$div(id = "loading-message",
                                    style = "position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);
                                                    z-index: 10000; background-color: rgba(255, 255, 255, 0.7);
                                                    padding: 20px; border-radius: 10px;",
                                    tags$p("Loading...", style = "font-size: 18px; font-weight: bold; text-align: center;"))),

          # downloadButton("downloadData", "Download meta-analysis results"),
          conditionalPanel(condition = "input.execute > 0",
                           h4('Meta-analysis results'),
                           withSpinner(DT::dataTableOutput("final_head")),
          ),


          # downloadButton("downloadData", "Download meta-analysis results"),
          conditionalPanel(condition = "input.execute > 0",
                           downloadButton("downloadData", "Download meta-analysis results")),

          # downloadButton("downloadData", "Download meta-analysis results"),
          conditionalPanel(condition = "input.ea_option =='YES' && input.execute > 0",
                           h4('Enrichment analysis results'),
                           withSpinner( DT::dataTableOutput("ea_results")),
          ),


          # downloadButton("downloadEAData", "Download Enrichment analysis results"),
          # downloadButton("downloadEAData", "Download Enrichment analysis results",
          #                conditionalPanel(condition = "input.execute > 0")),

          conditionalPanel(condition = "input.ea_option =='YES' && input.execute > 0",

                           downloadButton("downloadEAData", "Download Enrichment analysis results")),
          conditionalPanel(condition = "input.execute > 0",


                           h4('Plots'),
                           tabsetPanel(
                             tabPanel("Network",
                                      withSpinner(visNetworkOutput('network')),
                                      downloadButton("downloadPlot", "Download Network"),



                             ),
                             tabPanel("Volcano Plot",
                                      withSpinner(plotlyOutput("volc_plot"))
                             ),




                           ),

                           conditionalPanel(condition = "input.ea_option =='YES'",
                                            h4('Manhattan Plot (Enrichment Analysis)'),

                                            withSpinner(plotlyOutput("manhattan_plot")))
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
  ),
  tabPanel(
    "Help",
    fluidPage(
      # Your Help tab content goes here
      tags$h2("Instructions for using MetaPcor"),
      br(),
       br(),

      tags$p(
        "The Analysis tab contains several sidebar elements that allow you to customize your analysis. Below are explanations of each element:",
        tags$ol(
          tags$li(
            tags$b("Upload GSE Matrix file(s):"),
            "Use this to upload your GSE Matrix file(s) for analysis. You can upload multiple files if needed."
          ),
          br(),

          tags$li(
            tags$b("Choose an option:"),
            "Select one of the following options for your analysis:",
            tags$ul(
              tags$li(
                "Partial correlation meta-analysis: Conduct a meta-analysis of partial correlation coefficients."
              ),
              br(),

              tags$li(
                "Partial correlation meta-analysis with thresholds: Perform a meta-analysis of partial correlation coefficients while applying thresholds for correlation coefficients, p-values, or FDR."
              ),
              br(),

              tags$li(
                "DE Partial correlation meta-analysis (for case-control studies): First, perform Differential Expression analysis and then conduct a meta-analysis of partial correlation coefficients for the DEGs."
              )
            )
          ),
          br(),

          tags$li(
            tags$b("Choose meta-analysis model:"),
            "Select a meta-analysis model from the available options:",
            tags$ul(
              tags$li(
                "Random Effects model: Use a random effects model for the meta-analysis."
              ),

              tags$li(
                "Fixed Effects model: Use a fixed effects model for the meta-analysis."
              )
            )
          ),
          br(),

          tags$li(
            tags$b("Enrichment Analysis:"),
            "Toggle this option to enable or disable enrichment analysis. Enrichment analysis results will be available if enabled."
          ),
          br(),

          tags$li(
            tags$b("Volcano plot options:"),
            "When the Volcano plot is generated, you can customize it using the following options:",
            tags$ul(
              tags$li(
                "Coefficient threshold for Volcano plot: Set a threshold for coefficients in the Volcano plot."
              ),
              tags$li(
                "p-value threshold for Volcano plot: Define a threshold for p-values in the Volcano plot."
              )
            )
          ),
          br(),

          tags$li(
            tags$b("l1 (Norm Penalty):"),
            "Specify the 'l1' norm penalty value for your analysis."
          ),
          br(),

        )
      ),
      br(),

      tags$p(
        "Make sure to configure these settings in the sidebar panel before executing your analysis."
      )
    )
  )
)






