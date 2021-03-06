library(shiny)
library(shinyjs)
library(DT)
library(plotly)

appCSS <- "
#loading-content-container {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}

#loading-content-message{
  margin: 0;
  position: absolute;
  top: 50%;
  left: 50%;
  -ms-transform: translate(-50%, -50%);
  transform: translate(-50%, -50%);
  text-align: center;
  color: #FFFFFF;
  color: #FFFFFF;
}

#appIcon {
  height: 50px;
}

.no-screen-run{
  position: relative;
  left: 30px;
  margin-top: 30px;
}

.item-loading{
  margin: 0;
  position: relative;
  top: 50%;
  left: 50%;
  -ms-transform: translate(-50%, -50%);
  transform: translate(-50%, -50%);
  text-align: center;
  color: #FFFFFF;
  color: #FFFFFF;
}

.lds-facebook {
  display: inline-block;
  position: relative;
  width: 64px;
  height: 64px;
}

.black div:nth-child(1){
  background: #000;
}

.black div:nth-child(2){
  background: #000;
}

.black div:nth-child(3){
  background: #000;
}

.lds-facebook div {
  display: inline-block;
  position: absolute;
  left: 6px;
  width: 13px;
  background: #fff;
  animation: lds-facebook 1.2s cubic-bezier(0, 0.5, 0.5, 1) infinite;
}
.lds-facebook div:nth-child(1) {
  left: 6px;
  animation-delay: -0.24s;
}
.lds-facebook div:nth-child(2) {
  left: 26px;
  animation-delay: -0.12s;
}
.lds-facebook div:nth-child(3) {
  left: 45px;
  animation-delay: 0;
}
@keyframes lds-facebook {
  0% {
    top: 6px;
    height: 51px;
  }
  50%, 100% {
    top: 19px;
    height: 26px;
  }
}
.right{
  float: right;
}
"

# ================================== #
# R Shiny UI Object                  #
# ================================== #

# R-Shiny UI doesn't really display well for readability.
# TODO: clean up spacing so this is more readable
ui <- fluidPage(theme = "app.css",
  
  useShinyjs(),
  
  # inlineCSS(appCSS),
  
  # Loading Message 
  div(
    id = "loading-content-container",
    div(
      id = "loading-content-message",
      h2("Loading DepMap Data..."),
      h4("This can take about a minute..."),
      HTML('<div class="lds-facebook"><div></div><div></div><div></div></div>')
    )
  ),
  
  # Application title
  h1(
    img(src="appIcon.png", id="appIcon"),
    "DepMapScreener"
  ),
    
  # <-- Begin tabset panel -->
  tabsetPanel(
    tabPanel("Parameters",
             
       splitLayout(
         mainPanel(
               
           selectInput(
             "screenType", "Screen Type",
             c(`Knock Out` = "knockout",
               `Expression` = "expression",
               `Lineage` = "lineage",
               `Custom` = "custom")),
           
           # Knock Out Inputs
           conditionalPanel(
             condition = "input.screenType == 'knockout'",
             selectizeInput("myKnockoutGene",
                            "Gene(s) to Knock Out", 
                            choices = NULL,
                            selected = NULL,
                            multiple = TRUE,
                            options = list(
                              openOnFocus = FALSE,
                              maxOptions = 10,
                              placeholder = "Gene(s)"
                            )
             ),
             selectizeInput("myMutationAnnotation",
                            "Enter Mutation Type(s)", 
                            choices = NULL,
                            selected = NULL,
                            multiple = TRUE,
                            options = list(
                              openOnFocus = FALSE,
                              maxOptions = 10,
                              placeholder = "Mutation Type"
                            )
             )
           ),
           
           # Expression Inputs
           conditionalPanel(
             condition = "input.screenType == 'expression'",
             selectizeInput("myExpressionGene",
                            "Gene of Interest", 
                            choices = NULL,
                            selected = NULL,
                            multiple = FALSE,
                            options = list(
                              openOnFocus = FALSE,
                              maxOptions = 10,
                              placeholder = "Gene Expression"
                            )
             ),
             
             selectInput(
               "populationType", "Population Type",
               c(Top = "top",
                 Bottom = "bottom")
             ),
             
             conditionalPanel(
               condition = "input.populationType == 'top'",
               sliderInput("topPercentile", "Top X%", min=1, max=99, value=10)
             ),
             
             conditionalPanel(
               condition = "input.populationType == 'bottom'",
               sliderInput("botPercentile", "Bottom X%", min=1, max=99, value=10)
             )
           ),
           
           # Lineage Inputs
           conditionalPanel(
             condition = "input.screenType == 'lineage'",
             selectizeInput("myLineage",
                            "Lineage of Interest", 
                            choices = NULL,
                            selected = NULL,
                            multiple = FALSE,
                            options = list(
                              openOnFocus = FALSE,
                              maxOptions = 10
                            )
             )
           ),
           
           # Custom Inputs
           conditionalPanel(
             condition = "input.screenType == 'custom'",
             textAreaInput("custom_condition_IDs",
                            "Condition Group IDs",
                           resize = "vertical",
                           placeholder = "e.g.\nACH-000009\nACH-000015\nACH-000150\nACH-000235\nACH-000253\nACH-000258\n...",
                           height = "400px")
           ),
           
           actionButton("runAnalysis", "Run Virtual Screen"),
           style = "padding-bottom: 300px;",
           width = 12
        ),
         mainPanel(
           htmlOutput("stats1",
                      style = "margin-bottom: 50px"),
           div(
             class = "item-loading",
             HTML('<div class="black lds-facebook"><div></div><div></div><div></div></div>')
           ),
           plotlyOutput("cellLinePlot", height = "600px"),
           width = 12
         ),
        cellWidths = c("30%","70%"),
        
        # necessary for word wrap
        cellArgs = list(style='white-space: normal;') 
       )
    ),
    tabPanel("Enrichment Results", 
        tabsetPanel(
             tabPanel("Plot",
                      inputPanel(
                        selectizeInput("enrichmentPlot_yaxis",
                                       "Y Axis", 
                                       choices = c("-log10(p.value)", "-log10(q.value)"),
                                       selected = "-log10(p.value)",
                                       multiple = FALSE,
                                       options = list(
                                         openOnFocus = FALSE,
                                         maxOptions = 20
                                       )
                        )
                      ),
                      div(
                        class = "no-screen-run",
                        h4("No Current Screen")
                      ),
                      div(
                        class = "item-loading",
                        HTML('<div class="black lds-facebook"><div></div><div></div><div></div></div>')
                      ),
                      plotlyOutput("plot")
             ),
             tabPanel("Table",
               div(
                 class = "no-screen-run",
                 h4("No Current Screen")
               ),
               DT::dataTableOutput("mytable")
             ),
             tabPanel("Info",
               inputPanel(downloadButton("downloadEnrichment", "Download Enrichment Results")),
               strong("Plot Help:"),
               br(),
               em("X-Axis:"),
               p("Cell lines in the condition group have enriched dependence upon
                genes that are skewed left (negative mean difference)."),
               p("The CERES dependency score is scaled such that a negative value
               signifies a cell line is dependent upon that gene. The mean difference
               here is the mean dependency score across cell lines in  the conditioned
               group minus the mean dependency score across cell lines in the control
               group. Thus, a negative value means cell lines in the condition group
               were more dependent upon that particular gene."),
               em("Y-Axis:"),
               p("Typical choices for the y-axis will be -log10(p.value) or -log10(q.value).
               These choices give a standard volcano plot. The q-value is a FDR-corrected
               version of the p-value (Benjamini-Hochberg Correction).")
             )
        )
    ),
    
    tabPanel("Gene-Level",
        tabsetPanel(
          tabPanel("Plot",
             splitLayout(
               mainPanel(
                 
                 selectizeInput("myPlotGene",
                                "Gene to Plot", 
                                choices = NULL,
                                selected = NULL,
                                multiple = FALSE,
                                options = list(
                                  openOnFocus = FALSE,
                                  maxOptions = 10
                                )
                 ),
                 style = "padding-bottom: 450px;",
                 width = 12
               ),
               mainPanel(
                 div(
                   class = "no-screen-run",
                   h4("No Current Screen")
                 ),
                 plotlyOutput("genePlot", height = "600px"),
                 width = 12
               ),
               cellWidths = c("30%","70%"),
               
               # necessary for word wrap
               cellArgs = list(style='white-space: normal;') 
             )
          ),
          # tabPanel("Table",
          #      div(
          #        class = "no-screen-run",
          #        h4("No Current Screen")
          #      ),
          #      DT::dataTableOutput("geneLevelDep")
          # ),
          tabPanel("Info",
              inputPanel(
                downloadButton("downloadGeneDep", "Download Gene-Level Dependencies (~160 MB)")
              ),
              p(
                strong("Warning"),
                br(),
                p("This file is large and preparing it for download will take a few dozen seconds. Please be patient after clicking download.")
              ),
              p(
                strong("Note:"),
                br(),
                p("This file includes the gene-level dependency scores for each cell line. It also includes a column grouping the cell lines into the control or condition group."),
                p("The violin plots are generated by plotting the cell group indicator on the x-axis and a given gene's dependency scores on the y-axis.")
              )
          )
        )
    ),
    
    
    tabPanel("View Cell Lines",
        tabsetPanel(
          tabPanel("Condition Group",
             div(
               class = "no-screen-run",
               h4("No Current Screen")
             ),
             DT::dataTableOutput("conditionCellLines")
          ),
          tabPanel("Control Group",
             div(
               class = "no-screen-run",
               h4("No Current Screen")
             ),
             DT::dataTableOutput("controlCellLines")
          ),
          tabPanel("Info",
             inputPanel(
                downloadButton("downloadConditioned", "Download Condition Group"),
                downloadButton("downloadControl", "Download Control Group")
             )
          )
        )
    ),
    
    tabPanel("Mutation Lookup",
             
             inputPanel(
                   selectizeInput("mutationLookup_query",
                                  "Gene(s) or Cell Line", 
                                  choices = NULL,
                                  selected = NULL,
                                  multiple = TRUE,
                                  options = list(
                                    openOnFocus = FALSE,
                                    maxOptions = 10,
                                    placeholder = "e.g. APC or ACH-000004"
                                  )
                   ),
                   actionButton("lookupMutations", "Search Mutations",
                                style="margin-top: 25px")
             ),
             DT::dataTableOutput("mutationsTable"),
             downloadButton("downloadMutation", "Download Table")
    ),
    tabPanel("About", htmlOutput("about_HTML"))          
  ),
  # <-- end tabset panel -->
  width = 12
)
