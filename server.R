library(shiny)
library(DT)
library(plotly)

source("depMapAnalysis.R")

# Load in the DepMap data (into global vars)
loadFeatherFiles()

startUp <- function(session) {
  hideElement(selector = ".item-loading")
  hide(id = "loading-content-container", anim = TRUE, animType = "fade")
  
  updateSelectizeInput(session, 'myKnockoutGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myMutationAnnotation', choices = mutationAnnotations, server = TRUE, selected = "damaging")
  updateSelectizeInput(session, 'myExpressionGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myLineage', choices = lineages, server = TRUE)
  updateSelectizeInput(session, 'mutationLookup_gene', choices = geneList, server = TRUE)
}

createPlots <- function(input, output, analysisData){
  
  output$plot <- renderPlotly({
    
    plot_ly(analysisData[["Enrichment"]], x = ~EffectSize,
                 y = ~`-log10(p.value)`,
                 text = ~Gene,
                 hovertemplate = paste(
                   "<b>%{text}</b><br>",
                   "%{x}<br>",
                   "%{y}",
                   "<extra></extra>"
                 ),
                 color = ~`-log10(p.value)`,
                 colors = "OrRd",
                 size = ~`-log10(p.value)`,
                 type = "scattergl",
                 mode = "markers",
                 height = 500,
                 marker = list(
                   line = list(
                     color = 'rgb(0, 0, 0)',
                     width = 1
                   )
                 )) %>%
      hide_colorbar() %>%
      layout(xaxis = list(title = "Mean Depedency Score Difference"),
             yaxis = list(title = "-log(p.value)"),
             showlegend = FALSE)
  })
  
  output$cellLinePlot <- renderPlotly({
    
    numInCondition = dim(analysisData[["CellLineInfo"]][["condition"]])[1]
    numInControl = dim(analysisData[["CellLineInfo"]][["control"]])[1]
    
    df <- data.frame("Category" = c("Cell Lines in Condition","Cell Lines in Control"),
                     "Values" = c(numInCondition, numInControl))
    
    p <- plot_ly(df, labels = ~Category,
                 values = ~Values,
                 type = "pie",
                 textinfo = "percent+value")
  })
  
  output$genePlot <- renderPlotly({
    
    geneDep = analysisData[["GeneDependencies"]]
    
    print(colnames(geneDep))
    
    
    p <- plot_ly(geneDep,
                 x = ~`Cell.Line.Group`,
                 y = ~`CTNNB1 (1499)`,
                 text = ~`Cell.Line.Name`,
                 type = 'violin',
                 height = 500,
                 points="all",
                 jitter = 0,
                 # hoverinfo = 'text',
                 box = list(
                   visible = T
                 ),
                 meanline = list(
                   visible = T
                 )
    ) %>% 
      layout(
        yaxis = list(
          title = "DepMap Dependency Score",
          zeroline = T
        ),
        xaxis = list(
          title = "Cell Line Group"
        )
      )
  })
  
}

server <- function(input, output, session) {
  
  startUp(session)
  
  cellLineGroups <- eventReactive(input$runAnalysis, {return(getCellLineGroups(input))})
  analysisData <- eventReactive(cellLineGroups(), {return(performAnalysis(cellLineGroups()))})
  
  observeEvent(analysisData(), {
    updateSelectizeInput(session, 'myPlotGene', choices = analysisData()$Enrichment$Gene, server = TRUE)
    hideElement(selector = ".item-loading")
  })
  
  # ================================== #
  # HTML Page Renderers                #
  # ================================== #
  
  output$about_HTML <- renderUI({
    HTML(paste(readLines("www/src/about.html"), collapse=" "))
  })
  
  # ================================== #
  # PLOTS                              #
  # ================================== #
  
  createPlots(input, output, analysisData())
  
  # ================================== #
  # TABLES                             #
  # ================================== #

  # Generate an HTML table view of the data
  output$mytable = DT::renderDataTable(
    analysisData()$data,
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$conditionCellLines = DT::renderDataTable(
    analysisData()$conditionCellLines,
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$controlCellLines = DT::renderDataTable(
    analysisData()$controlCellLines,
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$mutationsTable = DT::renderDataTable(
    mutationLookupanalysisData(),
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(downloadData$name, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(downloadData$data, file)
    }
  )
  
  mutationLookupData <- eventReactive(input$lookupMutations, {
    
    myQuery <- input$mutationLookup_gene
    geneList <- lapply(myQuery, getGeneName)
    
    data <- mutsInDepMap[mutsInDepMap$Hugo_Symbol %in% geneList, ]
    
    return(data)
    
  })
  
  
  
}
