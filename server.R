library(shiny)
library(DT)
library(plotly)

source("depMapAnalysis.R")

# Load in the DepMap data (into global vars)
# loadFeatherFiles()

startUp <- function(session) {
  hideElement(selector = ".item-loading")
  hide(id = "loading-content-container", anim = TRUE, animType = "fade")
  
  updateSelectizeInput(session, 'myKnockoutGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myMutationAnnotation', choices = mutationAnnotations, server = TRUE, selected = "damaging")
  updateSelectizeInput(session, 'myExpressionGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myLineage', choices = lineages, server = TRUE)
  updateSelectizeInput(session, 'mutationLookup_query', choices = c(geneList, CCLE_mutations$DepMap_ID), server = TRUE)
}

createPlots <- function(input, output, analysisData){
  
  output$plot <- renderPlotly({
    
    plot_ly(analysisData[["Enrichment"]],
                 x = ~`EffectSize`,
                 y = ~get(input$enrichmentPlot_yaxis),
                 text = ~Gene,
                 hovertemplate = paste(
                   "<b>%{text}</b><br>",
                   "%{x}<br>",
                   "%{y}",
                   "<extra></extra>"
                 ),
                 color = ~get(input$enrichmentPlot_yaxis),
                 colors = "OrRd",
                 size = ~get(input$enrichmentPlot_yaxis),
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
      layout(xaxis = list(title = "Mean Dependency Score Difference"),
             yaxis = list(title = input$enrichmentPlot_yaxis),
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
    
    p <- plot_ly(geneDep,
                 x = ~`Cell.Line.Group`,
                 y = ~get(input$myPlotGene),
                 # split = ~`Cell.Line.Group`,
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
  
  cellLineGroups <- eventReactive(input$runAnalysis, {return(getCellLineGroups(input, output))})
  analysisData <- eventReactive(cellLineGroups(), {return(performAnalysis(cellLineGroups()))})
  
  observeEvent(analysisData(), {
    updateSelectizeInput(session, 'myPlotGene', choices = analysisData()$Enrichment$Gene, server = TRUE)
    
    plotAxisChoices = colnames(analysisData()$Enrichment)[-1] # without first `Gene` column
    
    updateSelectizeInput(session, 'enrichmentPlot_yaxis', selected = "-log10(p.value)", choices = plotAxisChoices, server = TRUE)
    hideElement(selector = ".item-loading")
    createPlots(input, output, analysisData())
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
  
  # createPlots(input, output, analysisData())
  
  # ================================== #
  # TABLES                             #
  # ================================== #

  # Generate an HTML table view of the data
  output$mytable = DT::renderDataTable(
    analysisData()[["Enrichment"]],
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$conditionCellLines = DT::renderDataTable(
    analysisData()[["CellLineInfo"]][["condition"]],
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$controlCellLines = DT::renderDataTable(
    analysisData()[["CellLineInfo"]][["control"]],
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$mutationsTable = DT::renderDataTable(
    mutationLookupData(),
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$geneLevelDep = DT::renderDataTable(
    analysisData()[["GeneDependencies"]],
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  # Downloadable csv of enrichment dataset ----
  output$downloadEnrichment <- downloadHandler(
    filename = function() {
      paste("enrichment_analysis", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["Enrichment"]], file)
    }
  )
  
  # Downloadable csv of conditioned cell lines dataset ----
  output$downloadConditioned <- downloadHandler(
    filename = function() {
      paste("condition_cell_lines", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["CellLineInfo"]][["condition"]], file)
    }
  )
  
  # Downloadable csv of control cell lines dataset ----
  output$downloadControl <- downloadHandler(
    filename = function() {
      paste("control_cell_lines", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["CellLineInfo"]][["control"]], file)
    }
  )
  
  # Downloadable csv of control cell lines dataset ----
  output$downloadGeneDep <- downloadHandler(
    filename = function() {
      paste("gene_dependency", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["GeneDependencies"]], file)
    }
  )
  
  output$downloadMutation <- downloadHandler(
    filename = function() {
      paste("mutation_annotations", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(mutationLookupData(), file)
    }
  )
  
  mutationLookupData <- eventReactive(input$lookupMutations, {
    
    query = input$mutationLookup_query
    
    # if the query is a cell line ID (depmap ID)...
    if (startsWith(query, "ACH-")){ 
      
      print("Cell Line")

      data <- mutsInDepMap[mutsInDepMap$DepMap_ID %in% query, ]
      
    # otherwise, it is a gene
    } else {
      
      print("Gene")
      
      geneList <- lapply(query, getGeneName)
      data <- mutsInDepMap[mutsInDepMap$Hugo_Symbol %in% geneList, ]

    }
    
    return(data)
    
  })
  
  
  
}
