library(shiny)
library(DT)
library(plotly)

# This file contains functions that perform the data analysis
source("depMapAnalysis.R")

# Load in the DepMap data (into global vars)
loadFeatherFiles()

# FUNCTION
# ========
# Hides the startup loading divs and updates the select input
# options with the data that loaded.
startUp <- function(session) {
  hideElement(selector = ".item-loading")
  hide(id = "loading-content-container", anim = TRUE, animType = "fade")
  
  updateSelectizeInput(session, 'myKnockoutGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myMutationAnnotation', choices = mutationAnnotations, server = TRUE, selected = "damaging")
  updateSelectizeInput(session, 'myExpressionGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myLineage', choices = lineages, server = TRUE)
  updateSelectizeInput(session, 'mutationLookup_query', choices = c(geneList, CCLE_mutations$DepMap_ID), server = TRUE)
}

# ========================================================#
# PLOTS                                                   #
#   - for some reason, these have to be before the server # 
#     object due to some promise default issue            #
# ======================================================= #

# This function just takes some of the code out of the server function
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

# ================================== #
# R Shiny Server Object              #
# ================================== #
server <- function(input, output, session) {
  
  # ================================== #
  # Main Server Logic                  #
  # ================================== #
  
  # Updates inputs with initial options
  startUp(session)
  
  # Note: promise structure so that no action is blocking
  
  # On runAnalysis button press, get cellLineGroups indicator vector
  cellLineGroups <- eventReactive(input$runAnalysis, {return(getCellLineGroups(input, output))})
  # Then run the analysis and return the analysisData object
  analysisData <- eventReactive(cellLineGroups(), {return(performAnalysis(cellLineGroups()))})
  
  # Once the analysisData object is obtained, update the select input options and update the plots
  observeEvent(analysisData(), {
    updateSelectizeInput(session, 'myPlotGene', choices = analysisData()$Enrichment$Gene, server = TRUE)
    
    plotAxisChoices = colnames(analysisData()$Enrichment)[-1] # without first `Gene` column
    
    updateSelectizeInput(session, 'enrichmentPlot_yaxis', selected = "-log10(p.value)", choices = plotAxisChoices, server = TRUE)
    hideElement(selector = ".item-loading")
    createPlots(input, output, analysisData())
  })
  
  # Handles the mutation data query. Returns data table of
  # annotations given the query (either depmap id or gene name)
  mutationLookupData <- eventReactive(input$lookupMutations, {
    query = input$mutationLookup_query
    return(getMutationData(query))
  })
  
  # ================================== #
  # HTML Page Renderers                #
  # ================================== #
  
  output$about_HTML <- renderUI({
    HTML(paste(readLines("www/src/about.html"), collapse=" "))
  })
  
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
  
  # ================================== #
  # FILE DOWNLOADS                     #
  # ================================== #
  # TODO: make this a single function
  
  # Downloads enrichment data table
  output$downloadEnrichment <- downloadHandler(
    filename = function() {
      paste("enrichment_analysis", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["Enrichment"]], file)
    }
  )
  
  # Downloads conditioned cell line annotations
  output$downloadConditioned <- downloadHandler(
    filename = function() {
      paste("condition_cell_lines", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["CellLineInfo"]][["condition"]], file)
    }
  )
  
  # Downloads control cell line annotations
  output$downloadControl <- downloadHandler(
    filename = function() {
      paste("control_cell_lines", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["CellLineInfo"]][["control"]], file)
    }
  )
  
  # Downloads gene-wise dependency (this is a big file and can take a while)
  output$downloadGeneDep <- downloadHandler(
    filename = function() {
      paste("gene_dependency", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(analysisData()[["GeneDependencies"]], file)
    }
  )
  
  # Downloads the mutation data table that the user is searching
  output$downloadMutation <- downloadHandler(
    filename = function() {
      paste("mutation_annotations", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(mutationLookupData(), file)
    }
  )
}
