library(shiny)
library(DT)
library(plotly)


source("depMapAnalysis.R")

geneList <<- scan("depmapData_feather/geneList.csv",
                  what="character",
                  sep=",")

lineages <<- scan("depmapData_feather/lineages.csv",
                  what="character",
                  sep=",")

depGeneList <<- scan("depmapData_feather/depGenes.csv",
                     what="character",
                     sep=",")

mutationAnnotations <<- scan("depmapData_feather/mutationAnnotations.csv",
                               what="character",
                               sep=",")

# Load in the DepMap data (into global vars)
# loadFeatherFiles()


server <- function(input, output, session) {
  
  hideElement(selector = ".item-loading")
  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content-container", anim = TRUE, animType = "fade")

  updateSelectizeInput(session, 'myKnockoutGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myMutationAnnotation', choices = mutationAnnotations, server = TRUE, selected = "damaging")
  updateSelectizeInput(session, 'myExpressionGene', choices = geneList, server = TRUE)
  updateSelectizeInput(session, 'myLineage', choices = lineages, server = TRUE)
  updateSelectizeInput(session, 'mutationLookup_gene', choices = geneList, server = TRUE)

  # Reactive expression to generate the requested distribution. This is
  # called whenever the inputs change. The renderers defined
  # below then all use the value computed from this expression
  data <- eventReactive(input$runAnalysis, {

    screenType <- input$screenType
    
    showElement(selector = ".item-loading")
    hideElement(selector = ".no-screen-run")
    
    if(screenType == 'knockout'){
      return(knockOut(input$myKnockoutGene,
                      input$myMutationAnnotation))
    }

    else if(screenType == 'expression'){
      
      if (input$populationType == "top"){
        percentile <- 1 - (input$topPercentile / 100)
      } else {
        percentile <- input$botPercentile / 100
      }
      
      print(percentile)
      
      return(expressionScreen(input$myExpressionGene,
                              percentile,
                              input$populationType))
    }

    else { # (screenType == 'lineage')
      return(compareLineage(input$myLineage))
    }

  })
  
  geneLevelData <- eventReactive(input$runGenePlot, {
                                 
        myGene <- input$myPlotGene
        print(myGene)
        effectVec <- data()$effectVec
        
        x <- ifelse(effectVec, "Condition", "Control")
        y <- Achilles_gene_effect[[myGene]]
        
        cellLine_ids <- Achilles_gene_effect$X1
        
        translation <- data.frame("DepMap_ID" = cellLine_ids)
        
        mini <- cellLines[,c("DepMap_ID", "CCLE_Name")]
        
        translation <- left_join(translation, mini, by = "DepMap_ID")
        
        name <- translation$CCLE_Name
        
        data <- data.frame("Cell.Line.Group" = x, "Dependency.Score" = y, "Cell.Line" = name)
        
        return(data)
        
  })
  
  mutationLookupData <- eventReactive(input$lookupMutations, {
    
    myQuery <- input$mutationLookup_gene
    geneList <- lapply(myQuery, getGeneName)
  
    data <- mutsInDepMap[mutsInDepMap$Hugo_Symbol %in% geneList, ]
    
    return(data)
    
  })
  
  observeEvent(data(), {
    updateSelectizeInput(session, 'myPlotGene', choices = data()$data$Gene, server = TRUE)
    hideElement(selector = ".item-loading")
  })

  # Generate a plot of the data. Also uses the inputs to build the
  # plot label. Note that the dependencies on both the inputs and
  # the 'data' reactive expression are both tracked, and all expressions
  # are called in the sequence implied by the dependency graph
  output$plot <- renderPlotly({
    
    df <- data.frame(`x.Effect.Size` = data()$data$EffectSize,
                     `y.neg.log10.p.value` = -log10(data()$data$p.value),
                     `text.Gene.Name` = data()$data$Gene)
    downloadData$data <<- df
    downloadData$name <<- "volcano_plot_data"
    
    p <- plot_ly(data()$data, x = ~EffectSize,
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
  

  # # Generate a summary of the data
  output$stats1 <- renderUI({
    HTML(paste("<b><p style='word-wrap: break-word; color: red'> ", sprintf(data()$status), "</p></b>", sep=""))
  })
  
  output$about_HTML <- renderUI({
    HTML(paste(readLines("www/src/about.html"), collapse=" "))
  })
  
  output$cellLinePlot <- renderPlotly({
    df <- data.frame("Category" = c("Cell Lines in Condition", "Cell Lines in Control"),
                     "Values" = c(data()$stats$trial, data()$stats$control))
    
    p <- plot_ly(df, labels = ~Category, 
                 values = ~Values,
                 type = "pie",
                 textinfo = "percent+value")
  })
  
  output$genePlot <- renderPlotly({
    
    downloadData$data <<- geneLevelData()
    downloadData$name <<- "gene_level_data"
    
    p <- plot_ly(geneLevelData(),
                 x = ~`Cell.Line.Group`,
                 y = ~`Dependency.Score`,
                 text = ~`Cell.Line`,
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

  # Generate an HTML table view of the data
  output$mytable = DT::renderDataTable(
    data()$data,
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$conditionCellLines = DT::renderDataTable(
    data()$conditionCellLines,
    extensions = 'FixedColumns',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE
    )
  )
  
  output$controlCellLines = DT::renderDataTable(
    data()$controlCellLines,
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
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(downloadData$name, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(downloadData$data, file)
    }
  )
  
  
  
}
