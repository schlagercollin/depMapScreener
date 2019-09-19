library(limma)
library(ggplot2)

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

getCellLineGroups <- function(input){
  
  showElement(selector = ".item-loading")
  hideElement(selector = ".no-screen-run")
  
  screenType = input$screenType
  print("Generating cell line groups")
  
  if(screenType == 'knockout'){
    return(knockOutScreen(input$myKnockoutGene, input$myMutationAnnotation))
  }
  else if(screenType == 'expression'){
    
    if (input$populationType == "top"){
      percentile <- 1 - (as.numeric(input$topPercentile) / 100)
    } else {
      percentile <- (as.numeric(input$botPercentile) / 100)
    }
    
    print(percentile)
    
    return(expressionScreen(input$myExpressionGene, input$populationType,  percentile))
  }
  else { # (screenType == 'lineage')
    return(lineageScreen(input$myLineage))
  }
};

doEnrichmentAnalysis <- function(depMatrix, cellLineGroups){
  
  enrichmentResults <- run_lm_stats_limma(depMatrix, cellLineGroups,
                                          covars = NULL, weights = NULL,
                                          target_type = 'Gene')
  enrichmentResults$`-log10(p.value)` = -log10(enrichmentResults$`p.value`)
  
  return(enrichmentResults)
}

getCellLineInfo <- function(cellLineGroups){
  conditionCellLines <- Achilles_gene_effect[cellLineGroups, ]$X1
  conditionCellLines_info <- cellLines[cellLines$DepMap_ID %in% conditionCellLines, ]
  
  controlCellLines <- Achilles_gene_effect[!cellLineGroups, ]$X1
  controlCellLines_info <- cellLines[cellLines$DepMap_ID %in% controlCellLines, ]
  
  return(list(condition = conditionCellLines_info, control = controlCellLines_info))
}

convertDepMapIDToCCLE <- function(depmapIDs){
  translation <- data.frame("DepMap_ID" = depmapIDs)
  mini <- cellLines[,c("DepMap_ID", "CCLE_Name")]
  translation <- left_join(translation, mini, by = "DepMap_ID")
  return(translation$CCLE_Name)
}

getGeneDependencies <- function(cellLineGroups){
  
  cellLineGroup <- ifelse(cellLineGroups, "Condition", "Control")
  geneDependencies <- copy(Achilles_gene_effect)
  
  oldColumns = colnames(geneDependencies)
  
  # convert DepMapIDs to CCLE_IDs
  depmapIDs <- Achilles_gene_effect$X1
  cellLineNames = convertDepMapIDToCCLE(depmapIDs)
  
  geneDependencies$Cell.Line.Group = cellLineGroup
  geneDependencies$Cell.Line.Name  = cellLineNames
  
  newColumns = c("Cell.Line.Name", "Cell.Line.Group", oldColumns)
  
  # reorder so that cell line labels are at the beginning
  geneDependencies = geneDependencies[newColumns]
  
  return(geneDependencies)
  
}

performAnalysis <- function(cellLineGroups){
  
  #TODO add error when a group has zero cell lines
  
  cellline_info = getCellLineInfo(cellLineGroups)
  enrichment_analysis = doEnrichmentAnalysis(depMatrix, cellLineGroups)
  gene_dependencies = getGeneDependencies(cellLineGroups)
  
  analysisData = list(CellLineInfo = cellline_info,
                      Enrichment = enrichment_analysis,
                      GeneDependencies = gene_dependencies)
  
  return(analysisData)
}

run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene', limma_trend = FALSE) {
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    #helper function for converting two-sided p-values to one-sided p-values
    one_sided_p <- two_sided_p / 2
    if (test_dir == 'right') {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
    } else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
    }
    return(one_sided_p)
  }
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                             p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                             q.left = p.adjust(p.left, method = 'BH'),
                             q.right = p.adjust(p.right, method = 'BH'))
  return(results)
}

getGeneName <- function(geneStr){
  return(strsplit(geneStr, " ")[[1]][1])
}

getGeneNames <- function(geneVector){
  return(lapply(geneVector, getGeneName))
}

getMutations <- function(geneName, mutationAnnotations = "damaging"){
  
  allMutationsOfGene <- CCLE_mutations[CCLE_mutations$Hugo_Symbol == geneName, ]
  subsetMutations <- allMutationsOfGene[allMutationsOfGene$Variant_annotation %in% mutationAnnotations, ]
  
  return(subsetMutations)
}

getCellLinesWithMutations <- function(geneNames, mutationTypes){
  
  geneMutations = lapply(geneNames, function(geneName){
    getMutations(geneName, mutationTypes)
  })
  
  cellLineIDs = lapply(geneMutations, function(x){
    x$DepMap_ID
  })
  
  cellLinesWithAllMutations = Reduce(intersect, cellLineIDs)
  return(cellLinesWithAllMutations)
}

# ================================================= #
# Screen Type Functions                             #
#   - these functions return a cellLineGroup vector #
#      for the given screen type                    #
# ================================================= #

knockOutScreen <- function(genes, mutationTypes){
  
  geneNames <- lapply(genes, getGeneName)
  geneMutations <- lapply(geneNames, function(x){
    return(getMutations(x, mutationTypes))
  })
  cellLineIDs <- lapply(geneMutations, function(x){
    return(x$DepMap_ID)
  })

  commonCellLinesWithMutations <- Reduce(intersect, cellLineIDs)

  effectVec <- Achilles_gene_effect$X1 %in% commonCellLinesWithMutations
  
  print(sum(effectVec))
  
  geneNames = getGeneNames(genes)
  cellLinesWithMutations = getCellLinesWithMutations(geneNames, mutationTypes)
  cellLineGroups = Achilles_gene_effect$X1 %in% cellLinesWithMutations
  
  print(sum(cellLineGroups))
  
  return(cellLineGroups)
}

expressionScreen <-function(gene, populationType, percentile){
  
  geneExpression = CCLE_expression[[gene]]
  
  percentileValue = quantile(geneExpression, c(percentile), names = FALSE)
  
  if (populationType == "top"){
    selectedCellLineIndices = geneExpression >= percentileValue
  } else {
    selectedCellLineIndices = geneExpression <= percentileValue
  }
  
  selectedCellLines = CCLE_expression[selectedCellLineIndices, ]$X1
  
  cellLineGroups = Achilles_gene_effect$X1 %in% selectedCellLines
  
  return(cellLineGroups)
}

lineageScreen <- function(lineage){
  
  cellLineNames = celllineinfo[celllineinfo$Lineage == lineage,]$Name
  cellLineIDs   = cellLines[cellLines$CCLE_Name %in% cellLineNames, ]$DepMap_ID
  
  cellLineGroups = Achilles_gene_effect$X1 %in% cellLineIDs
  
  return(cellLineGroups)
  
}

violinPlot <- function(effectVec, gene, trueLab = TRUE, falseLab = FALSE){
  
  x <- ifelse(colorectal$effectVec, trueLab, falseLab)
  y <- Achilles_gene_effect[[gene]]
  
  return(ggplot(data = Achilles_gene_effect, aes(x=x, y=`CTNNB1 (1499)`)) + 
           geom_violin())
}

createFeatherFiles <- function(){
  
  library(readr)
  
  celllineinfo <- read_csv("/Volumes/Rohatgi_CRISPR_Drive/DepMap Exploration/analysisScripts/depmapData/celllineinfo.csv",
                           progress = FALSE)
  CCLE_mutations <- read_csv("/Volumes/Rohatgi_CRISPR_Drive/DepMap Exploration/analysisScripts/depmapData/CCLE_mutations.csv",
                             progress = FALSE)
  Achilles_gene_effect <- read_csv("/Volumes/Rohatgi_CRISPR_Drive/DepMap Exploration/analysisScripts/depmapData/Achilles_gene_effect.csv",
                                   progress = FALSE)
  CCLE_expression <- read_csv("/Volumes/Rohatgi_CRISPR_Drive/DepMap Exploration/analysisScripts/depmapData/CCLE_expression_full.csv",
                              progress = FALSE)
  cellLines <- read_csv("/Volumes/Rohatgi_CRISPR_Drive/DepMap Exploration/analysisScripts/depmapData/DepMap-2019q1-celllines_v2.csv",
                        progress = FALSE)
  
  
  write_feather(celllineinfo, "celllineinfo.feather")
  write_feather(CCLE_mutations, "CCLE_mutations.feather")
  write_feather(CCLE_expression, "CCLE_expression.feather")
  write_feather(cellLines, "cellLines.feather")
  write_feather(Achilles_gene_effect, "Achilles_gene_effect.feather")
}

loadFeatherFiles <- function(){
  
  library(feather)
  
  print("STARTING LOAD")
  
  celllineinfo <<- read_feather("depmapData_feather/celllineinfo.feather")
  CCLE_mutations <<- read_feather("depmapData_feather/CCLE_mutations.feather")
  Achilles_gene_effect <<- read_feather("depmapData_feather/Achilles_gene_effect.feather")
  CCLE_expression <<- read_feather("depmapData_feather/CCLE_expression.feather")
  cellLines <<- read_feather("depmapData_feather/cellLines.feather")
  
  depMatrix <<- Achilles_gene_effect[, !names(Achilles_gene_effect) %in% c("X1")]
  mutsInDepMap <<- CCLE_mutations[CCLE_mutations$DepMap_ID %in% Achilles_gene_effect$X1, ]
  
  print("COMPLETE")
}

# OLD FUNCTIONS BELOW

# knockOut <- function(geneVec, mutationAnnotations, thresh = 4.6){
# 
#   geneNames <- lapply(geneVec, getGeneName)
#   geneMutations <- lapply(geneNames, function(x){
#     return(getMutations(x, mutationAnnotations))
#   })
#   cellLineIDs <- lapply(geneMutations, function(x){
#     return(x$DepMap_ID)
#   })
# 
#   commonCellLinesWithMutations <- Reduce(intersect, cellLineIDs)
# 
#   effectVec <- Achilles_gene_effect$X1 %in% commonCellLinesWithMutations
#   title <- paste(geneNames, collapse = "+", sep="")
#   title <- paste("knockOut_", title, "_", sep="")
#   mutationAnnotations_str <- paste(mutationAnnotations, collapse = "+", sep="")
#   title <- paste(title, "with_", mutationAnnotations_str, sep="")
#   results <- analyzeDifference(depMatrix, effectVec, title = title)
# 
#   return(results)
# }
# 
# expressionScreen <- function(gene, percentile, side){
#   
#   expression <- CCLE_expression[[gene]]
#   percentile_value <- quantile(expression, c(percentile), names = FALSE)
#   
#   if (side == "top"){
#     selectedIndices <- expression >= percentile_value
#   } else { # bottom
#     selectedIndices <- expression <= percentile_value
#   }
#   selectedCellLines <- CCLE_expression[selectedIndices, ]$X1
#   
#   effectVec <- Achilles_gene_effect$X1 %in% selectedCellLines
#   
#   title = paste("expression_", percentile, "_", gene, sep="")
#   
#   results <- analyzeDifference(depMatrix, effectVec, title = title)
#   
#   return(results)
# }
# 
# compareLineage <- function(lineage, thresh=4.6){
#   
#   cellLinesOfInterest_name <- celllineinfo[celllineinfo$Lineage == lineage,]$Name
#   cellLinesOfInterest_id <- cellLines[cellLines$CCLE_Name %in% cellLinesOfInterest_name,]
#   
#   effectVec <- Achilles_gene_effect$X1 %in% cellLinesOfInterest_id$DepMap_ID
#   
#   title <- paste("lineage_", lineage, sep="")
#   
#   results <- analyzeDifference(depMatrix, effectVec, title= title, thresh = thresh)
#   
#   return(results)
# }
