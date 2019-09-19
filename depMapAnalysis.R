library(limma)
library(ggplot2)

# ================================================= #
# Functions that Perform Data Analysis              #
# ================================================= #

# FUNCTION
# ========
# Group cell types for the eventual enrichment analysis.
# @param
#     - input: input object from server.R
#     - output: ouput object from server.R
# @return
#     - indicator vector denoting which group each cell type
#         in the gene dependency data (Achilles dataset) belongs
getCellLineGroups <- function(input, output){
  
  # UI stuff to update as this is working
  showElement(selector = ".item-loading")
  hideElement(selector = ".no-screen-run")
  
  screenType = input$screenType
  print("Generating cell line groups...")
  print(screenType)
  
  # Process the input data differently for each screen type
  if(screenType == 'knockout'){
    return(knockOutScreen(input$myKnockoutGene, input$myMutationAnnotation))
  }
  else if(screenType == 'expression'){
    
    if (input$populationType == "top"){
      percentile <- 1 - (as.numeric(input$topPercentile) / 100)
    } else {
      percentile <- (as.numeric(input$botPercentile) / 100)
    }
    return(expressionScreen(input$myExpressionGene, input$populationType,  percentile))
  }
  else if(screenType == 'lineage'){
    return(lineageScreen(input$myLineage))
  }
  else { # (screenType == 'custom')
    cellLineIDs = processCellLineList(input$custom_condition_IDs)
    unmatched = cellLineIDs[!(cellLineIDs %in% Achilles_gene_effect$X1)]
    cellLineGroups = customScreen(cellLineIDs)
    
    # Display any IDs that didn't match in the dataset
    if(length(unmatched) > 0){
      showModal(modalDialog(title = "Unmatched IDs", paste(unmatched ,collapse=" ")))
    }
    return(cellLineGroups)
  }
}

# FUNCTION
# ========
# Wrapper for the entirety of this analysis. This is what is called
# from server.R. The analysisData object has data produced in its
# parameters.
# @param
#   - cellLineGroups: indicator vector for cell line group
# @return
#   - list with the various analysis data bundled up:
#         - CellLineInfo <list with annotation information about cell lines by group>
#               - condition <condition cell lines only>
#               - control   <control cell lines only>
#         - Enrichment <data table returned from enrichment analysis>
#         - GeneDependencies <gene-level dependencies data table>
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

# FUNCTION
# ========
# Processes custom condition cell line string into a vector
# of cell line ids. Removes whitespace and then separates 
# by the newline character. 
# @param
#     - id_string: string with depmap ids separated by \n
# @return
#     - vector of cell line IDs split up
processCellLineList <- function(id_string){
  
  id_string_no_whitespace = gsub("[[:blank:]]", "", id_string)
  cellLineIDs = strsplit(id_string_no_whitespace, "\n")
  cellLineIDs = unlist(cellLineIDs)
  print(cellLineIDs)
  return(cellLineIDs)
  
}

# FUNCTION
# ========
# Wrapper function for the enrichment analysis (which is
# really performed in `run_lm_stats_limma`). Adds a few
# columns to the analysis (log-transformed for volcano plots).
# @param
#     - depMatrix (this is a global variable)
#     - cellLineGroups: indicator vector for cell line groups
# @return
#     - data table with enrichment results
doEnrichmentAnalysis <- function(depMatrix, cellLineGroups){
  
  enrichmentResults <- run_lm_stats_limma(depMatrix, cellLineGroups,
                                          covars = NULL, weights = NULL,
                                          target_type = 'Gene')
  enrichmentResults$`-log10(p.value)` = -log10(enrichmentResults$`p.value`)
  enrichmentResults$`-log10(q.value)` = -log10(enrichmentResults$`q.value`)
  
  return(enrichmentResults)
}

# FUNCTION
# ========
# Get cell line info and split the data into two data tables, one for
# the control group and the other for the condition group.
# @param
#     - cellLineGroups: indicator vector for cell groups
# @return
#     - list object with $condition and $control datatables
getCellLineInfo <- function(cellLineGroups){
  conditionCellLines <- Achilles_gene_effect[cellLineGroups, ]$X1
  conditionCellLines_info <- cellLines[cellLines$DepMap_ID %in% conditionCellLines, ]
  
  controlCellLines <- Achilles_gene_effect[!cellLineGroups, ]$X1
  controlCellLines_info <- cellLines[cellLines$DepMap_ID %in% controlCellLines, ]
  
  return(list(condition = conditionCellLines_info, control = controlCellLines_info))
}

# FUNCTION
# ========
# Converts vector of DepMapIDs (e.g. ACH-000009) to one with CCLE_Names
# The CCLE_Names are more human-readable (e.g. HEC251_ENDOMETRIUM)
convertDepMapIDToCCLE <- function(depmapIDs){
  translation <- data.frame("DepMap_ID" = depmapIDs)
  mini <- cellLines[,c("DepMap_ID", "CCLE_Name")]
  translation <- left_join(translation, mini, by = "DepMap_ID")
  return(translation$CCLE_Name)
}

# FUNCTION
# ========
# Returns a data table with gene-wise dependency for each cell line
# AND includes information about condition vs control group for each
# cell line. Basically, this just re-formats the raw Achilles_gene_effect
# data table into one that has columns for the cell line groups.
# @param
#   - cellLineGroups: indicator vector for cell line group
# @return
#   - data table with cell lines as rows, genes as cols, CERES dep score as data
getGeneDependencies <- function(cellLineGroups){
  
  cellLineGroup <- ifelse(cellLineGroups, "Condition", "Control")
  geneDependencies <- Achilles_gene_effect
  
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

# FUNCTION
# ========
# Perform enrichment analysis. Code from Chen et al, Nature 2019.
# @params
#   - mat: gene dependency matrix with cell lines as rows and genes as cols
#   - vec: indicator vector (1 or 0) for cell line group across cell lines
# @return
#   - data table with enrichment analysis statistics for each gene
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

# FUNCTION
# ========
# Extract just the gene name from a string
# Ex: "APC (ENSG00294)" ==> "APC"
getGeneName <- function(geneStr){
  return(strsplit(geneStr, " ")[[1]][1])
}

# FUNCTION
# ========
# Wrapper function to perform gene name extraction
# for each string in a vector of strings.
getGeneNames <- function(geneVector){
  return(lapply(geneVector, getGeneName))
}

# FUNCTION
# ========
# Get the mutations given a particular gene and type(s) of mutation
# annotations.
getMutations <- function(geneName, mutationAnnotations = "damaging"){
  
  allMutationsOfGene <- CCLE_mutations[CCLE_mutations$Hugo_Symbol == geneName, ]
  subsetMutations <- allMutationsOfGene[allMutationsOfGene$Variant_annotation %in% mutationAnnotations, ]
  
  return(subsetMutations)
}

# FUNCTION
# ========
# Returns cell line IDs (dep map IDs) for a set of genes and a set
# of mutations. Note: only cell lines that meet all conditions 
# will be returned (the intersection of the different gene mutations
# sets).
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
#      for the given screen type (indicator vec)    #
# ================================================= #

# FUNCTION
# ========
# Groups cell lines by gene(s) that have certain
# mutation type(s).Returns indicator vector for 
# cell line groups.
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

# FUNCTION
# ========
# Groups cell lines based on expression.
# E.g. cell lines with top 10% expression for APC
# Returns indicator vector of cell groups.
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

# FUNCTION
# ========
# Groups cell lines based on original lineage.
# E.g. Colorectal
lineageScreen <- function(lineage){
  
  cellLineNames = celllineinfo[celllineinfo$Lineage == lineage,]$Name
  cellLineIDs   = cellLines[cellLines$CCLE_Name %in% cellLineNames, ]$DepMap_ID
  
  cellLineGroups = Achilles_gene_effect$X1 %in% cellLineIDs
  
  return(cellLineGroups)
  
}

# FUNCTION
# ========
# User-defined condition group. 
# @params
#   - cellLineIDs: string with DepMapIDs separated by \n
# @return
#   - indicator vector for cell line groups
customScreen <- function(cellLineIDs){
  
  cellLineGroups = Achilles_gene_effect$X1 %in% cellLineIDs
  
  return(cellLineGroups)
  
}

# FUNCTION
# ========
# Retrieve mutation annotation data from a given query.
# Query can be a DepMap cell type ID (e.g. ACH-00009)
# or a gene name.
# @return
#   - data table with the right rows based on query
getMutationData <- function(query){
  # if the query is a cell line ID (depmap ID)...
  if (startsWith(query, "ACH-")){ 
    data <- mutsInDepMap[mutsInDepMap$DepMap_ID %in% query, ]
    
  # otherwise, it is a gene
  } else {
    geneList <- lapply(query, getGeneName)
    data <- mutsInDepMap[mutsInDepMap$Hugo_Symbol %in% geneList, ]
  }
  return(data)
}


# ================================================= #
# Load Some Data into Global                        #
#      - TODO: re-arrange so this isn't necessary   #
# ================================================= #
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

# ================================================= #
# DepMap Data Management Functions                  #
#   - these functions load the DepMap data into     #
#       global variables                            #
# ================================================= #

# FUNCTION
# ========
# Creates feather files from raw DepMap files. Use this to update the 
# database to a new dataset from DepMap. 
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

# FUNCTION
# ========
# Loads in DepMap data into global vars. Uses feather format because it's substantially
# faster than loading the raw .csv files.
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
