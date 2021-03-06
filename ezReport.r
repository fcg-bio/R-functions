## getCountData 
### returns : rawData with a child list of input$rawData with 4 extra elements :
###   1. x : log2 transformation of counts, filtering out non present genes
###   2. xNorm : normalized x using input method in normMethod parameter
###   3. sd : arrays containing the sd by gene
###   4. ord : gene order based on drecreasing values of sd
### input : rawData data
getCountData <- function (rawData, param, normMethod = param$normMethod) 
{
  dataset <- rawData$dataset
  seqAnno <- rawData$seqAnno
  
  if (is.null(rawData$signal)) {
    rawData$signal = ezNorm(rawData$counts, presentFlag = rawData$presentFlag, 
                            method = normMethod)
  }
  
  rawData$design = ezDesignFromDataset(dataset, param)
  rawData$samples = rownames(rawData$design)
  rawData$nSamples = length(rawData$samples)
  rawData$conds = ezConditionsFromDesign(rawData$design, maxFactors = 2)
  rawData$nConds = length(unique(rawData$conds))
  rawData$sampleColors = getSampleColors(rawData$conds)
  
  signal = shiftZeros(getSignal(rawData), param$minSignal)
  presentFlag = rawData$presentFlag
  signalRange = range(signal, na.rm = TRUE)
  log2Signal = log2(signal)
  isPresent = ezPresentFlags(signal, presentFlag = presentFlag, 
                             param = param, isLog = rawData$isLog)
  signalCond = 2^averageColumns(log2Signal, by = rawData$conds)
  isPresentCond = averageColumns(isPresent, by = rawData$conds) >= 0.5
  isPresentStudy = apply(isPresentCond, 1, mean) >= param$isPresentStudyThresh
  
  rawData$signal = signal
  isValid = isPresentStudy
  if (!is.null(seqAnno$IsControl)) {
    isValid = isValid & !seqAnno$IsControl
  }
  
  x = log2(2^log2Signal[isValid, ] + param$bgExpression)
  xNormed = sweep(x, 1, rowMeans(x))
  xSd = apply(x, 1, sd, na.rm = TRUE)
  ord = order(xSd, decreasing = TRUE)
  
  
  rawData$x <- x
  rawData$xNormed <- xNormed
  rawData$xSd <- xSd
  rawData$ord <- ord
  rawData$isPresent <- isPresent
  rawData$isPresentCond <- isPresentCond
  rawData$isPresentStudy <- isPresentStudy
  
  rawData
  
  
}



## plot_dendroAndColors 
### returns : annotated dendrogram plot using WGCNA::plotDendroAndColors function
### input : rawData object returned by getCountData()
plot_dendroAndColors <- function (rawData, topGenesSize = NULL, normalized = FALSE, 
                                  colorHeight = NULL, paletteList = NULL, multipalette = FALSE,
                                  addLegend = F, cex.legend = 1, add.title = "", ...) 
{
  
  x <- rawData$x
  xNormed <- rawData$xNormed
  xSd <- rawData$xSd
  ord <- rawData$ord
  sampleColors <- rawData$sampleColors
  
  # Annotation height
  if(is.null(colorHeight)) colorHeight <- 0.15
  
  # Get annotation from dataset
  annot <- rawData$dataset[,grep("\\[Factor\\]", colnames(rawData$dataset)), drop = FALSE]
  colnames(annot) <- gsub(" ", "",gsub("\\[Factor\\]", "", colnames(annot)))
  
  # Colors for annotation of dendrograms
  if(is.null(paletteList)) {
    if(!multipalette) paletteList <- list("grenYll" = c('#4db6ac','#aed581','#dce775','#ffd54f'))
    if(multipalette)  paletteList <- list("Cavalcanti" = wes_palette("Cavalcanti"),
                                          "Moonrise1" = wes_palette("Moonrise1"),
                                          "Darjeeling" = wes_palette("Darjeeling"),
                                          "Royal1" = wes_palette("Royal1"),
                                          "FantasticFox" = wes_palette("FantasticFox"),
                                          "Chevalier" = wes_palette("Chevalier"),
                                          'Moonrise2' = wes_palette("Moonrise2"))
  }
  
  colList <- list()
  for (j in 1:ncol(annot)) {
    gtab <- unique(annot[, j])
    colJ = length(paletteList) - (j %% length(paletteList))
    cols <- colorRampPalette(paletteList[[colJ]])(length(gtab))
    names(cols) = gtab
    colList[[colnames(annot)[j]]] = cols
  }
  colAnnot <- annot
  for (nm in colnames(colAnnot)){
    colAnnot[,nm] <-  colList[[nm]][annot[,nm]]
  }
  
  
  
  if(is.null(topGenesSize)) {
    title <- "all present genes"
  } else {
    topGenes <- ord[1:min(topGenesSize, length(ord))]
    title <- paste("top", topGenesSize, "genes")
  }
  
  if(normalized) {
    x <- xNormed
    title <- paste(title, "; gene-wise normalized")
  }
  
  if(add.title != "")
    title <- paste(title, add.title)
    
  if(is.null(topGenesSize)) {
    d <- as.dist(1 - cor(x, use = "complete.obs"))
  } else {
    d <- as.dist(1 - cor(x[topGenes, ], use = "complete.obs"))
  }
  
  # hcd <- as.dendrogram(hclust(d, method = "ward.D2"), hang = -0.1)
  hc <- hclust(d, method = "ward.D2")
  # hcd <- colorClusterLabels(hcd, sampleColors)
  
  
  if(!addLegend) {
    op <- par(no.readonly=TRUE)
    plotDendroAndColors(hc, colAnnot, autoColorHeight = F, hang = -0.1, main = title, ...)
    par(op)
  }
  if(addLegend) {
    opar <- par(no.readonly=TRUE)
    parMar0 <- par()$mar
    graphics::layout(matrix(c(1:4), 2, 2), 
                     heights = c(1 - colorHeight, colorHeight), 
                     widths = c(1 - 0.25, 0.25))
    plotDendroAndColors(hc, colAnnot, 
                        autoColorHeight = F,
                        marAll = c(1, 5, 3, 0),
                        setLayout = FALSE, hang = -0.1,
                        main = title, ...)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    lNames = gsub ("\\.", " ", names(unlist(colList)))
    legend("center", legend=lNames, fill = unlist(colList), bty = "n", cex = cex.legend)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    par(mar=parMar0)
    par(opar)
  }
  
  
}



## run_clusterResult 
### returns : list with the elements necessary (normalized counts and cluster results) to create a row-annotated heatmap with ezRun::clusterHeatmap
### input : rawData element of EzDeResult returned by getCountData(), contains the normalized counts

run_clusterResult <- function (rawData, topGenesSize = NULL, nClusters = 6, ...) 
{
  if(is.null(topGenesSize)) {
    # title <- "all present genes"
    signal <- rawData$x
    xNormed <- rawData$xNormed
  } else {
    # title <- paste("top", topGenesSize, "genes")
    topGenes <- rawData$ord[1:min(topGenesSize, length(rawData$ord))]
    signal <- rawData$x[topGenes, ]
    xNormed <- rawData$xNormed[topGenes, ]
  }
  
  # Row Clustering
  
  clusterColors <- c("red", "yellow", "orange", "green", 
                     "blue", "cyan", "#fb9a99", "#6a3d9a", "#33a02c", "#fb9a99")
  clusterResult <- clusterResults(xNormed, nClusters = nClusters, 
                                  clusterColors = clusterColors[1:nClusters])
  
  return(list(xNormed = xNormed, clusterResult = clusterResult))
  
}


## plot_heatmap 
### returns : heatmap plot using ezRun::clusterHeatmap (runs heatmap.2)
### input : list obtained with run_clusterResult (wrapper of ezRun::clusterResults), containing the normalized counts and the clusterResults ouput. rawData element of EzDeResult returned by getCountData(), contain the dataset used for sample annotation
plot_heatmap <- function (clusterResult, rawData, param, condition = NULL, paletteList = NULL, geneId = 'default', ...) 
{
  
  xNormed <- clusterResult$xNormed
  clusterResult <- clusterResult$clusterResult
  seqAnno <- rawData$seqAnno
  
  # Modify rownames
  if(geneId != 'default') {
    use <- names(clusterResult$clusterNumbers)
    use <- make.names(seqAnno[use, geneId], unique=TRUE) # Avoid duplicated names error
    names(clusterResult$clusterNumbers) <- use
    
    use <- rownames(xNormed)
    use <- make.names(seqAnno[use, geneId], unique=TRUE) # Avoid duplicated names error
    rownames(xNormed) <- use
    
  }
  
  # Get annotation from dataset
  annot <- rawData$dataset[,grep("\\[Factor\\]", colnames(rawData$dataset)), drop = FALSE]
  colnames(annot) <- gsub(" ", "",gsub("\\[Factor\\]", "", colnames(annot)))
  
  # Colors for column annotation
  if(is.null(paletteList)) paletteList <- list("Cavalcanti" = wes_palette("Cavalcanti"),
                                               "Moonrise1" = wes_palette("Moonrise1"),
                                               "Darjeeling" = wes_palette("Darjeeling"),
                                               "Royal1" = wes_palette("Royal1"),
                                               "FantasticFox" = wes_palette("FantasticFox"),
                                               "Chevalier" = wes_palette("Chevalier"),
                                               'Moonrise2' = wes_palette("Moonrise2"))
  
  # Assign colots to annotation
  colList <- list()
  for (j in 1:ncol(annot)) {
    gtab <- unique(annot[, j])
    colJ = length(paletteList) - (j %% length(paletteList))
    cols <- colorRampPalette(paletteList[[colJ]])(length(gtab))
    names(cols) = gtab
    colList[[colnames(annot)[j]]] = cols
  }
  colAnnot <- annot
  for (nm in colnames(colAnnot)){
    colAnnot[,nm] <-  colList[[nm]][annot[,nm]]
  }
  
  
  
  # Plot
  colColors <- if(is.null(condition)){ 
    rawData$sampleColors
  } else {
    colAnnot[, condition]
  }
  
  clusterHeatmap(xNormed, param, clusterResult,
                 colColors = colColors, ...)
}




## plot_mds 
### returns : annotated mds plot using condition as color annotation
### input : rawData object returned by getCountData()
plot_mds <- function (data, topGenesSize = NULL, condition = NULL, pointSize = 3, pointAlpha = 1, pointShape = 21, textVjust = 2, textAlpha = 1, textSize = 3, xlab = "Leading logFC dim1", ylab = "Leading logFC dim2") 
{
  if(is.null(topGenesSize)) {
    title <- "all present genes"
    signal <- data$x
  } else {
    title <- paste("top", topGenesSize, "genes")
    topGenes <- data$ord[1:min(topGenesSize, length(data$ord))]
    signal <- data$x[topGenes, ]
  }
  
  # Get annotation from dataset
  annot <- data$dataset[,grep("\\[Factor\\]", colnames(data$dataset)), drop = FALSE]
  colnames(annot) <- gsub(" ", "",gsub("\\[Factor\\]", "", colnames(annot)))
  
  
  # use ggplot
  sampleColors = data$sampleColors
  main = title
  requireNamespace("edgeR")
  y = DGEList(counts = signal, group = colnames(signal))
  pdf(file = NULL)
  mds = plotMDS(y, top = 500)
  dev.off()
  
  d <- data.frame(x = mds$x, y = mds$y, annot)
  
  if(is.null(condition)) {
    res <- ggplot(d, aes(x, y, label = rownames(d)))
    res <- res + geom_point(size = pointSize, alpha = pointAlpha, colour = "white", shape = pointShape, stroke = 0)
    res <- res + geom_text(check_overlap = FALSE, vjust = textVjust, size = textSize)
  } else {
    res <- ggplot(d, aes(x, y, label = rownames(d)), aes_string(color = condition))
    res <- res + geom_point(size = pointSize, alpha = pointAlpha, colour = "white", shape = pointShape, aes_string(fill = condition), stroke = 0)
    res <- res + geom_text(check_overlap = FALSE, vjust = textVjust, alpha = textAlpha, size = textSize, aes_string(color = condition))
  }
  res <- res + xlab(xlab)
  res <- res + ylab(ylab)
  res <- res + xlim(1.2 *  min(d$x), 1.2 * max(d$x))
  res <- res + ylim(1.2 * min(d$y), 1.2 * max(d$y))
  res <- res + ggtitle(title)
  res <- res + theme_bw()
  
  res
  
}


## filterRawData
### returns : list with filtered rawData list
### input : list from rawData as found in EzDeResult object and param containing the input filters from UI
# param$selectedSamples <- rownames(rawData$dataset)[1:9]
# rawData <- inputData$rawData
filterRawData <- function(rawData, param) {
  objects1 <- c('counts', 'presentFlag', 'rpkm' ,'tpm', 'x', 'xNormed','isPresent')
  objects2 <- c('dataset')
  objects3 <- c('seqAnno')
  objects4 <- c('xSd', 'ord')
  
  # Filter Samples
  if(!is.null(param$selectedSamples)) {
    for (i in objects1) {
      if(any(i %in% names(rawData)))
        rawData[[i]] <- rawData[[i]][, param$selectedSamples, drop = FALSE]
    }
    for (i in objects2) {
      if(any(i %in% names(rawData)))
        rawData[[i]] <- rawData[[i]][param$selectedSamples,]
    }
  }
  
  # Filter Genes
  if(!is.null(param$selectedGenes )) {
    if(param$selectedGenesId == 'Ensembl') {
      for (i in c(objects1, objects3)) {
        if(any(i %in% names(rawData))) {
          use <- param$selectedGenes[param$selectedGenes %in% rownames(rawData[[i]])]
          if(length(use) > 0)
            rawData[[i]] <- rawData[[i]][use, , drop = FALSE]
        }
      }
      for (i in c(objects4)) {
        if(any(i %in% names(rawData))) {
          use <- param$selectedGenes[param$selectedGenes %in% names(rawData[[i]])]
          if(length(use) > 0)
            rawData[[i]] <- rawData[[i]][use]
        }
      }
    } else if (param$selectedGenesId == 'GeneName') {
      use <- names(rawData$seqAnno$gene_name[rawData$seqAnno$gene_name %in% param$selectedGenes])
      for (i in c(objects1, objects3)) {
        if(length(use) > 0 & any(i %in% names(rawData)))
          rawData[[i]] <- rawData[[i]][use, , drop = FALSE]
      }
    }
  }
  rawData
  
}


## writeNgsTwoGroupReportSim
### Simplified version of ezRun::writeNgsTwoGroupReport
writeNgsTwoGroupReportSim <- function (dataset, deResult, output, htmlFile = "00index.html", 
                                       types = NULL) 
{
  param = deResult$param
  rawData = deResult$rawData
  result = deResult$result
  seqAnno = rawData$seqAnno
  titles = list()
  titles[["Analysis"]] = paste("Analysis:", param$name)
  doc = openBsdocReport(title = titles[[length(titles)]])
  addDataset(doc, dataset, param)
  addCountResultSummary(doc, param, result)
  titles[["Result Summary"]] = "Result Summary"
  addTitle(doc, titles[[length(titles)]], 2, id = titles[[length(titles)]])
  settings = character()
  settings["Number of features:"] = length(result$pValue)
  if (!is.null(result$isPresentProbe)) {
    settings["Number of features with counts above threshold:"] = sum(result$isPresentProbe)
  }
  addFlexTable(doc, ezGrid(settings, add.rownames = TRUE))
  titles[["Number of significants by p-value and fold-change"]] = "Number of significants by p-value and fold-change"
  addTitle(doc, titles[[length(titles)]], 3, id = titles[[length(titles)]])
  addSignificantCounts(doc, result)
  resultFile = addResultFile(doc, param, result, rawData)
  ezWrite.table(result$sf, file = "scalingfactors.txt", head = "Name", 
                digits = 4)
  resultObjFile = paste0("result--", param$comparison, "--", 
                         ezRandomString(length = 12), "--EzDeResult.RData")
  deResult$saveToFile(resultObjFile)
  # addParagraph(doc, ezLink(paste0("http://fgcz-176.uzh.ch/shiny/fgcz_rnaSeqInteractiveReport_app/?data=", 
  #                                 file.path(output$getColumn("Report"), resultObjFile)), 
  #                          "Explore result interactively", target = "_blank"))
  logSignal = log2(shiftZeros(result$xNorm, param$minSignal))
  result$groupMeans = cbind(rowMeans(logSignal[, param$grouping == 
                                                 param$sampleGroup, drop = FALSE]), rowMeans(logSignal[, 
                                                                                                       param$grouping == param$refGroup, drop = FALSE]))
  colnames(result$groupMeans) = c(param$sampleGroup, param$refGroup)
  if (param$writeScatterPlots) {
    testScatterTitles = addTestScatterPlots(doc, param, 
                                            logSignal, result, seqAnno, resultFile$resultFile, 
                                            types)
    titles = append(titles, testScatterTitles)
  }
  use = result$pValue < param$pValueHighlightThresh & abs(result$log2Ratio) > 
    param$log2RatioHighlightThresh & result$usedInTest
  use[is.na(use)] = FALSE
  if (sum(use) > param$maxGenesForClustering) {
    use[use] = rank(result$pValue[use], ties.method = "max") <= 
      param$maxGenesForClustering
  }
  if (sum(use) > param$minGenesForClustering) {
    xCentered = (logSignal - rowMeans(logSignal))[use, order(param$grouping)]
    sampleColors = getSampleColors(param$grouping)[order(param$grouping)]
    clusterPng = "cluster-heatmap.png"
    clusterColors = c("red", "yellow", "orange", "green", 
                      "blue", "cyan")
    clusterResult = clusterResults(xCentered, nClusters = 6, 
                                   clusterColors = clusterColors)
    plotCmd = expression({
      clusterHeatmap(xCentered, param, clusterResult, 
                     file = clusterPng, colColors = sampleColors, 
                     lim = c(-param$logColorRange, param$logColorRange))
    })
    clusterLink = ezImageFileLink(plotCmd, file = clusterPng, 
                                  width = max(800, 400 + 10 * ncol(xCentered)), height = 1000)
    if (doGo(param, seqAnno)) {
      clusterResult = goClusterResults(xCentered, param, 
                                       clusterResult, seqAnno = seqAnno, universeProbeIds = rownames(seqAnno)[result$isPresentProbe])
    }
    resultLoaded = ezRead.table(resultFile$resultFile)
    resultLoaded$Cluster = clusterResult$clusterColors[clusterResult$clusterNumbers[rownames(resultLoaded)]]
    ezWrite.table(resultLoaded, file = resultFile$resultFile)
    if (param$doZip) {
      zipFile(resultFile$resultFile)
    }
    titles[["Clustering of Significant Features"]] = "Clustering of Significant Features"
    addTitle(doc, titles[[length(titles)]], 2, id = titles[[length(titles)]])
    paragraphs = character()
    paragraphs["Significance threshold:"] = param$pValueHighlightThresh
    if (param$log2RatioHighlightThresh > 0) {
      paragraphs["log2 Ratio threshold:"] = param$log2RatioHighlightThresh
    }
    paragraphs["Number of significant features:"] = sum(use)
    addFlexTable(doc, ezGrid(paragraphs, add.rownames = TRUE))
    jsFile = system.file("extdata/enrichr.js", package = "ezRun", 
                         mustWork = TRUE)
    addJavascript(doc, jsFile)
    if (!is.null(clusterResult$GO)) {
      goTables = goClusterTable(param, clusterResult, 
                                seqAnno)
      if (any(c(grepl("Homo_", getOrganism(param$ezRef)), 
                grepl("Mus_", getOrganism(param$ezRef))))) {
        goAndEnrichr = cbind(goTables$linkTable, goTables$enrichrTable)
      }
      else {
        goAndEnrichr = goTables$linkTable
      }
      goAndEnrichrFt = ezFlexTable(goAndEnrichr, border = 2, 
                                   header.columns = TRUE, add.rownames = TRUE)
      bgColors = rep(gsub("FF$", "", clusterResult$clusterColors))
      goAndEnrichrFt = setFlexTableBackgroundColors(goAndEnrichrFt, 
                                                    j = 1, colors = bgColors)
      goAndEnrichrTableLink = as.html(ezGrid(rbind("Background color corresponds to the row colors in the heatmap plot.", 
                                                   as.html(goAndEnrichrFt))))
      goLink = as.html(ezGrid(rbind("Background color corresponds to the row colors in the heatmap plot.", 
                                    as.html(goTables$ft))))
    }
    else {
      goAndEnrichrTableLink = as.html(pot("No information available"))
      goLink = as.html(pot("No information available"))
    }
    goClusterTableDoc = openBsdocReport("GO Cluster tables")
    tbl = ezGrid(cbind(`Cluster Plot` = clusterLink, `GO categories of feature clusters` = goLink), 
                 header.columns = TRUE)
    addFlexTable(goClusterTableDoc, tbl)
    closeBsdocReport(goClusterTableDoc, "goClusterTable.html")
    addParagraph(doc, pot("GO cluster tables", hyperlink = "goClusterTable.html"))
    tbl = ezGrid(cbind(`Cluster Plot` = clusterLink, `GO categories of feature clusters` = goAndEnrichrTableLink), 
                 header.columns = TRUE)
    addFlexTable(doc, tbl)
  }
  if (doGo(param, seqAnno)) {
    goResult = twoGroupsGO(param, result, seqAnno, normalizedAvgSignal = rowMeans(result$groupMeans), 
                           method = param$goseqMethod)
    titles[["GO Enrichment Analysis"]] = "GO Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 2, id = titles[[length(titles)]])
    revigoTitle = addGoUpDownResult(doc, param, goResult)
    titles = append(titles, revigoTitle)
    if (any(c(grepl("Homo_", getOrganism(param$ezRef)), 
              grepl("Mus_", getOrganism(param$ezRef))))) {
      isSig = result$pValue < param$pValThreshGO & result$usedInTest
      isUp = result$log2Ratio > param$log2RatioThreshGO & 
        isSig
      isDown = result$log2Ratio < -param$log2RatioThreshGO & 
        isSig
      regulatedGenes = list()
      regulatedGenes$upGenes = na.omit(unique(seqAnno[isUp, 
                                                      "gene_name"]))
      regulatedGenes$downGenes = na.omit(unique(seqAnno[isDown, 
                                                        "gene_name"]))
      regulatedGenes$bothGenes = union(regulatedGenes$upGenes, 
                                       regulatedGenes$downGenes)
      enrichrLinks = ezMatrix("", rows = c("bothGenes", 
                                           "downGenes", "upGenes"), cols = 1)
      for (row in rownames(enrichrLinks)) {
        genesToUse = seqAnno$gene_name[which(seqAnno$gene_name %in% 
                                               regulatedGenes[[row]])]
        genesList = paste(genesToUse, collapse = "\\n")
        jsCall = paste0("enrich({list: \"", genesList, 
                        "\", popup: true});")
        enrichrLinks[row, 1] = as.html(pot(paste0("<a href='javascript:void(0)' onClick='", 
                                                  jsCall, "'>Enrichr</a>")))
      }
      titles[["Enrichr"]] = "Enrichr"
      addTitle(doc, titles[[length(titles)]], 3, id = titles[[length(titles)]])
      addFlexTable(doc, ezFlexTable(enrichrLinks, valign = "middle", 
                                    add.rownames = TRUE))
    }
  }
  if (param[["GAGEanalysis"]]) {
    gageRes = runGageAnalysis(result, param = param, output = output, 
                              rawData = rawData)
    titles[["GAGE Enrichment Analysis"]] = "GAGE Enrichment Analysis"
    addTitle(doc, titles[[length(titles)]], 3, id = titles[[length(titles)]])
    addGageTables(doc, param, gageRes)
  }
  addParagraph(doc, pot("advanced Plots", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}



