# Create HeatMap function
createHeatMap <- function(mat,
                          annot = NULL,
                          title = NULL,
                          krow = 4,
                          kcol = 2,
                          limAbs = ceiling(max(abs(mat))),
                          heatCol = colorRamp2(
                            quantile(-limAbs:limAbs,seq(0,1, 0.10)),
                            rev(brewer.pal(11, "RdBu"))),
                          rowCluster = TRUE,
                          colCluster = TRUE,
                          multipalette = TRUE,
                          palette.alpha = 1,
                          bottomAnnotation = FALSE,
                          show_legend_annotation = TRUE,
                          resCluster = FALSE,
                          row_dend_width =  unit(3, "char"),
                          column_dend_height =  unit(3, "char"),
                          clustering_distance_rows = "euclidean",
                          clustering_method_rows = "complete",
                          clustering_distance_columns = "euclidean",
                          clustering_method_columns = "complete",
                          ...) {
  # Description
  # Wrapper function for Heatmap {ComplexHeatmap}
  
  require(ComplexHeatmap, quietly = T)
  require(circlize, quietly = T)
  require(dendextend, quietly = T)
  require(RColorBrewer, quietly = T)
  require(plyr, quietly = T)
  require(wesanderson, quietly = T)
  require(scales, quietly = T)
  
  
  colAnnot <- list()
  
  ha = new("HeatmapAnnotation")
  hb = new("HeatmapAnnotation")
  if(!is.null(annot)) {
    # if(!multipalette) paletteList <- list("grenYll" = c('#4db6ac','#aed581','#dce775','#ffd54f'))
    if(!multipalette) paletteList <- list('Moonrise2' = wes_palette("Moonrise2")[c(1,3,2)])
    if(multipalette)  paletteList <- list(
                                          'Moonrise2' = wes_palette("Moonrise2")[c(1,3,2)],
                                          "Cavalcanti" = wes_palette("Cavalcanti"),
                                          "Moonrise1" = wes_palette("Moonrise1"),
                                          "Darjeeling" = wes_palette("Darjeeling"),
                                          "Royal1" = wes_palette("Royal1"),
                                          "FantasticFox" = wes_palette("FantasticFox"),
                                          "Chevalier" = wes_palette("Chevalier")
                                          )
    
    colList = list()
    for (j in 1:ncol(annot)) {
      gtab <- unique(annot[, j])
      # colJ <- length(paletteList) - (j %% length(paletteList))
      colJ <- j %% length(paletteList)
      if (colJ == 0)
        colJ <- length(paletteList)
      cols <- alpha(colorRampPalette(paletteList[[colJ]])(length(gtab)), palette.alpha)
      names(cols) <- gtab
      colList[[colnames(annot)[j]]] <- cols
    }
    
    ha = HeatmapAnnotation(df = annot,
                           col = colList,
                           show_legend = show_legend_annotation)
  }
  if(bottomAnnotation) hb = ha
  
  resClusterRes <- list()
  if(colCluster) {
    # dendC = as.dist(1-cor(mat, use="complete.obs"))
    # dendC = hclust(dendC, method="ward.2")
    dendC = dist(t(mat), method = clustering_distance_rows)
    dendC = hclust(dendC, method = clustering_method_rows)
    resClusterRes[['dendC']] <- dendC
    dendC = color_branches(dendC, k = krow)
    
  } else {
    dendC = FALSE
  }
  
  if(rowCluster) {
    # dendR = as.dist(1-cor(t(mat), use="complete.obs"))
    # dendR = hclust(dendR, method="ward.D2")
    dendR = dist(mat, method = clustering_distance_columns)
    dendR = hclust(dendR, method = clustering_method_columns)
    resClusterRes[['dendR']] <- dendR
    dendR = color_branches(dendR, k = krow)
    
  } else {
    dendR = FALSE
  }

  if(!rowCluster) dendR = FALSE
  

  res <- Heatmap(mat,
          col = heatCol,
          top_annotation = ha, 
          top_annotation_height = unit(2, "char"),
          bottom_annotation = hb, 
          bottom_annotation_height = unit(2, "char"),
          cluster_rows = dendR, 
          cluster_columns = dendC, 
          row_dend_width =  row_dend_width,
          column_dend_height =  column_dend_height,
          ...)
  
  if(resCluster) { 
    return(resClusterRes)
  } else {
    return(res)
  }
  
  
}


createClusteringReport <- function(x, annot, genes = row.names(x), prefix = "", doValidation = F, multipalette = F) {
  require(WGCNA, quietly = T)
  require(plyr, quietly = T)
  require(pvclust, quietly = T)
  require(RColorBrewer, quietly = T)
  require(cluster, quietly = T)
  require(fpc, quietly = T)
  require(clValid, quietly = T)
  require(RankAggreg, quietly = T)
  require(grid, quietly = T)
  require(wesanderson, quietly = T)
  
  res <- list()
  # Subset of genes
  x <- x[genes, ]
  #print(dim(x))
  
  # Colors for annotation of dendograms
  colAnnot <- apply(annot, 2, function(x) {
    gtab <- c(table(x))
    cols <- colorRampPalette(rev(brewer.pal(5, "RdBu")))(length(gtab))
    mapvalues(x, names(gtab), cols)
  })
  
  # Colors for annotation of dendograms with multipalette
  if(multipalette) {
    colorsPairs <- list("Zissou" = wes_palette("Zissou"),
                        "Chevalier" = wes_palette("Chevalier"),
                        "Moonrise1" = wes_palette("Moonrise1"),
                        'Moonrise2' = wes_palette("Moonrise2"),
                        "Royal1" = wes_palette("Royal1"), 
                        "Cavalcanti" = wes_palette("Cavalcanti")
    )
    
    i <- 1
    colAnnot <- annot
    for (j in 1:ncol(annot)) {
      gtab <- c(table(annot[, j]))
      cols <- colorRampPalette(colorsPairs[[i]])(length(gtab))
      colAnnot[, j] <- mapvalues(annot[, j], names(gtab), cols)
      i <- i + 1
      if(i > length(colorsPairs)) i <- 1
    }
    
  }
  
  # Bootstraping validation
  if(doValidation) {
    # Validation of Clustering Results Columns
    # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
    print("Validation of Clustering Mehods===================")
    xCl <- x
    if(nrow(xCl) > 100) xCl <- xCl[sample(1:nrow(xCl), 100),]  
    result <- clValid(xCl, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                      validation=c("internal","stability"))
    resCl <- getRanksWeights(result)
    # print(res$ranks[,1:3], quote=FALSE)
    CEWScol <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
    print(CEWScol)
    nClustOptColumn <- gsub('.*-','', CEWScol$top.list[1])
    # plot(CEWS)
    print("===================================================")
    
    # Validation of Clustering Results Rows
    # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
    #result <- clValid(t(x), 2:6, clMethods=c("hierarchical"),
    #                  validation=c("internal","stability"))
    #resCl<- getRanksWeights(result)
    #CEWSrow <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
    #nClustOptRow <- gsub('.*-','', CEWSrow$top.list[1])
  }
  
  # Dendograms
  d = as.dist(1-cor(x, use="complete.obs"));
  hc = hclust(d, method="ward.D2")
  hcd = as.dendrogram(hc, hang=-0.1)
  main = paste(prefix, "\nAnnotated dendogram", sep = "")
  op <- par(no.readonly=TRUE)
  plotDendroAndColors(hc, colAnnot, colorHeight = 0.05, 
                      cex.colorLabels = 0.5, cex.dendroLabels = 0.5, main=main)
  par(op)
  
  res$clusterAnnotSamples <- list()
  
  for (kii in 2:nClustOptColumn) {
    res$clusterAnnotSamples[[as.character(kii)]] <- cutree(hc, k = kii)
  }
  
  # Clustering with Bootstrapped p values
  #fit <- pvclust(x, method.hclust="ward.D2",
  #   method.dist="correlation", use.cor="complete.obs",
  #   nboot = 1000)
  #main = paste(prefix, "\nClustering with bootstrapped p-values (n=1000)", sep = "")
  #plot(fit, print.num=F, main = main, cex.pv=0.5, cex = 0.7) # dendogram with p values plot.pvclust
  #pvrect(fit, alpha=.95, type = "geq") 
  #pvsummary = pvpick(fit, alpha=.95, type = "geq")
  
  ## Plotting Cluster Solutions
  # fit <- kmeans(x, nClusters)
  # # Cluster Plot against 1st 2 principal components# vary parameters for most readable graph
  # clusplot(x, fit$cluster, color=TRUE, shade=TRUE, lines=0)
  # # Centroid Plot against 1st 2 discriminant functions
  # plotcluster(x, fit$cluster) 
  
  
  
  # Heatmap
  hm = createHeatMap(x, annot, kcol = nClustOptColumn)
  print(hm)
  
  return(res)
  
}

createDendogramReport <- function(x, annot, genes = row.names(x), doValidation = F, multipalette = F, addLegend = F, cex.legend = 1, paletteList = NULL, ...) {
  # Description
  # Wrapper function for plotDendroAndColors {WGCNA} -> plot(dendro)
  
  require(WGCNA, quietly = T)
  require(plyr, quietly = T)
  require(pvclust, quietly = T)
  require(RColorBrewer, quietly = T)
  require(cluster, quietly = T)
  require(fpc, quietly = T)
  require(clValid, quietly = T)
  require(RankAggreg, quietly = T)
  require(grid, quietly = T)
  require(wesanderson, quietly = T)
  
  # Setup different default parameters for plotDendroAndColors arguments if not specified in function call
  # if(!exists("cex.colorLabels")) cex.colorLabels = 1
  # if(!exists("cex.dendroLabels")) cex.dendroLabels = 0.7
  if(!exists("colorHeight")) colorHeight = 0.15
  
  res <- list()
  
  # Subset of genes
  x <- x[genes, ]
  #print(dim(x))
  
  
  
  # Colors for annotation of dendograms
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
  
  colList = list()
  for (j in 1:ncol(annot)) {
    gtab <- unique(annot[, j])
    colJ = length(paletteList) - (j %% length(paletteList))
    cols <- colorRampPalette(paletteList[[colJ]])(length(gtab))
    names(cols) = gtab
    colList[[colnames(annot)[j]]] = cols
  }

  res$colList = colList
  colAnnot = annot
  #for (nm in names(colAnnot)){
  #  colAnnot[[nm]] = colList[[nm]][annot[[nm]]]
  #}
  for (nm in colnames(colAnnot)){
    colAnnot[,nm] = colList[[nm]][annot[,nm]]
  }
  
  # Validation of Clustering Results Columns
  # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
  if(doValidation) {
    # print("Validation of Clustering Mehods===================")
    xCl <- x
    if(nrow(xCl) > 100) xCl <- xCl[sample(1:nrow(xCl), 100),]  
    result <- clValid(xCl, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                      validation=c("internal","stability"))
    resCl <- getRanksWeights(result)
    # print(res$ranks[,1:3], quote=FALSE)
    CEWScol <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
    # print(CEWScol)
    nClustOptColumn <- gsub('.*-','', CEWScol$top.list[1])
    # plot(CEWS)
    # print("===================================================")
  }
  # Validation of Clustering Results Rows
  # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
  #result <- clValid(t(x), 2:6, clMethods=c("hierarchical"),
  #                  validation=c("internal","stability"))
  #resCl<- getRanksWeights(result)
  #CEWSrow <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
  #nClustOptRow <- gsub('.*-','', CEWSrow$top.list[1])
  
  # Dendograms
  d = as.dist(1-cor(x, use="complete.obs"));
  hc = hclust(d, method="ward.D2")
  hcd = as.dendrogram(hc, hang=-0.1)
  
  
  if(!addLegend) {
    op <- par(no.readonly=TRUE)
    plotDendroAndColors(hc, colAnnot, autoColorHeight = F, ...)
    par(op)
  }
  if(addLegend) {
    opar <- par(no.readonly=TRUE)
    parMar0 <- par()$mar
    layout(matrix(c(1:4), 2, 2), heights = c(1 - colorHeight, colorHeight), widths = c(1 - 0.25, 0.25))
    plotDendroAndColors(hc, colAnnot, 
                        autoColorHeight = F,
                        marAll = c(1, 5, 3, 0),
                        setLayout = FALSE, ...)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    lNames = gsub ("\\.", " ", names(unlist(colList)))
    legend("center", legend=lNames, fill = unlist(colList), bty = "n", cex = cex.legend)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    par(mar=parMar0)
  }
  

  if(doValidation) {
    res$clusterAnnotSamples <- list()
    for (kii in 2:nClustOptColumn) {
      res$clusterAnnotSamples[[as.character(kii)]] <- cutree(hc, k = kii)
    }
    res$nClustOptColumn <- nClustOptColumn
  } else {
    if(!doValidation) res$nClustOptColumn <- 2
  }
  
  # Clustering with Bootstrapped p values
  #fit <- pvclust(x, method.hclust="ward.D2",
  #   method.dist="correlation", use.cor="complete.obs",
  #   nboot = 1000)
  #main = paste(prefix, "\nClustering with bootstrapped p-values (n=1000)", sep = "")
  #plot(fit, print.num=F, main = main, cex.pv=0.5, cex = 0.7) # dendogram with p values plot.pvclust
  #pvrect(fit, alpha=.95, type = "geq") 
  #pvsummary = pvpick(fit, alpha=.95, type = "geq")
  
  ## Plotting Cluster Solutions
  # fit <- kmeans(x, nClusters)
  # # Cluster Plot against 1st 2 principal components# vary parameters for most readable graph
  # clusplot(x, fit$cluster, color=TRUE, shade=TRUE, lines=0)
  # # Centroid Plot against 1st 2 discriminant functions
  # plotcluster(x, fit$cluster) 
  
  
  return(res)
  
}

createHeatmapReport <- function(x, annot=NULL, genes = row.names(x), prefix = "", k = 2, doValidation = T) {
  require(WGCNA, quietly = T)
  require(plyr, quietly = T)
  require(pvclust, quietly = T)
  require(RColorBrewer, quietly = T)
  require(cluster, quietly = T)
  require(fpc, quietly = T)
  require(clValid, quietly = T)
  require(RankAggreg, quietly = T)
  require(grid, quietly = T)
  require(wesanderson, quietly = T)
  
  res <- list()
  # Subset of genes
  x <- x[genes, ]
  #print(dim(x))
  
  # Colors for annotation of dendograms
  if(!is.null(annot)) {
    colAnnot <- apply(annot, 2, function(x) {
      gtab <- c(table(x))
      cols <- colorRampPalette(rev(brewer.pal(5, "RdBu")))(length(gtab))
      mapvalues(x, names(gtab), cols)
    })
  }
  
  # Validation of Clustering Results Columns
  # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
  if(doValidation) {
    print("Validation of Clustering Mehods===================")
    result <- clValid(x, 2:6, clMethods=c("hierarchical","kmeans","pam"),
                      validation=c("internal","stability"))
    resCl <- getRanksWeights(result)
    # print(res$ranks[,1:3], quote=FALSE)
    CEWScol <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
    print(CEWScol)
    nClustOptColumn <- gsub('.*-','', CEWScol$top.list[1])
    # plot(CEWS)
    print("===================================================")
    
    # Validation of Clustering Results Rows
    # http://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
    #result <- clValid(t(x), 2:6, clMethods=c("hierarchical"),
    #                  validation=c("internal","stability"))
    #resCl<- getRanksWeights(result)
    #CEWSrow <- RankAggreg(x=resCl$ranks, k=5, weights=resCl$weights, seed=123, verbose=FALSE)
    #nClustOptRow <- gsub('.*-','', CEWSrow$top.list[1])
  }
  
  # Dendograms
  d = as.dist(1-cor(x, use="complete.obs"));
  hc = hclust(d, method="ward.D2")
  hcd = as.dendrogram(hc, hang=-0.1)
  main = paste(prefix, "\nAnnotated dendogram", sep = "")
  op <- par(no.readonly=TRUE)
  plotDendroAndColors(hc, colAnnot, colorHeight = 0.05, 
                      cex.colorLabels = 0.5, cex.dendroLabels = 0.5, main=main)
  par(op)
  
  res$clusterAnnotSamples <- list()
  
  for (kii in 2:nClustOptColumn) {
    res$clusterAnnotSamples[[as.character(kii)]] <- cutree(hc, k = kii)
  }
  
  # Clustering with Bootstrapped p values
  #fit <- pvclust(x, method.hclust="ward.D2",
  #   method.dist="correlation", use.cor="complete.obs",
  #   nboot = 1000)
  #main = paste(prefix, "\nClustering with bootstrapped p-values (n=1000)", sep = "")
  #plot(fit, print.num=F, main = main, cex.pv=0.5, cex = 0.7) # dendogram with p values plot.pvclust
  #pvrect(fit, alpha=.95, type = "geq") 
  #pvsummary = pvpick(fit, alpha=.95, type = "geq")
  
  ## Plotting Cluster Solutions
  # fit <- kmeans(x, nClusters)
  # # Cluster Plot against 1st 2 principal components# vary parameters for most readable graph
  # clusplot(x, fit$cluster, color=TRUE, shade=TRUE, lines=0)
  # # Centroid Plot against 1st 2 discriminant functions
  # plotcluster(x, fit$cluster) 
  
  
  
  # Heatmap
  hm = createHeatMap(x, annot, kcol = nClustOptColumn)
  print(hm)
  
  return(res)
  
}
