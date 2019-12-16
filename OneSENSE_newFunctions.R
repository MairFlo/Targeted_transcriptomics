#All OneSENSE scripts prepared by Evan Newell and Timothy Bi
#Contact: Fred Hutchinson Cancer Research Center

# Evan Newell 2015, Ming 2016, Etienne Becht 2016, Evan 2017, Timothy Bi 2018




#### Main Functions (Algorithm Handlers and Plotters) --------------------------------------------------
FCStSNEone <- function (LoaderPATH ="fcs", #Folder with samples (even 1 sample ok)
                        useCSV = F,
                        TransformOutput = T,
                        ceil = 20, #number of events to take per sample (unless sample has less)
                        FNnames="names2.csv", #Parameters to analyze in csv file
                        OutputSuffix = "out2",
                        DotSNE = T,
                        tSNEperplexity = 30, #default = 30; increased perplexity => increased spread
                        DoOneSENSE = F,
                        DoPhenograph = F,
                        kValue = 10, #default = 30; k is low => more clusters and inversely
                        DoFlowSOM = F,
                        MaxClusters = 30,
                        DoIsomap = F,
                        DoDiffMap = F,
                        ManDist = F,
                        DoPCA = F,
                        Do3DtSNE =F,
                        DoUMAP = F,
                        Do3DUMAP = F,
                        DoOneSUMAP = F,
                        DoManualGrouping = F,
                        GroupingName = "TissueType",
                        GroupingOutput = "GroupingOutput",
                        markerCap = 0) {
  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = useCSV)
  
  FFdata<- prepData$FFdata              #getting data
  OrigNames <- prepData$OrigNames       #
  forNewFF <- prepData$forNewFF         #
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  
  keeptable <- read.csv(FNnames, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  data <- FFdata[,which (colnames(FFdata) %in% as.character(keeprows[,1]))]
  
  ## Names file directed transformation including autotransformation
  
  #transType <- as.character(keeptable[keeprowbool,2])
  
  nfTransOut <- nfTransform(keeprows,data,FFdata, doNormalization = doNormalization)
  ####  Run algorithms
  dataTransformed <- nfTransOut$dataA1
  FFdataTransformed <- nfTransOut$dataB1
  
  includeBool = sapply(keeptable[,3], function(x) !x=='x')
  if(!is.na(includeBool[1])) {
    includeNames <-  keeptable[includeBool,1]
    includeNames <- intersect(includeNames,keeprows[,1])
    dataTransformed <- dataTransformed[,includeNames]
    print("did this")
  }
  
  #setting all matrices to null
  score <- NULL
  tSNEmat <- NULL
  OneStSNEmat <- NULL
  usordata <- NULL
  DiffMapMat <- NULL
  TdtSNEmat <- NULL
  umapMat <- NULL
  TdumapMat <- NULL
  OneSUMAPmat <- NULL
  Rphenographmat <- NULL
  FlowSOMmat <- NULL
  groupingMat <- NULL
  sparseScore <- NULL
  
  #running algorithms
  if(usingSeuratData) {sparseScore <- pca(dataTransformed)}
  if(DoPCA) {
    if (is.null(sparseScore)) {
      score <- pca(dataTransformed)
    } else {
      score <- sparseScore
    }
  }
  if(DotSNE) {tSNEmat <- tSNE_2D(dataTransformed, tSNEperplexity, FFdata, sparseScore)}
  if(DoOneSENSE) {OneStSNEmat <- tSNE_oneSENSE(keeptable, FFdata, tSNEperplexity= 30, markerCap = markerCap)}
  if(DoIsomap) {usordata <- runIsomap(dataTransformed)}
  if(DoDiffMap) {DiffMapMat <- diffMap(ManDist, dataTransformed)}
  if(Do3DtSNE) {TdtSNEmat <- tSNE_3D(dataTransformed, tSNEperplexity, sparseScore)}
  
  #set UMAP parameters
  n <- 15
  mdist <- 0.2
  metric <- "euclidean"
  if(DoUMAP) {umapMat <- umap_2D(dataTransformed, n, mdist, metric, ncomp = 2, FFdata, sparseScore)}
  if(Do3DUMAP) {TdumapMat <- umap_3D(dataTransformed, n, mdist, metric, ncomp = 3, sparseScore)}
  if(DoOneSUMAP) {OneSUMAPmat <- umap_oneSENSE(FFdata, keeptable, markerCap = markerCap)}
  
  score2 <- cbind(score, tSNEmat, OneStSNEmat, usordata, TdtSNEmat, DiffMapMat, umapMat, TdumapMat, OneSUMAPmat)
  
  #prepare for export
  if(DoPhenograph) {Rphenographmat <- phenograph(dataTransformed, score2, kValue)}
  if(DoFlowSOM) {FlowSOMmat <- flowSOM(dataTransformed, MaxClusters)}
  if(DoManualGrouping) {
    for(i in 1:length(GroupingName)){
      groupingMat <- cbind(groupingMat, getGrouping(FFdata, FcsFileNames, GroupingName[i], GroupingOutput))
      colnames(groupingMat)[i] <- GroupingName[i]
    }
  }
  nontransformedData <- cbind(Rphenographmat,FlowSOMmat,groupingMat)
  ## 1 existing data - needs tranfrom
  ## 2 new data exported needs tranfrom
  ## 3 new data exported should be used as is
  
  normalizedScore <- reverseTransform(score2 = score2,
                                      nontransformedData = nontransformedData,
                                      FFdata = FFdata,
                                      TransformOutput = TransformOutput)
  NIscore <- normalizedScore[[1]]
  NewTransformBase <- normalizedScore[[2]]
  NewTransformNew <- normalizedScore[[3]]
  NewTransform <- normalizedScore[[4]]
  
  # Export Data to FCS file
  exportToFCS(useCSV = useCSV,
              numFiles = numFiles,
              forNewFF = forNewFF,
              OrigNames = OrigNames,
              FFdata = FFdata,
              NIscore = NIscore,
              NewTransform = NewTransform,
              NewTransformNew = NewTransformNew,
              LoaderPATH = LoaderPATH,
              OutputSuffix = OutputSuffix,
              FcsFileNames = FcsFileNames)
}

# Meaning Plots
meaningPlot <- function (LoaderPATH =sourceFcsForMeaningPlot,
                         useCSV = F, # T if using .csv files instead of fcs files - will make csv output also
                         ceil = 100,
                         FNnames= "names.csv",
                         TransformOutput = T,
                         MeaningTopPercentile = .99,
                         MeaningBotPercentile = .01,
                         PC1 = "tSNE1",
                         PC2 = "tSNE2",
                         DoOneForEach = F,
                         prefix = OutputSuffix,
                         palette=c("black","blue","green","yellow","red"),
                         color.scale.type=c("relative","absolute")[1],
                         resolution=72,
                         cex=0.5,
                         pch=16,
                         plotClusters = F,
                         plotInFile = F,
                         plotGroupings = F,
                         GroupingName = "TissueType",
                         GroupingOutput = "GroupingOutput",
                         MPNColors = 0,
                         MPColorList = NULL,
                         clusteringParameter = "Phenograph",
                         labelClusters = F,
                         DoMontagePlot = F,
                         plotDensity = F,
                         DoFileMontage = F,
                         parameterName = F,
                         DoOverlay = F,
                         baseFile = NULL){
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = F)
  
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  filenames <- sapply(strsplit(FcsFileNames, split = "\\."), "[", 1)
  
  keeptable <- read.csv(FNnames, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  suppressWarnings(dir.create("MeaningPlots"))
  LoaderPATH <- paste0("MeaningPlots/", LoaderPATH)
  
  data1u <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]  # existing data - needs transform
  dataForNormalize <- FFdata[,transList == "1"]
  data2u <- FFdata[,transList == "2"]                            # new data exported - needs transform
  data3 <- FFdata[,transList == "3"]                             # new data exported - use as is
  
  #transType <- as.character(keeptable[keeprowbool,2])
  #Transform data1u (using names file info) and data2u (using cytof parameters)
  nfTransOut <- nfTransform(keeprows, data1u, dataForNormalize, doNormalization = doNormalization)
  data1 <- nfTransOut$dataA1
  
  lgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  data2 <- apply(data2u,2,lgcl)
  
  keeptable <- read.csv(FNnames,stringsAsFactors=F, fileEncoding="UTF-8-BOM")
  data <- cbind(data1,data2,data3)
  print (plotDensity)
  if(plotDensity) {
    get_density <- function(x, y, n = 200) {
      dens <- MASS::kde2d(x = x, y = y, n = n)
      ix <- findInterval(x, dens$x)
      iy <- findInterval(y, dens$y)
      ii <- cbind(ix, iy)
      return(dens$z[ii])
    }
    density <- get_density(data[, PC1], data[, PC2])
    data <- cbind(data, Density = density)
  }
  data.range = range(data)
  discreteData <- NULL
  if (plotClusters) {
    clusterIndex <- grep(clusteringParameter, colnames(FFdata))
    clusteringParameter <- colnames(FFdata)[clusterIndex[length(clusterIndex)]]
    discreteData <- cbind(discreteData, FFdata[, clusteringParameter])
    colnames(discreteData) <- clusteringParameter
  }
  if (plotInFile) {
    discreteData <- cbind(discreteData, InFile = FFdata[, "InFile"])
  }
  if (plotGroupings) {
    for (name in GroupingName) {
      discreteData <- cbind(discreteData, Grouping = FFdata[, name])
      colnames(discreteData)[ncol(discreteData)] <- name
    }
  }
  #PCA has been run and appended to the FCS <-? not sure what this applies to
  x <- data[,PC1]
  y <- data[,PC2]
  
  bottoms <- apply(data,2,function(a) quantile(a, MeaningBotPercentile))
  
  bottoms[bottoms < 0] <- 0
  tops <- apply(data,2,function(a) quantile(a, MeaningTopPercentile))
  plotdata <- data[, which (colnames(data) %in% keeprows[,1])]
  if(plotDensity) {
    plotdata <- cbind(plotdata, Density = density)
  }
  
  print("Making concatenated meaning plots")
  rasters <- generateRasters(plotdata, x, y, color.scale.type, bottoms, tops, resolution, PC1, PC2, doAxes = T) # returns list of images
  if(!is.null(discreteData)){
    discreteRasters <- generateDiscreteRasters(discreteData, FFdata, clusteringParameter, MPNColors, PC1, PC2, resolution, MPColorList,
                                               labelClusters, FcsFileNames, GroupingOutput, pch, cex)
  }
  names_table = read.csv(FNnames, fileEncoding="UTF-8-BOM")
  pathName = ""
  last = ""
  count = 0
  for(i in rownames(names_table)) {
    current = names_table[i,2]
    if (last != current & i != rownames(names_table)[1]) {
      if (last != '')
      pathName = paste0(pathName,"_", names_table[as.numeric(i)-count,1], '-', names_table[as.numeric(i)-1,1], "(",last,')')
      count = 0
    } else {
      count = count + 1
    }
    last = current
  }
  if (last != '') {
    pathName = paste0(pathName,"_", names_table[as.numeric(i)-count,1], '-', names_table[as.numeric(i)-1,1], "(",last,')')
  }
  pathName = gsub('[/"]',"'",pathName)
  pathName = paste0("meaningplot")
  print(pathName)
  pltFN <- paste0(LoaderPATH, pathName, ".pdf")
  pdf(pltFN)
  if(color.scale.type=="absolute"){
    plot.new()
    grid.raster(rasters[[1]]$raster.scale)
  }
  
  
  sapply(rasters,function(x){
    par("mar"=c(0,0,0,0))
    if(color.scale.type=="relative"){
      grid.newpage()
      grid.raster(x$raster.main,y=0.6,height=0.8)
      grid.raster(x$raster.scale,y=0.1,height=0.2)
    }
    if(color.scale.type=="absolute"){
      par("mar"=c(0,0,0,0))
      grid.newpage()
      grid.raster(x$raster.main)
    }
  })
  if(!is.null(discreteData)) {
    sapply(discreteRasters,function(x){
      par("mar"=c(0,0,0,0))
      grid.newpage()
      grid.raster(x$raster.main,y=0.6,height=0.8)
      grid.raster(x$raster.legend,y=0.1,height=0.2)
    })
  }
  
  
  grb = tableGrob(names_table)
  
  mytheme <- ttheme_default(
    core = list(fg_params=list(cex = 0.3)),
    colhead = list(fg_params=list(cex = 0.3)),
    rowhead = list(fg_params=list(cex = 0.3)))
  
  d = names_table  ## all of this puts the names table in the pdf for future reference
  
  tg <- tableGrob(d, rows = seq_len(nrow(d)), theme = mytheme) 
  tg$heights = tg$heights
  fullheight <- convertHeight(sum(tg$heights), "cm", valueOnly = TRUE)
  margin <- unit(5,"in")
  margin_cm <- convertHeight(margin, "cm", valueOnly = TRUE)
  a4height <- 29.7 - margin_cm
  nrows <- nrow(tg)
  npages <- ceiling(fullheight / a4height)
  
  heights <- convertHeight(tg$heights, "cm", valueOnly = TRUE) 
  rows <- cut(cumsum(heights), include.lowest = FALSE,
              breaks = c(0, cumsum(rep(a4height, npages))))
  
  groups <- split(seq_len(nrows), rows)
  
  gl <- lapply(groups, function(id) tg[id,])
  blank_table = table(list(c(0,0),c(0,0)))
  
  
  gl = c(gl,rep(list(tableGrob(blank_table)),3-length(gl)%%3))
  
  print (length(gl))
  
  for(page in seq_len(npages/3+1)){
    #grid.newpage()
    grid.arrange(gl[[page*3-2]], gl[[page*3-1]], gl[[page*3]], ncol = 3)
  }
  
  #grid.draw(grb)
  print('done')
  dev.off()
  
  if(DoMontagePlot) { 
    print("Making montage plots")
    rasters <- generateRasters(plotdata, x, y, color.scale.type, bottoms, tops, resolution = resolution, PC1, PC2, doAxes = F)
    pltFN <- paste0(LoaderPATH, "_MontageMP.pdf")
    pdf(pltFN, paper="A4r", width = 11, height = 8.5)
    par("mar"=c(0,0,0,0))
    viewport()
    i <- c(1:length(rasters))
    sapply(i, function(x){
      if ((x-1)%%48 == 0){
        par("mar"=c(0,0,0,0))
        grid.newpage()
      }
      grid.raster(rasters[[x]]$raster.main, x = (((x-1)%%8)*0.125 + 0.05), y = (1 - ((((x-1)%/%8L)%%6)*0.17 + 0.1)), height = 0.17)
      
    })
    dev.off()
  }
  
  if(DoOneForEach) {
    print("Making meaning plots for each file:")
    for(FFs in 1:numFiles)
    {
      print(filenames[FFs])
      plotdataEach <- plotdata[FFdata[,"InFile"]==FFs, , drop = F]
      discreteDataEach <- discreteData[FFdata[, "InFile"] == FFs, ,drop = F]
      if(nrow(plotdataEach == 1) && ncol(discreteDataEach == 1)) {
        colnames(plotdataEach) <- colnames(plotdata)
        colnames(discreteDataEach) <- colnames(discreteData)
      }
      x <- data[FFdata[,"InFile"]==FFs, PC1]
      y <- data[FFdata[,"InFile"]==FFs, PC2]
      
      if(plotDensity && nrow(plotdataEach) != 1) {
        density <- get_density(x, y)
        plotdataEach[, "Density"] <- density
        plotdataEach.range = range(plotdataEach)
      }
      
      rasters <- generateRasters(plotdataEach, x, y, color.scale.type, bottoms, tops, resolution, PC1, PC2, doAxes = T)
      if(!is.null(discreteDataEach)){
        discreteRasters <- generateDiscreteRasters(discreteData, FFdata, clusteringParameter, MPNColors, PC1, PC2, resolution, MPColorList,
                                                   labelClusters, FcsFileNames, GroupingOutput, pch, cex)
      }
      BaseFN <- sapply(strsplit(filenames[FFs], split ="\\."), "[", 1)
      pltFN <- paste0("MeaningPlots/", BaseFN,prefix,".pdf")
      pdf(pltFN)
      if(color.scale.type=="absolute"){
        plot.new()
        grid.raster(rasters[[1]]$raster.scale)
      }
      sapply(rasters,function(x){
        par("mar"=c(0,0,0,0))
        if(color.scale.type=="relative"){
          grid.newpage()
          grid.raster(x$raster.main,y=0.6,height=0.8)
          grid.raster(x$raster.scale,y=0.1,height=0.2)
        }
        if(color.scale.type=="absolute"){
          par("mar"=c(0,0,0,0))
          grid.newpage()
          grid.raster(x$raster.main)
        }
      })
      if(!is.null(discreteData)) {
        sapply(discreteRasters,function(x){
          par("mar"=c(0,0,0,0))
          grid.newpage()
          grid.raster(x$raster.main,y=0.6,height=0.8)
          grid.raster(x$raster.legend,y=0.1,height=0.2)
        })
      }
      dev.off()
    }
  }
  
  if(DoFileMontage) {
    intensity <- data[, parameterName, drop = F]
    colnames(intensity) <- parameterName
    rasters <- sapply (c(1:numFiles), function(file) {
      x <- FFdata[FFdata[, "InFile"] == file, PC1]
      y <- FFdata[FFdata[, "InFile"] == file, PC2]
      if(plotDensity && parameterName == "Density") {
        density <- get_density(x, y)
        data <- as.matrix(density)
        colnames(data) <- "Density"
      } else {
        data <- intensity[FFdata[, "InFile"] == file, ,drop = F]
      }
      return(generateRasters(data, x, y, color.scale.type, bottoms, tops, resolution, PC1, PC2, doAxes = T, nameByFile = T, fileName = FcsFileNames[file]))
    })
    
    pltFN <- paste0(LoaderPATH, "_AllFiles_", parameterName, ".pdf")
    pdf(pltFN)
    i <- c(1:length(rasters))
    split <- ceiling(sqrt(length(rasters)))
    par <- c(0.5, 0.5, 0.5, 0.5)
    sapply(i, function(x){
      grid.raster(rasters[[x]]$raster.main, x = (((x-1)%%split)*(1/split) + 0.5/split),
                  y = (1 - (((x-1) %/% split)/split + 0.5/split)), height = (1/split - 0.05/split))
    })
    dev.off()
  }
  
  if(DoOverlay) {
    baseFileIndex <- grep(baseFile, FcsFileNames)
    baseX <- FFdata[FFdata[, "InFile"] == baseFileIndex, PC1]
    baseY <- FFdata[FFdata[, "InFile"] == baseFileIndex, PC2]
    overlayFileNames <- FcsFileNames[-baseFileIndex]
    
    overlayRasters <- sapply(1:length(overlayFileNames), function(x) {
      overlayX <- FFdata[FFdata[, "InFile"] == x, PC1]
      overlayY <- FFdata[FFdata[, "InFile"] == x, PC2]
      colors <- c("grey")
      
      mainplot=paste(tmpDir(),"/mainplot.png",sep="")
      png(mainplot,res=resolution,height=480*resolution/72,width=480*resolution/72)
      par("bty"="l")
      plot.new()
      plot(baseX, baseY, cex = cex, pch = pch, main = FcsFileNames[x], cex.main = 1,
           xlab = PC1, ylab = PC2, 
           col = "grey",
           axes = T  )
      points(overlayX, overlayY, cex = 1.2*cex, pch = pch, col = "red")
      dev.off()
      
      return(list(raster.main=readPNG(mainplot,native=T)))
    },simplify=F)
    
    pltFN <- paste0(LoaderPATH, "_overlayMP.pdf")
    pdf(pltFN)
    sapply(overlayRasters, function(x) {
      grid.newpage()
      grid.raster(x$raster.main)
    })
    dev.off()
  }
}

fileHeatplot <- function (LoaderPATH="fcs",
                          ceil=100,
                          TransformOutput=T,
                          FNnames="names.csv",
                          OutputSuffix="1",
                          palette=c("black","blue","lightblue","green","yellow","darkorange","darkred"),
                          DoClusters=F,
                          DoPositiveFrequency=F,
                          clustParam="Phenograph",
                          yparam="CD45RA",
                          yparam2="CD4"){  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = F)
  
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  pathForCoords <- LoaderPATH
  par(mar = c(7,4,4,2)+0.1)
  
  #transType <- as.character(keeptable[keeprowbool,2])
  sampID <- NULL
  hpmat <- NULL
  cellnum <- NULL
  
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  
  keeptable <- read.csv(FNnames, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  
  suppressWarnings(dir.create("Heatplots"))
  LoaderPATH <- paste0("Heatplots/", LoaderPATH)
  
  #data1u and data1 are from the names file
  data1u <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])] #Existing data - needs transform
  data2u <- FFdata[,transList == "2"]                           #New data - needs transform
  data3 <- FFdata[,transList == "3"]                            #New data - use as is (Phenograph, groupings, etc.)
  
  nfTransOut <- nfTransform(keeprows, data1u, data1u, doNormalization = doNormalization)
  data1 <- nfTransOut$dataA1
  lgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  data2 <- apply(data2u, 2, lgcl)
  
  data<-cbind(data1,data2,data3)
  for (i in 1:length(FcsFileNames)){
    filedata <- data[FFdata[, "InFile"] == i , which (colnames(data) %in% keeprows[,1]), drop = F]
    mtmdata <- apply(filedata, 2, median)
    cellnum <- cbind(cellnum, dim(data)[1])
    sampID <- c(sampID, rep(i,dim(data)[1]))
    hpmat <- cbind(hpmat, mtmdata)
  }
  
  #If only one file, need to duplicate for heatmap to work
  if (dim(hpmat)[2] == 1){
    hpmat <- cbind(hpmat, hpmat)
    colnames(hpmat) <- c(FcsFileNames, "")
    Rowv = F
    dendrogram = "column"
  } else {
    colnames(hpmat) <- FcsFileNames
    Rowv = T
    dendrogram = "both"
  }
  
  filenames <- strsplit(FcsFileNames, split = "\\.") 
  
  #row.names(keeprows) <- keeptable[keeprowbool,1]
  #hpmat<-hpmat[keeprowbool,]
  
  pdf(file=paste0(LoaderPATH,"Medians",OutputSuffix,".pdf"), width=14, height = 10)
  breaks = seq(0,max(hpmat),by=0.05)
  my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
  heatmap.2(t(hpmat), col=my_palette, 
            breaks = breaks,
            margins = c(20,28), 
            Colv = T,
            Rowv = Rowv,
            dendrogram = dendrogram,
            cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
            density.info=c("none"),
            keysize=1)    
  dev.off()
  
  write.csv(t(hpmat),paste0(LoaderPATH,"MedianValues",OutputSuffix,".csv"))
  
  if(DoClusters) {
    hpmat <- NULL
    for(cluster in 1:max(FFdata[,clustParam])) {
      mtmdata <- apply(data1[FFdata[,clustParam]==cluster,], 2, median)
      hpmat <- cbind(hpmat, mtmdata)
    }
    colnames(hpmat) <- 1:max(FFdata[,clustParam])
    
    pdf(file=paste0(LoaderPATH,"ClusterMedians",OutputSuffix,".pdf"), width=14, height =10)
    breaks = seq(0,max(hpmat),by=0.05)
    my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
    heatmap.2(t(hpmat), col=my_palette, 
              breaks = breaks,
              margins = c(20,10), 
              Colv = T,
              dendrogram = "both",
              cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
              density.info=c("none"),
              keysize=1)    
    dev.off()
    
    write.csv(t(hpmat),paste0(LoaderPATH,"ClusterMedianValues",OutputSuffix,".csv"))
    
    ####  Tabulate frequencies of each cluster across each sample and make heatplot:    
    mtmdata<- NULL
    hpmat <- NULL
    hpmatCnames <- NULL
    for(sample in 1:length(FcsFileNames)){
      clusterData <- data[FFdata[, "InFile"] == sample, ,drop = F]
      cellsInSample <- nrow(clusterData)
      clustFreq <- NULL
      for(cluster in 1:max(FFdata[,clustParam]))
      {
        cellsInSampleInCluster <- nrow(clusterData[clusterData[, clustParam] == cluster, , drop = F]) 
        clustFreq <- c(clustFreq, (cellsInSampleInCluster/cellsInSample)*100)
      }
      hpmat <- cbind(hpmat, clustFreq)
    }
    
    colnames(hpmat) <- c( sapply(filenames, "[", 1))
    row.names(hpmat) <- 1:max(FFdata[,clustParam])
    #hpmat<-hpmat[row.names(keeprows),]
    
    if (dim(hpmat)[2] == 1){
      hpmat <- cbind(hpmat, hpmat)
      colnames(hpmat) <- c(FcsFileNames, "")
    } else {
      colnames(hpmat) <- FcsFileNames
    }
    
    pdf(file=paste0(LoaderPATH,"ClusterFrequencies",OutputSuffix,".pdf"), width=14, height =15)
    breaks = seq(0,max(hpmat),by=max(hpmat)/50)
    my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
    heatmap.2(t(hpmat), col=my_palette, 
              breaks = breaks,
              margins = c(10,28), 
              Colv = T,
              Rowv = Rowv,
              dendrogram = dendrogram,
              cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
              density.info=c("none"),
              keysize=1)    
    dev.off()
    
    pdf(file=paste0(LoaderPATH,"ClusterLogFrequencies",OutputSuffix,".pdf"), width=14, height =15)
    breaks = seq(0,max(log10(hpmat+.001)),by=max(log10(hpmat+.001))/50)
    my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
    heatmap.2(t(log10(hpmat+.001)), col=my_palette, 
              breaks = breaks,
              margins = c(10,28), 
              Colv = T,
              Rowv = Rowv,
              dendrogram = dendrogram,
              cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
              density.info=c("none"),
              keysize=1)    
    dev.off()
    
    write.csv(t(hpmat),paste0(LoaderPATH,"ClusterFreqValues",OutputSuffix,".csv"))
    write.csv(t(log10(hpmat+.001)),paste0(LoaderPATH,"ClusterLogFrequencies",OutputSuffix,".csv"))
  }
  
  if(DoPositiveFrequency) {
    hpmat <- NULL
    cutoffs <- GetCoords(LoaderPATH = pathForCoords, FNnames = FNnames, yparam = yparam, yparam2 = yparam2,
                         CoordsFN = "Coords.csv")
    for(file in 1:length(FcsFileNames)){
      frequencyData <- data[FFdata[, "InFile"]== file, , drop = F]
      numCellsPerSample <- nrow(frequencyData)
      percentPositivePerParameter <- NULL
      for(parameter in rownames(cutoffs)) {
        numCellsPositive <- nrow(frequencyData[frequencyData[, parameter] > cutoffs[parameter, ], , drop = F])
        percentPositivePerParameter <- c(percentPositivePerParameter, (numCellsPositive / numCellsPerSample * 100))
      }
      hpmat <- cbind(hpmat, percentPositivePerParameter)
    }
    colnames(hpmat) <- FcsFileNames
    rownames(hpmat) <- rownames(cutoffs)
    
    pdf(file=paste0(LoaderPATH,"PercentPositive",OutputSuffix,".pdf"), width=14, height =15)
    breaks = seq(0,max(hpmat),by=max(hpmat)/50)
    my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
    heatmap.2(t(hpmat), col=my_palette, 
              breaks = breaks,
              margins = c(20,28), 
              Colv = T,
              Rowv = Rowv,
              dendrogram = dendrogram,
              cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
              density.info=c("none"),
              keysize=1)    
    dev.off()
    
    write.csv(t(hpmat),paste0(LoaderPATH,"PercentPositive",OutputSuffix,".csv"))
  }
  
  if(DoIndividualCell) {
    hpmat = data1
    pdf(file=paste0(LoaderPATH,"IndividualCells",OutputSuffix,".pdf"), width=14, height =15)
    breaks = seq(0,max(hpmat),by=max(hpmat)/50)
    my_palette <- colorRampPalette(palette)(n=length(breaks)-1)
    heatmap.2(t(hpmat), col=my_palette, 
              breaks = breaks,
              margins = c(10,20), 
              dendrogram = "none",
              cexCol = 1., cexRow =1., scale="none", key=TRUE,trace="none", 
              density.info=c("none"),
              keysize=1)    
    dev.off()
  }
}

oneSENSEplotter <- function(LoaderPATH ="fcs_out",
                            Bins = 250, # number of bins for annotation
                            FNnames="trafficknames2.csv",
                            resolution = 800,
                            category1 = "Function",
                            category2 = "Phenotype",
                            HeatmapProperty = "Medians",
                            algorithmUsed = "UMAP",
                            ceil = 1000,
                            useCSV = F,
                            PlotAllCategories = F,
                            markerCap = 0,
                            colorBy = NULL,
                            DoChisqTest = F,
                            chisqBins = 9,
                            biplotColors = NULL,
                            biplotColorList = NULL) {
  
  print ('Making OneSENSE plots')
  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = F)
  FFdata<- prepData$FFdata
  FcsFileNames <- prepData$FcsFileNames
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  
  BaseFNs <- sapply(strsplit(prepData$FcsFileNames, split ="\\."), "[", 1)
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  
  keeptable <- read.csv(FNnames, check.names = F, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  
  data1u <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
  data1 <- nfTransform(keeprows, data1u, data1u, doNormalization = doNormalization)$dataA1
  data2u <- FFdata[,transList == "2"] 
  data2 <- apply(data2u, 2, logicleTransform(w=0.25, t=16409, m=4.5, a=0))
  data3 <- FFdata[,transList == "3",drop=F]
  
  data <- cbind (data1, data2, data3)
  if (PlotAllCategories) {
    categories <- colnames(keeptable)[4:length(colnames(keeptable))]
    if (length(categories) < 2) {
      stop("Not enough categories.")
    }
  } else {
    categories <- c(category1, category2)
  }
  
  if(algorithmUsed == "both"){
    algorithmUsed = c("UMAP", "tSNE")
  }
  
  coloring = !is.null(colorBy)
  
  if (!coloring) {
    colorBy = 'this is never going to be used anyway'
  }
  
  colorBys = colorBy
  
  for(i in 1:length(algorithmUsed)) {
    for(category in categories) {
      assign(paste0("keeprowbool.", category), sapply(keeptable[, category], isTransform))
      assign(paste0("categorygroup.", category), keeptable$Parameters[get(paste0("keeprowbool.", category))])
      assign(paste0("dataX.", category), data[,which (colnames(data) %in% get(paste0("categorygroup.", category)))])
      assign(paste0("dataX.", category), get(paste0("dataX.", category))[, !duplicated(colnames(get(paste0("dataX.", category))))])
      assign(paste0("oneD", algorithmUsed[i], "mat.", category), data[, paste0(category, "_", algorithmUsed[i]), drop = F])
      assign(paste0(algorithmUsed[i],"Bins.",category), cut(get(paste0("oneD", algorithmUsed[i], "mat.", category)), breaks = Bins, labels = 1:Bins))
      
      if (ncol(get(paste0("dataX.", category))) > markerCap && markerCap != 0) {
        proteinIndices <- grepl("TotalSeqC", colnames(get(paste0("dataX.", category))))
        markerCapUpdated <- markerCap - length(which(proteinIndices))
        assign(paste0("dataX.", category, "trimmed"), get(paste0("dataX.", category))[, !proteinIndices])
        sums <- sapply(1:ncol(get(paste0("dataX.", category, "trimmed"))), function(x) var(get(paste0("dataX.", category, "trimmed"))[, x]))
        test <- sort.int(sums, index.return = TRUE, decreasing = T)[[2]][1:markerCapUpdated]
        assign(paste0("dataX.", category), get(paste0("dataX.", category))[, c(which(proteinIndices), test)])
        assign(paste0("label.", category), paste0("Top", markerCapUpdated, category))
      } else {
        assign(paste0("label.", category), category)
      }
      
      if(HeatmapProperty == "Medians") {
        plotdata <- apply(get(paste0("dataX.", category)), 2, function(x) (tapply(x, get(paste0(algorithmUsed[i],"Bins.",category)), FUN = median)))
      } else if (HeatmapProperty == "Frequencies") {
        Coords <- GetCoords(LoaderPATH = LoaderPATH, FNnames = FNnames, yparam = yparam, yparam2 = yparam2,
                            CoordsFN = "Coords.csv")
        plotdata <- matrix(nrow = Bins, ncol = ncol(get(paste0("dataX.", category))))
        colnames(plotdata) <- colnames(get(paste0("dataX.", category)))
        for(pname in colnames(get(paste0("dataX.", category)))){
          overthresh <- function(group){
            percpos <- (sum(group > Coords[pname,1]) / length(group))*100
            return(percpos)
          }
          plotdata[,pname] <- tapply(get(paste0("dataX.", category))[,pname], get(paste0("umapBins.", category)), FUN = overthresh)
        }
      } else {
        stop("Incorrect heatmap type (choose either \"Medians\" or \"Frequencies\"")
      }
      plotdata[is.na(plotdata)] <- 0
      assign(paste0("plotdata.", category), plotdata)
    }
    
    raster.concatenated <- vector("list", 1000)
    index <- 1
    
    
    
    
    
    
    
    for(category in 1:(length(categories) - 1)) {
      category1 <- categories[category]
      for(j in (category + 1):length(categories)) {
        category2 <- categories[j]
        print (paste0('making plot ',category1, ' vs ', category2, ' with ', algorithmUsed[i]))
        print (colorBys)
        for(colorBy in colorBys) {
          print (colorBy)
          
          if(coloring){
            # c <- FFdata[,colorBy]
            # if (Ncolors == 0){
            #   Ncolors <- max(FFdata[,colorBy])
            #   c <- cut(c, breaks=0:Ncolors)
            # } else {
            #   c<- ((c+1))/(4.5)*Ncolors 
            # }
            # # palette <- colorRampPalette(MPColorList)(n=Ncolors)
            # palette <- MPColorList
            # preIndex <- grep("pre", FcsFileNames)
            # palette[preIndex] <- "lightgrey"
            #col = palette[as.numeric(c)]
            c <- FFdata[, colorBy]
            top = max(c)
            bottom = min(c)
            palette=c("black","blue","green","yellow","red")
            color.scale=unique(colorRampPalette(palette)(1000))
            if(top == bottom) {
              top = top + 0.1
            }
            breaks=seq(bottom,top,length.out=length(color.scale)+1)
            bnum <- length(breaks)
            breaks[1] <- -Inf
            breaks[bnum] <- Inf
            col=as.character(cut(c,breaks=breaks,labels=color.scale))
          } else {
            col = "black"
          }

          par(mar = c(0,0,0,0))
          parOrig <- par()
          
          label1 <- get(paste0("label.",category1))
          label2 <- get(paste0("label.",category2))
          
          text1 = paste(tmpDir(),"/textbox1.png",sep="")
          png(text1, res=resolution,height=resolution*2,width=480*resolution/72)
          plot.new()
          text(0.5, 0.5, label1)
          dev.off()
          
          text2 = paste(tmpDir(),"/textbox2.png",sep="")
          png(text2, res=resolution,height=resolution*2,width=480*resolution/72)
          par(bg = NA)
          plot.new()
          text(0.5, 0.5, label2)
          dev.off()
          
          # print(get(paste0("plotdata.", category1)))
          # print(get(paste0("plotdata.", category2)))
          
          cexRow1 <- 30/(ncol(get(paste0("plotdata.", category1))))
          if(cexRow1 > 0.8){
            cexRow1 <- 0.8
          }
          cexRow2 <- 30/(ncol(get(paste0("plotdata.", category2))))
          if(cexRow2 > 0.8) {
            cexRow2 <- 0.8
          }
          heatmap1 = paste(tmpDir(),"/heatmap1.png",sep="")
          png(heatmap1,res=resolution,height=480*resolution/36,width=480*resolution/72)
          breaks = seq(0,max(get(paste0("plotdata.", category1))),by=0.005)  # you can put whatever colors you like here:
          #  my_palette <- colorRampPalette(c("blue","green","yellow"))(n=430)
          #my_palette <- colorRampPalette(c("black","orange","yellow","white"))(n=length(breaks)-1)
          #my_palette <- colorRampPalette(c("black","blue","green","yellow"))(n=length(breaks)-1)
          #my_palette <- colorRampPalette(c("blue","green","yellow"))(n=length(breaks)-1)
          my_palette <- colorRampPalette(c("black","blue","green","yellow","darkred"))(n=length(breaks)-1)
          #my_palette <- colorRampPalette(c("blue","lightblue","green","yellow","darkorange","darkred"))(n=length(breaks)-1)
          
          suppressWarnings(heatmap.2(t(get(paste0("plotdata.", category1))), col=my_palette, 
                                     breaks = breaks,
                                     margins = c(10,10), 
                                     Colv = F,
                                     dendrogram = "none",
                                     cexCol = 1., cexRow = cexRow1, scale="none", key = FALSE, trace="none", labCol = "",
                                     density.info=c("none"),lmat = rbind(1,4,3,2), lhei= c(3,2,2,2), lwid = 6,
                                     keysize=0.5))    
          dev.off()
          
          breaks = seq(0,max(get(paste0("plotdata.", category2))),by=0.005)
          my_palette <- colorRampPalette(c("black","blue","green","yellow","darkred"))(n=length(breaks)-1)
          
          heatmap2 = paste(tmpDir(),"/heatmap2.png",sep="")
          png(heatmap2,res=resolution,height=480*resolution/36,width=480*resolution/72)
          breaks = breaks
          suppressWarnings(heatmap.2(t(get(paste0("plotdata.", category2))), col=my_palette, 
                                     breaks = breaks,
                                     margins = c(10,10), 
                                     Colv = F,
                                     dendrogram = "none",
                                     cexCol = 1, cexRow = cexRow2, scale="none", key=F,trace="none", labCol = "",
                                     density.info=c("none"), lmat = rbind(1,4,3,2), lhei= c(3,2,2,2), lwid = 6,
                                     keysize=0.5))    
          dev.off()
          
          
          oneDmat1 <- get(paste0("oneD", algorithmUsed[i], "mat.", category1))
          oneDmat2 <- get(paste0("oneD", algorithmUsed[i], "mat.", category2))
          
          if(DoChisqTest) {
            chi = paste0(tmpDir(),'/chi.png')
            png(chi,res = resolution, height = 480*resolution/36, width = 480*resolution/36)
            plot.new()
            binsize = (max(oneDmat1)-min(oneDmat1)+max(oneDmat2)-min(oneDmat2))/(2*chisqBins)
            print(binsize)
            plot.new()
            rounded_x = round(oneDmat1/binsize)*binsize
            rounded_y = round(oneDmat2/binsize)*binsize
            chiresult = chisq.test(rounded_x, rounded_y)
            print(chiresult)
            pal <- colorRampPalette(c(rgb(0.5,1,0), rgb(0,1,1), rgb(1,1,1)))
            par(bg=rgb(0,0,0,0))
            heatmap(chiresult$stdres, NA, NA, margins = c(0,0), col = colorRampPalette(c(rgb(1,1,1),rgb(1,0.8,0.8)))(20))
            
            dev.off()
            
            background = readPNG(chi,native = T)
            
          }
          
          mainplot <- paste(tmpDir(),"/mainplot.png", sep="")
          png(mainplot, res = resolution, height = 480*resolution/72, width=480*resolution/72)
          plot.new()
          par(mar = c(0,0,0,0), bg = 'transparent')
          
          if (DoChisqTest) {
            rasterImage(background, -0.15, -0.20, 1.06,1.18 )#, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
            par(new = T)
          }
          
          
          plot(oneDmat2, oneDmat1, pch = ".", cex = 5, xaxt='n', yaxt='n', ann=FALSE, pin = c(4,4),
               xlim = c(min(oneDmat2) - max(oneDmat2)/40, max(oneDmat2) + max(oneDmat2/40)),
               ylim = c(min(oneDmat1) - max(oneDmat1)/40, max(oneDmat1) + max(oneDmat1/40)),
               col = col)
          box(lwd=5)
          
          
          
          dev.off()
          
         
          
          #tests = paste0(tmpDir(), "/tests.png")
          #plot.new()
          #text(0.5,0.5,chisq.test(oneDmat2,oneDmat1))
          #dev.off()
          

          if(coloring){
            # legend <- paste(tmpDir(), "/legend.png", sep="")
            # png(legend, res = resolution, height = 480*resolution/120, width=480*resolution/72)
            # plot.new()
            # legend(x = "top",
            #        inset = c(0, -0.1),
            #        legend = FcsFileNames,
            #        col = unique(col),
            #        pch = 16,
            #        bty = "n",
            #        pt.cex = 2,
            #        cex = 1,
            #        text.col = "black",
            #        horiz = F,
            #        ncol = max(FFdata[, "InFile"])%/%20 + 1,
            #        xpd = T)
            # dev.off()
            legend=paste(tmpDir(),"/legend.png",sep="")
            png(legend,res=resolution,height=480/2*resolution/72,width=480*resolution/72)
            plot.new()
            par("mar"=c(2,1,2,1))
            xlims=par("usr")[1:2]
            ylims=par("usr")[3:4]

            n=length(col)
            labels=signif((seq(bottom,top,length.out=length(color.scale)+1)),2)

            x_coords=seq(xlims[1],xlims[2],length.out=n+1)
            rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=NA)
            labels.x_coords=seq(x_coords[1],x_coords[length(x_coords)],length.out=5)
            labels=labels[round(seq(1,length(labels),length.out=5))]
            text(xpd=T,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
            text(xpd=T,y=ylims[2],pos=3,labels=paste(colorBy,"intensity"),x=mean(xlims))
            dev.off()
          } else {
            legend <- paste(tmpDir(), "/legend.png", sep="")
            png(legend, res = resolution, height = 480*resolution/120, width=480*resolution/72)
            plot.new()
            legend(x = "top",
                   inset = c(0, -0.1),
                   legend = FcsFileNames,
                   col = unique(col),
                   pch = 16,
                   bty = "n",
                   pt.cex = 2,
                   cex = 1,
                   text.col = "black",
                   horiz = F,
                   ncol = max(FFdata[, "InFile"])%/%5 + 1,
                   xpd = T)
            dev.off()
            
          }
          

          #hist(tSNEmat1, 100, main= paste("Histogram of",colnames(keeptable[factor] )))
          
          histogram1 = paste0(tmpDir(), '/histogram1.png')
          png(histogram1, res = resolution, height = 480*resolution/72, width = 480*resolution/72)
          plot.new()
          hist(oneDmat1, 100, density = 50, axes = F, labels = F, main = NULL, xlab = NULL, ylab = NULL)
          dev.off()
          
          histogram2 = paste0(tmpDir(), '/histogram2.png')
          png(histogram2, res = resolution, height = 480*resolution/72, width = 480*resolution/72)
          plot.new()
          hist(oneDmat2, 100, density = 50, axes = F, labels = F, main = NULL, xlab = NULL, ylab = NULL)
          dev.off()
          
          
          
          raster.textbox1 = readPNG(text1, native = T)
          raster.textbox2 = readPNG(text2, native = T)
          raster.heatplot1 = readPNG(heatmap1, native = T)
          raster.heatplot2 = readPNG(heatmap2, native = T)
          raster.mainplot = readPNG(mainplot, native = T)
          raster.legend = readPNG(legend, native = T)
          raster.hist1 = readPNG(histogram1, native = T)
          raster.hist2 = readPNG(histogram2, native = T)
          #raster.chi = readPNG(chi, native = T)
          
          
          

          concatenatedplot = paste(tmpDir(), "/concatenatedplot.png", sep="")
          png(concatenatedplot, res = resolution, height = 480*resolution/72, width=480*resolution/72)
          plot.new()
          pushViewport(viewport(angle = 90))
          grid.raster(raster.textbox1, y = 0.905, x = 0.64, height = 0.2)
          grid.raster(raster.heatplot1, y = 0.425, x = 0.7, height = 1.0)
          pushViewport(viewport(angle = 180))
          grid.raster(raster.hist1, y = 0.86, x = 0.365, height = -0.1, width = 0.53)
          pushViewport(viewport(angle = 90))
          grid.raster(raster.hist2, y = 0.917, x = 0.575, height = 0.1, width = 0.53)
          grid.raster(raster.legend, y = 0.2, x = 0.2, height = 0.17)
          grid.raster(raster.mainplot, y = 0.65, x = 0.59, height = 0.455)
          grid.raster(raster.heatplot2, y = -0.1, x = 0.64, height = 1.0)
          grid.raster(raster.textbox2, y = 0.135, x = 0.575, height = 0.2)
          #grid.raster(raster.chi, y = 0.2, x = 0.2, height = 0.4 )
          dev.off()
          
          raster.concatenated[[index]] <- readPNG(concatenatedplot, native = T)
          index <- index + 1
        }
      }
    }
    
    raster.concatenated <- raster.concatenated[-which(sapply(raster.concatenated, is.null))]
    
    suppressWarnings(dir.create("OneSENSE Plots"))
    
    pdf(paste0("OneSENSE Plots/OneSENSE_concatenated_",algorithmUsed[i],".pdf"))
    for (i in 1:length(raster.concatenated)) {
      grid.newpage()
      grid.raster(raster.concatenated[[i]])
    }
    dev.off()
  }
}

#### Data Processing Functions -------------------------------------------------------------------------

isTransform = function(tag) {
  return(!tag == '')
  #return( tag %in% c("y", "c", "l", "a","n","as5","as150","x","ln","0","+"))
}

#Transforms the data based on the data type and based on user specifications in the "names2.csv" file
nfTransform <- function(transTypeTable, dataA, dataB, doNormalization=F, factor=2){
  dataA1<-dataA
  dataB1<-dataB
  CyTOFlgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  ilgcl <- inverseLogicleTransform(trans = CyTOFlgcl)
  Fluorlgcl <- logicleTransform(w=0.1, t=500000, m=4.5, a=0)
  as5trans <- arcsinhTransform(a=0, b=(1/5))
  as150trans <- arcsinhTransform(a=0, b=(1/150))
  
  if(doNormalization) {
    rowAvg <- sapply(1:nrow(dataB1), function(x){
      return(sum(dataB1[x, -which("Cell_Index" %in% colnames(dataB1))]))
    })
    dataA <- lapply(1:nrow(dataA1), function(x) {
      return(dataA1[x,]/rowAvg[x]*10000)
    })
    dataB <- lapply(1:nrow(dataB1), function(x) {
      return(dataB1[x,]/rowAvg[x]*10000)
    })
    dataA <- do.call(rbind, dataA)
    dataB <- do.call(rbind, dataB)
  }
  # make groups
  
  groups = list()
  group_index = 0
  paramNames = as.character((transTypeTable[,1]))
  for (i in 1:length(paramNames)){
    ttParamNum <- which(as.character(transTypeTable[,1])==paramNames[i])
    ttParamType <- as.character(transTypeTable[ttParamNum,factor])
    if (ttParamType == '"') {
      groups[[group_index]] = c(groups[[group_index]],i)
    } else {
      group_index = group_index + 1
      groups[group_index] = i
    }
  }
  
  #print (groups)
  
  
  #perform transform
  last = 'a'
  for(i in groups)
  {
    paramName = paramNames[i]
    #ttParamNum <- which(paramName %in% as.character(transTypeTable[,1]))
    ttParamNum = i
    ttParamType <- as.character(transTypeTable[ttParamNum,factor])
    temp <- NULL
    dataNum <- which(colnames(dataA) %in% paramName)
    
    #print (paste0('name:',paramName,' ttpnuum:', ttParamNum, ' type:', ttParamType, ' datanum:',dataNum))
    ttParamType = ttParamType[1]
    
    if(grepl('y',ttParamType) || grepl('c',ttParamType)){   # Default CyTOF Parameters 
      y <- apply(dataA[,dataNum,drop=F],2, CyTOFlgcl)
      c = y
      #print(paste(paramName," CyTOFlgcl"))
    } else if( grepl('f',ttParamType) ) {                # Fluor Logicle Parameters (not recommended)
      f <- apply(dataA[,dataNum,drop=F],2, Fluorlgcl)
    } else if( grepl('ln',ttParamType) ) {
      ln <- log(dataA[,dataNum,drop=F] + 1)
    } else if (grepl('x',ttParamType)) { # no transform
      x = dataA[,dataNum]
    } else if (grepl('plus',ttParamType)) { # makes all val positive
      plus = dataA[,dataNum] - min(dataA[,dataNum])
    } else if( grepl('l',ttParamType) ) {                # Linear Normalized to max = 4.5
      l <- (dataA[,dataNum]/max(dataA[,dataNum]))*4.5
    } else if( grepl('n',ttParamType) ){                 # Linear Normalized to min = 0, max = 4.5
      n <- ((dataA[,dataNum]-min(dataA[,dataNum]))/(max(dataA[,dataNum])-min(dataA[,dataNum])))*4.5
    } else if( grepl('as5',ttParamType) ) {              # ArcSinh x/5 (for cytof)
      as5 <- apply(dataA[,dataNum,drop=F],2, as5trans)
    } else if( grepl('as150',ttParamType) ) {            # ArcSinh x/150 (for fluor)
      as150 <- apply(dataA[,dataNum,drop=F],2, as150trans)
    } else if (grepl('z', ttParamType)) {
      z <- (dataA[, dataNum] - mean(dataA[, dataNum])) / sd(dataA[, dataNum])
    } else if( grepl('a',ttParamType) ) {                 # Automatic determination of transform
      q <- 0.05
      m <- 4.5
      d <- dataA[,paramName]
      w <- 0
      t <- max(d)
      nd <- d[d < 0]
      nThres <- quantile(nd, 0.25) - 1.5 * IQR(nd)
      nd <- nd[nd >= nThres]
      #transId <- paste(p, "autolgclTransform", sep = "_")
      if (length(nd)) {
        r <- .Machine$double.eps + quantile(nd, q)
        if (10^m * abs(r) <= t) {
          w <- 0
        } else {
          w <- (m - log10(t/abs(r)))/2
          if (is.nan(w) || w > 2) {
            warning(paste0("autoLgcl failed for channel: ",
                           paramName, "; using default fluor logicle transformation!"))
            w <- 1
            t <- 500000
            m <- 4.5
          }
        }
      }
      templgcl <- logicleTransform(w=w, t=t, m=4.5, a=0)
      a <- apply(dataA[,dataNum,drop=F],2, templgcl)
      print(paste0(paramName, " w= ",w," t= ",t))
      #hist(temp, main=paramName,breaks=100)
    } else if(substring(ttParamType,0,3) == 'eq=') {
      equation = substring(ttParamType,4)
      x = dataA[,dataNum]
      temp = eval(parse(text = equation))
    }
    
    equation = ttParamType
    temp = eval(parse(text = equation))
    
    last = ttParamType
    
    dataA1[,dataNum] <- temp
    dataB1[,dataNum] <- temp
  }
  
  return(list(dataA1=dataA1, dataB1=dataB1))
}

# Prepares and concatenates a giant matrix of data based on all the input files for use in performing algorithms on.
prepFcsFolderData <- function(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = useCSV){
  if(!useCSV){
    FcsFileNames <- list.files(path = LoaderPATH, pattern = ".fcs")
    fs = list() # list of all fcs files in input directory
    for(FileNum in 1:length(FcsFileNames)){
      fs[[FileNum]] <- read.FCS(paste0(LoaderPATH,"/",FcsFileNames[FileNum]),transformation =FALSE,ignore.text.offset=T,truncate_max_range=FALSE,emptyValue = F)
    }
    numFiles <- length(FcsFileNames)
    FFdata <- NULL
    OrigNames <- fs[[1]]@parameters$name # getting names before they are removed. from first file bc they all sould be the same
    for (FFs in 1:numFiles){ 
      FFt <- exprs(fs[[FFs]]) 
      ## Downsample (remove samples above ceiling)
      if (nrow(FFt) <= ceil) { # ceiling is higher than number of samples in the file
        FFa <- FFt
      } else {                 # downsamples data from each input file
        FFa <- FFt[sample(nrow(FFt),ceil,replace=F),]
      }
      
      #Fixup column names (remove unused names)
      colnames(FFa) <- fs[[FFs]]@parameters$desc
      empties <- which(is.na(colnames(FFa)) | colnames(FFa)== " ")
      colnames(FFa)[empties] <- fs[[FFs]]@parameters$name[empties]
      fs[[FFs]]@parameters$desc <- colnames(FFa)
      
      #Add file label
      FFa <- cbind(FFa,rep(FFs,dim(FFa)[1]))
      colnames(FFa)[dim(FFa)[2]] <- "InFile"
      
      if(!is.null(FFdata)) {
        diffNames <- c(setdiff(colnames(FFa), colnames(FFdata)), setdiff(colnames(FFdata), colnames(FFa)))
        matchingNames <- intersect(colnames(FFa), colnames(FFdata))
        if(length(diffNames) != 0) {
          warning(paste("Parameter(s)", diffNames,"is/are missing from one or more of the .FCS files. Removing from FFdata now."), immediate. = T)
        }
        if(length(matchingNames == 0)) {
          stop("There are no matching parameters. Please double-check that you have put in the riht .FCS files and try again.")
        }
        FFdata <- FFdata[, matchingNames]
        FFa <- FFa[, matchingNames]
        OrigNames <- matchingNames[-length(matchingNames)]
      }
      
      #Concatenate
      FFdata <- rbind(FFdata,FFa)
    }
  } else {  #load csv files instead of fcs
    csvfilenames <- list.files(path = LoaderPATH, pattern=".csv")
    filenames <- csvfilenames
    csvdata <- lapply(paste0(LoaderPATH,"//",csvfilenames),function(x) read.csv(x, check.names = F,stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM"))
    FFdata<-NULL
    numFiles <- length(csvfilenames)
    
    for (FFs in 1:numFiles){
      FFt <- csvdata[[FFs]]
      FFt <- apply(FFt,2,function(x){as.numeric(gsub(",", "", x))})
      
      ## Downsample
      if (nrow(FFt)<=ceil) {
        FFa <- FFt
      } else {
        FFa <- FFt[sample(nrow(FFt),ceil,replace=F),]
      }
      
      FFa <- cbind(FFa,rep(FFs,dim(FFa)[1]))
      colnames(FFa)[dim(FFa)[2]] <- "InFile"
      
      #Concatenate
      FFdata <- rbind(FFdata,FFa)
    }
    OrigNames <- colnames(FFt[[1]])
    fs<-FFt
    FcsFileNames <- csvfilenames
  }
  return(list(FFdata=FFdata, OrigNames=OrigNames, forNewFF=fs[[1]], numFiles=numFiles, FcsFileNames = FcsFileNames))
  #           the data table    the col names before     first file?      number of files     file names
  #                             any were removed
}

reverseTransform <- function(score2 = NULL,
                             nontransformedData = NULL,
                             FFdata = NULL,
                             TransformOutput = F) {
  if(!is.null(score2))  lscore2 <- dim(score2)[2] else lscore2<-0
  if(!is.null(nontransformedData))  lnontransformedData <- dim(nontransformedData)[2] else lnontransformedData<-0
  
  NewTransformBase <- c (rep(1, dim(FFdata)[2]-1))
  if(TransformOutput){
    NewTransformNew <- c(rep(2,lscore2),
                         rep(3,lnontransformedData))
  }else  NewTransformNew <- c(rep(3,lscore2),
                              rep(3,lnontransformedData))
  NewTransform <- c(NewTransformBase,NewTransformNew)
  #
  if(TransformOutput){
    Nscore <- apply(score2,2,function(x) ((x-min(x))/(max(x)-min(x)))*3.7 )
    ilgcl <- inverseLogicleTransform(trans = logicleTransform(w=0.25, t=16409, m=4.5, a=0))
    NIscore <- apply(Nscore, 2, ilgcl)
    NIscore = cbind(NIscore,nontransformedData)
  } else if(!is.null(score2)){
    
    NIscore <- apply(score2,2,function(x) ((x-min(x))/(max(x)-min(x)))*10000 )
    #ilgcl <- inverseLogicleTransform(trans = lgcl)
    #NIscore <- apply(Nscore, 2, ilgcl)
    
    NIscore = cbind(NIscore,nontransformedData)
    #colnames(NIscore) <- c(colnames(score2),"phenograph")
  } else if(!is.null(nontransformedData)) {
    NIscore <- nontransformedData
  } else {
    print("NEED TO SELECT AT LEAST ONE ALGORITHM")
    NIscore <- nontransformedData
  }
  #Add * to duplicated parameter names
  preParams <- colnames(FFdata)
  newParams <- colnames(NIscore)
  
  while(length(intersect(preParams, newParams))>0)
  {
    dupNames <- intersect(preParams, newParams)
    for( dN in dupNames){
      newDN <- paste0(dN,"*")
      newParams[newParams == dN] <- newDN
    }
  }
  colnames(NIscore) <- newParams
  
  return(list(NIscore, NewTransformBase, NewTransformNew, NewTransform))
}

# Export data to FCS
exportToFCS <- function(useCSV = F,
                        numFiles = 1,
                        forNewFF = NULL,
                        OrigNames = NULL,
                        FFdata = NULL,
                        NIscore = NULL,
                        NewTransform = NULL,
                        NewTransformNew = NULL,
                        LoaderPATH = "fcs",
                        OutputSuffix = "out",
                        FcsFileNames = NULL) {
  if(!useCSV){
    #Splits into respective flowFrames.
    for (FFs in 1:numFiles){
      newFF <- forNewFF
      newBaseData <- FFdata[FFdata[,ncol(FFdata)]==FFs,-ncol(FFdata), drop = F]
      if (nrow(newBaseData) != 0) {
        newBaseDataDesc <- colnames(newBaseData)
        colnames(newBaseData) <- OrigNames
        exprs(newFF) <- newBaseData
        subsetNIscore <- NIscore[FFdata[,dim(FFdata)[2]]==FFs,,drop=FALSE]
        newFF <- cbind2(newFF, subsetNIscore)
        
        if(!is.null(newFF@description$CREATOR)) {
          if(newFF@description$CREATOR == "ENewellScript") { 
            preNewTransform <- newFF@description$TF
          }
        }
        
        newFF@parameters$desc <- c(newBaseDataDesc,colnames( subsetNIscore))
        suppressWarnings(dir.create(paste0(LoaderPATH,"_",OutputSuffix)))
        
        BaseFN <- sapply(strsplit(FcsFileNames[FFs], split ="\\."), "[", 1)
        FNresult <- paste0(LoaderPATH,"_",OutputSuffix,"/",BaseFN,"_",OutputSuffix,".fcs")
        newFF@description$'$FIL' <- paste0(BaseFN,"_",OutputSuffix,".fcs")
        
        if(!is.null(newFF@description$CREATOR)) {
          if( newFF@description$CREATOR == "ENewellScript") {
            newFF@description$TF <- paste0(preNewTransform,paste(unlist(NewTransformNew),collapse=''))
          } else {
            newFF@description$TF <- paste(unlist(NewTransform),collapse='')
          }
        } else {
          newFF@description$TF <- paste(unlist(NewTransform),collapse='')
        }
        newFF@description$CREATOR <- "ENewellScript"
        newFF@description$FILENAME <- paste0(BaseFN,"_",OutputSuffix,".fcs")
        identifier(newFF) <- paste0(BaseFN,"_",OutputSuffix)
        
        newFF@parameters@data$range <- rep(10000,length(newFF@parameters@data$range))
        
        for(r in 1:length(newFF@parameters@data$name))
        {
          row.names(newFF@parameters@data)[r] <- paste0("$P",r) 
        }
        
        print(FNresult)
        write.FCS(newFF, FNresult)
      }
    }
  } else {
    suppressWarnings(dir.create(paste0(LoaderPATH,"_",OutputSuffix)))
    for (FFs in 1:numFiles){
      newBaseData <- FFdata[FFdata[,ncol(FFdata)]==FFs,-ncol(FFdata)]
      subsetNIscore <- NIscore[FFdata[,dim(FFdata)[2]]==FFs,,drop=FALSE]
      newBaseData <- cbind(newBaseData, subsetNIscore)
      
      if (nrow(newBaseData) != 0) {
        meta <- data.frame(name=colnames(newBaseData),
                           desc=colnames(newBaseData))
        meta$range <- apply(apply(newBaseData,2,range), 2, diff)
        meta$minRange <- apply(newBaseData,2,min)
        meta$maxRange <- apply(newBaseData,2,max)
        
        newFF <- new("flowFrame",
                     exprs=newBaseData,
                     parameters=AnnotatedDataFrame(meta))
        
        if(!is.null(newFF@description$CREATOR)) {
          if(newFF@description$CREATOR == "ENewellScript") { 
            preNewTransform <- newFF@description$TF
          }
        }
        
        BaseFN <- sapply(strsplit(FcsFileNames[FFs], split ="\\."), "[", 1)
        FNresult <- paste0(LoaderPATH,"_",OutputSuffix,"/",BaseFN,"_",OutputSuffix,".fcs")
        newFF@description$'$FIL' <- paste0(BaseFN,"_",OutputSuffix,".fcs")
        
        if(!is.null(newFF@description$CREATOR)) {
          if( newFF@description$CREATOR == "ENewellScript") {
            newFF@description$TF <- paste0(preNewTransform,paste(unlist(NewTransformNew),collapse=''))
          } else {
            newFF@description$TF <- paste(unlist(NewTransform),collapse='')
          }
        } else {
          newFF@description$TF <- paste(unlist(NewTransform),collapse='')
        }
        newFF@description$CREATOR <- "ENewellScript"
        newFF@description$FILENAME <- paste0(BaseFN,"_",OutputSuffix,".fcs")
        identifier(newFF) <- paste0(BaseFN,"_",OutputSuffix)
        newFF@parameters@data$range <- rep(10000,length(newFF@parameters@data$range))
        
        for(r in 1:length(newFF@parameters@data$name)){
          row.names(newFF@parameters@data)[r] <- paste0("$P",r) 
        }
        
        print(FNresult)
        write.FCS(newFF, FNresult)
        write.csv(exprs(newFF),paste0(LoaderPATH,"_",OutputSuffix,"/",BaseFN,"_",OutputSuffix,".csv"),row.names=F)
      }
    }
  }
}

getGrouping <- function(FFdata = NULL,
                        FcsFileNames = NULL,
                        GroupingName = "TissueType",
                        GroupingOutput = "GroupingOutput"){ 
  if (file.exists(paste0(GroupingOutput, ".csv"))) {
    groupingMat <- read.csv(paste0(GroupingOutput, ".csv"), fileEncoding="UTF-8-BOM")
  } else {
    groupingMat <- NULL
  }
  groupingType <- list()
  groupingDict <- list()
  groupingCol <- NULL
  correspondingNumber <- NULL
  
  if(!is.null(groupingMat)) {
    if (GroupingName %in% colnames(groupingMat)){
      GroupNums <- groupingMat[, (grep(GroupingName, colnames(groupingMat)) + 1)]
      for (i in 1:length(GroupNums)){
        groupingCol <- c(groupingCol, rep(GroupNums[i], length(FFdata[FFdata[, "InFile"] == i, "InFile"])))
      }
    }
  } else {
    for(i in 1:max(FFdata[, 'InFile'])){
      groupingType[i] <- list(readline(prompt = paste0("What ",GroupingName," is the file \"", FcsFileNames[i], "\"? ")))
      if (is.na(match(groupingType[i][1], groupingDict))){
        groupingDict[length(groupingDict) + 1] <- groupingType[i][1]
      }
      correspondingNumber <- c(correspondingNumber, grep(groupingType[i][1], groupingDict))
      groupingCol <- c(groupingCol, rep(grep(groupingType[i][1], groupingDict), length(FFdata[FFdata[, "InFile"] == i, "InFile"])))
    }
    groupingMat <- cbind2(cbind(groupingType, groupingNumber = correspondingNumber), groupingMat)
    colnames(groupingMat)[1] <- GroupingName
    write.csv(groupingMat, file = paste0(GroupingOutput, ".csv"))
  }
  return(groupingCol)
}



# returns a llist of pngs each a plot for the pdf
generateRasters <- function(plotdata, x, y, color.scale.type, bottoms, tops, resolution, PC1, PC2, doAxes, nameByFile = F, fileName = NULL){
  rasters=sapply(colnames(plotdata),function(pname)
  {
    Exp <- plotdata[,pname]
    PltDat <- data.frame(x,y,Exp)
    
    palette=c("black","blue","green","yellow","red")
    color.scale=unique(colorRampPalette(palette)(1000))
    if(color.scale.type=="relative"){
      if(tops[pname] == bottoms[pname]) {
        tops[pname] = tops[pname] + 0.1
      }
      breaks=seq(bottoms[pname],tops[pname],length.out=length(color.scale)+1)
      bnum <- length(breaks)
      breaks[1] <- -Inf
      breaks[bnum] <- Inf
      points.colors=as.character(cut(Exp,breaks=breaks,labels=color.scale))
    }
    if(color.scale.type=="absolute"){
      breaks=seq(plotdata.range[1],plotdata.range[2],length.out=length(color.scale)+1)
      bnum <- length(breaks)
      breaks[1] <- -Inf
      breaks[bnum] <- Inf
      points.colors=as.character(cut(Exp,breaks=breaks,labels=color.scale))
    }
    
    mainplot=paste(tmpDir(),"/mainplot.png",sep="")
    png(mainplot,res=resolution,height=480*resolution/72,width=480*resolution/72)
    par("bty"="l")
    if(doAxes) {
      maincex <- 1
    } else {
      maincex <- 3
      PC1 = ""
      PC2 = ""
    }
    
    if(nameByFile) {
      main <- fileName
    } else {
      main <- pname
    }
    
    plot(x,y,col=points.colors,xlab=PC1,ylab=PC2,main=main,pch=pch,cex=cex,axes = doAxes, cex.main = maincex)#, xlim = c(0,10000), ylim = c(0,10000))
    dev.off()
    
    if(color.scale.type=="relative"){
      colorscale=paste(tmpDir(),"/colorscale.png",sep="")
      png(colorscale,res=resolution,height=480/2*resolution/72,width=480*resolution/72)
      plot.new()
      par("mar"=c(2,1,2,1))
      xlims=par("usr")[1:2]
      ylims=par("usr")[3:4]
      
      n=length(color.scale)
      labels=signif((seq(bottoms[pname],tops[pname],length.out=length(color.scale)+1)),2)
      
      x_coords=seq(xlims[1],xlims[2],length.out=n+1)
      rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
      labels.x_coords=seq(x_coords[1],x_coords[length(x_coords)],length.out=5)
      labels=labels[round(seq(1,length(labels),length.out=5))]
      text(xpd=T,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
      text(xpd=T,y=ylims[2],pos=3,labels=paste(pname,"intensity"),x=mean(xlims))
      dev.off()
    }
    if(color.scale.type=="absolute"){
      colorscale=paste(tmpDir(),"/colorscale.png",sep="")
      png(colorscale,res=resolution,height=480/2*resolution/72,width=480*resolution/72)
      plot.new()
      par("mar"=c(2,1,2,1))
      xlims=par("usr")[1:2]
      ylims=par("usr")[3:4]
      
      n=length(color.scale)
      labels=signif((seq(plotdata.range[1],plotdata.range[2],length.out=length(color.scale)+1)),2)
      
      x_coords=seq(xlims[1],xlims[2],length.out=n+1)
      rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
      labels.x_coords=seq(x_coords[1],x_coords[length(x_coords)],length.out=5)
      labels=labels[round(seq(1,length(labels),length.out=5))]
      text(xpd=T,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
      text(xpd=T,y=ylims[2],pos=3,labels=paste("Intensity"),x=mean(xlims))
      dev.off()
    }
    return(list(raster.main=readPNG(mainplot,native=T),raster.scale=readPNG(colorscale,native=T)))
  },simplify=F)
  return(rasters)
}

generateDiscreteRasters <- function(discreteData = NULL,
                                    FFdata = NULL,
                                    clustParam = "Phenograph",
                                    Ncolors = 0,
                                    PC1 = "UMAP1",
                                    PC2 = "UMAP2",
                                    resolution = 600,
                                    MPColorList = NULL,
                                    labelClusters = F,
                                    FcsFileNames = NULL,
                                    GroupingOutput = "GroupingOutput",
                                    pch = pch,
                                    cex = cex) {
  GroupingOutput = paste0(GroupingOutput, ".csv")
  if(file.exists(GroupingOutput)){
    groupings <- read.csv(GroupingOutput, fileEncoding="UTF-8-BOM")
  }
  rasters = sapply(colnames(discreteData),function(colparN) {
    c <- FFdata[,colparN]
    if (Ncolors == 0){
      Ncolors <- max(FFdata[,colparN])
      c <- cut(c, breaks=0:Ncolors)
    } else {
      c<- ((c+1))/(4.5)*Ncolors 
    }
    Ddata <- FFdata[,c(PC1, PC2, colparN)]
    palette <- colorRampPalette(MPColorList)(n=Ncolors)
    
    #scatterplot3d(Ddata[,-4], pch = ".", highlight.3d=TRUE, cex.symbols=1)
    x<-Ddata[,1]
    y<-Ddata[,2]
    
    # x1 <- (x - min(x))/(max(x) - min(x))
    # y1 <- (y - min(y))/(max(y) - min(y))
    
    mainplot=paste(tmpDir(),"/mainplot.png",sep="")
    png(mainplot,res=resolution,height=480*resolution/72,width=480*resolution/72)
    par("bty"="l")
    
    randOrder <- sample(length(x))
    col = palette[as.numeric(c)[randOrder]]
    
    plot(x=x[randOrder], y=y[randOrder], cex = cex, pch = pch, main = colparN, cex.main = 1,
         xlab = PC1, ylab = PC2, 
         col = col,
         axes = T)
    
    if (colparN == "InFile") {
      legendNames <- FcsFileNames
    } else if (colparN != clustParam) {
      legendNames <- unique(groupings[, colparN])
    }
    
    if(labelClusters && colparN == clustParam){
      for(cluster in 1:max(as.numeric(c))) {
        xclustlab <- median(x[as.numeric(c)==cluster])
        yclustlab <- median(y[as.numeric(c)==cluster])
        text(xclustlab,yclustlab, labels=paste0(cluster), cex=0.8, col="black")
      } 
    }
    dev.off()
    
    legendGraphic=paste(tmpDir(),"/legend.png",sep="")
    png(legendGraphic,res=resolution,height=480*resolution/72,width=480*5*resolution/72)
    ncol <- (max(discreteData[,colparN])%/%10 + 1)
    
    plot.new()
    if (colparN != clustParam) {
      legend("bottom",
             legend = legendNames,
             col = col,
             pch = 16,
             bty = "n",
             pt.cex = 2,
             cex = 2.5,
             text.col = "black",
             horiz = F,
             ncol = ncol,
             xpd = T)
    }
    dev.off()
    
    return(list(raster.main=readPNG(mainplot,native=T), raster.legend=readPNG(legendGraphic,native=T)))
  },simplify=F)
  return(rasters)
}

#### Algorithms -------------------------------------------------------------------------------------

#Runs a traditional 2D tSNE
tSNE_2D <- function(dataTransformed = NULL,
                    tSNEperplexity = 30,
                    FFdata = NULL,
                    sparseScore = NULL) {
  if(!is.null(sparseScore)) {
    dataTransformed <- sparseScore
  } else if(usingSeuratData) {dataTransformed <- pca(dataTransformed)}
  print(dim(dataTransformed))
  print("running tSNE")
  tSNEmat <- Rtsne(dataTransformed, dims=2, perplexity=tSNEperplexity, check_duplicates=F,verbose = T, pca = F)$Y
  
  columnMatch <- grep("tSNE", colnames(FFdata))
  if (length(columnMatch != 0)) {
    colnames(tSNEmat) <- c(colnames(FFdata)[columnMatch[length(columnMatch) - 1]], colnames(FFdata)[columnMatch[length(columnMatch)]])
  } else {
    colnames(tSNEmat) <- c("tSNE1","tSNE2")
  }
  plot(tSNEmat[, 1], tSNEmat[, 2], pch=".", xlab="tSNE1", ylab="tSNE2", cex=0.1)
  return(tSNEmat)
}

tSNE_oneSENSE <- function(keeptable = NULL,
                          FFdata = NULL,
                          tSNEperplexity= 30,
                          markerCap = 0) {
  print("running one-sense")
  OneStSNEmat <- NULL
  for (factor in 4:(dim(keeptable)[2])){
    print(factor)
    OneDtSNEname <- paste0(colnames(keeptable)[factor],"_tSNE")
    keeprowbool <- sapply(keeptable[,factor], isTransform)
    keeprows <- subset(keeptable, keeprowbool)
    dataX <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
    dataX <- nfTransform(keeprows,dataX,FFdata, doNormalization = doNormalization, factor=factor)$dataA1
    if(usingSeuratData) {
      dataX <- pca(dataX)
    } 
    else if (nrow(keeprows) > markerCap && markerCap != 0) {
      proteinIndices <- grepl("TotalSeqC", colnames(dataX))
      markerCapUpdated <- markerCap - length(which(proteinIndices))
      dataXtrimmed <- dataX[, !proteinIndices]
      sums <- sapply(1:ncol(dataXtrimmed), function(x) sum <- var(dataXtrimmed[, x]))
      test <- sort.int(sums, index.return = TRUE, decreasing = T)[[2]][1:markerCapUpdated]
      dataX <- dataX[, c(which(proteinIndices),test), drop = F]
    }
    tSNEdata3 <- Rtsne(dataX, dims=1, check_duplicates=F,verbose = T,perplexity = tSNEperplexity, pca = F) 
    tSNEmat1 <-  matrix(ncol=1,data=tSNEdata3$Y,dimnames=list(NULL,OneDtSNEname))
    colnames(tSNEmat1) <- OneDtSNEname
    hist(tSNEmat1, 100, main= paste("Histogram of",colnames(keeptable[factor] )))
    OneStSNEmat <- cbind(OneStSNEmat,tSNEmat1)
  }
  return(OneStSNEmat)
}

runIsomap <- function(dataTransformed = NULL) {
  print("running isomap")
  dis <- vegdist(dataTransformed, method="euclidean")
  ord <- isomap(dis,ndim=3,k=3, fragmentedOK = TRUE)
  usordata <-  matrix(ord$points[,1:3], ncol = 3)
  colnames(usordata) <- c("Isomap1","Isomap2","Isomap3")
  
  plot(usordata[, 1], usordata[, 2], pch=".", xlab="isomap1", ylab="isomap2", cex=0.1)
  return(usordata)
}

diffMap <- function(ManDist = F,
                    dataTransformed = NULL) {
  print("running diff map")
  if (ManDist){
    D <- D2.dist(dataTransformed, cov(dataTransformed))
  } else {
    D = dist(dataTransformed) # use Euclidean distance
  }
  
  dmap = diffuse(D, eps.val = epsilonCompute(D,  p = 0.01), neigen = 3, maxdim=3) # compute diffusion map & plot
  plot(dmap)
  DiffMapMat <- dmap$X[,1:3]
  plot(DiffMapMat[,1],DiffMapMat[,2],pch=".", xlab="DiffMap1", ylab="DiffMap2", cex=0.1)
  colnames(DiffMapMat) <- c("DiffMap1","DiffMap2","DiffMap3")
  return(DiffMapMat)
}

tSNE_3D <- function(dataTransformed = NULL,
                    tSNEperplexity = 30,
                    sparseScore = NULL) {
  print("running 3d tSNE")
  if(!is.null(sparseScore)) {
    dataTransformed <- sparseScore
  } else if(usingSeuratData) {dataTransformed <- pca(dataTransformed)}
  TDtSNEdata3 <- Rtsne(dataTransformed, dims=3, perplexity=tSNEperplexity, check_duplicates=F,verbose = T, pca = F)
  TDtSNEmat <- TDtSNEdata3$Y
  colnames(TDtSNEmat) <- c("ThreeDtSNE1","ThreeDtSNE2", "ThreeDtSNE3")
  plot(TDtSNEmat[, 1], TDtSNEmat[, 2], pch=".", xlab="3DtSNE1", ylab="3DtSNE2", cex=0.1)
  return(TDtSNEmat)
}

umap_2D <- function(dataTransformed = NULL,
                    n = 15,
                    mdist = 0.2,
                    metric = "euclidean",
                    ncomp = 2,
                    FFdata = NULL,
                    sparseScore = NULL){
  print("running UMAP")
  if(!is.null(sparseScore)) {
    dataTransformed <- sparseScore
  }else if(usingSeuratData) {dataTransformed <- pca(dataTransformed)}
  umapMat <- umap(X = dataTransformed, n_neighbors = n, n_components = ncomp, metric = metric, min_dist = mdist, verbose = T)
  
  columnMatch <- grep("UMAP", colnames(FFdata))
  if (length(columnMatch != 0)) {
    colnames(umapMat) <- c(colnames(FFdata)[columnMatch[length(columnMatch) - 1]], colnames(FFdata)[columnMatch[length(columnMatch)]])
  } else {
    colnames(umapMat) <- c("UMAP1","UMAP2")
  }
  plot(umapMat[, 1], umapMat[, 2], pch=".", xlab="UMAP1", ylab="UMAP2", cex=0.1)
  return(umapMat)
}

umap_3D <- function(dataTransformed = NULL,
                    n = 15,
                    mdist = 0.2,
                    metric = "euclidean",
                    ncomp = 3,
                    sparseScore = NULL) {
  print("running 3D UMAP")
  if(!is.null(sparseScore)) {
    dataTransformed <- sparseScore
  }else if(usingSeuratData) {dataTransformed <- pca(dataTransformed)}
  umapMat3D <- umap(X = dataTransformed, n_neighbors = n, n_components = ncomp, metric = metric, min_dist = mdist, verbose = T)
  colnames(umapMat3D) <- c("TD_UMAP1","TD_UMAP2","TD_UMAP3")
  plot(umapMat3D[, 1], umapMat3D[, 2], pch=".", xlab="TD_UMAP1", ylab="TD_UMAP2", cex=0.1)
  return(umapMat3D)
}

umap_oneSENSE <- function(FFdata = NULL,
                          keeptable = NULL,
                          markerCap = 0) {
  print("running one-SUMAP")
  OneSUMAPmat <- NULL
  for (factor in 4:ncol(keeptable)){
    OneDUMAPname <- paste0(colnames(keeptable)[factor],"_UMAP")
    keeprowbool <- sapply(keeptable[,factor], isTransform)
    keeprows <- subset(keeptable, keeprowbool)
    dataX <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
    dataX <- nfTransform(keeprows,dataX,FFdata,doNormalization = doNormalization, factor = factor)$dataA1
    if(usingSeuratData) {
      dataX <- pca(dataX)
    }
    else if (nrow(keeprows) > markerCap && markerCap != 0) {
      proteinIndices <- grepl("TotalSeqC", colnames(dataX))
      markerCapUpdated <- markerCap - length(which(proteinIndices))
      dataXtrimmed <- dataX[, !proteinIndices]
      sums <- sapply(1:ncol(dataXtrimmed), function(x) sum <- var(dataXtrimmed[, x]))
      test <- sort.int(sums, index.return = TRUE, decreasing = T)[[2]][1:markerCapUpdated]
      dataX <- dataX[, c(which(proteinIndices),test), drop = F]
    }
    
    n <- 15
    mdist <- 0.2
    metric <- "euclidean"
    ncomp <- 1
    ODumapMat <- umap(X = dataX, n_neighbors = n, n_components = ncomp, metric = metric, min_dist = mdist, verbose = T)
    ODumapMat1 <-  matrix(ncol=1,data=ODumapMat,dimnames=list(NULL,OneDUMAPname))
    colnames(ODumapMat1) <- OneDUMAPname
    hist(ODumapMat1, 100, main= paste("Histogram of",colnames(keeptable[factor] )))
    OneSUMAPmat <- cbind(OneSUMAPmat,ODumapMat1)
  }
  return(OneSUMAPmat)
}

pca <- function(dataTransformed = NULL) {
  print("running pca")
  mdpca <- prcomp(dataTransformed)
  
  score <- dataTransformed %*% mdpca$rotation
  eigs <- mdpca$sdev^2
  PCAsummary <- rbind(
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs))
  write.csv(PCAsummary, "PCA summary.csv")
  write.csv(mdpca$rotation,"PCAloading.csv")
  
  score <- score[,1:100]
  plot(score[, 1], score[, 2], pch=".", xlab="PC1", ylab="PC2", cex=0.1)
  
  return(score)
}

phenograph <- function(dataTransformed = NULL,
                       score2 = NULL,
                       kValue = 30){
  print("running phenograph")
  Rphenograph_out <- Rphenograph(dataTransformed, k = kValue)
  #return(Rphenograph_out)
  modularity(Rphenograph_out[[2]])
  membership(Rphenograph_out[[2]])
  #Rphenograph_cluster <- as.matrix(membership(Rphenograph_out[[2]]))
  Rpmat <- as.matrix(membership(Rphenograph_out[[2]]))
  #--------------------------------------------------
  # ggplot(data1, aes(x=Sepal.Length, y=Sepal.Width, shape=Rphenograph_cluster)) + geom_point(size = 3)+theme_bw()
  #--------------------------------------------------
  Rphenographmat <- matrix(ncol=1,data=Rpmat,dimnames=list(NULL,"Phenograph"))
  #colnames(Rphenographmat) <- "phenograph"
  
  #fixes missing names
  if(dim(Rphenographmat)[1] != dim(score2)[1]){
    rpnames <- as.numeric(row.names(Rpmat))
    newclustnum <- max(Rphenographmat)+1
    print(paste0((dim(score2)[1]-dim(Rphenographmat)[1]),
                 " events lost during running of Rphenograph - dunno why. Asigned as new cluster #",
                 newclustnum))
    losts <- setdiff(1:dim(score2)[1],rpnames)
    for(li in 1:length(losts)) {
      Rphenographmat <- append(Rphenographmat, newclustnum, after = (losts[li]-1) )
    }
    Rphenographmat <- matrix(ncol=1,data=Rphenographmat,dimnames=list(NULL,"Phenograph"))
  }
  return(Rphenographmat)
}

flowSOM <- function(dataTransformed = NULL,
                    MaxClusters = 0) {
  print("running FlowSOM")
  FlowSOMmat <- matrix(ncol=1,data=MetaClustering(dataTransformed,method="metaClustering_som",max=MaxClusters),
                       dimnames = list(NULL,"FlowSOM"))
  return(FlowSOMmat)
}

gatingDataPoints <- function(xAxis = 'UMAP1',
                             yAxis = 'UMAP2',
                             LoaderPATH = "fcs",
                             ceil = 10000,
                             useCSV = F,
                             cropped = T,
                             gatedPATH = "",
                             TransformOutput = F,
                             OutputSuffix = "out") {
  filePATH <- LoaderPATH
  if(!identical(gatedPATH, "")) {
    filePATH <- paste0(filePATH, "_", gatedPATH)
  }
  
  prepData <- prepFcsFolderData(LoaderPATH = filePATH, ceil = ceil, useCSV = F)
  
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  gatedModifier <- gatedPATH
  
  keeptable <- read.csv(FNnames)
  keeprowbool <- sapply(keeptable[,2], function(x) any(x=="y" | x=="c" | x=="l" | x=="a" | x=="n"|x=="as5"|x=="as150"|x=="x"|x=="ln"))
  keeprows <- subset(keeptable, keeprowbool)
  dataMat <- FFdata[,which (colnames(FFdata) %in% as.character(keeprows[,1]))]
  
  
  
  
  #if (useUMAP) {
  #  columnMatch <- grep("^UMAP", colnames(FFdata))
  #} else {
  #  columnMatch <- grep("^tSNE", colnames(FFdata))
  #}
  #labels <- c( colnames(FFdata)[ columnMatch[ length(columnMatch) - 1 ]], colnames(FFdata)[ columnMatch[ length(columnMatch) ]] )
  
  labels = c(xAxis, yAxis)
  
  dataMat <- cbind(dataMat, FFdata[,labels[1]], FFdata[,labels[2]])
  colnames(dataMat)[(ncol(dataMat) - 1):ncol(dataMat)] <- c(labels[1], labels[2])
  
  print("If using RStudio, hit 'Esc' to close the gate.")
  
  SAMPLE <- 20000
  if(SAMPLE > nrow(dataMat)) {SAMPLE <- nrow(dataMat)}
  
  gatedPoints <- gate_from_biplot(dataMat, labels[1], labels[2], sample = SAMPLE,
                                  axes = T)
  print(gatedPoints)
  print ('done gating')
  gatedDataBool <- sapply(gatedPoints, function(x) x==1)
  counter <- 0
  
  minX <- min(FFdata[gatedDataBool, ][,labels[1]])
  maxX <- max(FFdata[gatedDataBool, ][,labels[1]])
  minY <- min(FFdata[gatedDataBool, ][,labels[2]])
  maxY <- max(FFdata[gatedDataBool, ][,labels[2]])
  
  rescaledX <- ((FFdata[,labels[1]] - minX) / (maxX - minX)) * 9800 + 100
  rescaledY <- ((FFdata[,labels[2]] - minY) / (maxY - minY)) * 9800 + 100
  
  score2 <- cbind(rescaledX, rescaledY)
  
  if(cropped) {
    gatedModifier <- paste0(gatedModifier, "c")
    labels <- c(paste0(labels[1], "c"), paste0(labels[2], "c"))
    OutputSuffix <- "c"
  } else {
    gatedModifier <- paste0(gatedModifier, "s")
    labels <- c(paste0(labels[1], "s"), paste0(labels[2], "s"))
    OutputSuffix <- "s"
  }
  
  colnames(score2) <- c(labels[1], labels[2])
  print(colnames(score2))
  normalizedScore <- reverseTransform(FFdata = FFdata,
                                      nontransformedData = NULL,
                                      score2 = score2,
                                      TransformOutput = TransformOutput)
  
  NIscore <- normalizedScore[[1]]
  NewTransformBase <- normalizedScore[[2]]
  NewTransformNew <- normalizedScore[[3]]
  NewTransform <- normalizedScore[[4]]
  
  if (cropped) {
    FFdata <- FFdata[gatedDataBool, ]
    NIscore <- NIscore[gatedDataBool, ]
    score2 <- score2[gatedDataBool, ]
  }
  
  plot(score2[, 1], score2[, 2], xlab = labels[1], ylab = labels[2], xlim = c(0, 10000), ylim = c(0, 10000), pch=".")
  
  exportToFCS(useCSV = F,
              numFiles = numFiles,
              forNewFF = forNewFF,
              OrigNames = OrigNames,
              FFdata = FFdata,
              NIscore = NIscore,
              NewTransform = NewTransform,
              NewTransformNew = NewTransformNew,
              LoaderPATH = LoaderPATH,
              OutputSuffix = OutputSuffix,
              FcsFileNames = FcsFileNames)
}

GetCoords <- function(
  LoaderPATH ="All samples_out",  #Folder with samples (output folder from tSNE) 
  FNnames = "names.csv",
  yparam = "CD24",
  yparam2 = "Ly6GC",
  CoordsFN = "Coords.csv") {
  
  if(file.exists(CoordsFN)){
    Coords <- read.csv(CoordsFN, fileEncoding="UTF-8-BOM")
    rownames(Coords) <- Coords[,1]
    Coords <- Coords[,2, drop=F]
    return(Coords)
  }
  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = F)
  
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  keeptable <- read.csv(FNnames, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  keepnames <- keeptable[keeprowbool,1]
  
  data <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
  data <- cbind(data, FFdata[,which(colnames(FFdata) %in% colnames(keeptable[-1]))])
  data <- data[, !duplicated(colnames(data))]
  lgcl <- logicleTransform(w=0.25, t=16409, m=4.5, a=0)
  #lgcl <- logicleTransform(w=.1, t=4000, m=4.5, a=0)
  ilgcl <- inverseLogicleTransform(trans = lgcl)
  
  data1 <- apply(data, 2, lgcl)
  Params <- unique(keepnames)
  
  Cutoff <- numeric(length = length(Params))
  Coords <- data.frame(Cutoff)
  rownames(Coords) <- Params
  
  acceptLine <- function(prompt)
  {
    getGraphicsEvent(prompt = prompt, 
                     onMouseDown = NULL, onMouseMove = NULL,
                     onMouseUp = NULL, onKeybd = onKeybd, consolePrompt = "")
    Sys.sleep(0.01)
    return(keyPressed)
  }
  
  onKeybd <- function(key)
  {
    if (tolower(key) == "r"){
      keyPressed <<- F
    }
    if (tolower(key) == "y"){
      keyPressed <<- T
    }
    return(keyPressed)
  }
  
  print(keepnames)
  
  for(pname in keepnames){
    ExpX <- data1[,pname]
    ExpY <- data1[,yparam]
    
    if(pname == yparam) {
      ExpY <- data1[,yparam2]
    }
    
    print("[Press r to redraw, or y to move on.]")
    
    userSelectedCutoff <- F
    X11()
    while(!userSelectedCutoff) {
      plot(ExpX, ExpY, pch=".", xlab=pname, ylab=yparam, cex=0.1, ylim = c(-1,4.5), xlim = c(-1,4.5))
      tempXCoords <- locator(n=1, type = "l")$x
      abline(v = tempXCoords, col = "blue")
      userSelectedCutoff <- acceptLine(prompt = "Press r to redraw, or y to move on.")
    }
    Coords[pname, ] <- tempXCoords
    abline(v = Coords[pname, ], col = "blue")
    dev.off()
  }
  
  write.csv(Coords, CoordsFN)
  return(Coords)
}

developingGatingStrategies <- function(useCSV = F,
                                       sourceFcsFolder = "fcs",
                                       ceil = 10000,
                                       FNnames="trafficknames2.csv",
                                       gateBy = "Phenograph",
                                       categoryGateBy = "Genes",
                                       categoryToStrategize = "Proteins",
                                       biplotParameter = "UMAP",
                                       OneSENSE_parameter = "UMAP",
                                       category1 = "",
                                       category2 = "",
                                       kValue = 30) {
  prepData <- prepFcsFolderData(LoaderPATH = sourceFcsFolder, ceil = ceil, useCSV = F)
  
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  numFiles <- prepData$numFiles
  FcsFileNames <- prepData$FcsFileNames
  
  keeptable <- read.csv(FNnames, fileEncoding="UTF-8-BOM")
  keeprowbool <- sapply(keeptable[,2], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  data <- FFdata[, which(colnames(FFdata) %in% keeprows[,1])]
  
  nfTransOut <- nfTransform(keeprows,FFdata,FFdata, doNormalization = doNormalization)
  ####  Run algorithms
  dataTransformed <- nfTransOut$dataA1
  FFdataTransformed <- nfTransOut$dataB1
  
  #UMAP categories
  n <- 15
  mdist <- 0.2
  metric <- "euclidean"
  tSNEperplexity <- 30
  
  keeprowbool <- sapply(keeptable[, categoryGateBy], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  dataTransformed <- FFdataTransformed[, which(colnames(FFdata) %in% keeprows[,1])]
  
  if(biplotParameter == "UMAP") {
    if (paste0("UMAP_", categoryGateBy, 1) %in% FFdata && paste0("UMAP_", categoryGateBy, 2) %in% FFdata) {
      gatingData <- cbind(FFdata[, paste0("UMAP_", categoryGateBy, 1)], FFdata[, paste0("UMAP_", categoryGateBy, 1)])
    } else {
      gatingData <- umap_2D(dataTransformed, n, mdist, metric, ncomp = 2, FFdata)
    }
    colnames(gatingData) <- c(paste0("UMAP_", categoryGateBy, 1), paste0("UMAP_", categoryGateBy, 2))
  } else if (biplotParameter == "tSNE") {
    if (paste0("tSNE_", categoryGateBy, 1) %in% FFdata && paste0("tSNE_", categoryGateBy, 2) %in% FFdata) {
      gatingData <- cbind(FFdata[, paste0("tSNE_", categoryGateBy, 1)], FFdata[, paste0("tSNE_", categoryGateBy, 1)])
    } else {
      gatingData <- umap_2D(dataTransformed, n, mdist, metric, ncomp = 2, FFdata)
    }
    colnames(gatingData) <- c(paste0("tSNE_", categoryGateBy, 1), paste0("tSNE_", categoryGateBy, 2))
  } else if (biplotParameter == "OneSENSE") {
    gatingData <- cbind(FFdataTransformed[, paste0(category1, "_", OneSENSE_parameter)], FFdataTransformed[, paste0(category2, "_", OneSENSE_parameter)])
    colnames(gatingData) <- c(paste0(category1, "_", OneSENSE_parameter), paste0(category2, "_", OneSENSE_parameter))
  }
  
  nontransformedData <- NULL
  if (gateBy == "Phenograph") {
    if(paste0("Phenograph_", categoryGateBy) %in% FFdata) {
      nontransformedData <- FFdata[, paste0("Phenograph_", categoryGateBy)]
    } else {
      nontransformedData <- phenograph(dataTransformed = dataTransformed,
                                       score2 = gatingData,
                                       kValue = kValue)
    }
    colnames(nontransformedData) <- paste0("Phenograph_", categoryGateBy)
  }
  normalizedScore <- reverseTransform(FFdata = data,
                                      nontransformedData = nontransformedData,
                                      score2 = gatingData,
                                      TransformOutput = TransformOutput)
  
  NIscore <- normalizedScore[[1]]
  NewTransformBase <- normalizedScore[[2]]
  NewTransformNew <- normalizedScore[[3]]
  NewTransform <- normalizedScore[[4]]
  
  exportToFCS(useCSV = useCSV, numFiles = numFiles, forNewFF = forNewFF, OrigNames = OrigNames, FFdata = FFdata, NIscore = NIscore, NewTransform = NewTransform,
              NewTransformNew = NewTransformNew, LoaderPATH = sourceFcsFolder, OutputSuffix = "", FcsFileNames = FcsFileNames)
  
  if (gateBy == "Manual") {
    SAMPLE <- 20000
    if(SAMPLE > nrow(NIscore)) {SAMPLE <- nrow(NIscore)}
    gatedPoints <- gate_from_biplot(NIscore, colnames(NIscore)[1], colnames(NIscore)[2], xlim = c(0, 10000), ylim = c(0, 10000), sample = SAMPLE)
    gatedPoints[is.na(gatedPoints)] <- 0
    if (max(gatedPoints) > 1) {
      level <- readline("Which gate would you like to use? ")
    } else {
      level <- 1
    }
  } else if (gateBy == "Phenograph") {
    gatedPoints <- NIscore[, paste0("Phenograph_", categoryGateBy), drop = F]
    
    discreteRasters <- generateDiscreteRasters(discreteData = gatedPoints,
                                               FFdata = NIscore,
                                               clustParam = paste0("Phenograph_", categoryGateBy),
                                               PC1 = colnames(gatingData)[1],
                                               PC2 = colnames(gatingData)[2],
                                               MPColorList = MPColorList,
                                               labelClusters = T,
                                               pch = pch,
                                               cex = cex
    )
    X11()
    plot.new()
    grid.raster(discreteRasters[[1]]$raster.main)
    level <- readline("What cluster would you like to gate on? ")
  }
  
  keeprowbool <- sapply(keeptable[, categoryToStrategize], isTransform)
  keeprows <- subset(keeptable, keeprowbool)
  assign(paste0("data_", categoryToStrategize), FFdataTransformed[,which (colnames(FFdata) %in% keeprows[,1])])
  
  gatingStrategy <- hypergate(xp = get(paste0("data_", categoryToStrategize)), gate_vector = as.vector(gatedPoints), level = level, verbose = T)
  plot_gating_strategy(gate = gatingStrategy, xp = get(paste0("data_", categoryToStrategize)), gate_vector = gatedPoints, 
                       level = level, highlight = "firebrick3", path = "test")
  contributions <- channels_contributions(gate=gatingStrategy,xp=get(paste0("data_", categoryToStrategize)),level=level,gate_vector = gatedPoints)
  significant_channels=names(contributions)[contributions>0.01]
  hg_final<-reoptimize_strategy(gatingStrategy,gate_vector = gatedPoints,significant_channels,get(paste0("data_", categoryToStrategize)),level=level)
  
  plot_gating_strategy(gate = hg_final, xp = get(paste0("data_", categoryToStrategize)), gate_vector = gatedPoints, 
                       level = level, highlight = "firebrick3", path = "test_reoptimized")
}

ThreeDPlot <- function(LoaderPATH ="fcsfiles",  FNnames = "names.csv",
                       ceil = 300, OutputSuffix = "out", labelClusters = labelClusters,
                       parN1, parN2, parN3, colparN = "InFile", Ncolors = 4, ColorList =c("orange","red","green","blue") ){
  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = F)
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  
  BaseFNs <- sapply(strsplit(prepData$FcsFileNames, split ="\\."), "[", 1)
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  
  keeptable <- read.csv(FNnames)
  keeprowbool <- sapply(keeptable[,2], function(x) any(x=="y" | x=="c" | x=="l" | x=="a" | x=="n" | x=="as5"| x=="as150"|x=="x"|x=="ln"))
  keeprows <- subset(keeptable, keeprowbool)
  
  
  data1u <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
  data1 <- nfTransform(keeprows, data1u, data1u, doNormalization = doNormalization)$dataA1
  data2u <- FFdata[,transList == "2"] 
  data2 <- apply(data2u, 2, logicleTransform(w=0.25, t=16409, m=4.5, a=0))
  data3 <- FFdata[,transList == "3",drop=F]
  
  #backup plan to salvage colparN data for the user 
  if (!colparN %in% colnames(data)){
    data <- cbind(data, FFdata[,colparN])
    colnames(data)[dim(data)[2]] <- colparN
  }
  
  c <- data[,colparN]
  if (Ncolors == 0){
    Ncolors <- max(data[,colparN])
    c <- cut(c, breaks=0:Ncolors)
  } else c<- (c-min(c))/(max(c)-min(c))*Ncolors
  
  print(paste("Number of Colors",Ncolors))  
  Ddata <- data[,c(parN1, parN2, parN3, colparN)]
  palette <- colorRampPalette(ColorList)(n=Ncolors)
  
  # for cluster numbers only
  #if(colparN == "cluster") c<-sapply(c,ilgcl)
  
  #scatterplot3d(Ddata[,-4], pch = ".", highlight.3d=TRUE, cex.symbols=1)
  x<-Ddata[,1]
  y<-Ddata[,2]
  z<-Ddata[,3]
  x1 <- (x - min(x))/(max(x) - min(x))
  y1 <- (y - min(y))/(max(y) - min(y))
  z1 <- (z - min(z))/(max(z) - min(z))
  
  randOrder <- sample(length(x1))
  plot3d(x=x1[randOrder], y=y1[randOrder], z=z1[randOrder], 
         xlab = " ", ylab = " ", zlab = " ",
         col = palette[as.numeric(c)[randOrder]], type = 'p',  size=1, box =F,
         axes = F)
  
  if(labelClusters){
    
    for(cluster in 1:max(as.numeric(c)))
    {
      
      xclustlab <- median(x1[as.numeric(c)==cluster])
      yclustlab <- median(y1[as.numeric(c)==cluster])
      zclustlab <- median(z1[as.numeric(c)==cluster])
      text3d(xclustlab,yclustlab, zclustlab, text=paste0(cluster), cex=0.75, col="black")
    }
    
  }
  
  par3d(windowRect = c(0, 0, 1200, 1200)) # make the window large
  degrees <- seq(1,360, by = 2) # a sequence from 1 to 360
  
  suppressWarnings(dir.create(paste0(LoaderPATH,"_",OutputSuffix)))
  FNresult <- paste0(LoaderPATH,"_",OutputSuffix,"/tSNE3D","_",OutputSuffix,".fcs")
  
  for(i in 1:length(degrees)){
    view3d(degrees[i], phi = 0) # pick the angle of view
    rgl.snapshot(paste(paste(FNresult, "-", 
                             formatC(i, digits = 3, flag = "0"), sep = ""), "png", sep = "."))
  }
}

TwoDPlot <- function(LoaderPATH = sourceFcsForMeaningPlot,
                     FNnames = FNnames,
                     ceil = ceil,
                     OutputSuffix = TwoDOutputSuffix,
                     parN1 = TwoDparN1,
                     parN2 = TwoDparN2, 
                     Ncolors = TwoDNcolors,
                     labelClusters = TwoDlabelClusters,
                     colparN = TwoDColorBy,
                     pch = TwoDpch,
                     cex = TwoDcex,
                     resolution = TwoDresolution,
                     ColorList = TwoDColorList,
                     DoOneForEach2D = F){
  
  prepData <- prepFcsFolderData(LoaderPATH=LoaderPATH, ceil = ceil, useCSV = useCSV)
  
  #data <- prepData$data
  FFdata<- prepData$FFdata
  OrigNames <- prepData$OrigNames
  forNewFF <- prepData$forNewFF
  
  BaseFNs <- sapply(strsplit(prepData$FcsFileNames, split ="\\."), "[", 1)
  transList <- strsplit(unlist(forNewFF@description$TF),"")[[1]]
  
  keeptable <- read.csv(FNnames)
  keeprowbool <- sapply(keeptable[,2], function(x) any(x=="y" | x=="c" | x=="l" | x=="a" | x=="n"|x=="as5"|x=="as150"|x=="x"|x=="ln"))
  keeprows <- subset(keeptable, keeprowbool)
  
  data1u <- FFdata[,which (colnames(FFdata) %in% keeprows[,1])]
  data1 <- nfTransform(keeprows, data1u, data1u, doNormalization = doNormalization)$dataA1
  data2u <- FFdata[,transList == "2"] 
  data2 <- apply(data2u, 2, logicleTransform(w=0.25, t=16409, m=4.5, a=0))
  data3 <- FFdata[,transList == "3",drop=F]
  
  data <- cbind(data1, data2, data3, InFile = FFdata[,"InFile"])
  
  #backup plan to salvage colparN data for the user 
  if (!colparN %in% colnames(data)){
    data <- cbind(data, FFdata[,colparN])
    colnames(data)[dim(data)[2]] <- colparN
  }
  c <- data[,colparN]
  
  if (Ncolors == 0){
    Ncolors <- max(data[,colparN])
    c <- cut(c, breaks=0:Ncolors)
  } else {
    print("shading colors")
    print(paste("max",max(c),"min",min(c)))
    c<- ((c+1))/(4.5)*Ncolors 
  }
  print(paste("Number of Colors",Ncolors))  
  Ddata <- data[,c(parN1, parN2, colparN)]
  Ddata[, colparN] <- Ddata[,colparN]
  palette <- colorRampPalette(ColorList)(n=Ncolors)
  
  
  #scatterplot3d(Ddata[,-4], pch = ".", highlight.3d=TRUE, cex.symbols=1)
  x<-Ddata[,1]
  y<-Ddata[,2]
  
  x1 <- (x - min(x))/(max(x) - min(x))
  y1 <- (y - min(y))/(max(y) - min(y))
  
  FNresult <- paste0(LoaderPATH,"_",OutputSuffix,".svg")
  #png(FNresult,res=resolution,height=480*resolution/72,width=480*resolution/72)
  svg(FNresult)
  par("bty"="l")
  
  # randOrder <- sample(length(x1))
  # col = palette[as.numeric(c)[randOrder]]
  
  col = palette[as.numeric(c)]
  
  # plot(x=x1[randOrder], y=y1[randOrder], cex = cex, pch = pch,
  #      xlab = parN1, ylab = parN2, 
  #      col = col,
  #      axes = T)
  
  plot(x=x1, y=y1, cex = cex, pch = pch,
       xlab = parN1, ylab = parN2,
       col = col,
       axes = T)
  
  if (colparN == "InFile") {
    legendNames <- c(1:max(FFdata[, "InFile"]))
  } else if (colparN == "groupingCol") {
    legendNames <- groupingType$groupingDict
  }
  
  if(labelClusters && (colparN == "Phenograph" || colparN =="FlowSOM")){
    for(cluster in 1:max(as.numeric(c))) {
      xclustlab <- median(x1[as.numeric(c)==cluster])
      yclustlab <- median(y1[as.numeric(c)==cluster])
      text(xclustlab,yclustlab, labels=paste0(cluster), cex=0.8, col="black")
    } 
  }
  
  dev.off()
  
  if(DoOneForEach2D)
  {
    for(infileNum in 1:max(as.numeric(data[,"InFile"]) ))
    {
      FNresult <- paste0(LoaderPATH,"_",BaseFNs[infileNum],"_",OutputSuffix,".png")
      png(FNresult,res=resolution,height=480*resolution/72,width=480*resolution/72)
      par("bty"="l")
      
      pts <- data[,"InFile"] == infileNum
      x1a <- x1[pts]
      y1a <- y1[pts]
      ca <- c[pts]
      
      randOrder <- sample(length(x1a))
      
      plot(x=x1a[randOrder], y=y1a[randOrder], cex = cex, pch = pch,
           xlab = parN1, ylab = parN2, 
           col = palette[as.numeric(ca)[randOrder]], 
           axes = T)
      
      if(labelClusters){
        for(cluster in 1:max(as.numeric(c)))
        {
          xclustlab <- median(x1[as.numeric(c)==cluster])
          yclustlab <- median(y1[as.numeric(c)==cluster])
          text(xclustlab,yclustlab, labels=paste0(cluster), cex=0.8, col="black")
        }
      }
      
      dev.off()
    }
  }
}