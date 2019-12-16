#All OneSENSE scripts prepared by Evan Newell and Timothy Bi
#Contact: Fred Hutchinson Cancer Research Center

#Evan Newell 2015; Ming 2016; Etienne Becht 2016; Evan 2017, 2018

## New Names.csv file sytem:  using y to indicate channels to use still works but also indicates to use default cytof transform
# Now if "f" is indicated default (old version) fluor logicle params will be used (not recommended), "a" instructs auto-logicle transform (recommended for fluor data)
# in addition, "c" also means cytof transform, "l" linear normalized to max = 4.5, "n" linear normalized min=0, max=4.5, "as5" ArcSinh x/5 (for cytof), "as150" ArcSinh x/150 (for fluor)
# Additional transform options can be implemented with this sytem.

# Startup Settings:

if(!require(devtools)){
  install.packages("devtools")
}
library(devtools) # A workaround to Rtools incompatibility with more recent versions of R.

source("Scripts/package_install.R")





### USER INPUT STARTS HERE

# Parameters are split into basic and advanced.
## Basic Parameters:

##### Main Parameters --------------------------------------------------------------------------
sourceFcsFolder ="fcs"# Folder with samples (even 1 sample ok)
ceil =1000000# Number of events to take PER SAMPLE (unless sample has less)
FNnames="names.csv"# Parameters to analyze in csv file
OutputSuffix ="out"# Adds this to all new directories and file names

#Advanced: 
useCSV =F# T if using .csv files instead of fcs files - will make csv output also
TransformOutput =F# F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense
usingSeuratData =F# Using scaled data derived from the Seurat package
doNormalization =F
 
##### Choose functions to run-------------------------------------------------------------------
RunAlgorithms = T# change to F if you just want to make meaning plots or heatplots using output files (no need to rerun tSNE etc.)
RunMeaningPlot = T
RunFileHeatplots = T
RunOneSENSEPlot = F
RunGating =F #Gate (i.e. crop or scale) after 
RunHypergate =F

Run2DPlot =F
Run3DPlot =F

#((RunAlgorithms))##### Algorithms -------------------------------------------------------------------------------
DotSNE =F#
    tSNEperplexity =30# ((DotSNE)) default = 30; increased perplexity => increased spread
DoUMAP =T#
DoPhenograph =T# clusters cells using Rphenograph. 
    kValue =30#((DoPhenograph)) default = 30; k is low => more clusters and inversely
    
# For Categorical Dimensionality Reduction
DoOneSENSE =F# needs a modified names.csv file with column names for each category
DoOneSUMAP =F# umap version of One-SENSE
    markerCap =0#((DoOneSUMAP)) The max numbers of markers included when performing oneSENSE analysis (useful if all genes are empty).
                           # Set to 0 for unlimited

#Other Algorithms:
DoManualGrouping =F# Manually group FILES by a user-specified parameter (useful if multiple files are from same tissue type, etc.)
    GroupingOutput =""#((DoManualGrouping))Output name for grouping file
    GroupingName =c("TissueType")#((DoManualGrouping))
Do3DtSNE =F# May be useful for 3D visualization, but tends to be slow and clunky
Do3DUMAP =F
DoFlowSOM =F# Alternative clustering algorithm to Phenograph
MaxClusters =30
DoIsomap =F# Super slow, unuseable
DoDiffMap =F# same
DoPCA =F
RunOnServer =F

#((RunMeaningPlot))##### Meaningplot Parameters ------------------------------------------------------------------
Xaxis ="UMAP1"## examples: tSNE1 or UMAP1  Add * to color by an additional round of tSNE or UMAP
Yaxis ="UMAP2"## examples: tSNE2 or UMAP2
prefix =paste0("_UMAP")#edit that last part to give this plot a name
DoOneForEach =F#T => will do meaning plots for each individual file, F => for the concat file
DoMontagePlot =F
DoFileMontage =F
    plotParameter ="Density"
DoOverlay =F
    baseFile ="HD1Frac39-40_out.fcs"

#Discrete Plot Parameters
plotClusters =T#((DoPhenograph))
    clusteringParameter ="Phenograph"
    labelClusters =T
plotInFile =F
plotGroupings =F
plotDensity =T

#((RunFileHeatplots))### FileHeatplot Parameters -----------------------------------------------------------------------
DoPositiveFrequency=F
    yparam ="CD45RA"
    yparam2 ="CD4"
DoClusters =T# makes heatmaps to describe phenograph clusters - median marker intensity and frequencies (lin and log) within samples
    clustParam ="Phenograph"# e.g., "Phenograph" or "FlowSOM"

DoIndividualCell =F#For single-cell heatplots. Avoids medians.

#((RunOneSENSEPlot))### OneSensePlot Parameters ---------------------------------------------------------------
PlotAllCategories =F#Plots all the parameters against each other iteratively. If false, please specify two categories below.
    category1 ="TCR"
    category2 ="Phenotype"
HeatmapProperty ="Medians"#Can either be by medians or by positive frequencies (based on specified coordinates)
Bins =500#Number of "bins" for each OneSENSE axis heatplot
DoChisqTest =F
    chisqBins =9
OneSenseResolution =500#Can lower to speed up plotting speed if have a lot of categories
algorithmUsed ="tSNE"#The dimensionality reduction algorithm used to generate OneSENSE data (usually UMAP or tSNE, or both, capitalization matters!)
colorBy =c("clone_size")

#((RunGating))### Gating Parameters -------------------------------------------------------------------------
#### Gating Parameters -------------------------------------------------------------------------
cropped = T              #Whether you want to crop(T) or simply scale(F) the data points that you select with the gate.
gateXaxis = "Phenograph*"
gateYaxis = "CD8|Proteinc*"

#((RunHypergate))### Hypergate Strategies Parameters ---------------------------------------------------
gateBy ="Manual"# Manual, Phenograph, FlowSOM, etc.
categoryGateBy ="Genes"
categoryToStrategize ="Proteins"
biplotParameter ="UMAP"






#######  Advanced #######  #######  Advanced #######  

#((RunMeaningPlot))### Meaningplot Parameters ---------------------------------------------------------------------
#Plot Image Settings
resolution=150
cex=0.5
pch=16
palette=c("black","blue","green","yellow","red")
color.scale.type="relative"# choose "relative" or "absolute"

#Discrete plot color settings
MPNColors =0# Number of colors to use for the plots of discrete variables (set to 0 for max)
#qual_col_pals =brewer.pal.info[brewer.pal.info$category=='qual',]
#MPColorList =unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
MPColorList = brewer.pal(12, "Paired")

#Advanced:
MeaningTopPercentile =1
MeaningBotPercentile =0
sourceFcsForMeaningPlot =paste0(sourceFcsFolder,"_",OutputSuffix)

#((RunFileHeatplots))#### FileHeatplot Settings ---------------------------------------------------------------------
fileHPoutputsuffix =paste0(OutputSuffix,"HP")# Edit last part for unique suffix
HPpalette=c("black","blue","lightblue","green","yellow","darkorange","darkred")#adjust coloring of heatplot intensities

#((Run2DPlot))#### 2D Plot Settings --------------------------------------------------------------------------
TwoDNcolors =0# set to zero to have this automatically set to maximum number
TwoDpch ="."
TwoDcex =3
TwoDresolution =600

# code for making random colors, good for clustering
qual_col_pals =brewer.pal.info[brewer.pal.info$category=='qual',]
paletteFor2D =unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
paletteFor2D <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")
#same as for meaning plots
#paletteFor2D = c("black","blue","green","yellow","red")
TwoDColorList =paletteFor2D# can choose colors if numbers match otherwise it will interpolate
# ColorList = c("colorX", "colorY", ...) put n color names if n clusters in order to manually define the colors of the clusters
# Color palette in R: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)

#((Run3DPlot))### 3D Plot Settings --------------------------------------------------------------------------
Ncolors =0# set to zero to have this automatically set to maximum number
# code for making random colors, good for clustering
qual_col_pals =brewer.pal.info[brewer.pal.info$category=='qual',]
paletteFor3D =unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))
ColorList =paletteFor3D# can choose colors if numbers match otherwise it will interpolate
# ColorList = c("colorX", "colorY", ...) put n color names if n clusters in order to manually define the colors of the clusters
# Color palette in R: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf)


#### 2D Plot Parameters ------------------------------------------------------------------------
TwoDColorBy = "Phenograph"   # e.g. "FlowSOM", "Phenograph", "InFile", "GroupingType"
TwoDparN1 = "UMAP1"          # e.g., UMAP1 or tSNE1
TwoDparN2 = "UMAP2"          # e.g., UMAP2 or tSNE
TwoDlabelClusters = T
TwoDOutputSuffix = "UMAP"


#### 3D Plot Parameters ------------------------------------------------------------------------
Color3Dby = "Phenograph"     # e.g. "FlowSOM", "Phenograph", use "InFile" to color by sample file
parN1 = "TD_UMAP1"           # e.g., TD_UMAP1 or ThreeDtSNE1
parN2 = "TD_UMAP2"           # e.g., TD_UMAP2 or ThreeDtSNE2
parN3 = "TD_UMAP3"           # e.g., TD_UMAP3 or ThreeDtSNE3
labelClusters = T
TDOutputSuffix = "3D-tSNE"







###################################### EXECUTE FUNCTIONS ###################################################

##### ##### No edit below here ## No edit below here ## No edit below here #### ####

source('Scripts/newFunctions.R')




if (RunAlgorithms) 
  FCStSNEone(LoaderPATH = sourceFcsFolder, 
             useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output also
             TransformOutput = TransformOutput, # F Puts derived data on linear scale 1-10k - best for viewing in FJ10 otherwise T is setup for FJ9. Talk to Evan about implications for one-sense
             ceil = ceil, #number of events to take per sample (unless sample has less)
             FNnames=FNnames,#Parameters to analyze in csv file
             OutputSuffix = OutputSuffix, # Adds this to all new directories and file names
             DotSNE = DotSNE,
             tSNEperplexity = tSNEperplexity, #default = 30; increased perplexity => increased spread
             DoOneSENSE = DoOneSENSE, #needs a modified names.csv file with column names for each category
             DoPhenograph = DoPhenograph, #clusters cells using Rphenograpy
             kValue = kValue, #default = 30; k is low => more clusters and inversely
             DoFlowSOM = DoFlowSOM,
             MaxClusters = MaxClusters,
             DoIsomap = DoIsomap,
             DoDiffMap = DoDiffMap,
             DoPCA = DoPCA, 
             Do3DtSNE = Do3DtSNE, #T for running 3D tSNE
             DoUMAP = DoUMAP, #still in prep- will add many parameters for this
             Do3DUMAP = Do3DUMAP,
             DoOneSUMAP = DoOneSUMAP,
             DoManualGrouping = DoManualGrouping,
             GroupingName = GroupingName,
             GroupingOutput = GroupingOutput,
             markerCap = markerCap) 


## Still in prep:
if (RunMeaningPlot)
  meaningPlot(LoaderPATH =sourceFcsForMeaningPlot,
              useCSV = useCSV, # T if using .csv files instead of fcs files - will make csv output also
              FNnames = FNnames,
              ceil = ceil,
              TransformOutput = TransformOutput,
              MeaningTopPercentile = MeaningTopPercentile,
              MeaningBotPercentile = MeaningBotPercentile,
              PC1 = Xaxis,
              PC2 = Yaxis,
              DoOneForEach = DoOneForEach,
              prefix = prefix,
              palette=palette,
              color.scale.type=color.scale.type,
              resolution=resolution,
              cex=cex,
              pch=pch,
              plotClusters = plotClusters,
              plotInFile = plotInFile,
              plotGroupings = plotGroupings,
              MPNColors = MPNColors,
              MPColorList = MPColorList,
              GroupingName = GroupingName,
              GroupingOutput = GroupingOutput,
              labelClusters = labelClusters,
              clusteringParameter = clusteringParameter,
              DoMontagePlot = DoMontagePlot,
              DoFileMontage = DoFileMontage,
              plotDensity = plotDensity,
              parameterName = plotParameter,
              DoOverlay = DoOverlay,
              baseFile = baseFile)

if(RunFileHeatplots)
  fileHeatplot (LoaderPATH = sourceFcsForMeaningPlot,
                ceil = ceil,
                TransformOutput = TransformOutput,
                FNnames = FNnames,
                OutputSuffix = fileHPoutputsuffix,
                palette=HPpalette,
                DoClusters=DoClusters,
                DoPositiveFrequency=DoPositiveFrequency,
                clustParam=clustParam,
                yparam=yparam,
                yparam2=yparam2)


if(Run2DPlot){
  TwoDPlot (LoaderPATH =sourceFcsForMeaningPlot,
            FNnames = FNnames,
            ceil = ceil, 
            OutputSuffix = TwoDOutputSuffix,
            parN1 = TwoDparN1,
            parN2 = TwoDparN2, 
            Ncolors = TwoDNcolors,
            colparN = TwoDColorBy,
            labelClusters = TwoDlabelClusters,
            pch = TwoDpch,
            cex = TwoDcex,
            resolution = TwoDresolution,
            ColorList = TwoDColorList)  # can use "InFile" to color by sample file
}

if(Run3DPlot)
  ThreeDPlot (LoaderPATH =sourceFcsForMeaningPlot,
              FNnames = FNnames,
              ceil = ceil, 
              OutputSuffix = TDOutputSuffix,
              parN1 = parN1,
              parN2 = parN2, 
              parN3= parN3, 
              labelClusters = labelClusters,
              Ncolors = Ncolors,
              colparN = Color3Dby,
              ColorList = ColorList)  # can use "InFile" to color by sample file

if(RunGating) {
  gatingDataPoints(xAxis = gateXaxis,
                   yAxis = gateYaxis,
                   LoaderPATH = sourceFcsForMeaningPlot,
                   ceil = ceil,
                   useCSV = useCSV,
                   cropped = cropped,
                   gatedPATH = "",
                   TransformOutput = TransformOutput,
                   OutputSuffix = OutputSuffix)
}

if(RunOneSENSEPlot) {
  oneSENSEplotter(LoaderPATH = sourceFcsForMeaningPlot,
                  Bins = Bins,
                  FNnames = FNnames,
                  resolution = OneSenseResolution,
                  category1 = category1,
                  category2 = category2,
                  HeatmapProperty = HeatmapProperty,
                  ceil = ceil,
                  useCSV = useCSV,
                  algorithmUsed = algorithmUsed,
                  PlotAllCategories = PlotAllCategories,
                  markerCap = markerCap,
                  colorBy = colorBy,
                  DoChisqTest = DoChisqTest,
                  chisqBins = chisqBins,
                  biplotColors = MPNColors,
                  biplotColorList = MPColorList)
}

if(RunHypergate) {
  developingGatingStrategies(useCSV = useCSV,
                             sourceFcsFolder = sourceFcsFolder,
                             ceil = ceil,
                             FNnames = FNnames,
                             gateBy = gateBy,
                             categoryGateBy = categoryGateBy,
                             categoryToStrategize = categoryToStrategize,
                             biplotParameter = biplotParameter)
}

print("Completed without errors!")

##### Useful script for making names.csv file:
# if (!require(flowCore)) { 
#   source("http://bioconductor.org/biocLite.R")
#   biocLite("flowCore")
# } 
# 
# inFileN = "test.fcs"        #FCS file to make names file for 
# outNamesCsv = "names.csv"   #Output csv file name 
# 
# FF = read.FCS(inFileN)
# colNames = FF@parameters$desc
# empties = which(is.na(colNames) | colNames== " ")
# colNames[empties] = FF@parameters$name[empties]
# write.csv(colNames, outNamesCsv, row.names = F)