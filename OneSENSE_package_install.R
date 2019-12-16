if(!require(Rcpp)){
  install.packages("Rcpp")
  library(Rcpp)
}
if (!requireNamespace("BiocManager", quietly = TRUE))
{
  install.packages("BiocManager")
}
if (!require(Biobase)) { 
  BiocManager::install("Biobase")
  library(Biobase)
} 
if (!require(flowCore)) { 
  if(!require(robustbase)){
    install.packages("robustbase")
    library(robustbase)
  }
  BiocManager::install("flowCore")
  library(flowCore)
} 
if (!require(Rtsne)) {
  install.packages("Rtsne")
  library(Rtsne)
}
if (!require(reshape)) {
  install.packages("reshape")
  library(reshape)
}
if (!require(Rphenograph)) {
  if(!require(devtools)){
    install.packages("devtools")
  }
  options(buildtools.check = function(action) TRUE )
  install_github("ebecht/Rphenograph")
  library(Rphenograph)
  options(buildtools.check = function(action) FALSE )
}
if(!require(uwot)) {
  install.packages("uwot")
  library(uwot)
}
if(!require(hypergate)){
  install.packages("hypergate")
  library(hypergate)
}
if (!require(FlowSOM)) {
  BiocManager::install("FlowSOM")
  library(FlowSOM)
}
if (!require(vegan)) { 
  install.packages("vegan")    
  library(vegan)
} 
if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}
if (!require(diffusionMap)) { 
  install.packages("diffusionMap")    
  library(diffusionMap)
} 
if (!require(gplots)) {
  install.packages("gplots", dependencies=TRUE, 
                   repos="http://cran.cnr.Berkeley.edu")
}
if (!require(gridExtra)) {
  install.packages("gridExtra")
  library(gridExtra)
}
if (!require(grid)) {
  install.packages("grid")
  library(grid)
}
if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
} 
if (!require(scatterplot3d)) {
  install.packages("scatterplot3d")
  library(scatterplot3d)
}
if (!require(rgl)) {
  BiocManager::install("rgl")
  library(rgl)
}
sapply(c("png","raster"),function(package){
  if(!require(package,character.only=T)){
    install.packages(pkgs=package)
    library(package,character.only=T)
  }
})