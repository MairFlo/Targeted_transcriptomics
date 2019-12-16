#All OneSENSE scripts prepared by Evan Newell and Timothy Bi
#Contact: Fred Hutchinson Cancer Research Center

#Put your FCS file in the "fcs" folder first, then type in the filename of your FCS file below, and the output you want.
#Don't forget to set your working directory.

if (!require(flowCore)) { 
  source("http://bioconductor.org/biocLite.R")
  biocLite("flowCore")
} 


## Parameters (modify)
folderName = "fcs"                        # Folder containing FCS files
outNamesCsv <- "names.csv"       # Output name of the folder


## Actual code (DO NOT MODIFY UNLESS YOU KNOW WHAT YOU ARE DOING)
inFileN <- list.files(paste0(folderName, "/"))[1]       #Input FCS file to make names file for 
inFileExtension <- strsplit(inFileN, "\\.")[[1]][2]
if (inFileExtension == "fcs") {
  FF <- read.FCS(paste0(folderName, "/", inFileN), emptyValue = F)
  colNames <- FF@parameters$desc
  empties <- which(is.na(colNames) | colNames== " ")
  colNames[empties] <- FF@parameters$name[empties]
  namesFile <- cbind(Parameters = colNames, Transform = "", CategoryForOneSense = "")
} else if (inFileExtension == "csv"){
  FF <- read.csv(paste0(folderName, "/", inFileN), check.names = F)
  colNames <- colnames(FF)
  namesFile <- cbind(Parameters = colNames, Transform = "", CategoryForOneSense = "")
}

#Fixup column names

write.csv(namesFile, outNamesCsv, row.names = F)




# After you have created the names CSV, open it and type in which columns of data you want to select.
# To do a transform of the data, see below:

# Using y to indicate channels to use still works but also indicates to use default cytof transform
# Now if "f" is indicated default (old version) fluor logicle params will be used (not recommended),
# "a" instructs auto-logicle transform (recommended for fluor data). In addition, "c" also means cytof transform,
# "l" linear normalized to max = 4.5, "n" linear normalized min=0, max=4.5, "as5" ArcSinh x/5 (for cytof),
# "as150" ArcSinh x/150 (for fluor)