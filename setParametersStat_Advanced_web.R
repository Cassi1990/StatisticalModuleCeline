# for compatibility with R local script, set every variable to a boolean depending on whether it exists
# only to be run when the code is called from the webportal or GenePattern

if(!exists("refName")) refName <- "" 
if(!exists("normDataTable")) normDataTable <- "" 
if(!exists("descriptionFile")) descriptionFile <- ""
keepAnnotation <- TRUE
html <- TRUE
defaultContr <- FALSE # default contrasts are written in the matfile
if(!exists("matfileName")) matfileName <- NULL
cutOffTable <- exists("cutOffTable")
if(!exists("cutOffPval")) cutOffPval <- 0.5
if(!exists("cutOfflogFC")) cutOfflogFC <- 2
if(!exists("cutOffAveExpr")) cutOffAveExpr <- 7
plotPvalHist <- exists("plotPvalHist")
plotFCHist <- exists("plotFCHist")
plotVennDiagram <- exists("plotVennDiagram")
plotVolcanoPlot <- exists("plotVolcanoPlot")
summaryTable <- exists("summaryTable")
if(!exists("pvaluelist")) {
  pvaluelist <- c(0.001,0.01,0.05, 0.1)
}else{
  pvaluelist <- as.numeric(unlist(strsplit(gsub(" ", "", pvaluelist), ",")))
}
if(!exists("adjpvaluelist")) {
  adjpvaluelist <- c(0.05)
}else{
  adjpvaluelist <- as.numeric(unlist(strsplit(gsub(" ", "", adjpvaluelist), ",")))
}	
if(!exists("foldchangelist")) {
  foldchangelist <- c(1.1,1.5,3)
}else{
  foldchangelist <- as.numeric(unlist(strsplit(gsub(" ", "", foldchangelist), ",")))
}	

print ("Parameters have been registered")

# PAREMETERS DESCRIPTION : see arrayAnalysisStat.R