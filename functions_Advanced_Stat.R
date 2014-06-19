
source("C:/Users/Celine Adriaans/Documents/University/Major Internship/Extending Stat Module/functions_Stat.R")

normDataTable <- read.delim("C:/Users/Celine Adriaans/Documents/University/Major Internship/Extending Stat Module/Real documents/RMANormData_.txt", sep="\t", as.is=T) 
descriptionFile <- read.delim("C:/Users/Celine Adriaans/Documents/University/Major Internship/Extending Stat Module/Real documents/description_set1_extended.txt",sep="\t",as.is=T)

require("limma", quietly = TRUE)
require("gdata", quietly = TRUE) #trim function
print ("Libraries has been loaded")

#################################
## compute Advanced Statistics ##
#################################

Covariates_String <- "Time,n; Tissue,f; Treatment,f"
Interaction_String <-"Time*Tissue, Tissue*Treatment, Treatment*Time"
Covariates <- TRUE
Interaction <-TRUE
Paired_data <-"Person"
defaultContr<-TRUE


computeAdvancedStatistics <-function(normDataTable, descriptionFile, Covariates=TRUE, Covariates_String, Interaction=TRUE, Interaction_String, 
                                     Paired_data=TRUE, Paired_String, defaultContr=TRUE, CreateVolcanoPlot=TRUE, CreateVennPlot=TRUE, 
                                     matfileName=NULL, keepAnnotation=FALSE) {
  
  # Load the data	
  if(is.null(normDataTable)) stop("normDataTable has to be provided")
  if(is.null(descriptionFile)) stop("descriptionFile has to be provided")
  if(is.null(dim(normDataTable))) { # this is not a data.frame or matrix
    extension<-strsplit(normDataTable,"\\.")
    extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")
    try(ndt<-readExtFile(normDataTable, extension))
    if(!exists("ndt")) stop("normDataTable format was not recognized")
    normDataTable <- ndt
    rm(ndt)
  }
  
  descdata<-descriptionFile
  rownames(descdata) <- descdata$Description
  size <-length(normDataTable)
  print(arrayNames<-descdata[,1])
 

  #in R, proper names should fulfill criteria, which already have been automatically applied to the column (array) names when reading the data
  arrayNames <- make.names(arrayNames)
  #if(sum(arrayNames != descdata[,1])>0) warning("One or more array names have been adapted to be valid R names")
  #line above: warning not needed as arraynames are in non of the outcome files; if this changes, uncomment the line
  
  experimentFactor<-descdata[,2]
  
  #in R, proper names cannot start with a number, so add a prefix if they do
  experimentFactor <- make.names(experimentFactor)
  #if(sum(experimentFactor != descdata[,2])>0) warning("One or more experimental group names have been adapted to be valid R names")
  #for now the code still does not work with unvalid R names as group names
  if(sum(experimentFactor != descdata[,2])>0) {
    stop("One or more experimental group names are invalid R names\nplease don't use special characters other than . and _ and don't start names with a number")
  }
  
  experimentFactor<-as.factor(experimentFactor)
  print(experimentFactor)
  
  # Annotation  
  
  # if the first column header of the normalized data table is not in the 
  # list of arrayNames, this will be defined as first column of annotation
  # in the result tables (usually it contains the gene/probeset IDs)
  firstColumn = NULL
  if(sum(arrayNames==colnames(normDataTable)[1])==0) {
    firstColumn <- as.matrix(normDataTable[,1])
    colnames(firstColumn) <- colnames(normDataTable)[1]
  }
  
  # other annotation columns are kept only if keepAnnotation is TRUE, in this 
  # case, all other columns of normDataTable that were not recognized as
  # array names are set to annotation. This annotation will be added at the  
  # end of the result tables
  annotation = NULL
  if(keepAnnotation) {
    headers <- colnames(normDataTable)
    notArrayNames <- apply(as.matrix(headers[2:length(headers)]),1,function(x) {
      sum(arrayNames==x)})
    annotation <- normDataTable[,(which(notArrayNames==0)+1)]
  }
  
  # Normalized data
  # recreate the normDataTable based on the desc file to be sure that all  
  # array groups are defined and the columns are correctly ordered:
  normDataTable <- apply(as.matrix(arrayNames),1,function(x) {
    as.numeric(normDataTable[,which(colnames(normDataTable)==x)])})
  normDataTable <- as.data.frame(normDataTable)
  #colnames(normDataTable) <- arrayNames
  
  if(is.null(dim(normDataTable))) {
    stop("could not match array names from description file to normalized data
         file")
  }
  
  # Compute Advanced statistical analysis	
  
  #Covariates
    if(Covariates==TRUE){
      #read string  #for loop
      Cov <- unlist(strsplit(Covariates_String, ";"))
      Design_Cov = ""
      for (i in 1:length(Cov)){
        Cov1 <- unlist(strsplit(Cov[i], ","))
        print(Cov1)
        
        #print(Cov1[1])
        #print(Cov1[2])
        
        if(Cov1[2]=="emp"){
          print(Cov1[1])
        }else{
          if(Cov1[2]=="n"){
            suppressWarnings(as.numeric(Cov1[1]))
          }
          if(Cov1[2]=="f"){
            as.factor(Cov1[1])
          }
          Design_Cov<-paste(Design_Cov, Cov1[1], sep="+")
        }
      }
      
      print(Design_Cov)
  
     #Interaction model
       if(Interaction==TRUE){  
         Cov <- unlist(strsplit(Interaction_String, ","))
         Design_Int = ""
         for (i in 1:length(Cov)){
           Cov1 <- unlist(strsplit(Cov[i], ","))
           print(Cov1)
           
           #(Cov1[1])
           #print(Cov1[2])
           
           Design_Int<-paste(Design_Int, Cov1[1], sep="+")
         }
         
         print(Design_Int)
        #make design covariation + interaction
        design <- model.matrix(~ 0 + Design_Cov + Design_Int )
        } else {
        #make design covariation
        design <- model.matrix(~ 0 + Design_Cov)
        }
      }
    rownames(design) <- arrayNames
    colnames(design) <- gsub("Cov","",colnames(design))  
    print(design)
  
  #Paired vs unpaired data 
  if(Paired_data==TRUE){
    #Paired data
    fit <- lmFit(normDataTable[,2:(size-1)], design, ndups=1, block=Paired_String, correlation=corfit$consensus)
    fit <- eBayes(fit)
  } else {
    #unpaired data 
    fit <- lmFit(normDataTable[,2:(size-1)],design)
    fit <- eBayes(fit)
  }
  #Pass the Fit
  createVolcanoPlot(fit)
  createVennPlot(fit)
  
    # Use contrast matrix to generate group comparisons  
    filesNew<-NULL
    if(defaultContr) { # always FALSE when coming from the web form.  
      if(length(levels(experimentFactor))<=4){
        #make contrast
            cont.matrix <- defaultMatrix(design)
        
        defFiles<-saveComparison(cont.matrix,fit,normDataTable,
                                 annotation,firstColumn)
        filesNew<-c(filesNew,defFiles)
      }
      else
        print ("[[----Cant compute default matrices for more than four 
  		groups----]]")
    }
    if(!is.null(matfileName)) {
      #be careful, since the saved contrast matrix still uses the not corrected names, also pass the not corrected names
      cont.matrix <- enterMatrix(matfileName,descdata[,2])
      advFiles<-saveComparison(cont.matrix,fit,normDataTable,
                               annotation,firstColumn)
      filesNew<-c(filesNew,advFiles)
    }
    
    return(filesNew)

    #Create Volcano Plot
    if(CreateVolcanoPlot){
    createVolcanoPlot <- function(fit) {
      #extract data from the comparison tables
     
        cat (paste("--[[ Saving volcano_plot.png ]]--\n"))
        png("C:/Users/Celine Adriaans/Documents/University/Major Internship/Extending Stat Module/Real documents/Volcano.png",width=1000,
            height=1000, bg="white")
        
        volcanoplot(fit1,coef=2,highlight=0, las=1)
        dev.off()
        }
      }
  
  
    #Create Venn Diagram
    if(CreateVennPlot){
    createVennPlot <- function(fit) {
      #extract data from the comparison tables
      
        cat (paste("--[[ Saving","venn_plot.png ]]--\n"))
        png("C:/Users/Celine Adriaans/Documents/University/Major Internship/Extending Stat Module/Real documents/Venn.png",width=1000,
            height=1000, bg="white")
        
        VennCounts<-vennCounts(fit, include="both")
        vennDiagram(VennCounts, include="both", names=NULL, mar=rep(1,4), cex=c(1.5,1,0.7), lwd=1,
                    circle.col=NULL, counts.col=NULL, show.include=NULL)
        print("VennDiagram")
        dev.off()
        }
      }
    }


computeAdvancedStatistics (normDataTable, descriptionFile, Covariates=TRUE, Covariates_String, Interaction=TRUE, Interaction_String, 
                                     Paired_data=TRUE, Paired_String, defaultContr=TRUE, CreateVolcanoPlot=TRUE, CreateVennPlot=TRUE, 
                                     matfileName=NULL, keepAnnotation=FALSE)

