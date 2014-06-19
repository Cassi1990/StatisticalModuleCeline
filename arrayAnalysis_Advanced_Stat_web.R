#=============================================================================#
# ArrayAnalysis - affyAnalysisStat                                            #
# a tool for statistical analysis of Affymetrix expression data               #
#                                                                             #
# Copyright 2010-2011 BiGCaT Bioinformatics                                   #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
#=============================================================================#

arrayAnalysisStat<-function(...) {
  
  ###############################################################################
  # Get input parameters                            			              #
  ###############################################################################
  
  args <- list(...)
  for(i in 1:length(args)) {
    flag <- substring(args[[i]], 0, 2)
    value <- substring(args[[i]], 3, nchar(args[[i]]))
    if(flag=='-r'){
      refName = value
      print(refName)
    }			
    if(flag=='-n'){
      normDataTable = value
      print(normDataTable)
    }	
    if(flag=='-d'){
      descriptionFile = value
      print(descriptionFile)
    }				
    if(flag=='-k'){
      keepAnnotation = value
    }				
    if(flag=='-m'){
      matfileName = value
    }
    if(flag=='-C'){
      cutOffTable = value
    }
    if(flag=='-p'){
      cutOffPval = value
    }
    if(flag=='-f'){
      cutOfflogFC = value
    }
    if(flag=='-e'){
      cutOffAveExpr = value
    }
    if(flag=='-H'){
      plotPvalHist = value
    }
    if(flag=='-h'){
      plotFCHist = value
    }
    if(flag=='-S'){
      summaryTable = value
    }		  
    if(flag=='-P'){
      pvaluelist = value
    }
    if(flag=='-A'){
      adjpvaluelist = value
    }
    if(flag=='-F'){
      foldchangelist = value
    }
    ################################Added parameters
    if(flag=='-v'){
     plotVennDiagram= value
    }
    if(flag=='-V'){
      plotVolcanoPlot= value
    }
    if(flag=='-a'){
      Paired_data= value
    }
    if(flag=='-x'){
      Covariates= value
    }
    if(flag=='-X'){
      Interaction= value
    }
  }
  
  ###############################################################################
  # Set directories 		                          				              
  ###############################################################################
  
  SCRIPT.DIR <- getwd()
  WORK.DIR <- paste("../temp/",refName,"_Stat",sep="") # directory where the results are computed
  
  correctDIR <- function(d) { 
    lastChar <- substr(d,nchar(d),nchar(d))
    if((lastChar != "/") && (lastChar != "/")) d <- paste(d,"/",sep="")
    return(d)
  }
  SCRIPT.DIR <- correctDIR(SCRIPT.DIR)
  WORK.DIR <- correctDIR(WORK.DIR)
  
  if(file.exists(WORK.DIR)) {
    setwd(WORK.DIR)
    listfile <- list.files()
    unlink(listfile)
    unlink(WORK.DIR) 
  }
  dir.create(WORK.DIR) 
  
  ###############################################################################
  # Move input files to WORK.DIR	                          				              
  ###############################################################################
  
  setwd(SCRIPT.DIR)
  listfile<-paste("../temp/",c(normDataTable,descriptionFile,matfileName),sep="")
  file.copy(listfile, WORK.DIR)
  unlink(listfile)      
  
  ###############################################################################
  # Run functions 		                          				              
  ###############################################################################
  
  source(paste(SCRIPT.DIR,"setParameters_Advanced_Stat_web.R",sep=""),local=TRUE)
  source(paste(SCRIPT.DIR,"run_arrayAnalysis_Advanced_Stat.R",sep=""),local=TRUE)
  
  # PARAMETERS DESCRIPTION : see arrayAnalysisStat.R
  
}