####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.1.0
####
####  Copyright 2012, University of Miami
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Project Partially Funded By: 
####  ----------------------------------------------------------------
####  Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
####  National Science Foundation and grant R01 CA163739 from the National
####  Cancer Institute.
####
####  Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
####  National Cancer Institute.
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


rf2rfz <- function(object,
                   forestName = NULL,
                   ...)
{

  ## Ensure the forest object is coherent.
  rfsrcForest <- checkForestObject(object)

  if (is.null(forestName)) {
    stop("RFSRC forest name is NULL.  Please provide a valid name for the forest .rfz file.")
  }

  ## If the user has already provided the .rfz extension, remove it.
  if (nchar(forestName) > 4) {
    if (substr(forestName, nchar(forestName)-3, nchar(forestName)) == ".rfz") {
      forestName <- substr(forestName, 1, nchar(forestName)-4)
    }
  }
  
  ## Initialize the local variables extracted from the forest object.
  nativeArray <- rfsrcForest$nativeArray
  time.interest <- rfsrcForest$time.interest
  formula <- rfsrcForest$formula
  forestSeed <- rfsrcForest$seed
  xvar.names <- rfsrcForest$xvar.names

  ## Extract the xvar types.

  get.factor <- extract.factor(rfsrcForest$xvar, xvar.names)
  xvar.type <- get.factor$generic.types

  ## This may be null in the absence of factors.
  nativeFactorArray <- rfsrcForest$nativeFactorArray
        
  ## Count the number of trees in the forest.
  numTrees <- length(as.vector(unique(nativeArray$treeID)))

  ## Define the root elements of the PMML file.  This is a quick work-around for an issue
  ## with this version of the XML package, and the inablility to add namespace information
  ## and attributes concurrently. 
  rootString <- getRootString()

  ## Define the document and the root node.
  pmmlDoc <- xmlTreeParse(rootString, asText=TRUE)
  pmmlRoot <- xmlRoot(pmmlDoc)

  ## Add the DataDictionary to the root node.
  pmmlRoot <- append.XMLNode(pmmlRoot, getDataDictNode(xvar.names=xvar.names, xvar.type=xvar.type))

  ## Write the native array information.
  write.table(nativeArray, 
              paste(forestName, ".txt", sep=""), quote = FALSE)

  ## Write the native factor array information if it exists.
  ## *** WARNING ***  *** WARNING ***  *** WARNING *** 
  ## In 32-bit and 64-bit systems, the integer value 0x8000000 is
  ## interpreted as NA by R as it output from the native code SEXP
  ## object.  Thus write.table will contain NA's.  These need to be
  ## handled specially in the JUNG code that parses the .rfz file.
  ## *** WARNING ***  *** WARNING ***  *** WARNING ***   
  write.table(nativeFactorArray, 
                paste(forestName, ".factor.txt", sep=""), col.names=FALSE, quote = FALSE)


  ## Write the xvar names and types.
  xmlFile <- file(paste(forestName, ".xml", sep=""), open="w")
  saveXML(pmmlRoot, xmlFile)
  close(xmlFile)

  zipCommand <- paste("zip", sep=" ",
    paste(forestName, ".rfz", sep=""),
    paste(forestName, ".txt", sep=""),
    paste(forestName, ".factor.txt", sep=""),
    paste(forestName, ".xml", sep="")) 

  system(command = zipCommand)

  unlink(paste(forestName, ".txt", sep=""))
  unlink(paste(forestName, ".factor.txt", sep=""))
  unlink(paste(forestName, ".xml", sep=""))

}

## Coherency check of the forest object.  Failure stops execution, hence there is no return value.
checkForestObject <- function(object) {

  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }

  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    }
    rfForest <- object$forest
  }
  else {
    rfForest <- object
  }

  if (is.null(rfForest$nativeArray)) {
    stop("RFsrc nativeArray content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rfForest$xvar.names)) {
    stop("RFsrc xvar.names content is NULL.  Please ensure the object is valid.")
  }

  if (is.null(rfForest$xvar)) {
    stop("RFsrc xvar content is NULL.  Please ensure the object is valid.")
  }

  return (rfForest)
}


##  The root string is accessed here by rf2rfz().
getRootString <- function() {
  rootString <- 
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
         <PMML version=\"3.1\" xmlns=\"http://www.dmg.org/PMML-3_1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">
           <Header copyright=\"Copyright 2008, Cleveland Clinic\" description=\"Random Survival Forest Tree Model\">
              <Application name=\"Random Survival Forest\" version=\"3.0\"/>
           </Header>
         </PMML>
       "
  return (rootString)
}

## Form the data dictionary node.  This is also accessed by rf2rfz().
getDataDictNode <-  function(xvar.names, xvar.type) {

  ## Define the DataDictionary node for the document.
  dataDictNode <- xmlNode("DataDictionary", attrs=c(numberOfFields=length(xvar.names)))

  ## Add the xvar names to the DataDictionary.
  for (k in 1:length(xvar.names)) {
    if (xvar.type[k] == "C") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=xvar.names[k], optype="categorical", dataType="string")))
    }
    if (xvar.type[k] == "I") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=xvar.names[k], optype="ordinal", dataType="integer")))
    }
    if (xvar.type[k] == "R") {
      dataDictNode <- append.XMLNode(dataDictNode, xmlNode("DataField", attrs=c(name=xvar.names[k], optype="continuous", dataType="double")))
    }
  }

  return (dataDictNode)
  
}
