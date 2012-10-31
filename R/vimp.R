####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.0.0
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
####    email:  kogalurshear@gmail.com
####    URL:    http://www.kogalur-shear.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


vimp.rfsrc <- function(object,
                 xvar.names,
                 importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
                 joint = FALSE,
                 newdata,
                 subset,
                 na.action = c("na.omit", "na.impute"),
                 seed = NULL,
                 do.trace = FALSE,
                 num.threads = -1,
                 papply = if (require("multicore")) mclapply else lapply,
                 ...)

{

  ## Incoming parameter checks.  All are fatal.
  if (missing(object)) stop("Object is missing")

  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'.")

  ## Process the importance specification
  if (!is.logical(joint)) {
    stop("joint must be a logical value")
  }
  importance <- importance[1]
  if (joint & importance != "none") {
    i.str <- unlist(strsplit(importance, "\\."))
    if (length(i.str) == 1) {
      importance <- paste(i.str[1], ".joint", sep = "")
    }
    else if (length(i.str) == 2) {
      importance <- paste(i.str[1], ".joint.", i.str[2], sep = "")
    }
  }
  importance <- match.arg(importance,
      c("permute", "random", "permute.ensemble", "random.ensemble", "none", 
        "permute.joint", "random.joint", "permute.joint.ensemble", "random.joint.ensemble"))

  
  ### absence of newdata indicates grow scenario
  ### this has implications for whether peformance measures must be present in the object
  ### and whether outcome = train or outcome = test for the generic.predict call

  if (missing(newdata)) {
    
    ## grow objects under non-standard bootstrapping are devoid of
    ## performance values
    if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) {
        stop("The forest is empty.  Re-run rfsrc (grow) call with forest=TRUE")
      }
      else {
        bootstrap <- object$forest$bootstrap
      }
    }
    else {
      bootstrap <- object$bootstrap
    }
    if (bootstrap != "by.root") {
      stop("grow objects under non-standard bootstrapping are devoid of performance values")
    }
    
    ## combine the xvar and yvar data into a "test" data
    ## call generic predict with outcome = "test"
    object$yvar <- as.data.frame(object$yvar)
    colnames(object$yvar) <- object$yvar.names
    newdata <- cbind(object$yvar, object$xvar)

    ## set the outcome specification 
    outcome.type <- "test"
    
  }
  
  else {

    ## nothing to do except confirm that newdata is a data frame
    if (!is.data.frame(newdata)) {
      stop("newdata must be a data frame")
    }

    ## set the outcome specification 
    outcome.type <- "train"
    
  }
  
  ## Initialize sample size.
  n <- nrow(newdata)

  ## Process the subsetted index 
  ## Assumes the entire data set is to be used if not specified
  if (missing(subset)) {
    subset <- 1:n
  }
  else {
    ## convert the user specified subset into a usable form
    if (is.logical(subset)) subset <- which(subset)
    subset <- unique(subset[subset >= 1 & subset <= n])

    if (length(subset) == 0) {
      stop("'subset' not set properly.")
    }
  }

  ##subset the data
  newdata <- newdata[subset,, drop = FALSE]

  result <- generic.predict.rfsrc(object,
                                  newdata = newdata,
                                  na.action = na.action,
                                  importance = importance,
                                  importance.xvar = xvar.names,
                                  outcome = outcome.type,
                                  proximity = FALSE,
                                  var.used = FALSE,
                                  split.depth = FALSE,
                                  perf = TRUE,
                                  seed = seed,
                                  do.trace = do.trace,
                                  num.threads = num.threads,
                                  membership = FALSE,
                                  papply = papply,
                                  restore.only = FALSE,
                                  ...)


 return(result)

}

vimp <- vimp.rfsrc
