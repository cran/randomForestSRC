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


########################################################################
# 
# impute data using rfsrc
#
#
########################################################################

impute.rfsrc <- function(formula,
                         data,
                         ntree = 1000,
                         mtry = NULL,
                         nodesize = NULL,
                         splitrule = NULL,
                         nsplit = 0,
                         nimpute = 1,
                         xvar.wt = NULL,
                         seed = NULL,
                         do.trace = FALSE,
                         ...)
{

  ## set parameters accordingly
  importance <- "none"
  na.action <- "na.impute"
  forest <- FALSE
  proximity <- FALSE
  var.used <- FALSE
  split.depth <- FALSE
  membership <- FALSE

  ## terminate if there is no data
  if (missing(data)) {
    stop("data is missing")
  }

  ## save the rownames
  row.names.data <- rownames(data)

  ## handle the case where there is no formula
  ## basically this amounts to random splitting
  ## we mimic this by adding a fake y-outcome

  if (missing(formula) | (!missing(formula) && is.null(formula))) {
    formula <- as.formula("junkFAKEYjunkjunkFAKEYjunk ~ .")
    data$junkFAKEYjunkjunkFAKEYjunk <- rnorm(nrow(data))
    missing.formula <- TRUE
    splitrule <- "random"
  }
  else {
    missing.formula <- FALSE
  }
  
  ## rfsrc grow call
  object <- rfsrc(formula = formula,
                  data = data,
                  ntree = ntree,
                  mtry = mtry,
                  nodesize = nodesize,
                  splitrule = splitrule,
                  nsplit = nsplit,
                  nimpute = nimpute,
                  xvar.wt = xvar.wt,
                  seed = seed,
                  do.trace = do.trace,
                  importance = importance,
                  na.action = na.action,
                  forest = forest,
                  proximity = proximity,
                  var.used = var.used,
                  split.depth = split.depth,
                  membership = membership,
                  impute.only = TRUE)

  
  ## preliminary results of imputation
  imputed.result <- cbind(object$yvar, object$xvar)
  colnames(imputed.result) <- c(object$yvar.names, object$xvar.names)

  ##Overlay the data (only necessary when nimpute = 1)
  if (nimpute == 1) {
    imputed.result[object$imputed.indv, ] <- object$imputed.data
  }

  ## strip away the fake y-outcome in the missing formula scenario
  if (missing.formula) {
    imputed.result$junkFAKEYjunkjunkFAKEYjunk <- NULL
  }

  ## restore the row names
  rownames(imputed.result) <- row.names.data
  
  
  ## return the goodies 
  invisible(imputed.result) 

}


