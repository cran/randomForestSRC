####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.6.1
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
####    email:  commerce@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


impute.rfsrc <- function(formula,
                         data,
                         ntree = 250,
                         mtry = NULL,
                         nodesize = 1,
                         splitrule = NULL,
                         nsplit = 1,
                         na.action = c("na.impute", "na.random"),
                         nimpute = 1, 
                         xvar.wt = NULL, 
                         do.trace = FALSE,
                         ...)
{
  if (missing(data)) {
    stop("data is missing")
  }
  which.na <- is.na(data)
  if (!any(which.na) || all(which.na)) {
    return(invisible(data))
  }
  p <- ncol(data)
  n <- nrow(data)
  all.r.na <- rowSums(which.na) == p
  all.c.na <- colSums(which.na) == n
  data <- data[!all.r.na, !all.c.na, drop = FALSE]
  which.na <- which.na[!all.r.na, !all.c.na, drop = FALSE]
  if (!any(which.na)) {
    return(data)
  }
  p <- ncol(data)
  n <- nrow(data)
  all.var.names <- colnames(data)
  mforest <- FALSE
  blocks <- list(1:nrow(data))
  if (!mforest) {
    if (missing(formula)) {
      ytry <- min(p - 1, max(25, ceiling(sqrt(p))))
      formula <- as.formula(paste("Unsupervised(", ytry, ") ~ ."))
    }
    nullBlocks <- lapply(blocks, function(blk) {
      dta <- data[blk,, drop = FALSE]
      retO <- tryCatch({generic.impute.rfsrc(formula = formula,
                                   data = dta,
                                   ntree = ntree,
                                   nimpute = nimpute,
                                   mtry = mtry,
                                   nodesize = nodesize,
                                   splitrule = splitrule,
                                   nsplit = nsplit,
                                   na.action = na.action,
                                   xvar.wt = xvar.wt,
                                   do.trace = do.trace)}, error = function(e) {NULL})
      if (!is.null(retO)) {
        if (!is.null(retO$missing$row)) {
          blk <- blk[-retO$missing$row]
        }
        if (!is.null(retO$missing$col)) {
          ynames <- all.var.names[-retO$missing$col]
        }
        else {
          ynames <- all.var.names
        }
        data[blk, ynames] <<- retO$data[, ynames, drop = FALSE]
      }
      NULL
    })
    rm(nullBlocks)
  }
  invisible(data)
}
