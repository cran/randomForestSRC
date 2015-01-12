####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.6.0
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


rfsrcSyn.rfsrc <-
  function(formula,          
           data,
           object,
           newdata,
           ntree = 1000,
           mtry = NULL,
           mtrySeq = NULL,
           nodesize = 5,
           nodesizeSeq = c(1:10,20,30,50,100),
           nsplit = 0,
           min.node = 3,
           use.org.features = TRUE,
           na.action = c("na.omit", "na.impute"),
           verbose = TRUE,
           ...
           )
{
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  if (!missing(object)) {
    if (sum(inherits(object, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) != 2) {
      stop("this function only works for objects of class `(rfsrc, synthetic)'")
    }
    M <- length(object$rfMachines)
    fmly <- object$rfMachines[[1]]$family
    xvar.names <- object$rfMachines[[1]]$xvar.names
    yvar.names <- object$rfMachines[[1]]$yvar.names
    J <- nlevels(object$rfMachines[[1]]$yvar)
    rfMachines <- object$rfMachines
    rfSyn <- object$rfSyn
    synthetic <- object$synthetic
    opt.machine <- object$opt.machine
    list.names <- unlist(lapply(synthetic, function(ss) {colnames(ss)}))
  }
  else {
    if (missing(formula) || missing(data)) {
      stop("need to specify 'formula' and 'data' or provide a grow forest object")
    }
    f <- as.formula(formula)
    if (na.action == "na.impute" && any(is.na(data))) {
      if (verbose) {
        cat("\t imputing the data\n")
      }
      data <- impute.rfsrc(data = data, ntree = ntree, nodesize = nodesize, nsplit = nsplit)
    }
    preObj <- rfsrc(f, data, ntree = 1, importance = "none",
                    nodesize = nrow(data), splitrule = "random")
    fmly <- preObj$family
    if (!(fmly == "regr" | fmly == "class")) {
      stop("this function only applies to regression and classification families")
    }
    xvar.names <- preObj$xvar.names
    yvar.names <- preObj$yvar.names
    p <- length(xvar.names)
    if (is.null(mtrySeq)) {
      mtrySeq <- ceiling(p/3)
    }
    else {
      mtrySeq <- unique(ceiling(mtrySeq))
      mtrySeq <- mtrySeq[mtrySeq>=1 & mtrySeq <= p]
      if (length(mtrySeq) == 0) {
        stop("invalid choice for mtrySeq:", mtrySeq)
      }
    }
    nodesizeSeq <- sort(nodesizeSeq)
  }
  na.action <- match.arg(na.action, c("na.omit", "na.impute"))
  if (missing(object)) {
    N <- length(nodesizeSeq)
    J <- nlevels(preObj$yvar)
    rfMachines <- lapply(nodesizeSeq, function(nn) {
      lapply(mtrySeq, function(mm) {
      if (verbose) {
        cat("\t RF nodesize:", nn, "mtry:", mm, "\r")
      }
      rfsrc(f, data, ntree = ntree, mtry = mm, nodesize = nn, 
            nsplit = nsplit, importance = "none")
    })})
    rfMachines <- unlist(rfMachines, recursive = FALSE)
    list.names <- paste(rep(nodesizeSeq, each = length(mtrySeq)), mtrySeq, sep = ".")
    M <- length(rfMachines)                         
    if (is.numeric(min.node) && min.node > 0) {
      good.machines <- which(sapply(1:M, function(m) {
        mean(rfMachines[[m]]$leaf.count, na.rm = TRUE) > min.node}))
      if (length(good.machines) == 0) {
        good.machines <- 1
      }
      list.names <- list.names[good.machines]
      rfMachines <- lapply(good.machines, function(m) {rfMachines[[m]]})
      M <- length(rfMachines)
    }
    opt.machine <- rf.opt(rfMachines)
    if (verbose) {
      cat("\t making the synthetic features\n")
    }
    if (fmly == "class") {
      synthetic <- lapply(1:M, function(m) {
        prb <- rfMachines[[m]]$predicted.oob[ ,1:(J-1), drop = FALSE]
        if (J > 2) {
          colnames(prb) <- paste(1:(J-1), list.names[m], sep = ".")
        }
        else {
          colnames(prb) <- list.names[m]
        }
        prb
      })
    }
    else {
      synthetic <- lapply(1:M, function(m) {
        yhat <- cbind(rfMachines[[m]]$predicted.oob)
        colnames(yhat) <- list.names[m]
        yhat
      })
    }
    x.s <- do.call("cbind", synthetic)
    if (verbose) {
      cat("\t making the synthetic forest\n")
    }
    if (use.org.features) {
      data <- data.frame(y = preObj$yvar, preObj$xvar, x.s = x.s)
    }
    else {
      data <- data.frame(y = preObj$yvar, x.s = x.s)
    }
    rfSyn <- rfsrc(y ~ ., data, ntree = ntree, mtry = mtry, nodesize = nodesize,
                     nsplit = nsplit, ... )
  }
  if (!missing(newdata)) {
    if (na.action == "na.impute" && any(is.na(newdata))) {
      if (verbose) {
        cat("\t imputing the test data\n")
      }
      newdata <- impute.rfsrc(data = newdata, ntree = ntree, nodesize = nodesize, nsplit = nsplit)
    }
    if (verbose) {
      cat("\t making the test set synthetic features\n")
    }
    xtest <- newdata[, xvar.names, drop = FALSE]
    if (fmly == "class") {
      synthetic <- lapply(1:M, function(m) {
        prb <- predict(rfMachines[[m]], xtest, importance = "none")$predicted[, 1:(J-1), drop = FALSE]
        if (J > 2) {
          colnames(prb) <- paste(1:(J-1), list.names[m], sep = ".")
        }
        else {
          colnames(prb) <- list.names[m]
        }
        prb
      })
    }
    else {
      synthetic <- lapply(1:M, function(m) {
        yhat <- cbind(predict(rfMachines[[m]], xtest, importance = "none")$predicted)
        colnames(yhat) <- list.names[m]
        yhat
      })
    }
    xtest.s <- do.call("cbind", synthetic)
    data.test <- data.frame(x.s = xtest.s)
    if (use.org.features) {
      data.test <- data.frame(data.test, xtest)
    }
    yhat <- predict(rfSyn, data.test, importance = "none")$predicted
    if (!is.null(data.test$yvar.names)) {
      if (fmly == "class") {
        err.rate <- brier(data.test$yvar.names, yhat)
      }
      else {
        err.rate <- mean((data.test$yvar.names - yhat)^2, na.rm = TRUE)
      }
    }
    else {
       err.rate <- NULL
    }
  }
  else {
    yhat <- rfSyn$predicted
    if (fmly == "class") {
      err.rate <- brier(rfSyn$yvar, rfSyn$predicted.oob)
    }
    else {
      err.rate <- mean((rfSyn$yvar - rfSyn$predicted.oob)^2, na.rm = TRUE)
    }
  }
  retObj <- list(rfMachines = rfMachines,
                 rfSyn = rfSyn,
                 synthetic = synthetic,
                 predicted = yhat,
                 err.rate = err.rate,
                 opt.machine = opt.machine)
  class(retObj) <- c("rfsrc", "synthetic")
  retObj
}
rf.opt <- function(obj)
{
  if (obj[[1]]$family == "regr") {
    ntree <- length(obj[[1]]$err.rate)
    which.min(unlist(lapply(1:length(obj), function(ll) {obj[[ll]]$err.rate[ntree]})))[1]
  }
  else {
    which.min(unlist(lapply(1:length(obj), function(ll) {
    brier(obj[[ll]]$yvar, obj[[ll]]$predicted.oob)})))[1]
  }
}
rfsrcSyn <- rfsrcSyn.rfsrc
