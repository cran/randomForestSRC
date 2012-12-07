####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.0.2
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
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


adrop <- function(x, d) {

  #this function is for arrays only
  if (!is.array(x)) {
    x
  }

  else {
    
    if (d > 1) {
      x[,,1:d, drop = FALSE]
    }
    else {
      if (dim(x)[1] == 1) {
        rbind(x[,,1, drop = TRUE])
      }
      else {
        if (dim(x)[2] == 1) {
          cbind(x[,,1, drop = TRUE])
        }
        else {
          x[,,1, drop = TRUE]
        }
      }
    }

  }
    
}

amatrix <- function(x, d, names) {

  x <- matrix(x, d, dimnames = names)
  
  if (ncol(x) > 1) {
    x
   }
  else {
    c(x)
  }
  
}


amatrix.remove.names <- function(x) {

  if (!is.null(dim(x)) && ncol(x) == 1) {
    unlist(c(x), use.names = FALSE)
  }
  else {
    x
  }
  
}


atmatrix <- function(x, d, names, keep.names = FALSE) {

  x <- t(matrix(x, ncol = d, dimnames = names))
  
  if (ncol(x) > 1) {
    x
  }
  else {
    if (keep.names == FALSE) {
      c(x)
    }
    else {
      x.names <- rownames(x) 
      x <- c(x)
      names(x) <- x.names
      x
    }
  }
  
}

avector <- function(x, name = FALSE) {
  if (!is.null(dim(x)) && nrow(x) > 1 && ncol(x) == 1) {
    x.names <- rownames(x)
    x <- unlist(c(x))
    if (name) names(x) <- x.names else names(x) <- NULL
    x
  }
  else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) > 1) {
    x.names <- colnames(x)
    x <- unlist(c(x))
    if (name) names(x) <- x.names else names(x) <- NULL
    x
  }
  else if (!is.null(dim(x)) && nrow(x) == 1 && ncol(x) == 1) {
    unlist(c(x))
  }
  else {
    x
  }
}

available <- function (package, lib.loc = NULL, quietly = TRUE) 
{
  package <- as.character(substitute(package))
  installed <- package %in% installed.packages()
  if (installed) {
    require(package, quietly = TRUE, character.only = TRUE)
  }
  else {
    return(invisible(FALSE))
  }
}

bayes.rule <- function(prob) {
  levels.class <- colnames(prob)
  factor(levels.class[apply(prob, 1, function(x) {
    if (!all(is.na(x))) {
      resample(which(x == max(x, na.rm = TRUE)), 1)
    }
    else {
      NA
    }
  })], levels = levels.class)
}
                          
data.matrix <- function(x) {
  as.data.frame(lapply(x, function(xi) {
     if (is.integer(xi) || is.numeric(xi)) {
       xi
     }
     else if (is.logical(xi) || is.factor(xi)) {
       as.integer(xi)
     }
     else {
       as.numeric(xi)
     }
  }))
}

extract.pred <- function(obj, type, subset, percentile, which.outcome, oob = FALSE) {
  # first decide if OOB or in-bag values are requested
  if (oob == FALSE) {
    pred <- obj$predicted
    surv <- obj$survival
    chf <- obj$chf
    cif <- obj$cif
  }
  else {
    pred <- obj$predicted.oob
    surv <- obj$survival.oob
    chf <- obj$chf.oob
    cif <- obj$cif.oob
  }
  #now extract the information dependent on the family
  if (obj$family == "surv") {
    n <- length(pred)
    if (missing(subset)) subset <- 1:n
    surv.type <- match.arg(type, c("mort", "rel.freq", "surv"))
    time.idx <-  max(which(obj$time.interest <=
                quantile(obj$yvar[, 1], probs = percentile, na.rm = TRUE)))
    return(switch(surv.type,
             "mort" = pred[subset],
             "rel.freq" = pred[subset]/max(n, na.omit(pred)),
             "surv" =  100 * surv[subset, time.idx]
    ))
  }
  #now extract the information dependent on the family
  else if (obj$family == "surv-CR") {
    n <- length(pred)
    if (missing(subset)) subset <- 1:n
    if (missing(which.outcome)) which.outcome <- 1#default is first event type
    cr.type <- match.arg(type, c("years.lost", "cif", "chf"))
    time.idx <-  max(which(obj$time.interest <=
                quantile(obj$yvar[, 1], probs = percentile, na.rm = TRUE)))
    return(switch(cr.type,
             "years.lost" = pred[subset, which.outcome],
             "cif" = cif[subset, time.idx, which.outcome],
             "chf" = chf[subset, time.idx, which.outcome]
    ))
  }
  else if (obj$family == "class") {
    class.type <- match.arg(type, c("response", "prob"))
    if (missing(subset)) subset <- 1:nrow(pred)
    if (missing(which.outcome)) which.outcome <- 1#default is first class label
    prob <- pred[subset,, drop = FALSE]
    return(switch(class.type,
                  "prob" = prob[, which.outcome],
                  "response" =  bayes.rule(prob)))
  }
  else {
    if (missing(subset)) subset <- 1:length(pred)
    return(pred[subset])
  }
}

family.pretty <- function(fmly) {
  switch(fmly,
         "surv"     = "Random Forests [S]RC",
         "surv-CR"  = "Random Forests [S]RC",
         "regr"     = "Random Forests S[R]C",
         "class"    = "Random Forests SR[C]")
}

finalizeFormula <- function(formula.obj, data) {

  ## parse the formula object
  yvar.names <- formula.obj$yvar.names
  index <- length(yvar.names)
  fNames <- formula.obj$fNames
  fmly <- formula.obj$family

  ## extract the xvar names
  if (fNames[index + 1] == ".") {
    xvar.names <- names(data)[!is.element(names(data), fNames[1:index])]
  }
  else {
    xvar.names <- fNames[-c(1:index)]
    not.specified <- !is.element(xvar.names, names(data))
    if (sum(not.specified) > 0) {
      stop("formula appears misspecified, object ", xvar.names[not.specified], " not found")
    }
  }

  # return the goodies
  return (list(family=fmly, yvar.names=yvar.names, xvar.names=xvar.names))
  
}

finalizeData <- function(fnames, data, na.action) {

  ## Data conversion to numeric mode
  data <- data[ , is.element(names(data), fnames), drop = FALSE]
  data <- as.data.frame(data.matrix(data))

  if (na.action == "na.omit") {
    data <- na.omit(data)
  }

  ## is anything left?
  if (nrow(data) == 0) {
    stop("No records in the NA-processed data.  Consider na.impute.")
  }

  return (data)
  
}

get.importance.xvar <- function(importance.xvar, importance, object) {

  ## Check that importance has been requested
  if (!is.null(importance)) {
    ## Map vimp names to columns of GROW x-matrix
    ## Ensure names are coherent
    if (missing(importance.xvar) || is.null(importance.xvar)) {
      importance.xvar <- object$xvar.names
    }
    else {
      importance.xvar <- unique(importance.xvar)
      importance.xvar <- intersect(importance.xvar, object$xvar.names)
    }

    if (length(importance.xvar) == 0) {
      stop("xvar names do not match object xvar matrix")
    }
  }
  else {
    importance.xvar <- 0
  }
  return (importance.xvar)
}

get.nmiss <- function(xvar, yvar = NULL) {
  if (!is.null(yvar)) {
    sum(apply(yvar, 1, function(x){any(is.na(x))}) | apply(xvar, 1, function(x){any(is.na(x))}))
  }
  else {
    sum(apply(xvar, 1, function(x){any(is.na(x))}))
  }
}

    
get.grow.nodesize <- function(fmly, nodesize) {

  ## Default node size for right-censored survival
  if (fmly == "surv"){
    if (is.null(nodesize)) {
      nodesize <- 3
    }
    
  }

  ## Default node size for competing risks
  else if (fmly == "surv-CR"){
   if (is.null(nodesize)) {
     nodesize <- 6
   }
 }

  ## Default node size for classification
  else if (fmly == "class") {
    if (is.null(nodesize)) {
      nodesize <- 1
    }
  }

  ## Default node size for regression
  else {
    if (is.null(nodesize)) {
      nodesize <- 5
    }
  }

  ##nodesize should be rounded if non-integer
  nodesize <- round(nodesize)
  
}  

get.sexp.dim <- function(fmly, event.type, yfactor, splitrule = NULL) {
  
  if (grepl("surv", fmly)) {
    
    if (!is.null(splitrule)) {
      if ((length(event.type) > 1) & (splitrule != "logrank") & (splitrule != "logrankscore")) {
        sexp.dim <- length(event.type)
      }
      else {
        sexp.dim <- 1
      }
    }
    else {
      if (length(event.type) > 1) {
        sexp.dim <- length(event.type)
      }
      else {
        sexp.dim <- 1
      }
    }
  }
  
  else {

    sexp.dim <- 0
    
    if (fmly == "class") {
      ## TBD TBD TBD
      ## This hack assumes only one response to determine the numbers of
      ## levels in the factor.  This needs to be generalized to handle
      ## multiple mixed responses. We can use extract.factor() for this
      ## purpose and create a generalized response matrix with
      ## associated levels and an immutable mapping.

      if (!is.null(yfactor$levels)) {
        sexp.dim <- length(yfactor$levels[[1]]) + 1
      }
      else {
        if (!is.null(yfactor$order.levels)) {
          sexp.dim <- length(yfactor$order.levels[[1]]) + 1
        }
      }
      if (sexp.dim == 0) {
        stop("The classification outcome is not a factor")
      }

    }
    
    else {
      ## Regression family.
      sexp.dim <- 1
    }
    
  }

  sexp.dim
}

get.event.info <- function(obj, subset = NULL) {


  ##survival case
  if (grepl("surv", obj$family)) {

    if (!is.null(obj$yvar)) {

      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
  
      r.dim <- 2
      time <- obj$yvar[subset, 1]
      cens <- obj$yvar[subset, 2]
      
      ## Extract the unique event types.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- unique(event)
      
    }

    ##everything else
    else {
      
      r.dim <- 0
      event <- event.type <- cens <- cens <- time <- NULL
      
    }
    
    ## Set grid of time points.
    time.interest <- obj$time.interest
    
 
  }
  
  else {#NULL for other families
    
    ## TBD TBD TBD
    ## Currently all REGR and CLAS analysis is restricted to a single outcome.
    r.dim <- 1
    event <- event.type <- cens <- time.interest <- cens <- time <- NULL
    
  }
  
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
  
}

get.grow.event.info <- function(yvar, fmly, need.deaths = TRUE, ntime) {

  
  if (grepl("surv", fmly)) {

    r.dim <- 2
    time <- yvar[, 1]
    cens <- yvar[, 2]
    if (need.deaths & all(na.omit(cens) == 0)) {
      stop("no deaths in data!")
    }
    
    ## Check for event time consistency.
    if (!all(na.omit(time) >= 0)) {
      stop("time must be  positive")
    }
    
    ## Extract the unique event types.    
    event.type <- unique(na.omit(cens))
    ## Ensure they are all greater than or equal to zero.
    if (sum(event.type >= 0) != length(event.type)) {
      stop("censoring variable must be coded as NA, 0, or greater than 0.")    
    }
    
    ## Discard the censored state, if it exists.
    event <- na.omit(cens)[na.omit(cens) > 0]
    event.type <- unique(event)
    
    ## Set grid of time points.
    nonMissingOutcome <- which(!is.na(cens) & !is.na(time))
    nonMissingDeathFlag <- (cens[nonMissingOutcome] != 0)
    time.interest <- sort(unique(time[nonMissingOutcome[nonMissingDeathFlag]]))
    ## trim the time points if the user has requested it
    if (!missing(ntime) && length(time.interest) > ntime) {
      time.interest <- time.interest[
         unique(round(seq.int(1, length(time.interest), length.out = ntime)))]
    }
    
  }
  
  else {#NULL for other families
    ## TBD TBD TBD
    ## Currently all REGR and CLAS analysis is restricted to a single outcome.
    r.dim <- 1
    event <- event.type <- cens <- time.interest <- cens <- time <- NULL
  }
  
  return(list(event = event, event.type = event.type, cens = cens, 
              time.interest = time.interest,
              time = time, r.dim = r.dim))
  
}


get.yvar.type <- function(fmly) {
  switch(fmly, "surv"  = c("T", "S"), "surv-CR"  = c("T", "S"), "regr"  = "R", "class" = "C")
}



get.grow.splitinfo <- function (fmly, splitrule, nsplit, event.type) {

  ##                                     ## HARD CODED ON NATIVE SIDE  
  splitrule.names <- c("logrank",        ## 1  
                       "",               ## 2
                       "logrankscore",   ## 3
                       "logrankCR",      ## 4
                       " ",              ## 5
                       " ",              ## 6
                       "regression",     ## 7
                       "classification", ## 8
                       "random")         ## 9

  ## Preliminary check for consistency.
  nsplit <- round(nsplit)
  if (nsplit < 0) {
    stop("Invalid nsplit value specified.")    
  }
  
  if (grepl("surv", fmly)) {    
    if (is.null(splitrule)) {
      ## No split rule specified, use default.
      if (length(event.type) ==  1) {
        splitrule.idx <- which(splitrule.names == "logrank")
      }
      else {
        splitrule.idx <- which(splitrule.names == "logrankCR")
      }
      splitrule <- splitrule.names[splitrule.idx]
    }
    else {
      ## User split rule specified.
      splitrule.idx <- which(splitrule.names == splitrule)
      if (length(splitrule.idx) != 1) {
        stop("Invalid split rule specified:  ", splitrule)
      }
      if ((length(event.type) ==  1) & (splitrule.idx == which(splitrule.names == "logrankCR"))) {
        stop("Cannot specify logrankCR splitting for right-censored data")
      }
    }
  }
  if (fmly == "class") {    
    if (is.null(splitrule)) {
      ## No split rule specified, use default.
      splitrule.idx <- which(splitrule.names == "classification")
      splitrule <- splitrule.names[splitrule.idx]      
    }
    else {
      ## User specified split rule.
      if ((splitrule != "classification") & (splitrule != "random")) {
        stop("Invalid split rule specified:  ", splitrule)
      }
      splitrule.idx <- which(splitrule.names == splitrule)
    }
  }
  if (fmly == "regr") {    
    if (is.null(splitrule)) {
      ## No split rule specified, use default.
      splitrule.idx <- which(splitrule.names == "regression")
      splitrule <- splitrule.names[splitrule.idx]      
    }
    else {
      ## User specified split rule.
      if ((splitrule != "regression") & (splitrule != "random")) {
        stop("Invalid split rule specified:  ", splitrule)
      }
      splitrule.idx <- which(splitrule.names == splitrule)
    }
  }

  ## Override an nsplit value of zero (0) in the case of pure random
  ## splitting.  It must be set to one (1).
  if ((splitrule == "random") & (nsplit == 0)) {
    nsplit <- 1
  }

  splitinfo <- list(name = splitrule, index = splitrule.idx, nsplit = nsplit)

  return (splitinfo)

}


get.grow.xvar.wt <- function(weight, n.xvar) {

  ## Set the xvar weight.
  if (is.null(weight)) {
    weight <- rep(1/n.xvar, n.xvar)
  }
  else {
    if (any(weight < 0) | length(weight) != n.xvar | all(weight == 0)) {
      weight <- rep(1/n.xvar, n.xvar)
    }
    else {
      weight <-weight/sum(weight)
    }
  }

  return (weight)
}


get.grow.mtry <- function (mtry = NULL, n.xvar, fmly) {
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > n.xvar) mtry <- max(1, min(mtry, n.xvar))
  }
  else {
    if (grepl("surv", fmly) | (fmly == "class")) {
      mtry <- max(ceiling(sqrt(n.xvar)), 1)
    }
    else {
      mtry <- max(ceiling(n.xvar/3), 1)
    }
  }

  return (mtry)
}


parseFormula <- function (formula, data) {

  if (!inherits(formula, "formula")) {
    stop("'formula' is not a formula object.")
  }

  if (is.null(data)) {
    stop("'data' is missing.")
  }
  
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }

  fNames <- all.vars(formula, max.names = 1e7)

  if ((all.names(formula)[2] == "Surv")) {
    if (sum(is.element(names(data), fNames[1:2])) != 2) {
        stop("Survival formula incorrectly specified.")
    }
    family <- "surv"
    yvar.names <- fNames[1:2]
  }
  else {
    ## Currently all REGR and CLAS analysis is restricted to a single outcome.
    if (sum(is.element(names(data), fNames[1])) != 1) {
        stop("Regression (or) classification formula incorrectly specified.")
    }
    yvar.names <- fNames[1]
    Y <- data[, yvar.names]
    if ( !(is.factor(Y) | is.real(Y) | is.integer(Y))) {
      stop("the y-outcome must be either real or a factor.")
    }
    if (is.factor(Y)) {
      family <- "class"
    }
    else {
      family <- "regr"
    }
  }

  return (list(fNames=fNames, family=family, yvar.names=yvar.names))

}

parseMissingData <- function(formula.obj, data) {

  ## parse the formula object
  yvar.names <- formula.obj$yvar.names
  resp <- data[, yvar.names, drop = FALSE]
  
  ## determine whether all the outcomes are missing
  ## works for any dimension
  col.resp.na <- apply(data[, yvar.names, drop = FALSE], 2, function(x) {all(is.na(x))})
  if (any(col.resp.na)) {
    stop("All records are missing for the yvar(s)")
  }

  ## remove all x variables with missing values for all records
  colPt <- apply(data, 2, function(x){all(is.na(x))})
  ## terminate if all columns have all missing data 
  if (sum(colPt) >= (ncol(data) - length(yvar.names))) {
    stop("All x-variables have all missing data:  analysis not meaningful.")
  }
  data <- data[, !colPt, drop = FALSE]

  ## remove all records with missing values for all outcomes(s) and xvar
  rowPt <- apply(data, 1, function(x){all(is.na(x))})
  if (sum(rowPt) == nrow(data)) {
    stop("Rows of the data have all missing data:  analysis not meaningful.")
  }
  data <- data[!rowPt,, drop = FALSE]

  ## return the NA processed data
  return(data)
  
}

resample <- function(x, size, ...) {
    if (length(x) <= 1) {
      if (!missing(size) && size == 0) x[FALSE] else x
    }
    else {
      sample(x, size, ...)
    }  
}




