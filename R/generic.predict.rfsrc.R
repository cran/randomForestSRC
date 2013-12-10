####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.4
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


generic.predict.rfsrc <-
  function(object,
           newdata, 
           importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
           importance.xvar,
           na.action = c("na.omit", "na.impute"),
           outcome = c("train", "test"),
           proximity = FALSE,
           var.used = c(FALSE, "all.trees", "by.tree"),
           split.depth = c(FALSE, "all.trees", "by.tree"),
           seed = NULL,
           do.trace = FALSE,
           membership = TRUE,
           statistics = FALSE,
           ptn.count = 0,
           ...)
{
  if (missing(object)) {
    stop("object is missing!")
  }
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2  &
      sum(inherits(object, c("rfsrc", "partial"), TRUE) == c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'")
  }
  if (sum(inherits(object, c("rfsrc", "partial"), TRUE) == c(1, 2)) != 2) {
    partial.class <- FALSE
  }
  else {
   partial.class <- TRUE
  }
  importance <- match.arg(importance[1],
      c("permute", "random", "permute.ensemble", "random.ensemble", "none", 
        "permute.joint", "random.joint", "permute.joint.ensemble", "random.joint.ensemble"))
  if (grepl("joint", importance)) {
    vimp.joint <- TRUE
  }
  else {
    vimp.joint <- FALSE
  }
  importance.xvar <- get.importance.xvar(importance.xvar, importance, object)
  importance.xvar.idx <- match(importance.xvar, object$xvar.names)
  na.action <- match.arg(na.action, c("na.omit", "na.impute")) 
  outcome <- match.arg(outcome, c("train", "test"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (split.depth == "FALSE") split.depth <- FALSE
  if (missing(newdata)) {
    outcome <- "train"
    restore.only <- TRUE
  }
  else {
    restore.only <- FALSE
  }
  restore.only.bits <- get.restore.only(restore.only)
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
    if (is.null(object$forest)) {
      stop("The forest is empty.  Re-run rfsrc (grow) call with forest=TRUE")
    }
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
    else {
      big.data <- FALSE
    }
    object <- object$forest
  }
  else {
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
    else {
      big.data <- FALSE
    }
  }
  outcome.target <- get.outcome.target(object$family, outcome.target)
  yvar.types <- get.yvar.type(object$family, object$yvar)
  yvar.target <- get.yvar.target(object$family, yvar.types, outcome.target)
  object$yvar <- as.data.frame(object$yvar)
  colnames(object$yvar) <- object$yvar.names
  yfactor <- extract.factor(object$yvar)
  event.info <- get.event.info(object)
  sexp.dim <- get.sexp.dim(object$family, event.info$event.type, yfactor)
  cr.bits <- get.cr.bits(object$family)
  xfactor <- extract.factor(object$xvar)
  if (object$family == "unsupv") {
    outcome <- "train"
    perf.flag <- FALSE 
    importance <- "none"
  }
  if (missing(newdata)) {
    if (object$family != "unsupv") {
      newdata <- cbind(object$yvar, object$xvar)
    }
    else {
      newdata <- object$xvar
    }
  }
  if (!is.data.frame(newdata)) {
    stop("test data must be a data frame.")
  }
  if (!restore.only) {
    if (!partial.class) {
      newdata <- newdata[, is.element(names(newdata),
                       c(object$yvar.names, object$xvar.names)), drop = FALSE]
      newdata <- rm.na.levels(newdata, object$xvar.names)
      newdata.xfactor <- extract.factor(newdata, object$xvar.names)
      if (!setequal(xfactor$factor, newdata.xfactor$factor)) {
        stop("x-variable factors from test data do not match original training data")
      }
      if (!setequal(xfactor$order, newdata.xfactor$order)) {
        stop("(ordered) x-variable factors from test data do not match original training data")
      }
      if (object$family == "class") {    
        if (sum(is.element(names(newdata), object$yvar.names)) > 0) {
          newdata <- rm.na.levels(newdata, object$yvar.names)
          newdata.yfactor <- extract.factor(newdata, object$yvar.names)
          if (!setequal(yfactor$factor, newdata.yfactor$factor)) {
            stop("class outcome from test data does not match original training data")
          }
          if (!setequal(yfactor$order, newdata.yfactor$order)) {
            stop("(ordered) class outcome from test data does not match original training data")
          }
        }
      }
      if (length(object$xvar.names) != sum(is.element(object$xvar.names, names(newdata)))) {
        stop("x-variables in test data do not match original training data")
      }
      newdata[, object$xvar.names] <- check.factor(object$xvar, newdata[, object$xvar.names, drop = FALSE], xfactor)
      if (object$family == "class") {
          if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
            newdata[, object$yvar.names] <- check.factor(object$yvar, newdata[, object$yvar.names, drop = FALSE], yfactor)
        }
      }
      if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
        fnames <- c(object$yvar.names, object$xvar.names)
      }
      else {
        fnames <- object$xvar.names
      }
      newdata <- finalizeData(fnames, newdata, na.action) 
      xvar.newdata  <- as.matrix(newdata[, object$xvar.names, drop = FALSE])
      n.newdata <- nrow(newdata)
      if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
        yvar.newdata <- as.matrix(newdata[, object$yvar.names, drop = FALSE])
        if(dim(yvar.newdata)[2] == 0) {
          yvar.newdata = NULL
        }
        event.info.newdata <- get.grow.event.info(yvar.newdata, object$family, need.deaths = FALSE)
        r.dim.newdata <- event.info.newdata$r.dim
        if (r.dim.newdata > 0) {
          perf.flag <- TRUE
        }
        else {
          perf.flag <- FALSE
        }
        if (grepl("surv", object$family) &&
            length(setdiff(na.omit(event.info.newdata$cens), na.omit(event.info$cens))) > 1) {
          stop("survival events in test data do not match training data")
        }
      }       
      else {
        if (outcome == "test") {
          stop("outcome=TEST, but the test data has no y values, which is not permitted")
        }
        r.dim.newdata <- 0
        yvar.newdata <-  NULL
        perf.flag <- FALSE
        importance <- "none"
      }
    } 
    else {
      newdata <- as.data.frame(data.matrix(newdata))
      xvar.newdata <- as.matrix(newdata)
      n.newdata <- nrow(xvar.newdata)
      r.dim.newdata <- 0
      yvar.newdata <- NULL 
      perf.flag <- FALSE
      importance <- "none"
    }
    if (outcome != "test") {
      rownames(xvar.newdata) <- colnames(xvar.newdata) <- NULL
    }
    remove(newdata)
  }
  else {
    n.newdata <- 0
    r.dim.newdata <- 0
    xvar.newdata <- NULL
    yvar.newdata <-  NULL
    outcome <- "train"
    if (object$bootstrap != "by.root" | object$family == "unsupv") {
      importance <- "none"
      perf.flag <- FALSE
    }
    else {
      perf.flag <- TRUE
    }
  } 
  if (outcome == "train") {
    xvar <- as.matrix(data.matrix(object$xvar))    
    yvar <- as.matrix(data.matrix(object$yvar))
  }
  else {
    xvar <- xvar.newdata
    yvar <- yvar.newdata
    restore.only <- TRUE 
    perf.flag <- TRUE
    n.newdata <- 0
    r.dim.newdata <- 0
  }
  r.dim <- ncol(cbind(yvar))
  rownames(xvar) <- colnames(xvar) <- NULL
  n.xvar <- ncol(xvar)
  n <- nrow(xvar)
  split.null <- object$split.null
  ntree <- object$ntree
  importance.bits <- get.importance(importance)
  var.used.bits <- get.var.used(var.used)
  split.null.bits <- get.split.null(split.null)
  split.depth.bits <- get.split.depth(split.depth)
  membership.bits <-  get.membership(membership)  
  seed <- get.seed(seed)
  proximity.bits <- get.proximity(restore.only, proximity)
  outcome.bits <- get.outcome(outcome)
  perf.bits <-  get.perf.bits(perf.flag)
  statistics.bits <- get.statistics(statistics)
  do.trace <- get.trace(do.trace)
  nativeOutput <- .Call("rfsrcPredict",
                        as.integer(do.trace),
                        as.integer(seed),
                        as.integer(restore.only.bits +
                                   importance.bits +
                                   object$bootstrap.bits + 
                                   proximity.bits +
                                   split.null.bits +                                   
                                   split.depth.bits +
                                   var.used.bits +
                                   outcome.bits +
                                   perf.bits +
                                   membership.bits +
                                   cr.bits +
                                   statistics.bits),
                        as.integer(ntree),
                        as.integer(n),
                        as.integer(r.dim),
                        as.character(yvar.types),
                        as.integer(yvar.target),
                        as.integer(yfactor$nlevels),
                        as.double(as.vector(yvar)),
                        as.integer(ncol(xvar)),
                        as.character(xfactor$generic.types),
                        as.integer(xfactor$nlevels),
                        as.double(xvar),
                        as.integer(n.newdata),
                        as.integer(r.dim.newdata),
                        as.double(if (outcome != "test") yvar.newdata else NULL),                        
                        as.double(if (outcome != "test") xvar.newdata else NULL),
                        as.integer(length(event.info$time.interest)),
                        as.double(event.info$time.interest),
                        as.integer((object$nativeArray)$treeID),
                        as.integer((object$nativeArray)$nodeID),
                        as.integer((object$nativeArray)$parmID),
                        as.double((object$nativeArray)$contPT),
                        as.integer((object$nativeArray)$mwcpSZ),                        
                        as.integer(object$nativeFactorArray),
                        as.integer(object$totalNodeCount),
                        as.integer(object$seed),
                        as.integer(length(importance.xvar.idx)),
                        as.integer(importance.xvar.idx),
                        as.integer(ptn.count),
                        as.integer(get.rf.cores()))
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
  if (restore.only) {
    n.miss <- get.nmiss(xvar, yvar)
  }
  else {
    n.miss <- get.nmiss(xvar.newdata, yvar.newdata)
  }
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss)
    nativeOutput$imputation <- NULL
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.data.frame(imputed.data[, -1, drop = FALSE])
    if (r.dim.newdata > 0 | perf.flag) {
      colnames(imputed.data) <- c(object$yvar.names, object$xvar.names)
    }
    else {
      colnames(imputed.data) <- object$xvar.names      
    }
  }
  if ( (!partial.class) & (!restore.only | outcome == "test") ) {
    xvar.newdata <- as.data.frame(xvar.newdata)
    colnames(xvar.newdata) <- object$xvar.names
    xvar.newdata <- map.factor(xvar.newdata, xfactor)
    if (perf.flag) {
      yvar.newdata <- as.data.frame(yvar.newdata)
      colnames(yvar.newdata) <- object$yvar.names
      yvar.newdata <- map.factor(yvar.newdata, yfactor)
    }
  }
  if (n.miss > 0) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (perf.flag) {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  if (grepl("surv", object$family)) {
    if (object$family == "surv-CR") {
      ens.names <- list(NULL, NULL, paste("condCHF.", 1:(sexp.dim), sep = ""))
      cif.names <- list(NULL, NULL, paste("CIF.", 1:(sexp.dim), sep = ""))
      err.names <- list(paste("event.", 1:(sexp.dim), sep = ""), NULL)
      vimp.names <- list(paste("event.", 1:(sexp.dim), sep = ""),
                           if(vimp.joint) "joint" else importance.xvar)
      mortality.names <- list(NULL, paste("event.", 1:(sexp.dim), sep = ""))
    }
    else {
      survival.names <- ens.names <- list(NULL, NULL, NULL)
      mortality.names <- err.names <- list(NULL, NULL)
      vimp.names <- list(NULL, if (vimp.joint) "joint" else importance.xvar)
    }
  }
  else {
    if (object$family == "class") {
      err.names <- list(c("all", yfactor$levels[[1]]), NULL)
      vimp.names <- list(c("all", yfactor$levels[[1]]), if(vimp.joint) "joint" else importance.xvar)
      ens.names <- list(NULL, yfactor$levels[[1]], NULL)
    }
    else {
      err.names <- list(NULL, NULL)
      ens.names <- list(NULL, NULL, NULL)
      vimp.names <- list(NULL, if(vimp.joint) "joint" else importance.xvar)
    }
  }
  if (proximity != FALSE) {
    if (restore.only) {
      prox.n <- n
    }
    else {
      prox.n <- n.newdata
    }
    proximity.out <- matrix(0, prox.n, prox.n)
    count <- 0
    for (k in 1:prox.n) {
      proximity.out[k,1:k] <- nativeOutput$proximity[(count+1):(count+k)]
      proximity.out[1:k,k] <- proximity.out[k,1:k]
      count <- count + k
    }
    proximity.out <- proximity.out / diag(proximity.out)
    nativeOutput$proximity <- NULL
  }
  else {
    proximity.out <- NULL
  }
  if (membership) {
    if (restore.only) {
      membership.n <- n
    }
    else {
      membership.n <- n.newdata
    }
    membership.out <- matrix(nativeOutput$nodeMembership, c(membership.n, ntree))
    if (ptn.count > 0) {
      ptn.membership.out <- matrix(nativeOutput$ptnMembership, c(membership.n, ntree))
    }
    else {
      ptn.membership.out <- NULL
    }
    inbag.out <- matrix(nativeOutput$bootMembership, c(membership.n, ntree))
    nativeOutput$nodeMembership <- NULL
    nativeOutput$bootMembership <- NULL
    if (ptn.count > 0) {
      nativeOutput$ptnMembership <- NULL
    }
  }
  else {
    membership.out <- NULL
    ptn.membership.out <- NULL    
    inbag.out <- NULL
  }
  if (var.used != FALSE) {
    if (var.used == "all.trees") {
      var.used.out <- nativeOutput$varUsed
      names(var.used.out) <- object$xvar.names
    }
    else {
      var.used.out <- matrix(nativeOutput$varUsed, nrow = ntree, byrow = TRUE)
      colnames(var.used.out) <- object$xvar.names
    }
    nativeOutput$varUsed <- NULL
  }
  else {
    var.used.out <-  NULL
  }
  if (split.depth != FALSE) {
    if (split.depth == "all.trees") {
      split.depth.out <- array(nativeOutput$splitDepth, c(n, n.xvar))
    }
    else {
      split.depth.out <- array(nativeOutput$splitDepth, c(n, n.xvar, ntree))
    }
    nativeOutput$splitDepth <- NULL
  }
  else {
    split.depth.out <-  NULL
  }
  if (statistics == TRUE) {
    node.stats <- as.data.frame(cbind(nativeOutput$spltST))
    names(node.stats) <- c("spltST")
  }
  else {
    node.stats <- NULL
  }
  if (perf.flag) {
    ERR <- atmatrix(nativeOutput$performance, ntree, err.names)
    nativeOutput$performance <- NULL
  }
  else {
    ERR <- NULL
  }
  if ((importance != "none") & perf.flag) {
    VIMP <- atmatrix(nativeOutput$importance, (if (vimp.joint) 1 else length(importance.xvar)),
                      vimp.names, keep.names = TRUE)
    nativeOutput$importance <- NULL
  }
  else {
    VIMP <- NULL
  }
  rfsrcOutput <- list(
    call = match.call(),
    family = object$family,
    n = (if (restore.only) n else n.newdata),
    ntree = ntree,
    yvar = (if (perf.flag | partial.class) {
      if (outcome != "test" & (restore.only | partial.class))
        amatrix.remove.names(object$yvar) else amatrix.remove.names(yvar.newdata)} else NULL),
    yvar.names = object$yvar.names,
    xvar = (if(outcome != "test" & restore.only) object$xvar else xvar.newdata),                      
    xvar.names = object$xvar.names,
    leaf.count = nativeOutput$leafCount,
    forest = object,
    proximity = proximity.out,
    membership = membership.out,
    ptn.membership = ptn.membership.out,                      
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss>0) imputed.indv else NULL),
    imputed.data = (if (n.miss>0) imputed.data else NULL),
    split.depth  = split.depth.out,
    err.rate = ERR,
    importance = VIMP
  )
  nativeOutput$leafCount <- NULL
  object.family <- object$family
  remove(object)
  remove(proximity.out)
  remove(membership.out)
  remove(ptn.membership.out)  
  remove(inbag.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(var.used.out)
  remove(split.depth.out)
  remove(node.stats)
  remove(ERR)
  remove(VIMP)
  if (grepl("surv", object.family)) {
    predicted <- (if (!is.null(nativeOutput$fullMortality))
             amatrix(nativeOutput$fullMortality, c(if (restore.only) n else n.newdata, sexp.dim),
                     mortality.names) else NULL)
    nativeOutput$fullMortality
    survOutput <- list(predicted = predicted)
    remove(predicted)
    predicted.oob <-  (if (!is.null(nativeOutput$oobMortality))
             amatrix(nativeOutput$oobMortality, c(n, sexp.dim), mortality.names) else NULL)
    nativeOutput$oobMortality <- NULL
    survOutput <- c(survOutput, predicted.oob = list(predicted.oob))
    remove(predicted.oob)
    survival = (if (!is.null(nativeOutput$fullSurvival))
             adrop(array(nativeOutput$fullSurvival, c(if (restore.only) n else n.newdata,
                 length(event.info$time.interest), sexp.dim), dimnames=survival.names), sexp.dim) else NULL)
    nativeOutput$fullSurvival <- NULL
    survOutput <- c(survOutput, survival = list(survival))
    remove(survival)
    survival.oob = (if (!is.null(nativeOutput$oobSurvival))
             adrop(array(nativeOutput$oobSurvival, c(n, length(event.info$time.interest), sexp.dim),
                 dimnames=survival.names), sexp.dim) else NULL)
    nativeOutput$oobSurvival <- NULL
    survOutput <- c(survOutput, survival.oob = list(survival.oob))
    remove(survival.oob)
    chf = (if (!is.null(nativeOutput$fullEnsemble))
             adrop(array(nativeOutput$fullEnsemble, c(if (restore.only) n else n.newdata,
                 length(event.info$time.interest), sexp.dim),
                 dimnames = ens.names), sexp.dim) else NULL)
    nativeOutput$fullEnsemble <- NULL
    survOutput <- c(survOutput, chf = list(chf))
    remove(chf)
    chf.oob = (if (!is.null(nativeOutput$oobEnsemble))
             adrop(array(nativeOutput$oobEnsemble,  c(n, length(event.info$time.interest), sexp.dim),
                 dimnames=ens.names), sexp.dim) else NULL)
    nativeOutput$oobEnsemble <- NULL
    survOutput = c(survOutput, chf.oob = list(chf.oob))
    remove(chf.oob)
    cif = (if (!is.null(nativeOutput$fullCIF))
             adrop(array(nativeOutput$fullCIF, c(if (restore.only) n else n.newdata,
                 length(event.info$time.interest), sexp.dim),
                 dimnames = cif.names), sexp.dim) else NULL)
    nativeOutput$fullCIF <- NULL
    survOutput <- c(survOutput, cif = list(cif))
    remove(cif)
    cif.oob = (if (!is.null(nativeOutput$oobCIF))
             adrop(array(nativeOutput$oobCIF,  c(n, length(event.info$time.interest), sexp.dim),
                 dimnames=cif.names), sexp.dim) else NULL)
    nativeOutput$oobCIF <- NULL
    survOutput = c(survOutput, cif.oob = list(cif.oob))
    remove(cif.oob)
    survOutput = c(
      survOutput, list(
      time.interest = event.info$time.interest,
      ndead = (if (perf.flag) sum((if (restore.only) yvar[, 2] else yvar.newdata[, 2]) !=0 , na.rm=TRUE) else NULL))
    )
    rfsrcOutput = c(rfsrcOutput, survOutput)
  }
  else {
    if (object.family == "class") {        
      predicted <- adrop(array(nativeOutput$fullEnsemble, c(if (restore.only) n else n.newdata, sexp.dim - 1, 1),
            dimnames=ens.names), 1, TRUE)
      nativeOutput$fullEnsemble <- NULL
      clasOutput <- list(predicted = predicted)
      classResp <- (if (!is.null(predicted)) bayes.rule(predicted) else NULL)
      remove(predicted)
      predicted.oob = (if (!is.null(nativeOutput$oobEnsemble))
             adrop(array(nativeOutput$oobEnsemble, c(n, sexp.dim - 1, 1), dimnames=ens.names), 1, TRUE) else NULL)
      nativeOutput$oobEnsemble <- NULL
      clasOutput <- c(clasOutput, predicted.oob = list(predicted.oob))
      classResp.oob <- (if (!is.null(predicted.oob)) bayes.rule(predicted.oob) else NULL)
      remove(predicted.oob)
      clasOutput <- c(clasOutput, class = list(classResp), class.oob = list(classResp.oob))
      remove(classResp)
      remove(classResp.oob)
      rfsrcOutput = c(rfsrcOutput, clasOutput)
    }
    else {
      predicted <- c(nativeOutput$fullEnsemble)
      nativeOutput$fullEnsemble <- NULL
      regrOutput <- list(predicted = predicted)
      remove(predicted)
      predicted.oob <-  (if (!is.null(nativeOutput$oobEnsemble)) c(nativeOutput$oobEnsemble) else NULL)
      nativeOutput$oobEnsemble <- NULL
      regrOutput <- c(regrOutput, predicted.oob = list(predicted.oob))
      remove(predicted.oob)
      rfsrcOutput = c(rfsrcOutput, regrOutput)
    }
  }
  class(rfsrcOutput) <- c("rfsrc", "predict",   object.family)
  return(rfsrcOutput)
}
