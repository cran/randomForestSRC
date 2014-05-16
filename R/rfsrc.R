####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.5.1
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


rfsrc <- function(formula,
                  data,
                  ntree = 1000,
                  bootstrap = c("by.root", "by.node", "none"),
                  mtry = NULL,
                  nodesize = NULL,
                  nodedepth = NULL,
                  splitrule = NULL,
                  nsplit = 0,
                  split.null = FALSE,
                  importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
                  na.action = c("na.omit", "na.impute", "na.random"),
                  nimpute = 1,
                  ntime,
                  cause,
                  xvar.wt = NULL,
                  proximity = FALSE,
                  forest = TRUE,
                  var.used = c(FALSE, "all.trees", "by.tree"),
                  split.depth = c(FALSE, "all.trees", "by.tree"),
                  seed = NULL,
                  do.trace = FALSE,
                  membership = TRUE,
                  ...)
{
  user.option <- match.call(expand.dots = TRUE)
  impute.only <- is.hidden.impute.only(user.option)
  miss.tree <- is.hidden.impute.only(user.option)
  bootstrap <- match.arg(bootstrap, c("by.root", "by.node", "none"))
  importance <- match.arg(importance, c("permute", "random", "permute.ensemble", "random.ensemble", "none"))
  na.action <- match.arg(na.action, c("na.omit", "na.impute", "na.random"))
  proximity <- match.arg(as.character(proximity), c(FALSE, TRUE, "inbag", "oob", "all"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  if (split.depth == "FALSE") split.depth <- FALSE
  if (missing(formula) | (!missing(formula) && is.null(formula))) {
    formula <- as.formula("Unsupervised() ~ .")
  }
  if (missing(data)) stop("data is missing")
  formulaPrelim <- parseFormula(formula, data)
  if (any(is.na(data))) {
    data <- parseMissingData(formulaPrelim, data)
    miss.flag <- TRUE
  }
  else {
    miss.flag <- FALSE
  }
  formulaDetail <- finalizeFormula(formulaPrelim, data)
  ntree <- round(ntree)
  if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
  if (!is.null(nodesize) && nodesize < 1) stop("Invalid choice of 'nodesize'. Cannot be less than 1.")
  if (!is.null(nodedepth)) nodedepth = round(nodedepth) else nodedepth = -1
  nimpute <- round(nimpute)
  if (nimpute < 1) stop("Invalid choice of 'nimpute'.  Cannot be less than 1.")
  seed <- get.seed(seed)
  fmly <- formulaDetail$family
  xvar.names <- formulaDetail$xvar.names
  yvar.names <- formulaDetail$yvar.names
  # .. are there any x-variables?  (can happen when for multivariate formula)
  if (length(xvar.names) == 0) {
    stop("something seems wrong: your formula did not define any x-variables")
  }
  # .. are there any y-variables?  (do not test for the unsupervised case)
  if (fmly != "unsupv" && length(yvar.names) == 0) {
    stop("something seems wrong: your formula did not define any y-variables")
  }
  if (fmly == "class") {
    if (length(setdiff(levels(data[, yvar.names]), unique(data[, yvar.names]))) > 0) {
      warning("empty classes found when implementing classification\n")
    }
  }
  data <- rm.na.levels(data, xvar.names)
  data <- rm.na.levels(data, yvar.names)
  yfactor <- extract.factor(data, yvar.names)
  xfactor <- extract.factor(data, xvar.names)
  yvar.types <- get.yvar.type(fmly, data[, yvar.names, drop = FALSE])
  data <- finalizeData(c(yvar.names, xvar.names), data, na.action, miss.flag)
  xvar <- as.matrix(data[, xvar.names])
  rownames(xvar) <- colnames(xvar) <- NULL
  split.wt <- NULL
  statistics <- FALSE
  forest.wt <- FALSE
  n <- nrow(xvar)
  n.xvar <- ncol(xvar)
  mtry <- get.grow.mtry(mtry, n.xvar, fmly)
  xvar.wt  <- get.grow.x.wt(xvar.wt, n.xvar)
  split.wt <- get.grow.x.wt(split.wt, n.xvar)
  forest.wt <- match.arg(as.character(forest.wt), c(FALSE, TRUE, "inbag", "oob", "all"))
  yvar <- as.matrix(data[, yvar.names, drop = FALSE])
  if(dim(yvar)[2] == 0) {
    yvar <- NULL
  }
  if (miss.flag) {
    n.miss <- get.nmiss(xvar, yvar)
  }
  else {
    n.miss <- 0
  }
  if (impute.only && n.miss == 0) {
    return(data)
  }
  remove(data)
  big.data <- FALSE
  event.info <- get.grow.event.info(yvar, fmly, ntime = ntime)
  splitinfo <- get.grow.splitinfo(formulaDetail, splitrule, nsplit, event.info$event.type)
  sexp.dim <- get.sexp.dim(fmly, event.info$event.type, yfactor, splitinfo$name)
  if (fmly == "surv") {
    if (length(event.info$event.type) > 1) {
      if (missing(cause)) {
        cause.wt <- rep(1, length(event.info$event.type))
      }
      else {
        if (length(cause) == 1) {
          if (cause >= 1 && cause <= length(event.info$event.type)) {
            cause.wt <- rep(0, length(event.info$event.type))
            cause.wt[cause] <- 1
          }
          else {
            cause.wt <- rep(1, length(event.info$event.type))
          }
        }
        else {
          if (length(cause) == length(event.info$event.type) && all(cause >= 0) && !all(cause == 0)) {
            cause.wt <- cause / sum(cause)
          }
          else {
            cause.wt <- rep(1, length(event.info$event.type))
          }
        }
      }
    }
    else {
      cause.wt = 1
    }
    if (sexp.dim > 1) {
      fmly = "surv-CR"
    }
  }
  else {
    cause.wt <- NULL
  }
  nodesize <- get.grow.nodesize(fmly, nodesize)
  perf <- NULL
  if (bootstrap != "by.root") {
    importance <- "none"
    perf <- FALSE
  }
  if (fmly == "class+" | fmly == "regr+" | fmly == "mix+" | fmly == "unsupv") {
    importance <- "none"
    perf <- FALSE
  }
  if (impute.only) {
    na.action    <- "na.impute"
    forest       <- FALSE
    proximity    <- FALSE
    forest.wt    <- FALSE
    var.used     <- FALSE
    split.depth  <- FALSE
    membership   <- FALSE
    perf         <- FALSE
    importance   <- "none"
  }
  impute.only.bits <- get.impute.only(impute.only, n.miss)
  var.used.bits <- get.var.used(var.used)
  split.depth.bits <- get.split.depth(split.depth)
  importance.bits <- get.importance(importance)
  bootstrap.bits <- get.bootstrap(bootstrap)
  forest.bits <- get.forest(forest)
  proximity.bits <- get.proximity(TRUE, proximity)
  split.null.bits <- get.split.null(split.null)
  membership.bits <-  get.membership(membership)
  statistics.bits <- get.statistics(statistics)
  # This is dependent on impute.only and the family being initialized.
  perf.flag <- get.perf(perf, impute.only, fmly)
  perf.bits <-  get.perf.bits(perf.flag)
  forest.wt.bits <- get.forest.wt(TRUE, bootstrap, forest.wt)
  na.action.bits <- get.na.action(na.action)
  do.trace <- get.trace(do.trace)
  nativeOutput <- .Call("rfsrcGrow",
                        as.integer(do.trace),
                        as.integer(seed),
                        as.integer(impute.only.bits +
                                   var.used.bits +
                                   split.depth.bits +
                                   importance.bits +
                                   bootstrap.bits +
                                   forest.bits +
                                   proximity.bits +
                                   split.null.bits +
                                   perf.bits +
                                   membership.bits +
                                   statistics.bits),
                        as.integer(0 +
                                   na.action.bits),
                        as.integer(splitinfo$index),
                        as.integer(splitinfo$nsplit),
                        as.integer(mtry),
                        as.integer(if(is.na(formulaDetail$ytry)) 0 else formulaDetail$ytry),
                        as.integer(nodesize),
                        as.integer(nodedepth),
                        as.double(cause.wt),
                        as.integer(ntree),
                        as.integer(n),
                        as.integer(length(yvar.types)),
                        as.character(yvar.types),
                        as.integer(yfactor$nlevels),
                        as.double(as.vector(yvar)),
                        as.integer(n.xvar),
                        as.character(xfactor$generic.types),
                        as.integer(xfactor$nlevels),
                        as.double(xvar.wt),
                        as.double(split.wt),
                        as.double(xvar),
                        as.integer(length(event.info$time.interest)),
                        as.double(event.info$time.interest),
                        as.double(miss.tree),
                        as.integer(nimpute),
                        as.integer(get.rf.cores()))
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }
  if (n.miss > 0) {
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss, byrow = FALSE)
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.matrix(imputed.data[, -1])
    if (n.miss == 1) imputed.data <- t(imputed.data)
    nativeOutput$imputation <- NULL
    if (nimpute > 1) {
      if (grepl("surv", fmly)) {
        yvar[imputed.indv, 1] <- imputed.data[, 1]
        yvar[imputed.indv, 2] <- imputed.data[, 2]
        xvar[imputed.indv, ] <- imputed.data[, -c(1:2)]
      }
      else {
        if (!is.null(yvar.types)) {
          yvar[imputed.indv, ] <- imputed.data[, length(yvar.types)]
          xvar[imputed.indv, ] <- imputed.data[, -c(1:length(yvar.types))]
        }
        else {
          xvar[imputed.indv, ] <- imputed.data
        }
      }
      imputed.indv <- NULL
      imputed.data <- NULL
      imputedOOBData <- NULL
    }
    else {
      colnames(imputed.data) <- c(yvar.names, xvar.names)
      imputed.data <- as.data.frame(imputed.data)
    }
  }
  xvar <- as.data.frame(xvar)
  colnames(xvar) <- xvar.names
  xvar <- map.factor(xvar, xfactor)
  if (fmly != "unsupv") {
    yvar <- as.data.frame(yvar)
    colnames(yvar) <- yvar.names
  }
  else {
    yvar <- NULL
  }
  if (fmly != "unsupv") {
    yvar <- amatrix.remove.names(map.factor(yvar, yfactor))
  }
  if ((n.miss > 0) & (nimpute < 2)) {
    imputed.data <- map.factor(imputed.data, xfactor)
    if (fmly != "unsupv") {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
  }
  if (forest) {
    nativeArray <- as.data.frame(cbind(nativeOutput$treeID,
                                       nativeOutput$nodeID,
                                       nativeOutput$parmID,
                                       nativeOutput$contPT,
                                       nativeOutput$mwcpSZ))
    names(nativeArray) <- c("treeID", "nodeID", "parmID", "contPT", "mwcpSZ")
    nativeFactorArray <- nativeOutput$mwcpPT
    forest.out <- list(nativeArray = nativeArray,
                       nativeFactorArray = nativeFactorArray,
                       totalNodeCount = dim(nativeArray)[1],
                       nodesize = nodesize,
                       nodedepth = nodedepth,
                       split.null = split.null,
                       ntree = ntree,
                       family = fmly,
                       yvar = yvar,
                       yvar.names = yvar.names,
                       xvar = xvar,
                       xvar.names = xvar.names,
                       seed = nativeOutput$seed,
                       bootstrap.bits = bootstrap.bits,
                       bootstrap = bootstrap)
    if (grepl("surv", fmly)) {
      forest.out$time.interest <- event.info$time.interest
    }
    class(forest.out) <- c("rfsrc", "forest", fmly)
    if (big.data) {
      class(forest.out) <- c(class(forest.out), "bigdata")
    }
  }
  else {
    forest.out <- NULL
  }
  if (grepl("surv", fmly)) {
    if (fmly == "surv-CR") {
      ens.names <- list(NULL, NULL, c(paste("condCHF.", 1:(sexp.dim), sep = "")))
      cif.names <- list(NULL, NULL, c(paste("CIF.", 1:(sexp.dim), sep = "")))
      err.names <- list(paste("event.", 1:(sexp.dim), sep = ""), NULL)
      vimp.names <- list(paste("event.", 1:(sexp.dim), sep = ""), xvar.names)
      mortality.names <- list(NULL, paste("event.", 1:(sexp.dim), sep = ""))
    }
    else {
      survival.names <- ens.names <- list(NULL, NULL, NULL)
      mortality.names <- err.names <- list(NULL, NULL)
      vimp.names <- list(NULL, xvar.names)
    }
  }
  else {
    if (fmly == "class") {
      err.names <- list(c("all", yfactor$levels[[1]]), NULL)
      vimp.names <- list(c("all", yfactor$levels[[1]]), xvar.names)
      ens.names <- list(NULL, yfactor$levels[[1]], NULL)
    }
    else {
      err.names <- list(NULL, NULL)
      ens.names <- list(NULL, NULL, NULL)
      vimp.names <- list(NULL, xvar.names)
    }
  }
  if (proximity != FALSE) {
    proximity.out <- matrix(0, n, n)
    count <- 0
    for (k in 1:n) {
      proximity.out[k,1:k] <- nativeOutput$proximity[(count+1):(count+k)]
      proximity.out[1:k,k] <- proximity.out[k,1:k]
      count <- count + k
    }
    nativeOutput$proximity <- NULL
  }
  else {
    proximity.out <- NULL
  }
  if (forest.wt != FALSE) {
    forest.wt.out <- matrix(nativeOutput$weight, c(n, n), byrow = TRUE)
    nativeOutput$weight <- NULL
  }
  else {
    forest.wt.out <- NULL
  }
  if (membership) {
    membership.out <- matrix(nativeOutput$nodeMembership, c(n, ntree))
    inbag.out <- matrix(nativeOutput$bootMembership, c(n, ntree))
    nativeOutput$nodeMembership <- NULL
    nativeOutput$bootMembership <- NULL
  }
  else {
    membership.out <- NULL
    inbag.out <- NULL
  }
  if (var.used != FALSE) {
    if (var.used == "all.trees") {
      var.used.out <- nativeOutput$varUsed
      names(var.used.out) <- xvar.names
    }
    else {
      var.used.out <- matrix(nativeOutput$varUsed, nrow = ntree, byrow = TRUE)
      colnames(var.used.out) <- xvar.names
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
  if (statistics) {
    node.stats <- as.data.frame(cbind(nativeOutput$spltST))
    names(node.stats) <- c("spltST")
  }
  else {
    node.stats <- NULL
  }
  if (importance != "none") {
    VIMP <- atmatrix(nativeOutput$importance, n.xvar, vimp.names, keep.names = TRUE)
    nativeOutput$importance <- NULL
  }
  else {
    VIMP <- NULL
  }
  if (!is.null(nativeOutput$performance)) {
    ERR <- atmatrix(nativeOutput$performance, ntree, err.names)
    nativeOutput$performance <- NULL
  }
  else {
    ERR <- NULL
  }
  rfsrcOutput <- list(
    call = match.call(),
    family = fmly,
    n = n,
    ntree = ntree,
    nimpute = nimpute,
    mtry = mtry,
    nodesize = nodesize,
    nodedepth = nodedepth,
    splitrule = splitinfo$name,
    nsplit = splitinfo$nsplit,
    yvar = yvar,
    yvar.names = yvar.names,
    xvar = xvar,
    xvar.names = xvar.names,
    xvar.wt = xvar.wt,
    split.wt = split.wt,
    leaf.count = nativeOutput$leafCount,
    proximity = proximity.out,
    forest = forest.out,
    membership = membership.out,
    inbag = inbag.out,
    var.used = var.used.out,
    imputed.indv = (if (n.miss > 0) imputed.indv else NULL),
    imputed.data = (if (n.miss > 0) imputed.data else NULL),
    split.depth  = split.depth.out,
    err.rate = ERR,
    importance = VIMP
  )
  remove(yvar)
  remove(xvar)
  nativeOutput$leafCount <- NULL
  remove(proximity.out)
  remove(forest.out)
  remove(forest.wt.out)
  remove(membership.out)
  remove(inbag.out)
  remove(var.used.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(split.depth.out)
  remove(ERR)
  remove(VIMP)
  if (grepl("surv", fmly)) {
    predicted <- (if (!is.null(nativeOutput$fullMortality))
             amatrix(nativeOutput$fullMortality, c(n, sexp.dim), mortality.names) else NULL)
    nativeOutput$fullMortality <- NULL
    survOutput <- list(predicted = predicted)
    remove(predicted)
    predicted.oob <- (if (!is.null(nativeOutput$oobMortality))
             amatrix(nativeOutput$oobMortality, c(n, sexp.dim), mortality.names) else NULL)
    nativeOutput$oobMortality <- NULL
    survOutput <- c(survOutput, predicted.oob = list(predicted.oob))
    survival <-  (if (!is.null(nativeOutput$fullSurvival))
             adrop(array(nativeOutput$fullSurvival, c(n, length(event.info$time.interest), sexp.dim),
                dimnames=survival.names), sexp.dim) else NULL)
    nativeOutput$fullSurvival <- NULL
    survOutput <- c(survOutput, survival = list(survival))
    remove(survival)
    survival.oob = (if (!is.null(nativeOutput$oobSurvival))
             adrop(array(nativeOutput$oobSurvival, c(n, length(event.info$time.interest), sexp.dim),
                dimnames=survival.names), sexp.dim) else NULL)
    nativeOutput$oobSurvival <- NULL
    survOutput <- c(survOutput, survival.oob = list(survival.oob))
    remove(survival.oob)
    chf <- (if (!is.null(nativeOutput$fullEnsemble))
           adrop(array(nativeOutput$fullEnsemble, c(n, length(event.info$time.interest), sexp.dim),
                dimnames=ens.names), sexp.dim) else NULL)
    nativeOutput$fullEnsemble <- NULL
    survOutput <- c(survOutput, chf = list(chf))
    remove(chf)
    chf.oob <- (if (!is.null(nativeOutput$oobEnsemble))
          adrop(array(nativeOutput$oobEnsemble,  c(n, length(event.info$time.interest), sexp.dim),
               dimnames=ens.names), sexp.dim) else NULL)
    nativeOutput$oobEnsemble <- NULL
    survOutput = c(survOutput, chf.oob = list(chf.oob))
    remove(chf.oob)
    cif <- (if (!is.null(nativeOutput$fullCIF))
           adrop(array(nativeOutput$fullCIF, c(n, length(event.info$time.interest), sexp.dim),
               dimnames=cif.names), sexp.dim) else NULL)
    nativeOutput$fullCIF <- NULL
    survOutput <- c(survOutput, cif = list(cif))
    remove(cif)
    cif.oob <- (if (!is.null(nativeOutput$oobCIF))
           adrop(array(nativeOutput$oobCIF,  c(n, length(event.info$time.interest), sexp.dim),
               dimnames=cif.names), sexp.dim) else NULL)
    nativeOutput$oobCIF <- NULL
    survOutput = c(survOutput, cif.oob = list(cif.oob))
    remove(cif.oob)
    survOutput = c(
      survOutput, list(
       time.interest = event.info$time.interest, ndead = sum(na.omit(event.info$cens) != 0))
    )
    rfsrcOutput <- c(rfsrcOutput, survOutput)
  }
  else {
    if (fmly == "class") {
      predicted <- (if (!is.null(nativeOutput$fullEnsemble))
          adrop(array(nativeOutput$fullEnsemble, c(n, sexp.dim - 1, 1), dimnames=ens.names), 1, TRUE) else NULL)
      nativeOutput$fullEnsemble <- NULL
      clasOutput <- list(predicted = predicted)
      classResp <- (if (!is.null(predicted)) bayes.rule(predicted) else NULL)
      remove(predicted)
      predicted.oob <- (if (!is.null(nativeOutput$oobEnsemble))
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
      predicted <- (if (!is.null(nativeOutput$fullEnsemble)) c(nativeOutput$fullEnsemble) else NULL)
      nativeOutput$fullEnsemble <- NULL
      regrOutput <- list(predicted = predicted)
      remove(predicted)
      predicted.oob <- (if (!is.null(nativeOutput$oobEnsemble)) c(nativeOutput$oobEnsemble) else NULL)
      nativeOutput$oobEnsemble <- NULL
      regrOutput <- c(regrOutput, predicted.oob = list(predicted.oob))
      remove(predicted.oob)
      rfsrcOutput = c(rfsrcOutput, regrOutput)
    }
  }
  class(rfsrcOutput) <- c("rfsrc", "grow", fmly)
  if (big.data) {
    class(rfsrcOutput) <- c(class(rfsrcOutput), "bigdata")
  }
  return(rfsrcOutput)
}
