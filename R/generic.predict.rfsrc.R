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
    restore.only = FALSE,
    ...)
{

  ## Incoming parameter checks.  All are fatal.
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
  
  ## Verify the importance option
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

  ## Verify key options
  na.action <- match.arg(na.action, c("na.omit", "na.impute")) 
  outcome <- match.arg(outcome, c("train", "test"))
  var.used <- match.arg(as.character(var.used), c("FALSE", "all.trees", "by.tree"))
  if (var.used == "FALSE") var.used <- FALSE
  split.depth <- match.arg(as.character(split.depth),  c("FALSE", "all.trees", "by.tree"))
  if (split.depth == "FALSE") split.depth <- FALSE

  ## Assignment for restore.only
  restore.only.bits <- get.restore.only(restore.only)
  #if outcome=test it will be convenient for processing to assume restore.only = FALSE
  #later after pre-processing we change this to TRUE
  #thus outcome=test is always restore.only mode for the native code
  if (outcome == "test") {
    restore.only <- FALSE
  }
   
  ## acquire the forest.
  ## memory management "big.data" not currently implemented: TBD TBD TBD 
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
    ## object is already a forest
    if (inherits(object, "bigdata")) {
      big.data <- TRUE
    }
    else {
      big.data <- FALSE
    }
  }

  ## Initialize the yvar types
  yvar.types <- get.yvar.type(object$family)

  ## Determine the immutable yvar factor map which is needed for classification sexp dimensioning
  ## First convert object$yvar to a data frame which is required for factor processing
  object$yvar <- as.data.frame(object$yvar)
  colnames(object$yvar) <- object$yvar.names
  yfactor <- extract.factor(object$yvar)

  ## Get event information for survival families
  event.info <- get.event.info(object)

  ## SEXP dimensioning
  sexp.dim <- get.sexp.dim(object$family, event.info$event.type, yfactor)

  ## CR.bits assignment
  cr.bits <- get.cr.bits(object$family)
 
  ## Determine the immutable xvar factor map.
  xfactor <- extract.factor(object$xvar)

  ## Coherence check for test data
  ## Use the training data when test data is not present
  if (missing(newdata)) {
    newdata <- cbind(object$yvar, object$xvar)
  }
  if (!is.data.frame(newdata)) {
    stop("test data must be a data frame.")
  }

  #####################################################################
  ## From the native code's perspective, PRED mode can process one (1)
  ## or two (2) data sets.  If one (1) data set is sent in, we assume
  ## we wish to restore the forest with original-training or
  ## pseudo-training data.  If two (2) data sets are sent in, we
  ## assume we are sending in the original-training data set, and a
  ## test data set.  When one (1) data set is sent in, we call this
  ## REST mode, even though it is a sub-mode of PRED in the native
  ## code.
  ###################################################################

  ############## NON-RESTORE MODE ################
  if (!restore.only) {

    ## Automatic pass given if object is of internal class "partial".
    ## Send void outcomes to the native code. Factors must not have NA
    ## in their levels (but not enforced here)
    if (!partial.class) {

      ## Filter the test data based on the formula.
      newdata <- newdata[, is.element(names(newdata),
                       c(object$yvar.names, object$xvar.names)), drop = FALSE]

      ## Check that test/train factors are the same.  If factor has an
      ## NA in its levels, remove it.  Confirm factor labels overlap.
      newdata <- rm.na.levels(newdata, object$xvar.names)
      newdata.xfactor <- extract.factor(newdata, object$xvar.names)
      
      if (!setequal(xfactor$factor, newdata.xfactor$factor)) {
        stop("x-variable factors from test data do not match original training data")
      }
      if (!setequal(xfactor$order, newdata.xfactor$order)) {
        stop("(ordered) x-variable factors from test data do not match original training data")
      }

      # In classification we need to check that train/test y-outcomes are compatible.
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
      
      ## Confirm coherence of train/test xvar. 
      if (length(object$xvar.names) != sum(is.element(object$xvar.names, names(newdata)))) {
        stop("x-variables in test data do not match original training data")
      }

      ## Force test factor levels to equal grow factor levels (this last step
      ## is crucial to ensuring an immutable map).
      newdata <- check.factor(object$xvar, newdata, xfactor)

      # In classification we need to check y-outcome factor coherence. 
      if (object$family == "class") {
          if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
          newdata <- check.factor(object$yvar, newdata, yfactor)
        }
      }

      ## Extract test yvar names (if any) and xvar names.
      if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
        fnames <- c(object$yvar.names, object$xvar.names)
      }
      else {
        fnames <- object$xvar.names
      }

      ## Data conversion to numeric mode
      newdata <- finalizeData(fnames, newdata, na.action) 
                           
      ## Extract the test x-matrix and sort the columns.  This
      ## accomodates for an incoming test x-matrix that is in a
      ## different order than the training x-matrix.
      xvar.newdata  <- newdata[, is.element(names(newdata), object$xvar.names), drop = FALSE]
      xvar.newdata <- as.matrix(xvar.newdata[,
                           na.omit(match(colnames(xvar.newdata), object$xvar.names)), drop = FALSE])
      n.newdata <- nrow(newdata)

      ## Process the y-outcomes.  Set their dimension.
      if (sum(is.element(object$yvar.names, names(newdata))) == length(object$yvar.names)) {
        yvar.newdata <- as.matrix(newdata[, object$yvar.names, drop = FALSE])
        event.info.newdata <- get.grow.event.info(yvar.newdata, object$family, need.deaths = FALSE)
        r.dim.newdata <- event.info.newdata$r.dim
        perf.flag <- TRUE

        ## Event type consistency
        if (grepl("surv", object$family) &&
            length(setdiff(na.omit(event.info.newdata$cens), na.omit(event.info$cens))) > 1) {
          stop("survival events in test data do not match training data")
        }
        
      }       
      else {

        ##disallow outcome=TEST without y-outcomes
        if (outcome == "test") {
          stop("outcome=TEST, but the test data has no y values, which is not permitted")
        }

        ## there is no test data
        r.dim.newdata <- 0
        yvar.newdata <-  NULL
        perf.flag <- FALSE
        importance <- "none"
        
      }
                         
    }#end !partial.class

    ##partial class
    else {

      ## this is an internal mode, so we do not perform any checks

      ##extract the data
      newdata <- as.data.frame(data.matrix(newdata))
      xvar.newdata <- as.matrix(newdata)

      ## Initialize sample size.
      n.newdata <- nrow(xvar.newdata)

      ## there is no test data
      r.dim.newdata <- 0
      yvar.newdata <- NULL 
      perf.flag <- FALSE
      importance <- "none"


    }

    ## remove xvar row and column names for proper processing by the native library
    ## does not apply when outcome = TEST because the xvar TEST data has been made NULL
    if (outcome != "test") {
      rownames(xvar.newdata) <- colnames(xvar.newdata) <- NULL
    }

    ## we don't need the test data anymore
    remove(newdata)

  }
  
  else {

    ## ############ RESTORE MODE ONLY ############## ##
    ## We are processing the restore only mode

    ## there cannot be test data in restore mode
    n.newdata <- 0
    r.dim.newdata <- 0
    xvar.newdata <- NULL
    yvar.newdata <-  NULL

    ## Outcome is set to train for the native code
    ## Determine whether performance values are requested
    outcome <- "train"
    if (object$bootstrap != "by.root") {
      importance <- "none"
      perf.flag <- FALSE
    }
    else {
      perf.flag <- TRUE
    }
  }#ends restore.mode check

  
  ##############################################################
  #### we have completed the restore/non-restore mode processing
  ##############################################################
  
  ## final processing of xvar and yvar test data
  ##  depends on "outcome"

  # outcome=train
  if (outcome == "train") {
    
    ## data conversion for training data
    xvar <- as.matrix(data.matrix(object$xvar))    
    yvar <- as.matrix(data.matrix(object$yvar))
    
  }
  # outcome=test
  else {

    ##we have assumed non-restore mode but now we
    ##swap the training data out with the test data
    ##thus we convert outcome = test to restore mode
    ##performance is always requested for this setting
    xvar <- xvar.newdata
    yvar <- yvar.newdata
    restore.only <- TRUE 
    perf.flag <- TRUE
    
    ##pretend there is no test data, but do *not* get rid of it
    ##we need (and use) this data *after* the native code call
    n.newdata <- 0
    r.dim.newdata <- 0
    
  }

  ## Set the y dimension
  r.dim <- ncol(cbind(yvar))

  ## remove row and column names for proper processing by the native library
  ## set the dimensions
  rownames(xvar) <- colnames(xvar) <- NULL
  n.xvar <- ncol(xvar)
  n <- nrow(xvar)

  ## Initialize the number of trees in the forest.
  ntree <- object$ntree

  ## Initialize the bits
  importance.bits <- get.importance(importance)
  var.used.bits <- get.var.used(var.used)
  split.depth.bits <- get.split.depth(split.depth)
  membership.bits <-  get.membership(membership)  
  seed <- get.seed(seed)
  proximity.bits <- get.proximity(proximity)
  outcome.bits <- get.outcome(outcome)
  perf.bits <-  get.perf.bits(perf.flag)

  do.trace <- get.trace(do.trace)
  nativeOutput <- .Call("rfsrcPredict",
                        as.integer(do.trace),
                        as.integer(seed),
                        as.integer(restore.only.bits +
                                   importance.bits +
                                   object$bootstrap.bits + 
                                   proximity.bits +
                                   split.depth.bits +
                                   var.used.bits +
                                   outcome.bits +
                                   perf.bits +
                                   membership.bits +
                                   cr.bits),
                        as.integer(ntree),
                        as.integer(n),
                        as.integer(r.dim),
                        as.character(yvar.types),
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
                        as.integer(object$nodesize),
                        as.integer(object$seed),
                        as.integer(length(importance.xvar.idx)),
                        as.integer(importance.xvar.idx),
                        as.integer(get.rf.cores()))

  ## Check for error return condition in the native code.
  if(is.null(nativeOutput)) {
    stop("Error occurred in algorithm.  Please turn trace on for further analysis.")
  }

  ## Determine missingness according to mode (REST vs. PRED)
  if (restore.only) {
    n.miss <- get.nmiss(xvar, yvar)
  }
  else {
    n.miss <- get.nmiss(xvar.newdata, yvar.newdata)
  }

  ## Extract the imputed data if there was missingness
  if (n.miss > 0) {
    
    imputed.data <- matrix(nativeOutput$imputation, nrow = n.miss)
    nativeOutput$Imputation <- NULL
    imputed.indv <- imputed.data[, 1]
    imputed.data <- as.data.frame(imputed.data[, -1, drop = FALSE])

    if (r.dim.newdata > 0 | perf.flag) {
      colnames(imputed.data) <- c(object$yvar.names, object$xvar.names)
    }
    
    else {
      colnames(imputed.data) <- object$xvar.names      
    }
    
  }
  
  ## Post-process the test data
  ## For restore mode there is no test data *except* when outcome=TEST
  ## For partial class mode there is no y-outcome test values and
  ## anyways post-processing is not necessary
  if ( (!partial.class) & (!restore.only | outcome == "test") ) {

    ## Add column names to test xvar matrix. 
    xvar.newdata <- as.data.frame(xvar.newdata)
    colnames(xvar.newdata) <- object$xvar.names

    ## Map xvar factors back to original values
    xvar.newdata <- map.factor(xvar.newdata, xfactor)

    if (perf.flag) {

      ## Add column names to test response matrix
      yvar.newdata <- as.data.frame(yvar.newdata)
      colnames(yvar.newdata) <- object$yvar.names

      ## Map response factors back to original values
      yvar.newdata <- map.factor(yvar.newdata, yfactor)

    }
    
  }


  ## Map imputed data factors back to original values
  if (n.miss > 0) {

    imputed.data <- map.factor(imputed.data, xfactor)
    
    if (perf.flag) {
      imputed.data <- map.factor(imputed.data, yfactor)
    }
    
  }


  ## pretty names + dimensioning of various objects
  
  if (grepl("surv", object$family)) {
    ## Determine whether the grow forest is CR or right-censoring
    if (object$family == "surv-CR") {
      ## CR analysis
      ens.names <- list(NULL, NULL, paste("condCHF.", 1:(sexp.dim), sep = ""))
      cif.names <- list(NULL, NULL, paste("CIF.", 1:(sexp.dim), sep = ""))
      err.names <- list(paste("event.", 1:(sexp.dim), sep = ""), NULL)
      vimp.names <- list(paste("event.", 1:(sexp.dim), sep = ""),
                           if(vimp.joint) "joint" else importance.xvar)
      mortality.names <- list(NULL, paste("event.", 1:(sexp.dim), sep = ""))
    }
    else {
      ## right-censoring.
      survival.names <- ens.names <- list(NULL, NULL, NULL)
      mortality.names <- err.names <- list(NULL, NULL)
      vimp.names <- list(NULL, if (vimp.joint) "joint" else importance.xvar)
    }
  }
  
  else {
    ## all other families.
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


  ## proximity
  if (proximity) {
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
    proximity.out <- proximity.out/diag(proximity.out)
    nativeOutput$proximity <- NULL
  }
  else {
    proximity.out <- NULL
  }

  ## variable importance (VIMP)
  if ((importance != "none") & perf.flag) {
    VIMP <- atmatrix(nativeOutput$importance, (if (vimp.joint) 1 else length(importance.xvar)),
                      vimp.names, keep.names = TRUE)
    nativeOutput$importance <- NULL
  }
  else {
    VIMP <- NULL
  }

  
  ## split depth
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

  ## error Rate
  if (perf.flag) {
    ERR <- atmatrix(nativeOutput$performance, ntree, err.names)
    nativeOutput$performance <- NULL
  }
  else {
    ERR <- NULL
  }

  ## membership
  if (membership) {
    if (restore.only) {
      membership.n <- n
    }
    else {
      membership.n <- n.newdata
    }
    membership.out <- matrix(nativeOutput$nodeMembership, c(membership.n, ntree))
    inbag.out <- matrix(nativeOutput$bootMembership, c(membership.n, ntree))
    nativeOutput$nodeMembership <- NULL
    nativeOutput$bootMembership <- NULL
  }
  else {
    membership.out <- NULL
    inbag.out <- NULL
  }

  ## make the output object
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
    inbag = inbag.out,
    imputed.indv = (if (n.miss>0) imputed.indv else NULL),
    imputed.data = (if (n.miss>0) imputed.data else NULL),
    split.depth  = split.depth.out,
    err.rate = ERR,
    importance = VIMP
  )
  #memory management
  nativeOutput$leafCount <- NULL
  object.family <- object$family
  remove(object)
  remove(proximity.out)
  remove(membership.out)
  remove(inbag.out)
  if (n.miss > 0) remove(imputed.indv)
  if (n.miss > 0) remove(imputed.data)
  remove(split.depth.out)
  remove(ERR)
  remove(VIMP)


  ## family specific additions to the predict object
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
            dimnames=ens.names), 1)
      nativeOutput$fullEnsemble <- NULL
      clasOutput <- list(predicted = predicted)
      classResp <- (if (!is.null(predicted)) bayes.rule(predicted) else NULL)
      remove(predicted)
      
      predicted.oob = (if (!is.null(nativeOutput$oobEnsemble))
             adrop(array(nativeOutput$oobEnsemble, c(n, sexp.dim - 1, 1), dimnames=ens.names), 1) else NULL)
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


