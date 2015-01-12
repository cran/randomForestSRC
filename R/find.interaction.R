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


find.interaction.rfsrc <- function(
    object, 
    xvar.names,
    cause,
    importance = c("permute", "random", "permute.ensemble", "random.ensemble"),
    method = c("maxsubtree", "vimp"),
    sorted = TRUE,
    nvar = NULL, 
    nrep  = 1,
    subset,
    na.action = c("na.omit", "na.impute", "na.random"),
    seed = NULL,
    do.trace = FALSE,
    verbose = TRUE,
    ...)
{
  if (is.null(object)) stop("Object is empty!")
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2    &
      sum(inherits(object, c("rfsrc", "forest"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rfsrc, grow)' or '(rfsrc, forest)'.")
    if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) 
        stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    }
  if (missing(cause)) {
    cause <- 1
  }
  outcome.target <- get.outcome.target(object$family, outcome.target)
  if (object$family == "unsupv") {
    method <- "maxsubtree"
  }
  method <- match.arg(method,  c("maxsubtree", "vimp"))
  importance <- match.arg(importance, c("permute", "random", "permute.ensemble", "random.ensemble"))
  event.info <- get.event.info(object)
  n.event <- max(1, length(event.info$event.type))
  if (n.event > 1) {
    object.imp <- NULL
  }
  else {
    object.imp <- cbind(object$importance)[, cause, drop = FALSE]
  }
  cov.names <- object$xvar.names
  n.interact <- n.cov <- length(cov.names)
  if (!missing(xvar.names)) {
    if (sum(is.element(cov.names, xvar.names)) == 0) {
      stop("Variables do not match original analysis:", cov.names)
    }
    xvar.names <- unique(xvar.names[is.element(xvar.names, cov.names)])
    if (length(xvar.names) == 1)
      stop("Pairwise comparisons require more than one candidate variable.")
    cov.names <- xvar.names
    sorted <- FALSE
    n.interact <- length(xvar.names)
  }
  if (sorted) {
    if (!is.null(object.imp)) {
      o.r <- order(object.imp, decreasing = TRUE)
      cov.names <- cov.names[o.r]
    }
  }
  if (!missing(nvar)) {
    n.interact <- min(n.cov, max(round(nvar), 1))
  }
  if (n.interact == 1) {
    stop("Pairwise comparisons require more than one candidate variable.")
  }
  if (method == "vimp") {
    if (n.event > 1) {
      interact.imp.list <- vector("list", n.event)
      names(interact.imp.list) <- colnames(object$importance)
    }
    for (j in 1:n.event) {
      rownames.interact.imp <- interact.imp <- NULL
      target.dim <- ifelse(n.event > 1, j, 1)
      if (verbose && n.event > 1) {
        cat("--> event", j, "\n")
      }
      for (k in 1:(n.interact-1)) {
        n.joint.cov <- n.interact - k
        imp <- rep(0 , 1 + n.joint.cov)
        imp.joint <- rep(0, n.joint.cov)
        for (l in (k+1):n.interact) {
          if (verbose) {
            cat("Pairing",cov.names[k],"with",cov.names[l],"\n")
          }
          for (m in 1:nrep) {
            imp.indv.m <- c(cbind(vimp(object, cov.names[c(k,l)], outcome.target = outcome.target,
                                       importance = importance, joint = FALSE, 
                                       na.action = na.action, subset = subset, seed = seed, 
                                       do.trace = do.trace)$importance)[, target.dim])
            imp.joint.m <- vimp(object, cov.names[c(k,l)], outcome.target = outcome.target,
                                importance = importance, joint = TRUE,
                                na.action = na.action, subset = subset, seed = seed,
                                do.trace = do.trace)$importance[target.dim]
            imp[1] <- imp[1] + imp.indv.m[1]
            imp[l-k+1] <- imp[l-k+1] + imp.indv.m[2]
            imp.joint[l-k] <- imp.joint[l-k] + imp.joint.m
          }
        }
        imp[1] <- imp[1] / n.joint.cov
        imp <- imp/nrep
        imp.joint <- imp.joint/nrep
        interact.imp <- rbind(interact.imp,
                              cbind(imp[1], imp[-1], imp.joint, (imp[1] + imp)[-1], imp.joint - (imp[1] + imp)[-1]))
        rownames.interact.imp <- c(rownames.interact.imp,
                                   paste(cov.names[k],":",cov.names[(k+1):n.interact],
                                         sep=""))
      }
      colnames(interact.imp) <- c("Var 1", "Var 2","Paired","Additive","Difference")
      rownames(interact.imp) <- rownames.interact.imp
      if (n.event > 1) {
        interact.imp.list[[j]] <- interact.imp
      }
    }
    if (n.event > 1) {
      interact.imp <- interact.imp.list
    }
    if (verbose) {
      cat("\n")
      cat("                              Method: ", method,                       "\n", sep="")
      cat("                    No. of variables: ", n.cov,                        "\n", sep="")
      cat("           Variables sorted by VIMP?: ", sorted,                       "\n", sep="")
      cat("   No. of variables used for pairing: ", n.interact,                   "\n", sep="")
      cat("    Total no. of paired interactions: ", length(rownames.interact.imp),"\n", sep="")
      cat("            Monte Carlo replications: ", nrep,                         "\n", sep="")
      cat("    Type of noising up used for VIMP: ", importance,                   "\n", sep="")
      cat("\n")
      if (n.event == 1) print(round(interact.imp, 4)) else print(interact.imp)
    }
    invisible(interact.imp)
  }
  else {
    max.obj <- max.subtree(object, sub.order = TRUE, max.order = 1)
    sub.order <- max.obj$sub.order
    if (sorted) {
      o.r <- order(diag(sub.order), decreasing = FALSE)
      sub.order <- sub.order[o.r, o.r]
    }
    cov.pt <- is.element(colnames(sub.order), cov.names[1:n.interact])
    sub.order <- sub.order[cov.pt, cov.pt]
    if (verbose) {
      cat("\n")
      cat("                              Method: ", method,              "\n", sep="")
      cat("                    No. of variables: ", n.cov,               "\n", sep="")
      cat("  Variables sorted by minimal depth?: ", sorted,              "\n", sep="")
      cat("\n")
      print(round(sub.order, 2))
    }
    invisible(sub.order)
  }
}
find.interaction <- find.interaction.rfsrc
