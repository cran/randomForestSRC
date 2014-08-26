####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.5.5
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


plot.variable.rfsrc <- function(
  x,
  xvar.names,
  which.outcome, 
  time,
  surv.type = c("mort", "rel.freq", "surv", "years.lost", "cif", "chf"),
  partial = FALSE,
  show.plots = TRUE,
  plots.per.page = 4,
  granule = 5,
  sorted = TRUE,
  nvar,
  npts = 25,
  smooth.lines = FALSE,
  subset,
  ...)
{
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 2)) == 2) {
    x <- x$rfSyn
  }
  object <- x
  remove(x)
  if (sum(inherits(object, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rfsrc", "predict"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(object, c("rfsrc", "plot.variable"), TRUE) == c(1,2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)', '(rfsrc, predict)' or '(rfsrc, plot.variable)'")
  }
  if (object$family == "unsupv") {
    stop("this function does not apply to unsupervised forests")
  }
  if (partial && is.null(object$forest)) {
    stop("forest is empty:  re-run rfsrc (grow) call with forest=TRUE")
  }
  xvar <- object$xvar
  if (!is.null(object$imputed.indv)) {
    xvar[object$imputed.indv, ] <- object$imputed.data[, object$xvar.names]
  }
  n <- nrow(xvar)
  if (!inherits(object, "plot.variable")) {
    if (missing(subset)) {
      subset <- 1:n
    }
    else {
      if (is.logical(subset)) subset <- which(subset)
      subset <- unique(subset[subset >= 1 & subset <= n])
      if (length(subset) == 0) {
        stop("'subset' not set properly")
      }
    }
    xvar <- xvar[subset,, drop = FALSE]
    n <- nrow(xvar)
    fmly <- object$family
    outcome.target <- get.outcome.target(object$family, outcome.target)
    if (grepl("surv", fmly)) {
      event.info <- get.event.info(object, subset)
      cens <- event.info$cens
      event.type <- event.info$event.type
      if (missing(time)) {
        time <- median(event.info$time.interest, na.rm = TRUE)
      }
      if (fmly == "surv-CR") {
        if (missing(which.outcome)) {
          which.outcome <- 1
        }
        else {
          if (which.outcome < 1 || which.outcome > max(event.type, na.rm = TRUE)) {
            stop("'which.outcome' is specified incorrectly")
          }
        }
        VIMP <- object$importance[, which.outcome]
        surv.type <- setdiff(surv.type, c("mort", "rel.freq", "surv"))
        pred.type <- match.arg(surv.type, c("years.lost", "cif", "chf"))
        ylabel <- switch(pred.type,
                         "years.lost" = paste("Years lost for event ", which.outcome),
                         "cif" = paste("CIF for event ", which.outcome, " (time=", time, ")", sep = ""),
                         "chf" = paste("CHF for event ", which.outcome, " (time=", time, ")", sep = ""))
      }
      else {
        which.outcome <- 1
        VIMP <- object$importance
        surv.type <- setdiff(surv.type, c("years.lost", "cif", "chf"))
        pred.type <- match.arg(surv.type, c("mort", "rel.freq", "surv"))
        ylabel <- switch(pred.type,
                         "mort"      = "mortality",
                         "rel.freq"  = "standardized mortality",
                         "surv"      = paste("predicted survival (time=", time, ")", sep = ""))
      }
    }
    else {
      event.info <- time <- NULL
      if (fmly == "class" || fmly == "class+" || (fmly ==  "mix+" && is.factor(object$yvar[, outcome.target]))) {
        object.yvar <- data.frame(object$yvar)[, outcome.target]
        if (missing(which.outcome)) {
          which.outcome <- 1
        }
        else if (is.character(which.outcome)) {
          which.outcome <- match(match.arg(which.outcome, levels(object.yvar)), levels(object.yvar))
        }
        else {
          if (which.outcome > length(levels(object.yvar)) | which.outcome < 1) {
            stop("which.outcome is specified incorrectly:", which.outcome)
          }
        }
        pred.type <- "prob"
        VIMP <- object$importance[, 1 + which.outcome]
        ylabel <- paste("probability", levels(object.yvar)[which.outcome])
        remove(object.yvar)
      }
      else {
        pred.type <- "y"
        which.outcome <- NULL
        VIMP <- object$importance
        ylabel <- expression(hat(y))
      }
    }
    if (missing(xvar.names)) {
      xvar.names <- object$xvar.names
    }
    else {
      xvar.names <- intersect(xvar.names, object$xvar.names)
      if (length(xvar.names) ==  0){
        stop("none of the x-variable supplied match available ones:\n", object$xvar.names)
      }
    }
    if (sorted & !is.null(VIMP)) {
      xvar.names <- xvar.names[rev(order(VIMP[xvar.names]))]
    }
    if (!missing(nvar)) {
      nvar <- max(round(nvar), 1)
      xvar.names <- xvar.names[1:min(length(xvar.names), nvar)]
    }
    nvar <- length(xvar.names)
    if (!partial) {
      yhat <- extract.pred(object, pred.type, subset, time, which.outcome)
    }
    else {
      class(object$forest) <- c("rfsrc", "partial", class(object)[3])
      if (npts < 1) npts <- 1 else npts <- round(npts)
      prtl <- lapply(1:nvar, function(k) {        
        x <- na.omit(object$xvar[, object$xvar.names == xvar.names[k]])
        if (is.factor(x)) x <- factor(x, exclude = NULL)          
        n.x <- length(unique(x))
        if (!is.factor(x) & n.x > npts) {
          x.uniq <- sort(unique(x))[unique(as.integer(seq(1, n.x, length = min(npts, n.x))))]
        }
        else {
          x.uniq <- sort(unique(x))
        }
        n.x <- length(x.uniq)
        yhat <- yhat.se <- NULL
        newdata.x <- xvar
        factor.x <- !(!is.factor(x) & (n.x > granule))
        for (l in 1:n.x) {        
          newdata.x[, object$xvar.names == xvar.names[k]] <- rep(x.uniq[l], n)
          pred.temp <- extract.pred(predict.rfsrc(object$forest, newdata.x, importance = "none"),
                                    pred.type, 1:n, time, which.outcome)
          mean.temp <- mean(pred.temp , na.rm = TRUE)
          if (!factor.x) {
            yhat <- c(yhat, mean.temp)
            if (fmly == "class") {
              yhat.se <- c(yhat.se, mean.temp * (1 - mean.temp) / sqrt(n))
            }
            else {
              yhat.se <- c(yhat.se, sd(pred.temp / sqrt(n) , na.rm = TRUE))
            }
          }
          else {
            pred.temp <- mean.temp + (pred.temp - mean.temp) / sqrt(n)
            yhat <- c(yhat, pred.temp)
          }
        }
        list(xvar.name = xvar.names[k], yhat = yhat, yhat.se = yhat.se, n.x = n.x, x.uniq = x.uniq, x = x)
      })
    }
    plots.per.page <- max(round(min(plots.per.page,nvar)), 1)
    granule <- max(round(granule), 1)
    plot.variable.obj <- list(family = fmly,
                    partial = partial,
                    event.info = event.info,
                    which.outcome = which.outcome,
                    ylabel = ylabel,
                    n = n,
                    xvar.names = xvar.names,
                    nvar = nvar, 
                    plots.per.page = plots.per.page,
                    granule = granule,
                    smooth.lines = smooth.lines)
    if (partial) {
      plot.variable.obj$pData <- prtl
    }
    else {
      plot.variable.obj$yhat <- yhat
      plot.variable.obj$xvar <- xvar
    }
    class(plot.variable.obj) <- c("rfsrc", "plot.variable", fmly)
  }
  else {
    plot.variable.obj <- object
    remove(object)
    fmly <- plot.variable.obj$family
    partial <- plot.variable.obj$partial
    event.info <- plot.variable.obj$event.info
    which.outcome <- plot.variable.obj$which.outcome
    ylabel <- plot.variable.obj$ylabel
    n <- plot.variable.obj$n
    xvar.names <- plot.variable.obj$xvar.names
    nvar <- plot.variable.obj$nvar 
    plots.per.page <- plot.variable.obj$plots.per.page
    granule <- plot.variable.obj$granule
    smooth.lines <- plot.variable.obj$smooth.lines
    if (partial) {
      prtl <- plot.variable.obj$pData
    }
    else {
      yhat <- plot.variable.obj$yhat
      xvar <- plot.variable.obj$xvar
    }
    if (!is.null(event.info)){
      cens <- event.info$cens
      event.type <- event.info$event.type
    }
  }
  if (show.plots) {
    old.par <- par(no.readonly = TRUE)
  }
  if (!partial && show.plots) {
    par(mfrow = c(min(plots.per.page, ceiling(nvar / plots.per.page)), plots.per.page))
    if (n > 500) cex.pt <- 0.5 else cex.pt <- 0.75
    for (k in 1:nvar) {
      x <- xvar[, which(colnames(xvar) == xvar.names[k])]
      x.uniq <- unique(x)
      n.x <- length(x.uniq)
      if (!is.factor(x) & n.x > granule) {
        plot(x,
             yhat,
             xlab = xvar.names[k],
             ylab = ylabel,
             type = "n", ...) 
        if (grepl("surv", fmly)) {
          points(x[cens == which.outcome], yhat[cens == which.outcome], pch = 16, col = 4, cex = cex.pt)
          points(x[cens == 0], yhat[cens == 0], pch = 16, cex = cex.pt)
        }
        lines(lowess(x[!is.na(x)], yhat[!is.na(x)]), col = 2, lwd=3)
      }
      else {
        if (is.factor(x)) x <- factor(x, exclude = NULL)          
        boxplot(yhat ~ x, na.action = "na.omit",
                xlab = xvar.names[k],
                ylab = ylabel,
                notch = TRUE,
                outline = FALSE,
                col = "bisque",
                names = rep("", n.x),
                xaxt = "n", ...)
        at.pretty <- unique(round(pretty(1:n.x, min(30, n.x))))
        at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
        axis(1,
             at = at.pretty,
             labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
             tick = TRUE)
      }
    }
  }
  if (partial && show.plots) {
    plots.per.page <- max(round(min(plots.per.page,nvar)), 1)
    granule <- max(round(granule),1)
    par(mfrow = c(min(plots.per.page, ceiling(nvar/plots.per.page)), plots.per.page))
    for (k in 1:nvar) {
      x <- prtl[[k]]$x
      if (is.factor(x)) x <- factor(x, exclude = NULL)          
      x.uniq <- prtl[[k]]$x.uniq
      n.x <- prtl[[k]]$n.x
      if (n.x > 25) cex.pt <- 0.5 else cex.pt <- 0.75
      yhat <- prtl[[k]]$yhat
      yhat.se <- prtl[[k]]$yhat.se
      factor.x <- !(!is.factor(x) & (n.x > granule))
      if (!factor.x) {
        plot(c(min(x), x.uniq, max(x), x.uniq, x.uniq),
             c(NA, yhat, NA, yhat + 2 * yhat.se, yhat - 2 * yhat.se),
             xlab = prtl[[k]]$xvar.name,
             ylab = ylabel,
             type = "n", ...)
        points(x.uniq, yhat, pch = 16, cex = cex.pt, col = 2)
        if (!is.na(yhat.se) && any(yhat.se > 0)) {
          if (smooth.lines) {
            lines(lowess(x.uniq, yhat + 2 * yhat.se), lty = 3, col = 2)
            lines(lowess(x.uniq, yhat - 2 * yhat.se), lty = 3, col = 2)
          }
          else {
            lines(x.uniq, yhat + 2 * yhat.se, lty = 3, col = 2)
            lines(x.uniq, yhat - 2 * yhat.se, lty = 3, col = 2)
          }
        }
        if (smooth.lines) {
          lines(lowess(x.uniq, yhat), lty = 2, lwd=2)
        }
        else {
          lines(x.uniq, yhat, lty = 2, lwd=2)
        }
        rug(x, ticksize=0.03)
      }
      else {
        y.se <- 0.005
        bxp.call <- boxplot(yhat ~ rep(x.uniq, rep(n, n.x)), range = 2, plot = FALSE)
        boxplot(yhat ~ rep(x.uniq, rep(n, n.x)),
                xlab = xvar.names[k],
                ylab = ylabel,
                notch = TRUE,
                outline = FALSE,
                range = 2,
                ylim = c(min(bxp.call$stats[1,], na.rm=TRUE) * ( 1 - y.se ),
                         max(bxp.call$stats[5,], na.rm=TRUE) * ( 1 + y.se )),
                col = "bisque",
                names = rep("",n.x),
                xaxt = "n", ...)
        at.pretty <- unique(round(pretty(1:n.x, min(30,n.x))))
        at.pretty <- at.pretty[at.pretty >= 1 & at.pretty <= n.x]
        axis(1,
             at = at.pretty,
             labels = format(sort(x.uniq)[at.pretty], trim = TRUE, digits = 4),
             tick = TRUE)
      }
    }
  }
  if (show.plots) {
    par(old.par)
  }
  invisible(plot.variable.obj)
}
plot.variable <- plot.variable.rfsrc
