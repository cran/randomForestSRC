plot.subsample.rfsrc <- function(x, alpha = .01, xvar.names,
                          standardize = TRUE, normal = TRUE, jknife = FALSE,
                          target, m.target = NULL, pmax = 75, main = "",
                          sorted = TRUE, show.plots = TRUE,
                          ...)
{
  ##--------------------------------------------------------------
  ##
  ## was subsampling or double-bootstrap used?
  ##
  ##--------------------------------------------------------------
  if (sum(c(grepl("rfsrc", class(x))), grepl("subsample", class(x))) == 2) {
    subsample <- TRUE
  }
  else if (sum(c(grepl("rfsrc", class(x))), grepl("bootsample", class(x))) == 2) {
    subsample <- FALSE
  }
  else {
    stop("object must be obtained from call to 'subsample' function")
  }
  ##--------------------------------------------------------------
  ##
  ## coerce the (potentially) multivariate rf object
  ##
  ##--------------------------------------------------------------
  m.target <- get.univariate.target(x$rf, m.target)
  x$rf <- coerce.multivariate(x$rf, m.target)
  ##--------------------------------------------------------------
  ##
  ## family specific details
  ## - set the target if not supplied
  ## - assign pretty labels for the horizontal axis
  ##
  ##--------------------------------------------------------------
  fmly <- x$rf$family
  ## labels
  if (standardize) {
    if (fmly == "regr") {
      xlab <- "standardized vimp"
    }
    else {
      xlab <- "vimp"
    }
  }
  else {
    xlab <- "vimp"
  }
  ## over-ride x label when only generalization error is requested
  if (x$performance.only) {
    xlab <- "generalization error"
  }
  ## extract vimp and vimp column names - depends on family and univariate/multivariate
  ## special mod needed for generalization error
  if (is.null(m.target)) {  
    vmp <- x$vmp[[1]]
    if (!x$performance.only) {
      vmp.col.names <- colnames(x$vmp[[1]])
    }
    else {
      vmp.col.names <- names(x$vmp)
    }
  }
  else {
    vmp <- x$vmp[[m.target]]
    if (!x$performance.only) {
      vmp.col.names <- colnames(x$vmp[[m.target]])
    }
    else {
      vmp.col.names <- m.target
    }
  }
  if (fmly == "regr" || fmly == "surv") {
    target <- 0
    vmp <- c(vmp)
    xlab <- paste(xlab, " (", vmp.col.names, ")", sep = "")
  }
  else if (fmly == "class") {
    if (missing(target)) {
      target <- 0
    }
    else {
      yvar.levels <- levels(x$rf$yvar)
      if (is.character(target)) {
        target <- match(match.arg(target, yvar.levels), yvar.levels)
      }
      else {
        if (target < 0 || target > length(yvar.levels)) {
          stop("'target' is specified incorrectly")
        }
      }
    }
    vmp <- c(vmp[, 1 + target])
    xlab <- paste(xlab, " (", vmp.col.names[1 + target], ")", sep = "")
  }
  else if (fmly == "surv-CR") {    
    if (missing(target)) {
      target <- 0
    }
    else {
      n.event <- length(get.event.info(x$rf)$event.type)
      if (target < 1 || target > n.event) {
        stop("'target' is specified incorrectly")
      }
      target <- target - 1
    }
    vmp <- c(vmp[, 1 + target])
    xlab <- paste(xlab, " (", vmp.col.names[1 + target], ")", sep = "")
  }
  ##--------------------------------------------------------------
  ##
  ## over-ride x label if the user has supplied their own
  ##
  ##--------------------------------------------------------------
  if (!is.null(list(...)$xlab)) {
    xlab <- list(...)$xlab
  }
  ##--------------------------------------------------------------
  ##
  ## extract necessary objects
  ##
  ##--------------------------------------------------------------
  if (subsample) {
    oo <- extract.subsample(x, alpha = alpha, target = target, m.target = m.target,
                            standardize = standardize, raw = TRUE)
    boxplot.dta <- oo$boxplot.dta
  }
  else {
    oo <- extract.bootsample(x, alpha = alpha, target = target, m.target = m.target,
                             standardize = standardize, raw = TRUE)
    boxplot.dta <- oo[[1]]
  }
  if (normal) {
    if (subsample && jknife) {
      ci <- oo$ci.jk.Z
    }
    else {
      ci <- oo$ci.Z
    }
  }
  else {
    ci <- oo$ci
  }
  ##--------------------------------------------------------------
  ##
  ## final processing
  ##
  ##--------------------------------------------------------------
  ## (a) user has provided a list of target variables
  if (!missing(xvar.names)) {
    trim.pt <- colnames(boxplot.dta) %in% xvar.names
    if (sum(trim.pt) > 0) {
      boxplot.dta <- boxplot.dta[, trim.pt, drop = FALSE]
      ci <- ci[, trim.pt, drop = FALSE]
      vmp <- vmp[trim.pt]
    }
    else {
      stop("user supplied variable names do not match any of the variables from analysis")
    }
  }
  ## (b) sort by vimp?
  p <- ncol(boxplot.dta)
  if (sorted || p > pmax) {
    o.pt <- order(vmp, decreasing = TRUE)
  }
  else {
    o.pt <- 1:p
  }
  ## (c) trim if too many variables
  if (p > pmax) {
    o.pt <- o.pt[1:pmax]
  }
  ## redefine the data to match the processing
  boxplot.dta <- boxplot.dta[, o.pt, drop = FALSE]
  ci <- ci[, o.pt, drop = FALSE]
  ##--------------------------------------------------------------
  ##
  ## skeleton boxplot 
  ##
  ##--------------------------------------------------------------
  bp <- boxplot(boxplot.dta,
                yaxt="n",
                outline = FALSE,
                horizontal = TRUE,
                plot = FALSE)
  bp$stats <- ci
  col.pt <- bp$stats[1, ] > 0 
  col.pt[is.na(col.pt)] <- FALSE
  colr <- c("blue", "red")[1 + col.pt]
  ##--------------------------------------------------------------
  ##
  ## finesse ... unnamed options to be passed to bxp and axis
  ##
  ##--------------------------------------------------------------
  ## pull the unnamed options
  dots <- list(...)
  bxp.names <- c(names(formals(bxp)),
                 "xaxt", "yaxt", "las", "cex.axis", 
                 "col.axis", "cex.main",
                 "col.main", "sub", "cex.sub", "col.sub", 
                 "ylab", "cex.lab", "col.lab")
  axis.names <- formals(axis)
  axis.names$tick <- axis.names$las <- axis.names$labels <- NULL
  axis.names <- c(names(axis.names), "cex.axis") 
  ## override xlab
  if (!is.null(dots$xlab)) {
    xlab <- dots$xlab
  }
  ## overlay ylab when user mistakenly uses it  
  if (!is.null(dots$ylab) && length(dots$ylab) == p) {
    bp$names <- dots$ylab[o.pt]
    bxp.names <- bxp.names[bxp.names != "ylab"]
  }
  ## overlay names
  if (!is.null(dots$names)) {
    bp$names <- dots$names[o.pt]
  }
  ## special mod needed for generalization error
  if (x$performance.only) {
    bp$names <- ""
  }
  ##--------------------------------------------------------------
  ##
  ## draw the core bxp plot
  ##
  ##--------------------------------------------------------------
  if (show.plots) {
    do.call("bxp", c(list(z = bp, main = main, xlab = xlab, 
                          boxfill = colr, xaxt = "n", yaxt = "n",
                          outline = FALSE, horizontal = TRUE),
                     dots[names(dots) %in% bxp.names]))
    do.call("axis", c(list(side = 1, at = pretty(c(bp$stats)), tick = .02), dots[names(dots) %in% axis.names]))
    do.call("axis", c(list(side = 2, at = 1:length(bp$names), labels = bp$names,
                           las = 2, tick = FALSE), dots[names(dots) %in% axis.names]))
    abline(h = 1:length(bp$names), col = gray(.9), lty = 1)
    abline(v = 0, lty = 1, lwd = 1.5, col = gray(.8))
    bxp(bp, boxfill=colr,xaxt="n",yaxt="n",
        outline=FALSE,horizontal=TRUE,add=TRUE,
        whisklty=1,whisklwd=2)
  }
  ##--------------------------------------------------------------
  ##
  ## return the invisible boxplot data
  ##
  ##--------------------------------------------------------------
  invisible(bp)
}
plot.subsample <- plot.subsample.rfsrc
