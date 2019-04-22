sgreedy.rfsrc <- function(formula, data,
                          ntree = 500,
                          hdim = 5,
                          treesize = function(x){min(50, x * .25)},
                          tune = TRUE, lag = 8, strikeout = 5,
                          mtry = NULL,
                          nodesize = 1,
                          nsplit = 5,
                          bootstrap = "by.root",
                          sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x},
                          samptype = "swor",
                          samp = NULL,
                          ...)
{
  ## --------------------------------------------------------------
  ##   
  ##   preliminary processing
  ##
  ## --------------------------------------------------------------
  ## verify key options
  if (!is.function(treesize) && !is.numeric(treesize)) {
    stop("treesize must be a function or number specifying size of tree")
  }
  if (is.function(treesize)) {
    treesize <- treesize(nrow(data))
  }
  ##--------------------------------------------------------------
  ##
  ## extract additional options specified by user
  ## we lock this down to allowed types
  ##
  ##--------------------------------------------------------------
  ## list of forest parameters
  rfnames <- names(formals(rfsrc))
  ## add key hidden parameters
  rfnames <- c(rfnames, "rfq", "perf.type", "gk.quantile", "prob", "prob.epsilon", "vtry")
  rfnames <- c(rfnames, "lot", "empirical.risk")
  ## restrict to allowed values
  rfnames <- rfnames[rfnames != "ntree"     &
                     rfnames != "mtry"      &
                     rfnames != "nodesize"  &
                     rfnames != "nsplit"    &
                     rfnames != "bootstrap" &
                     rfnames != "sampsize"  &
                     rfnames != "samptype"  ]
  ## get the permissible hidden options
  dots <- list(...)
  dots <- dots[names(dots) %in% rfnames]
  ## set the splitrule to the default super greedy splitting unless user overrides
  if (is.null(dots$splitrule)) {
    dots$splitrule = "sg.regr"
  }
  ## add formula if present --> for future unsupervised sgreedy
  if (!missing(formula)) {
    dots$formula <- formula
  }
  ## set bootstrap accordingly if the user has provided their own sampling
  ## ntree and sampsize are handled in rfsrc
  if (!is.null(samp)) {
    bootstrap <- "by.user"
  }
  ## set the lot (Laterally Optimized Tree) parameters
  dots$lot <- get.lot(hdim, treesize, lag, strikeout)
  ## over-ride auto termination (tune=FALSE) we set lag=0 and safe strikeout.
  ## !!must be integer!!
  if (!tune) {
    dots$lot$lag = dots$lot$strikeout <- as.integer(0)
  }
  ##--------------------------------------------------------------
  ##
  ## make the grow call and return the object
  ##
  ##--------------------------------------------------------------
  ## TBD TBD currently hard-coded for regression only TBD TBD
  return(do.call("rfsrc",
                 c(list(data = data,
                 ntree = ntree,
                 mtry = mtry,
                 nodesize = nodesize,
                 nsplit = nsplit,
                 bootstrap = bootstrap,
                 sampsize = sampsize,
                 samptype = samptype,
                 samp = samp), dots)))
}
sgreedy <- sgreedy.rfsrc
