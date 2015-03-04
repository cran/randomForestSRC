\name{impute.rfsrc}
\alias{impute.rfsrc}
\title{Impute Only Mode}
\description{
  Fast imputation mode.  A random forest is grown and used to impute
  missing data.  No ensemble estimates or error rates are calculated.
}
\usage{\method{impute}{rfsrc}(formula, data, ntree = 250, mtry = NULL,
  nodesize = 1, splitrule = NULL, nsplit = 1,
  na.action = c("na.impute", "na.random"),
  nimpute = 1, 
  xvar.wt = NULL, 
  do.trace = FALSE, ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.  Can be
      left unspecified if there are no outcomes or we don't care to
      distinguish between y-outcomes and x-variables in the imputation.}
    \item{data}{Data frame containing the data to be imputed.}
    \item{ntree}{Number of trees to grow.}
    \item{mtry}{Number of variables randomly sampled at each split.}
    \item{nodesize}{Minimum terminal node size.}
    \item{splitrule}{Splitting rule used to grow trees.}
    \item{nsplit}{Non-negative integer value used to specify random splitting.}
    \item{na.action}{Missing value action. See details below.}
    
    \item{nimpute}{Number of iterations of missing data algorithm.}
    
    
    \item{xvar.wt}{Weights for selecting variables for splitting on.}

    \item{do.trace}{Logical. Should trace output be enabled on each
      iteration? Default is \code{FALSE}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \enumerate{
  \item Grow a forest and use this to impute data.  All external
  calculations such as ensemble calculations, error rates, etc. are
  turned off.  Use this function if your only interest is imputing the
  data.

  \item Data is imputed using the default missing data algorithm,
  however users can select \option{na.action="na.random"} for a cruder,
  but faster imputation.  Unlike \option{na.action="na.impute"}, data is not imputed
  as the tree is grown, instead tree nodes are split using non-missing
  in-bag data.  Following the split to a node, data points with
  missing values on the variable used to split the node are randomly
  assigned to daughter nodes. 

  
  \item If no formula is specified, unsupervised splitting
  is implemented which treats the data as if there are no y-outcomes.
  
  
  

  \item Prior to imputation, the data is processed and records with all
  values missing are removed, as are variables having all missing values.

  \item If there is no missing data, either before or after
  processing of the data, the algorithm returns the processed
  data and no imputation is performed.
    
  \item All options are the same as \command{rfsrc} and the user should
  consult the \command{rfsrc} help file for details.
  }
}   
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlayed.
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
      Statist.}, 2:841-860.

    Stekhoven D.J. and Buhlmann P. (2012). MissFores--non-parametric
    missing value imputation for mixed-type data.
    \emph{Bioinformatics}, 28(1):112-118.
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## example of survival imputation
## ------------------------------------------------------------

#imputation using outcome splitting
data(pbc, package = "randomForestSRC")
pbc.d <- impute.rfsrc(Surv(days, status) ~ ., data = pbc, nsplit = 3)

#when no formula is given we default to unsupervised splitting
pbc2.d <- impute.rfsrc(data = pbc, nodesize = 1, nsplit = 10, nimpute = 5)

#random splitting can be reasonably good
pbc3.d <- impute.rfsrc(Surv(days, status) ~ ., data = pbc,
          splitrule = "random", nodesize = 1, nimpute = 5)

## ------------------------------------------------------------
## example of regression imputation
## ------------------------------------------------------------

air.d <- impute.rfsrc(Ozone ~ ., data = airquality, nimpute = 5)
air2.d <- impute.rfsrc(data = airquality, nimpute = 5, nodesize = 1)
air3.d <- impute.rfsrc(Ozone ~ ., data = airquality, nimpute = 5,
           splitrule = "random", nodesize = 1)


}
}
\keyword{missing data}

