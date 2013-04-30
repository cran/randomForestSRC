\name{impute.rfsrc}
\alias{impute.rfsrc}
\title{Impute Only Mode}
\description{
  Fast imputation mode.  A random forest is grown and used to impute
  missing data.  No ensemble estimates or error rates are calculated.
}
\usage{\method{impute}{rfsrc}(formula, data, ntree = 1000, mtry = NULL,
  nodesize = NULL, splitrule = NULL, nsplit = 0, nimpute = 1,
  xvar.wt = NULL, seed = NULL, do.trace = FALSE, ...)
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
    \item{nimpute}{Number of iterations of missing data algorithm.}
    \item{xvar.wt}{Weights for selecting variables for splitting on.}
    \item{seed}{Seed for random number generator.}
    \item{do.trace}{Should trace output be enabled?} 
    \item{...}{Further arguments passed to or from other methods.}
  }
\details{
  Grow a forest and use this to impute data.  All external calculations
  such as ensemble calculations, error rates, etc. are turned off.  Use
  this function if your only interest is imputing the data.

  All options are the same as \command{rfsrc} and the user should consult the
  help file for \command{rfsrc} for details.
}   
\value{
  Invisibly, the data frame containing the orginal data with imputed
  data overlayed.
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## ------------------------------------------------------------
## survival example
## ------------------------------------------------------------

#default split rule
data(pbc, package = "randomForestSRC")
pbc.d <- impute.rfsrc(Surv(days, status) ~ ., data = pbc, nsplit = 3)

#random splitting is fast and works well
pbc2.d <- impute.rfsrc(Surv(days, status) ~ ., data = pbc, splitrule = "random")
summary(pbc.d - pbc2.d)

## ------------------------------------------------------------
## regression.example
## ------------------------------------------------------------

air.d <- impute.rfsrc(Ozone ~ ., data = airquality, nimpute = 5)
air2.d <- impute.rfsrc(Ozone ~ ., data = airquality, nimpute = 5, splitrule = "random")

## ------------------------------------------------------------
## unsupervised example
## we impute without distinction between y-outcomes and x-variables
## we are invoking random splitting
## all variables are imputed
## ------------------------------------------------------------

data(pbc, package = "randomForestSRC")
pbcR.d <- impute.rfsrc(data = pbc, nimpute = 5)
airR.d <- impute.rfsrc(data = airquality, nimpute = 5)
}
}
\keyword{missing data}

