\name{rfsrcSyn}
\alias{rfsrcSyn}
\alias{rfsrcSyn.rfsrc}
\title{Synthetic Random Forests}
\description{
  Grows a synthetic random forest (RF) using RF machines as synthetic
  features.  Applies only to regression and classification settings.
}
\usage{\method{rfsrcSyn}{rfsrc}(formula, data, object, newdata,
  ntree = 1000,
  mtry = NULL,
  mtrySeq = NULL,
  nodesize = 5,
  nodesizeSeq = c(1:10,20,30,50,100),
  nsplit = 0,
  min.node = 3,
  use.org.features = TRUE,
  na.action = c("na.omit", "na.impute"),
  verbose = TRUE,
  ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fit.
       Must be specified unless \code{object} is given.}
  \item{data}{Data frame containing the y-outcome and x-variables in
      the model. Must be specified unless \code{object} is given.}
  \item{object}{An object of class \code{(rfsrc, synthetic)}.
    Not required when \code{formula} and \code{data} are supplied.}
  \item{newdata}{Test data used for prediction.  If missing, predicted
    values for the training data are returned.}
  \item{ntree}{Number of trees.}
  \item{mtry}{mtry value for synthetic forest.}
  \item{mtrySeq}{Sequence of mtry values used for fitting the
    collection of RF machines.  If \code{NULL}, set to the default
   value \code{p}/3.}
  \item{nodesize}{Nodesize value for the synthetic forest.}
  \item{nodesizeSeq}{Sequence of nodesize values used for the fitting the
    collection of RF machines.}
  \item{nsplit}{If non-zero, nsplit-randomized splitting is used which can
    significantly increase speed.}
  \item{min.node}{Minimum forest averaged number of nodes a RF machine
    must exceed in order to be used as a synthetic feature.}
  \item{use.org.features}{In addition to synthetic features, should
    the original features be used when fitting synthetic forests?}
  \item{na.action}{Missing value action. The default \code{na.omit}
    removes the entire record if even one of its entries is \code{NA}.
    The action \code{na.impute} pre-imputes the data using fast
    imputation via \command{impute.rfsrc}.}
  \item{verbose}{Set to \code{TRUE} for verbose output.}
  \item{...}{Further arguments to be passed to the \code{rfsrc}
    function used for fitting the synthetic forest.}
}
\details{
  A collection of random forests are fit using different nodesize
  values.  The predicted values from these machines are then used as
  synthetic features (called RF machines) to fit a synthetic random
  forest (the original features are also used when fitting the synthetic
  forest).  Applies only to regression and classification settings.

  Note that synthetic features are constructed using out-of-bag (OOB)
  data in order to avoid double dipping into training data.
  Nevertheless, the internal OOB error rate for the synthetic forest may
  be biased and thus it is recommended that cross-validation be used for
  determining test set performance.

  If \code{mtrySeq} is set, RF machines are constructed for each
  combination of nodesize and mtry values specified by
  \code{nodesizeSeq} \code{mtrySeq}.  However, a sequence of values for
  \code{mtrySeq} generall does not work as well as using a fixed value.
  Generally, performance gains are observed when one of the two
  sequences is fixed: \code{mtrySeq} is fixed and \code{nodesizeSeq} is
  varied, or \code{nodesizeSeq} is fixed and \code{mtrySeq} is varied.
  However, see the examples below.
}
\value{
     A list with the following components:
     \item{rfMachines}{RF machines used to construct the synthetic
       features.}
     \item{rfSyn}{The synthetic RF.}
     \item{synthetic}{List containing the synthetic features.}
     \item{predicted}{The synthetic forest predicted value.}
     \item{err.rate}{OOB error rate or test set error rate.  For regression,
       MSE values are reported.  For classification, the Brier score.}
     \item{opt.machine}{Optimal machine: RF machine with smallest
       OOB error rate}
}
\author{
    Hemant Ishwaran and Udaya B. Kogalur
}
\references{
  Ishwaran H. and Malley J.D. (2014).  Synthetic learning machines.
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{impute.rfsrc}}
}
\examples{
\donttest{
## ------------------------------------------------------------
## compare synthetic forests to regular forest (classification)
## ------------------------------------------------------------

## rfsrc and rfsrcSyn calls
if (library("mlbench", logical.return = TRUE)) {

  ## simulate the data 
  ring <- data.frame(mlbench.ringnorm(250, 20))

  ## classification forests
  ringRF <- rfsrc(classes ~., data = ring)

  ## synthetic forests:
  ## 1 = nodesize varied
  ## 2 = nodesize/mtry varied
  ringSyn1 <- rfsrcSyn(classes ~., data = ring)
  ringSyn2 <- rfsrcSyn(classes ~., data = ring, mtrySeq = c(1, 10, 20))

  ## test-set performance
  ring.test <- data.frame(mlbench.ringnorm(500, 20))
  print(predict(ringRF, newdata = ring.test))
  print(rfsrcSyn(object = ringSyn1, newdata = ring.test))
  print(rfsrcSyn(object = ringSyn2, newdata = ring.test))

}

## ------------------------------------------------------------
## compare synthetic forest to regular forest (regression)
## ------------------------------------------------------------

## simulate the data
n <- 250
ntest <- 1000
N <- n + ntest
d <- 50
std <- 0.1
x <- matrix(runif(N * d, -1, 1), ncol = d)
y <- 1 * (x[,1] + x[,4]^3 + x[,9] + sin(x[,12]*x[,18]) + rnorm(n, sd = std)>.38)
dat <- data.frame(x = x, y = y)
test <- (n+1):N

## regression forests
regF <- rfsrc(y ~ ., data = dat[-test, ], )
pred.regF <- predict(regF, dat[test, ])$predicted

## synthetic forests
## we pass both the training and testing data
## but this can be split into separate commands as in the
## previous classification example
synF1 <- rfsrcSyn(y ~ ., data = dat[-test, ],
  newdata = dat[test, ])
synF2 <- rfsrcSyn(y ~ ., data = dat[-test, ],
  newdata = dat[test, ], mtrySeq = c(1, 10, 20, 30, 40, 50))

## standardized MSE performance
mse <- c(mean((y[test] - pred.regF)^2, na.rm = TRUE),
         mean((y[test] - synF1$predicted)^2, na.rm = TRUE),
         mean((y[test] - synF2$predicted)^2, na.rm = TRUE)) / var(y[-test])
names(mse) <- c("forest", "synthetic1", "synthetic2")
print(mse)
}
}
\keyword{forest}
\keyword{predict}
