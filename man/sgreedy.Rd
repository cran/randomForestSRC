\name{sgreedy.rfsrc}
\alias{sgreedy.rfsrc}
\alias{sgreedy}
\alias{sgreedy.cart}
\title{Super Greedy Forests}

\description{
  Random forests using super greedy trees.
}

\usage{\method{sgreedy}{rfsrc}(formula, data,
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

## convenient interface for growing a super greedy tree
sgreedy.cart(formula, data, ntree = 1, bootstrap = "none", ...)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit.}

  \item{data}{Data frame containing the y-outcome and x-variables.} 

  \item{ntree}{Number of trees.}

  \item{hdim}{Depth (dimension) of the hypercube subtree.  When
    \code{hdim}=0 reverts to usual CART.  When \code{hdim}=1 reverts to
    CART but in place of recursive partitioning (RPART) using lateral
    optimized trees (LOT).  Otherwise uses LOT with subtree hypercubes
    grown to a depth of \code{hdim}>1.}

  \item{treesize}{Function specifying the size of the tree (number of
    tree splits) relative to the sample size. Can also be a number.}

  \item{tune}{Automatically determine the optimal tree size using
  out-of-bag empirical risk?  Uses a smoothed OOB empirical risk and
  stops the first time the slope is positive exactly \code{strikeout}
  number of times. }

  \item{lag}{Running average lag used to smooth OOB empirical risk.}
  
  \item{strikeout}{Strikeout used to determine optimal tree size.}
  
  \item{mtry}{Number of variables randomly selected for splitting a
  subtree node.}

  \item{nodesize}{Number of cases needed to split a subtree node.  When
    \code{nodesize}=1 the subtree depth is entirely controlled by
    \code{hdim}.}

  \item{nsplit}{Non-negative integer value specifying number of 
    random split points used to split a subtree node.}

  \item{bootstrap}{Bootstrap protocol used in growing a tree.}

  \item{sampsize}{Function specifying size of subsampled data. Can also be a number.}

  \item{samptype}{Type of bootstrap used.}

  \item{samp}{Bootstrap specification when \code{"by.user"} is used.}

  \item{...}{Further arguments to be passed to \code{\link{rfsrc}}.}
  
}

\details{
  
  While CART has been a central tenet of machine learning methods, its
  binary tree structure resulting from left versus right splitting
  severly hinders local adaptivity that ultimately affects prediction
  performance.  Super Greedy Trees (SGT's) dispenses with CART splitting
  and uses instead a higher level of abstraction of splitting that
  yields a richer partition structure for greater adaptivity.

  \code{sgreedy} calls \code{\link{rfsrc}} with tuning parameters set
  appropriately for computing super greedy trees and forests.  The use
  of lateral optimized trees (LOT) allows tree growing to be terminated
  on the fly when a tree becomes overly greedy.  Optimal \code{treesize}
  is determined using out-of-bag (OOB) empirical smoothed by a moving
  average with tree size determined when risk flattens and/or starts to
  increase.

}

\value{
  An object of class \code{(rfsrc, grow)}.
}

\author{
  Hemant Ishwaran and Udaya B. Kogalur
}

\references{

  Devroye, L., Gyorfi, L., and Lugosi, G. (1996).  \emph{A Probabilistic
  Theory of Pattern Recognition}, Volume 31. Springer.
    
  Ishwaran H. (2019).  Super greedy trees.
}

\seealso{
  \command{\link{rfsrc}}
}

\examples{
\donttest{
## ------------------------------------------------------------
## simple regression example
## ------------------------------------------------------------

n <- 1e4
x <- matrix(runif(10 * n), ncol = 10)
fx <- 20 * (x[, 1] - 0.5)^2 + 10 * x[, 2] 
y <- fx + rnorm(n)
simple <- data.frame(y = y / var(y), x)
simple.hp <- data.frame(simple, x.12 = x[, 1] + x[, 2])

## compare CART-LOT to SGT
o1 <- sgreedy.cart(y ~ ., simple, hdim = 1)
o2 <- sgreedy.cart(y ~ ., simple)
o3 <- sgreedy.cart(y ~ ., simple.hp)

## plot the two predicted surfaces side by side
plot.yhat  <- function(o, nd = 1, main = NULL) {
  x <- o$xvar[, 1]
  y <- o$xvar[, 2]
  plot(x, y, cex = o$predicted, pch = 16, col = 2,
       xlab = expression(x[1]), ylab = expression(x[2]))
}

par(mfrow = c(2,2))
plot.yhat(o1)
mtext("CART-LOT")
plot.yhat(o2)
mtext("SGT")
plot.yhat(o3)
mtext("SGT+hyperplane")
plot(x[, 1:2], cex = fx / var(y), pch = 16,
    xlab = expression(x[1]), ylab = expression(x[2]))
mtext("Truth")

## ------------------------------------------------------------
## "ice cream sandwich" regression example
## ------------------------------------------------------------

## set the graphical windows
par(mfrow = c(2,2))

## simulate the oblique decision boundary
std <- .1
n  <- 500
q  <- 10
p <- q + 2
x1 <- runif(n, -2, 2)
x2 <- runif(n, -2, 2)
noise <- matrix(runif(n * q, -2, 2), nrow = n)
fx <- 1 * (((x1 + x2) < 1) & ((x1 + x2) > -1))
y  <- fx + rnorm(n, mean = 0, sd = std)
ice <- data.frame(y = y, x1 = x1, x2 = x2, noise)


## standard random forests
o <- rfsrc(y ~ ., ice, ntree = 500, nodesize = 1)
print(o)

## super greedy random forests
so <- sgreedy(y ~ ., ice, mtry = p)
print(so)

## compute OOB empirical risk
ntree <- 100
so <- sgreedy(y ~ ., ice, mtry = p,
   empirical.risk = TRUE, ntree = ntree, block.size = 1)

## compute running average of OOB empirical risk
runavg <- function(x, lag = 8) {
  
  x <- c(na.omit(x))
  lag <- min(lag, length(x))
  cx <- c(0,cumsum(x))
  rx <- cx[2:lag] / (1:(lag-1))
  c(rx, (cx[(lag+1):length(cx)] - cx[1:(length(cx) - lag)]) / lag)

}
oob.risk <- apply(so$oob.empr.risk / var(so$yvar), 2, runavg)

## compare OOB empirical tree risk to OOB forest error
plot(c(1,max(so$leaf.count)), range(c(oob.risk)), type="n",
   xlab = "Tree size", ylab = "OOB empirical risk")
lo <- lapply(oob.risk, function(rsk){lines(rsk,col=grey(0.8))})
plot(1:ntree, so$err.rate / var(so$yvar), type = "s", xlab = "Trees", ylab = "OOB error")

## super greedy random forests with manual tree size
so2 <- sgreedy(y ~ ., ice, treesize = 10, tune = FALSE, mtry = p)
print(so2)

## compare single CART tree to SGT
o1 <- rfsrc.cart(y ~ ., ice, mtry = p)
o2 <- sgreedy.cart(y ~ ., ice, mtry = p)

plot.yhat(o1)
mtext("CART")
plot.yhat(o2)
mtext("SGT")



## ------------------------------------------------------------
## peak data regression example
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

## load the data
peak <- data.frame(mlbench.peak(500, 25))

## standard random forests
o <- rfsrc(y ~ ., peak, ntree = 500, nodesize = 1)
print(o)

## super greedy random forests
so <- sgreedy(y ~ ., peak, hdim = 20)
print(so)

## same as above but using mse splitting
## not nearly as good as super greedy splitting
so <- sgreedy(y ~ ., peak, treesize = 100, tune = FALSE, hdim = 20, splitrule = "mse")
print(so)

}

## ------------------------------------------------------------
## friedman 1 (illustrates synthetic super greedy)
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

## load the data
fr1 <- data.frame(mlbench.friedman1(500))

## standard random forests
o <- rfsrc(y ~ ., fr1, ntree = 500, nodesize = 1)
print(o)

## super greedy random forests
so <- sgreedy(y ~ ., fr1, hdim = 20)
print(so)

## super greedy synthetic forests  - slow but good
sso <- synthetic(y ~ ., fr1, hdim = c(20, 5), nsplit = c(5, 10))
print(sso)

}

## ------------------------------------------------------------
## friedman 1 (illustrates hyperplane and virtual twin splits)
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

## augment data with hyperplane variables
make.hp <- function(f, dta, subset = NULL) {
  ynm <- all.vars(f)[1]
  y <- dta[, ynm]
  dta <- data.matrix(dta[, colnames(dta) != ynm])
  if (is.null(subset)) {
    subset <- 1:ncol(dta)
  }
  nm <- c(ynm, colnames(dta))
  hp <- do.call(cbind, lapply(subset, function(j) {
    nm <<- c(nm, paste("hp.", j, ".", setdiff(subset, j), sep = ""))
    dta[, j] + dta[, setdiff(subset, j), drop = FALSE]
  }))
  dta <- data.frame(y = y, dta, hp = hp)
  colnames(dta) <- nm
  dta
}

## augment data with virtual twin variables
make.vt <- function(f, dta, subset = NULL) {
  ynm <- all.vars(f)[1]
  y <- dta[, ynm]
  dta <- data.matrix(dta[, colnames(dta) != ynm])
  if (is.null(subset)) {
    subset <- 1:ncol(dta)
  }
  nm <- c(ynm, colnames(dta))
  vt <- do.call(cbind, lapply(subset, function(j) {
    nm <<- c(nm, paste("vt.", j, ".", setdiff(subset, j), sep = ""))
    dta[, j] * dta[, setdiff(subset, j), drop = FALSE]    
  }))
  dta <- data.frame(y = y, dta, vt = vt)
  colnames(dta) <- nm
  dta
}

## load the data 
fr1 <- data.frame(mlbench.friedman1(500))

## regular super greedy call
so <- sgreedy(y ~ ., fr1)
print(so)

## super greedy call with hyper splits
fr1hp <- make.hp(y ~ ., fr1, 1:5)
sohp <- sgreedy(y ~ ., fr1hp, mtry = 30)
print(sohp)

## super greedy call with virtual twin splits
fr1vt <- make.vt(y ~ ., fr1, 1:5)
sovt <- sgreedy(y ~ ., fr1vt, mtry = 30)
print(sovt)

}


## ------------------------------------------------------------
## boston housing regression example
## ------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {

## load the data
data(BostonHousing)

## standard random forests
o <- rfsrc(medv ~ ., BostonHousing, ntree = 500, nodesize = 1)
print(o)

## super greedy random forests
## use a larger strikeout to encourage slightly larger trees
so <- sgreedy(medv ~ ., BostonHousing, strikeout = 10)
print(so)

}

}}
\keyword{forest}
