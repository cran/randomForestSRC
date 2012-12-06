\name{predict.rfsrc}
\alias{predict.rfsrc}
\title{Prediction for Random Forests for Survival, Regression, and Classification}
\description{
    Prediction on test data using a fitted RF-SRC object.
}
\usage{\method{predict}{rfsrc}(object, newdata,
  importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
  na.action = c("na.omit", "na.impute"), outcome = c("train", "test"),
  proximity = FALSE, var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"), seed = NULL,
  do.trace = FALSE, membership = TRUE, ...)
}
\arguments{
    \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	forest)}. Requires \option{forest=TRUE} in the original \code{rfsrc}
        call.}
    \item{newdata}{Test data.  If missing, the original grow
      (training) data is used.}    
    \item{importance}{Method for computing variable importance (VIMP).
      See \command{rfsrc} for details. Only applies when the test data
      contains y-outcome values.}
    \item{na.action}{Missing value action. The default
         \code{na.omit} removes the entire record if
        even one of its entries is \code{NA}.  Use \option{na.impute} to
        impute the test data.}
    \item{outcome}{Determines which data should be used to calculate the
        grow predictor.  The default and natural choice \code{train}
	uses the training data, but see details below.}
    \item{proximity}{Should proximity measure between
        test observations be calculated? Can be large.}
    \item{var.used}{Record the number of times a variable is split?}
    \item{split.depth}{Return minimal depth for each variable for each case?}
    \item{seed}{Negative integer specifying seed for the random number
      generator.} 
    \item{do.trace}{Should trace output be enabled? 
        Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{membership}{Should terminal node membership and inbag
      information be returned?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  Predicted values are obtained by dropping test data down the grow
  forest (the forest grown using the the training data).  If the test
  data contains y-outcome values, the overall error rate and VIMP are
  also returned.  Single as well as joint VIMP measures can be
  requested.  Note that calculating VIMP can be computationally
  expensive, especially when the dimension is high.  If VIMP is not
  needed, computational times can be significantly improved by setting
  \option{importance="none"} which turns VIMP calculations off entirely.
  
  Setting \option{na.action="na.impute"} imputes missing
  test data (x-variables and/or y-outcomes).  Imputation uses the
  grow-forest such that only training data is used when imputing test
  data to avoid biasing error rates and VIMP (Ishwaran et al. 2008).
   
  If \option{outcome="test"}, the predictor is calculated by using
  y-outcomes from the test data (outcome information must be present).
  In this case, the terminal nodes from the grow-forest are recalculated
  using the y-outcomes from the test set.  This yields a modified
  predictor in which the topology of the forest is based solely on the
  training data, but where the predicted value is based on the test
  data.  Error rates and VIMP are calculated by bootstrapping the test
  data and using out-of-bagging to ensure unbiased estimates.  See the
  examples for illustration.
}    
\value{
     An object of class \code{(rfsrc, predict)}, which is a list with the
     following components:
     \item{call}{The original grow call to \command{rfsrc}.}
     \item{family}{The family used in the analysis.}
     \item{n}{Sample size of test data (depends upon \code{NA} values).}
     \item{ntree}{Number of trees in the grow forest.}
     \item{yvar}{Test set y-outcomes or original grow y-outcomes if none.}
     \item{yvar.names}{A character vector of the y-outcome names.}
     \item{xvar}{Data frame of test set x-variables.}
     \item{xvar.names}{A character vector of the x-variable names.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
       grow forest.  Vector of length \code{ntree}.}
     \item{forest}{The grow forest.}
     \item{proximity}{Symmetric proximity matrix of the test data.}
     \item{membership}{Matrix of dimension \code{n}x\code{ntree}
       recording terminal node membership for the
       test data where each column contains the node number that a
       case falls in for that tree.}
     \item{inbag}{Matrix of dimension \code{n}x\code{ntree} recording
       inbag membership for the test data where each column contains the
       number of times that a case appears in the bootstrap sample for
       that tree.}
     \item{imputed.indv}{Vector of indices of records in test data with
         missing values.}
     \item{imputed.data}{Data frame comprising imputed test data.  First
       columns are the y-outcomes.}
     \item{split.depth}{Matrix of size \code{n}x\code{p} where entry
            [i][j] is the minimal depth for variable [j] for test case [i].}
     \item{err.rate}{Cumulative OOB error rate for the test data if
       y-outcomes are present.}
     \item{importance}{Test set variable importance (VIMP).  Can be \code{NULL}.}
     \item{predicted}{Test set predicted value.}
     \item{predicted.oob}{OOB predicted value (\code{NULL} unless
       \option{outcome="test"}).}\cr
     
\item{...... class}{for classification settings, additionally the following ......} \cr     
      \item{class}{In-bag predicted class labels.}
      \item{class.oob}{OOB predicted class labels (\code{NULL} unless \option{outcome="test"}).}\cr

\item{...... surv}{for survival settings, additionally the following ......} \cr
     
     \item{chf}{Cumulative hazard function (CHF).}
     \item{chf.oob}{OOB CHF (\code{NULL} unless \option{outcome="test"}).}
     \item{survival}{Survival function.}
     \item{survival.oob}{OOB survival function (\code{NULL} unless \option{outcome="test"}).}
     \item{time.interest}{Ordered unique death times.}
     \item{ndead}{Number of deaths.}\cr

\item{...... surv-CR}{for competing risks, additionally the following ......} \cr

     \item{chf}{Cause-specific cumulative hazard function (CSCHF)
     for each event.}
     \item{chf.oob}{OOB CSCHF for each event (\code{NULL} unless \option{outcome="test"}).}
     \item{cif}{Cumulative incidence function (CIF) for each event.}
     \item{cif.oob}{OOB CIF for each event (\code{NULL} unless \option{outcome="test"}).}
     \item{time.interest}{Ordered unique event times.}
     \item{ndead}{Number of events.}
}
\note{
  The dimensions and values of returned objects depend heavily on the
  underlying family and whether y-outcomes are present in the test data.
  In particular, items related to performance will be \code{NULL} when
  y-outcomes are not present.

  For detailed definitions of returned values (such as \code{predicted})
  see the help file for \command{rfsrc}.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
  (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
  
  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}	
\seealso{
  \command{\link{rfsrc}}
}
\examples{
### typical train/testing scenario

data(veteran, package = "randomForestSRC")
train <- sample(1:nrow(veteran), round(nrow(veteran) * 0.80))
veteran.grow <- rfsrc(Surv(time, status) ~ ., veteran[train, ], ntree = 100) 
veteran.pred <- predict(veteran.grow, veteran[-train , ])
print(veteran.grow)
print(veteran.pred)


\dontrun{
### unique feature of rfsrc
### cross-validation can be used when factor labels differ over
### training and test data

# first we convert all x-variables to factors
veteran.factor <- data.frame(lapply(veteran, factor))
veteran.factor$time <- veteran$time
veteran.factor$status <- veteran$status

# split the data into unbalanced train/test data (5/95)
# the train/test data have the same levels, but different labels
train <- sample(1:nrow(veteran), round(nrow(veteran) * .05))
summary(veteran.factor[train,])
summary(veteran.factor[-train,])

# grow the forest on the training data and predict on the test data
veteran.f.grow <- rfsrc(Surv(time, status) ~ ., veteran.factor[train, ]) 
veteran.f.pred <- predict(veteran.f.grow, veteran.factor[-train , ])
print(veteran.f.grow)
print(veteran.f.pred)


### example illustrating the flexibility of outcome = "test"
### use this to get the out-of-bag error rate

# first, we make the grow call
data(pbc, package = "randomForestSRC")
pbc.grow <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)

# now use predict with outcome = TEST
pbc.pred <- predict(pbc.grow, pbc, outcome = "test")

# notice that error rates are the same!!
print(pbc.grow)
print(pbc.pred)

# similar example, but with na.action = "na.impute"
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
print(airq.obj)
print(predict(airq.obj, outcome = "test", na.action = "na.impute"))

# similar example for classification
iris.obj <- rfsrc(Species ~., data = iris)
print(iris.obj)
print(predict.rfsrc(iris.obj, outcome = "test"))


### another example illustrating outcome = "test"
### unique way to check reproducibility of the forest

# primary call
set.seed(542899)
data(pbc, package = "randomForestSRC")
train <- sample(1:nrow(pbc), round(nrow(pbc) * 0.50))
pbc.out <- rfsrc(Surv(days, status) ~ .,  data=pbc[train, ],
        nsplit = 10)

# standard predict call
pbc.train <- predict(pbc.out, pbc[-train, ], outcome = "train")
#non-standard predict call: overlays the test data on the grow forest
pbc.test <- predict(pbc.out, pbc[-train, ], outcome = "test")

# check forest reproducibilility by comparing "test" predicted survival
# curves to "train" predicted survival curves for the first 3 individuals
Time <- pbc.out$time.interest
matplot(Time, t(exp(-pbc.train$chf)[1:3,]), ylab = "Survival", col = 1, type = "l")
matlines(Time, t(exp(-pbc.test$chf)[1:3,]), col = 2)
}
}
\keyword{predict}
\keyword{forest}

