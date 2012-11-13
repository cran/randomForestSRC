\name{plot.rfsrc}
\alias{plot.rfsrc}
\title{Plot Error Rate and Variable Importance from a RF-SRC analysis}
\description{
  Plot out-of-bag (OOB) error rates and variable importance (VIMP)
  from a RF-SRC analysis. This is the default plot method for the package.
}
\usage{\method{plot}{rfsrc}(x, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE,  ...)}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc, predict)}.}
  \item{plots.one.page}{Should plots be placed on one page?}
  \item{sorted}{Should variables be sorted by importance values?}
  \item{verbose}{Should VIMP be printed?}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Plot cumulative OOB error rates as a function of number of trees.
  Plot variable importance (VIMP) if available.
}
\author{
  Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
  Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
  Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

  Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
  \emph{Rnews}, 7(2):25-31.
}
\seealso{
  \command{\link{predict.rfsrc}},
  \command{\link{rfsrc}}
}
\examples{
\dontrun{
## classification example

iris.obj <- rfsrc(Species ~ ., data = iris)
plot(iris.obj)


## competing risk example
## use the pbc data from the survival package
## events are transplant (1) and death (2)

if (library("survival", logical.return = TRUE)) {
  data(pbc, package = "survival")
  pbc$id <- NULL
  plot(rfsrc(Surv(time, status) ~ ., pbc, nsplit = 10))
}
}}
\keyword{plot}
