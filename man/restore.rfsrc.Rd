\name{restore}
\alias{restore}
\alias{restore.rfsrc}
\title{Restoration for Random Forests for Survival, Regression, and Classification}
\description{  
  Restoration mode for RF-SRC.
}
\usage{\method{restore}{rfsrc}(object,
  importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
  proximity = FALSE, var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"), seed = NULL, do.trace = FALSE,
  membership = FALSE, ...)
}
\arguments{
    \item{object}{An object of class \code{(rfsrc, grow)} or \code{(rfsrc,
	forest)}. Requires \option{forest=TRUE} in the original \code{rfsrc}
        call.}
    \item{importance}{Method used to compute variable importance.}
      \item{proximity}{Should the proximity between observations be
        calculated?}
    \item{var.used}{Analyzes which variables are split on.}
    \item{split.depth}{Return minimal depth for each variable for each case.}
    \item{seed}{Negative integer specifying the random number
    generator seed.}
    \item{do.trace}{Should trace output be enabled? Default is
      \code{FALSE}. Integer values can also be passed.  A positive value
      causes output to be printed each \code{do.trace} iteration.}
    \item{membership}{Should terminal node membership and inbag
      information be returned?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \command{restore} takes the forest grow object and gives the user the
    ability to restore all grow outputs.
}   
\value{
     An object of class \code{(rfsrc, predict)}.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\seealso{
  \command{\link{rfsrc}}
}
\examples{
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute", ntree=100)
restore(airq.obj)
print(airq.obj)
}
\keyword{predict}
\keyword{forest}
