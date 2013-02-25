\name{plot.variable}
\alias{plot.variable}
\alias{plot.variable.rfsrc}
\title{Plot Marginal Effect of Variables}
\description{
  Plot the marginal effect of an x-variable on the class probability
  (classification), response (regression), mortality (survival), or
  the expected years lost (competing risk) from a RF-SRC analysis.  Users
  can select between marginal (unadjusted, but fast) and partial plots
  (adjusted, but slow).
}
\usage{\method{plot.variable}{rfsrc}(x, xvar.names, 
  surv.type = c("mort", "rel.freq", "surv", "years.lost", "cif", "chf"),
  time, which.outcome, partial = FALSE, plots.per.page = 4, granule = 5,
  sorted = TRUE, nvar, npts = 25, smooth.lines = FALSE, subset,
  ...)
}
\arguments{
  \item{x}{An object of class \code{(rfsrc, grow)} or
    \code{(rfsrc, predict)}.}
  \item{xvar.names}{Names of the x-variables to be used.} 
  \item{surv.type}{For survival families, specifies the predicted value.}
  \item{time}{For survival families, the value of time used for plotting
    survival quantities.  Default is the median follow up time.}
  \item{which.outcome}{For classification families, an integer or
    character value specifying the class to focus on (defaults to the
    first class).  For competing risk families, an integer value between
    1 and \code{J} indicating the event of interest, where \code{J} is
    the number of event types.  The default is to use the first event
    type.}
  \item{partial}{Should partial plots be used?}
  \item{plots.per.page}{Integer value controlling page layout.}
  \item{granule}{Integer value controlling whether a plot for a
    specific variable should be given as a boxplot or scatter 
    plot. Larger values coerce boxplots.}
  \item{sorted}{Should variables be sorted by importance values.}
  \item{nvar}{Number of variables to be plotted. Default is all.}
  \item{npts}{Maximum number of points used when generating partial
    plots for continuous variables.}
  \item{smooth.lines}{Use lowess to smooth partial plots.}
  \item{subset}{Vector indicating which rows of the x-variable matrix
    \code{x$xvar} to use. All rows are used if not specified.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
    The vertical axis displays the ensemble predicted value, while
    x-variables are plotted on the horizontal axis.

    \enumerate{
      \item For regression, the predicted response is used.

      \item For classification, it is the predicted class probability
      specified by \option{which.outcome}.

      \item For survival, the choices are:
      \itemize{
	\item  Mortality (\code{mort}).
	\item  Relative frequency of mortality (\code{rel.freq}).
	\item Predicted survival (\code{surv}), where the predicted
        survival is for the time point specified using
        \code{time} (the default is the median follow up time).
      }
      \item For competing risks, the choices are:
      \itemize{
	\item  The expected number of life years lost (\code{years.lost}).
	\item  The cumulative incidence function (\code{cif}).
	\item  The cumulative hazard function (\code{chf}).
      }
      In all three cases, the predicted value is for the event type
      specified by \option{which.outcome}.  For \code{cif} and
      \code{chf} the quantity is evaluated at the time point specified
      by \code{time}.
    }
    
    For partial plots use \option{partial=TRUE}.  Their interpretation are
    different than marginal plots.  The y-value for a variable \eqn{X},
    evaluated at \eqn{X=x}, is   
    \deqn{
      \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,o}),
    }
    where \eqn{x_{i,o}} represents the value for all other variables
    other than \eqn{X} for individual \eqn{i} and \eqn{\hat{f}} is the
    predicted value. Generating partial plots can be very slow.
    Choosing a small value for \code{npts} can speed up computational
    times as this restricts the number of distinct \eqn{x} values used
    in computing \eqn{\tilde{f}}.
    
    For continuous variables, red points are used to indicate partial
    values and dashed red lines indicate a smoothed error bar of +/- two
    standard errors.  Black dashed line are the partial values.  Set
    \option{smooth.lines=TRUE} for lowess smoothed lines.  For discrete
    variables, partial values are indicated using boxplots with whiskers
    extending out approximately two standard errors from the mean.
    Standard errors are meant only to be a guide and should be
    interpreted with caution.
    
    Partial plots can be slow.  Setting \option{npts} to a smaller number
    can help.
}
\author{
    Hemant Ishwaran  \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
  Friedman J.H. (2001). Greedy function approximation: a gradient
    boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.

    Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.
    
    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
      Statist.}, 2:841-860.

    Ishwaran H., Gerds, T.A. Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2012). Random survival forests for competing risks.
}
\seealso{
  \command{\link{rfsrc}},
  \command{\link{predict.rfsrc}}
}
\examples{
\dontrun{
### survival/CR examples

# survival
data(veteran, package = "randomForestSRC") 
v.obj <- rfsrc(Surv(time,status)~., veteran, nsplit = 10, ntree = 100)
plot.variable(v.obj, plots.per.page = 3)
plot.variable(v.obj, plots.per.page = 2, xvar.names = c("trt", "karno", "age"))
plot.variable(v.obj, surv.type = "surv", nvar = 1, time = 200)
plot.variable(v.obj, surv.type = "surv", partial = TRUE, smooth.lines = TRUE)
plot.variable(v.obj, surv.type = "rel.freq", partial = TRUE, nvar = 2)

# competing risks
data(follic, package = "randomForestSRC")
follic.obj <- rfsrc(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 100)
plot.variable(follic.obj, which.outcome = 2)

### regression examples

# airquality
airq.obj <- rfsrc(Ozone ~ ., data = airquality)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

# motor trend cars
mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)
plot.variable(mtcars.obj, partial = TRUE, smooth.lines = TRUE)

### classification example

# iris
iris.obj <- rfsrc(Species ~., data = iris)
plot.variable(iris.obj, partial = TRUE)

# motor trend cars: predict number of carburetors
mtcars2 <- mtcars
mtcars2$carb <- factor(mtcars2$carb,
   labels = paste("carb", sort(unique(mtcars$carb))))
mtcars2.obj <- rfsrc(carb ~ ., data = mtcars2)
plot.variable(mtcars2.obj, partial = TRUE)
}
}
\keyword{plot}
