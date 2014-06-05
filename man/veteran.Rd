\name{veteran}
\docType{data}
\alias{veteran}
\title{Veteran's Administration Lung Cancer Trial}
\description{
    Randomized trial of two treatment regimens for lung cancer.
    This is a standard survival analysis data set.
}
\format{A data frame containing:
    \tabular{ll}{
        trt     \tab treatment: 1=standard 2=test\cr
        celltype\tab cell-type: 1=squamous,  2=smallcell,  3=adeno,  4=large\cr
        time    \tab survival time\cr
        status  \tab censoring status\cr
        karno   \tab Karnofsky performance score (100=good)\cr
        diagtime\tab months from diagnosis to randomisation\cr
	age     \tab age in years\cr
        prior   \tab prior therapy 0=no, 1=yes
    }
  }
\source{Kalbfleisch and Prentice, \emph{The Statistical
    Analysis of Failure Time Data.}}
\references{
    Kalbfleisch J. and Prentice R, (1980) \emph{The Statistical
    Analysis of Failure Time Data.}  New York: Wiley.
}
\examples{data(veteran, package = "randomForestSRC")}
\keyword{datasets}
