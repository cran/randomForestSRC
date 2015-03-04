\name{pbc}
\docType{data}
\alias{pbc}
\title{Primary Biliary Cirrhosis (PBC) Data}
\description{
  Data from the Mayo Clinic trial in primary biliary cirrhosis (PBC) of
  the liver conducted between 1974 and 1984.  A total of 424 PBC
  patients, referred to Mayo Clinic during that ten-year interval, met
  eligibility criteria for the randomized placebo controlled trial of
  the drug D-penicillamine.  The first 312 cases in the data set
  participated in the randomized trial and contain largely complete
  data.  
}
\format{A data frame containing:
    \tabular{ll}{
           days       \tab survival time in days \cr
           status     \tab censoring indicator \cr
           drugs      \tab 1=D-penicillamine, 2=placebo \cr
           age        \tab age in days \cr
           sex        \tab 0=female, 1=male \cr
           asictes    \tab presence of asictes, 0=no 1=yes \cr
           hepatom    \tab presence of hepatomegaly, 0=no 1=yes \cr
           spiders    \tab presence of spiders, 0=no 1=yes \cr
	   edema      \tab presence of edema (0, 0.5, 1) \cr
	   bili       \tab serum bilirubin in mg/dl \cr
	   chol       \tab serum cholesterol in mg/dl \cr
	   albumin    \tab albumin in gm/dl \cr
	   copper     \tab urine copper in ug/day \cr
	   alk        \tab alkaline phosphatase in U/liter \cr
           sgot       \tab SGOT in U/ml \cr
	   trig       \tab triglicerides in mg/dl \cr
           platelet   \tab platelets per cubic ml/1000 \cr
	   prothrombin\tab prothrombin time in seconds \cr
	   stage      \tab histologic stage of disease
      }
}
\source{Flemming and Harrington, 1991, Appendix D.1.}
\references{
    Flemming T.R and Harrington D.P., (1991) \emph{Counting Processes
    and Survival Analysis.} New York: Wiley.
}
\examples{
\dontrun{
data(pbc, package = "randomForestSRC")
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 3)
}}
\keyword{datasets}
