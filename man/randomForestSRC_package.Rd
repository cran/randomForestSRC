\name{randomForestSRC-package}
\alias{randomForestSRC-package}
\docType{package}
\title{
  Random Forests for Survival, Regression, and Classification (RF-SRC)
}
\description{
  This package provides a unified treatment of Breiman's random forests
  (Breiman 2001) for survival, regression, and classification problems.
  The underlying code is based on Ishwaran and Kogalur's now retired 
  \pkg{randomSurvivalForest} package (Ishwaran and Kogalur 2007), and
  has been significantly refactored for improved computational speed.  }
\details{
  Implements Breiman's random forests for a variety of data settings.
  Numeric or categorical (factor) responses yield regression and
  classification forests.  Survival and competing risk forests are grown
  when the response is right-censored (Ishwaran et al. 2008, 2012).
  Different splitting rules invoked under deterministic or random
  splitting are available for all families.  Variable predictiveness can
  be assessed using variable importance (VIMP) measures for single, as
  well as grouped variables.  Variable selection is implemented using
  minimal depth variable selection (Ishwaran et al. 2010). Missing data
  (for x-variables and y-outcomes) can be imputed on both training and
  test data.

  This package implements OpenMP shared-memory parallel programming.
  However, the default installation will only execute serially.  To
  utilize OpenMP, the target system must first support it.  To install
  the package with OpenMP compiler options turned on:
  
  \enumerate{
    \item Download the source code for the package.

    \item From the root directory of the package source run the command
    \option{autoconf}.

    \item Run \option{R CMD INSTALL randomForestSRC} on the modified
    package directory.
  }
  
  This package contains many useful functions and users should read the
  help file in its entirety for details.  However, we briefly mention
  several key functions that may make it easier to navigate and
  understand the layout of the package.

  \enumerate{
    \item \command{\link{rfsrc}}
    
    This is the main entry point to the package.  It grows a random forest
    using user supplied training data.  We refer to the resulting object
    as a RF-SRC grow object.  Formally, the resulting object has class
    \code{(rfsrc, grow)}.

    \item \command{\link{predict.rfsrc}} (\command{predict})
    
    Used for prediction.  Predicted values are obtained by dropping the
    user supplied test data down the grow forest.  The resulting object
    has class \code{(rfsrc, predict)}.


    \item \command{\link{max.subtree}}, \command{\link{var.select}}
    
    Used for variable selection.  The function \command{max.subtree}
    extracts maximal subtree information from a RF-SRC object which is
    used for selecting variables by making use of minimal depth variable
    selection.  The function \command{var.select} provides
    an extensive set of variable selection options and is a wrapper to
    \command{max.subtree}.

    \item \command{\link{impute.rfsrc}}
  
    Fast imputation mode for RF-SRC.  Both \command{rfsrc} and
    \command{predict.rfsrc} are capable of imputing missing data.
    However, for users whose only interest is imputing data, this function
    provides an efficient and fast interface for doing so.
  }
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
    
    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

    Ishwaran H., Gerds, T.A. Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2012). Random survival forests for competing risks.
}
\seealso{
  \command{\link{find.interaction}},
  \command{\link{impute.rfsrc}},
  \command{\link{max.subtree}},
  \command{\link{plot.competing.risk}},
  \command{\link{plot.rfsrc}},
  \command{\link{plot.survival}},
  \command{\link{plot.variable}},
  \command{\link{predict.rfsrc}},
  \command{\link{print.rfsrc}},
  \command{\link{rfsrc}},
  \command{\link{restore}},
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\keyword{package}

  
