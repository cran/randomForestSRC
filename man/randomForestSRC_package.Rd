\name{randomForestSRC-package}
\alias{randomForestSRC-package}
\docType{package}
\title{
  Random Forests for Survival, Regression and Classification (RF-SRC)
}
\description{
  This package provides a unified treatment of Breiman's random forests
  (Breiman 2001) for a variety of data settings.  Regression and
  classification forests are grown when the response is numeric or
  categorical (factor) while survival and competing risk forests
  (Ishwaran et al. 2008, 2012) are grown for right-censored survival
  data.  Different splitting rules invoked under deterministic or random
  splitting are available for all families.  Variable predictiveness can
  be assessed using variable importance (VIMP) measures for single, as
  well as grouped variables.  Variable selection is implemented using
  minimal depth variable selection (Ishwaran et al. 2010). Missing data
  (for x-variables and y-outcomes) can be imputed on both training and
  test data.  The underlying code is based on Ishwaran and Kogalur's now
  retired \pkg{randomSurvivalForest} package (Ishwaran and Kogalur
  2007), and has been significantly refactored for improved
  computational speed.
}
\section{OpenMP Parallel Processing -- Installation}{
 This package implements OpenMP shared-memory parallel programming.
  However, the default installation will only execute serially.  To
  utilize OpenMP, the target architecture and operating system must
  first support it.  

  To install the package with OpenMP parallel processing enabled, on most
  non-Windows systems, do the following:
  
  \enumerate{
    \item Download the package source code from CRAN.

    \item Open a console, navigate to the directory containing the
    tarball, and untar it using the command
    \command{tar -xvf randomForestSRC_X.x.x.tar.gz}

    \item This will create a directory structure with the root directory
    of the package named \command{randomForestSRC}.  Change into the root
    directory of the package using the command \command{cd randomForestSRC}

    \item Run autoconf using the command \command{autoconf}

    \item Change back to your working directory using the command
    \command{cd ..}

    \item Run \command{R CMD INSTALL randomForestSRC} on the modified
    package.  Ensure that you do not target the unmodified tarball, but
    instead act on the directory structure you just modified.
  }

  To install the package with OpenMP parallel processing enabled, on
  most Windows systems, do the following:
  
  \enumerate{
    \item Download the Windows specific custom binary build from
    \url{http://www.ccs.miami.edu/~hishwaran/rfsrc.html}

    \item If you are using the R GUI, start the GUI.  From the menu
    click on

    \code{Packages > Install package(s) from local zip files}

    Then navigate to the directory where you downloaded the zip file and
    click on it. 
  }
}
\section{OpenMP Parallel Processing -- Setting the Number of CPUs}{
  There are several ways to control the number of CPU cores that the
  package accesses during OpenMP parallel execution.  First, you will
  need to determine the number of cores on your local machine.  Do this
  by starting an R session and issuing the command
  \command{detectCores()}.

  Then you can do the following:

  At the start of every R session, you can set the number of cores
  accessed during OpenMP parallel execution by issuing the command
  \command{options(rf.cores = x)}, where \command{x} is the number of
  cores.  If \command{x} is a negative number, the package will access
  the maximum number of cores on your machine.  The options command can
  also be placed in the users .Rprofile file for convenience.  You can,
  alternatively, initialize the environment variable \command{RF_CORES}
  in your shell environment.

  The default value for rf.cores is -1 (-1L), if left unspecified, which
  uses all but one available cores, with a minimum of two.
}  
\section{R-side Parallel Processing -- Setting the Number of CPUs}{
  The package also implements R-side parallel processing by replacing
  the R function \command{lapply} with \command{mclapply} found in the
  \pkg{parallel} package.  You can set the number of cores accessed by
  \command{mclapply} by issuing the command \command{options(mc.cores =
  x)}, where \command{x} is the number of cores.  The options command
  can also be placed in the users .Rprofile file for convenience.  You
  can, alternatively, initialize the environment variable
  \command{MC_CORES} in your shell environment.  See the help files in
  \pkg{parallel} for more information.

  The default value for \command{mclapply} on non-Windows systems is
  two (2L) cores.  On Windows systems, the default value is one (1L)
  core.
}
\section{Example: Setting the Number of CPUs}{
  As an example, issuing the following options command uses all
  available cores minus one for both OpenMP and R-side processing:

  \command{options(rf.cores=detectCores()-1, mc.cores=detectCores()-1)}

  As stated above, this option command can be placed in the users
  .Rprofile file.
}
\section{Package Overview}{
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

  
