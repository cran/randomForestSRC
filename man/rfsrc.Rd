\name{rfsrc}
\alias{rfsrc}
\alias{randomForestSRC}
\title{Random Forests for Survival, Regression, and Classification (RF-SRC)}
\description{
  A random forest (Breiman, 2001) is grown using user supplied training
  data.  Applies when the response is numeric, categorical (factor), or
  right-censored (including competing risk), and yields regression,
  classification, and survival forests, respectively.  The resulting
  forest, informally referred to as a RF-SRC object, contains many
  useful values which can be directly extracted by the user and/or parsed
  using additional functions (see the examples below).  This is the main
  entry point to the \pkg{randomForestSRC} package.

  The package implements OpenMP shared-memory parallel programming.
  However, the default installation will only execute serially.  Users
  should consult the randomForestSRC-package help file for details on
  installing the OpenMP version of the package.  The help file is
  readily accessible via the command
  \code{package?randomForestSRC}.

}
\usage{
rfsrc(formula, data, ntree = 1000,
  bootstrap = c("by.root", "by.node", "none"), mtry = NULL, nodesize = NULL,
  splitrule = NULL, nsplit = 0, split.fast = FALSE,
  importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none"),
  big.data = FALSE, na.action = c("na.omit", "na.impute"), nimpute = 1,
  cause, xvar.wt = NULL, forest = TRUE, proximity = FALSE,
  var.used = c(FALSE, "all.trees", "by.tree"),
  split.depth = c(FALSE, "all.trees", "by.tree"), seed = NULL,
  do.trace = FALSE, membership = TRUE, ...)
}
\arguments{
    \item{formula}{A symbolic description of the model to be fit.}
    \item{data}{Data frame containing the y-outcome and x-variables in
      the model.}
    \item{ntree}{Number of trees in the forest.}
    \item{bootstrap}{Bootstrap protocol.  The default is \code{by.root}
      which bootstraps the data by sampling with replacement at the root
      node before growing the tree.  If \code{by.node} is choosen, the
      data is bootstrapped at each node during the grow process.  If
      \code{none} is chosen, the data is not bootstrapped at all.  See
      the details below on prediction error when the default choice
      is not in effect.}
    \item{mtry}{Number of variables randomly selected as candidates for
      each split. For survival and classification the default is
      sqrt(\code{p}), where \code{p} equals the number of variables.
      For regression it is \code{p/3}.  Values are rounded up.}
    \item{nodesize}{Minimum number of unique cases (data points) in a
      terminal node.  The defaults are: survival (3), competing
      risk (6), regression (5), classification (1).}
    \item{splitrule}{Splitting rule used to grow trees.  Available rules
      are as follows:
      \describe{
	\item{\kbd{Regression:}}{The default rule is weighted mean-squared error
	  splitting (Breiman et al. 1984, Chapter 8.4).}
	
	\item{\kbd{Classification:}}{The default rule is Gini index splitting (Breiman
	  et al. 1984, Chapter 4.3).}
	
	\item{\kbd{Survival:}}{Two rules are available. (1) The default rule is
	  \code{logrank} which implements log-rank splitting (Segal, 1988;
	  Leblanc and Crowley, 1993); (2) \code{logrankscore} implements
	  log-rank score splitting (Hothorn and Lausen, 2003).}
	
	\item{\kbd{Competing risks:}}{The default rule is \code{logrankCR} which
	  implements a modified log-rank splitting rule modeled after Gray's
	  test (Gray, 1988).}
      }
    }
    \item{nsplit}{Non-negative integer value.  If non-zero, the
      specified tree splitting rule is randomized which can
      significantly increase speed.}
    \item{split.fast}{Set this value to TRUE when speed more than
      accuracy is important in growing a forest.}
    \item{importance}{Method for computing variable importance (VIMP).}
    \item{big.data}{Improved memory performance for large sample sizes.}
    \item{na.action}{Action to be taken if the data contains \code{NA}'s.
        Possible values are \code{na.omit} and \code{na.impute}.
        Default is \code{na.omit}, which removes the entire record if
        even one of its entries is \code{NA} (for x-variables this applies only
        to those specifically listed in 'formula').  The action
        \code{na.impute} implements a sophisticated tree imputation
        technique.}
    \item{nimpute}{Number of iterations of the missing data algorithm.
	  Performance measures such as out-of-bag (OOB) error rates tend
	  to become optimistic if \code{nimpute} is greater than 1.}
    \item{cause}{An integer value between 1 and \code{J} indicating
	  the event of interest for competing risks, where \code{J} is
	  the number of event types (this option applies only to
	  competing risks and is ignored otherwise).  While growing a
	  tree, the splitting of a node is restricted to the event type
	  specified by \code{cause}.  If not specified, the default is
	  to use a composite splitting rule which is an average over the
	  entire set of event types (a democratic approach).  Users can
	  also pass a vector of non-negative weights of length \code{J}
	  if they wish to use a customized composite split statistic
	  (for example, passing a vector of ones reverts to the default
	  composite split statistic).  In all instances when
	  \code{cause} is set incorrectly, splitting reverts to the
	  default.  Finally, note that regardless of how \code{cause} is
	  specified, the returned forest object always provides
	  estimates for all event types.}
    \item{xvar.wt}{Vector of non-negative weights where entry
        \code{k}, after normalizing, is the probability of selecting
        variable \code{k} as a candidate for splitting a node.  Default is to
        use uniform weights.  Vector must be of dimension \code{p},
        where \code{p} equals the number of variables, otherwise the
	  default is invoked.}
    \item{forest}{Should the forest object be returned?  Used for
        prediction on new data and required by many of the functions
	  used to parse the RF-SRC object.}
    \item{proximity}{Should the proximity between observations be
        calculated?  Creates an \code{n}x\code{n} matrix, which can be
        large.}
    \item{var.used}{Should a record of which variables were used for
        splitting be kept?  Default is \code{FALSE}.  Possible values
        are \code{all.trees} and \code{by.tree}.}
   \item{split.depth}{Records the minimal depth for each variable.
      Default is \code{FALSE}.  Possible values are \code{all.trees} and
      \code{by.tree}. Used for variable selection.}
    \item{seed}{Negative integer specifying seed for the random number
      generator.} 
    \item{do.trace}{Should trace output be enabled? Default is
        \code{FALSE}. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{membership}{Should terminal node membership and inbag
      information be returned?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  \enumerate{
  \item \emph{Families}

  There are four families of random forests: \code{regr}, \code{class},
  \code{surv}, and \code{surv-CR}.  Regression forests (\code{regr}) are
  constructed when the response is continuous, classification forests
  (\code{class}) when the response is a factor, survival forests
  (\code{surv}) when the response corresponds to right-censored survival
  settings, and competing risk survival forests (\code{surv-CR}) for
  right-censored competing risk scenarios.  See below for how to code
  the response in the two different survival scenarios.

    
  \item \emph{Randomized Splitting Rules}
 
  A random version of a splitting rule can be invoked using
  \option{nsplit}.  If \option{nsplit} is set to a non-zero positive
  integer, then a maximum of \option{nsplit} split points are chosen
  randomly for each of the \option{mtry} variables within a node (this is
  in contrast to non-random splitting, i.e. \option{nsplit=0}, where
  all possible split points for each of the \option{mtry} variables are
  considered).  The splitting rule is applied to the random split
  points and the node is split on that variable and random split point
  yielding the best value (as measured by the splitting rule).

  Pure random splitting can be invoked by setting
  \option{splitrule="random"}.  For each node, a variable is randomly
  selected and the node is split using a random split point (Cutler and
  Zhao, 2001; Lin and Jeon, 2006). 
  
  Trees tend to favor splits on continuous variables (Loh and Shih,
  1997), so it is good practice to use the \option{nsplit} option when
  the data contains a mix of continuous and discrete variables.  Using a
  reasonably small value mitigates bias and may not compromise accuracy.

  \item \emph{Fast Splitting}
  
  The value of \option{nsplit} has a significant impact on the time taken
  to grow a forest.  When non-random splitting is in effect
  (\option{nsplit=0}), iterating over each split point can sometimes be
  CPU intensive.  However, when \code{nsplit} > 0, or when pure random
  splitting is in effect, CPU times are drastically reduced.  Further,
  under such settings, the time to sort covariates is a significant
  portion of the overall grow process. Setting \option{split.fast=TRUE}
  in such cases omits the sorting/uniquifying of covariates and may
  further decrease execution times. This option should only be used when
  the data set is considered to contain unique values along
  covariates. In such situations sorting/uniquifying becomes redundant.

  \item \emph{Variable Importance}

  The option \option{importance} allows four distinct ways to calculate
  variable importance (VIMP).  The default \code{permute} returns
  Breiman-Cutler permutation VIMP as described in Breiman (2001).  For
  each tree, the prediction error on the out-of-bag (OOB) data is
  recorded.  Then for a given variable \emph{x}, out-of-bag (OOB) cases
  are randomly permuted in \emph{x} and the prediction error is
  recorded.  The VIMP for \emph{x} is defined as the difference between
  the perturbed and unperturbed error rate averaged over all trees.  If
  \code{random} is used, then \emph{x} is not permuted, but rather an
  OOB case is assigned a daughter node randomly whenever a split on
  \emph{x} is encountered in the in-bag tree.  The OOB error rate is
  compared to the OOB error rate without randomly splitting on \emph{x}.
  The VIMP is the difference averaged over the forest.  If the options
  \code{permute.ensemble} or \code{random.ensemble} are used, then VIMP is
  calculated by comparing the error rate for the perturbed OOB forest
  ensemble to the unperturbed OOB forest ensemble where the perturbed
  ensemble is obtained by either permuting \emph{x} or by random
  daughter node assignments for splits on \emph{x}.  Thus, unlike the
  Breiman-Cutler measure, here VIMP does not measure the tree average
  effect of \emph{x}, but rather its overall forest effect.  See
  Ishwaran et al. (2008) for further details.  Finally, the option
  \code{none} turns VIMP off entirely.

  Note that the function \command{vimp} provides a friendly user
  interface for extracting VIMP.
  
  \item \emph{Predition Error}

  Prediction error is calculated using OOB data.  Performance is
  measured in terms of mean-squared-error for regression and
  misclassification error for classification.  

  For survival, prediction error is measured by 1-C, where C is
  Harrell's (Harrell et al., 1982) concordance index.  Prediction error
  is between 0 and 1, and measures how well the predictor correctly
  ranks (classifies) two random individuals in terms of survival.  A
  value of 0.5 is no better than random guessing. A value of 0 is
  perfect.

  When bootstrapping is \code{by.node} or \code{none}, a coherent OOB
  subset is not available to assess prediction error.  Thus, all outputs
  dependent on this are suppressed.  In such cases, prediction error is
  only available via classical cross-validation (the user will need to
  use \command{predict.rfsrc}, for example).

  
  \item \emph{Survival, Competing Risks}

  \itemize{
  \item Survival settings require a time and censoring variable which
  should be identifed in the formula as the response using the standard
  \code{Surv} formula specification (see examples below).
    
  \item For survival forests (Ishwaran et al. 2008), the censoring
  variable must be coded as a non-negative integer with 0 reserved for
  censoring and (usually) 1=death (event).  The event time must be
  non-negative.

  \item For competing risk forests (Ishwaran et al., 2012), the
  implementation is similar to survival, but with the following
  caveats:

  \itemize{
    \item Censoring must be coded as a non-negative integer, where 0
    indicates right-censoring, and non-zero values indicate different
    event types.  While {0,1,2,..,J} is standard, and recommended,
    events can be coded non-sequentially, although 0 must always be used
    for censoring.
    
    \item Over-riding the default splitting rule \code{logrankCR} by
    manually selecting any split rule other than \code{logrankCR} or
    \code{random} will result in a survival analysis in which all
    (non-censored) events are treated as if they are the same type
    (indeed, they will coerced as such).  Note that \option{nsplit} works as
    in survival.

    \item Generally, \code{nodesize} should be set larger for competing
    risks than survival settings.
  }
  }
  
  \item \emph{Imputation}

  Setting \option{na.action="na.impute"} implements a tree imputation
  method whereby missing data (x-variables or y-outcomes) are imputed
  dynamically as a tree is grown by randomly sampling from the
  distribution within the current node (Ishwaran et al. 2008).  OOB data
  is not used in imputation to avoid biasing prediction error and VIMP
  estimates.  Final imputation for integer valued variables and
  censoring indicators use a maximal class rule, whereas continuous
  variables and survival time use a mean rule.  Records in which all
  outcome and x-variable information are missing are removed.  Variables
  having all missing values are removed.  The algorithm can be iterated
  by setting \code{nimpute} to a positive integer greater than 1.  A few
  iterations should be used in heavy missing data settings to improve
  accuracy of imputed values (see Ishwaran et al., 2008).  Note if the
  algorithm is iterated, a side effect is that missing values in the
  returned objects \code{xvar}, \code{yvar} are replaced by imputed
  values.  Further, imputed objects such as \code{imputed.data} are set
  to \code{NULL}.

  See the function \command{impute.rfsrc} for a fast impute interface.

  \item \emph{Large n}

  For increased efficiency, users should consider setting
  \option{big.data=TRUE} when the sample size (number of rows of the data)
  is large.  This invokes efficient memory management.  For survival
  families a side-effect is that survival functions are constrained to a
  restricted grid of time points.
  
  \item \emph{Factors}
  
  Variables encoded as factors are treated as such.  If the factor is
  ordered, then splits are similar to real valued variables.  If the
  factor is unordered, a split will move a subset of the levels in the
  parent node to the left daughter, and the complementary subset to the
  right daughter.  All possible complementary pairs are considered and
  apply to factors with an unlimited number of levels.  However, there
  is an optimization check to ensure that the number of splits attempted
  is not greater than the number of cases in a node (this internal check
  will override the \code{nsplit} value in random splitting mode if
  \code{nsplit} is large enough).

  All x-variables other than factors are coerced and treated as real
  valued.
  }
  \emph{Miscellanea}
  \enumerate{
    \item Setting \option{var.used="all.trees"} returns a vector of size
    \code{p} where each element is a count of the number of times a split
    occurred on a variable.  If \option{var.used="by.tree"}, a matrix of
    size \code{ntree}x\code{p} is returned.  Each element [i][j] is the
    count of the number of times a split occurred on variable [j] in tree
    [i].
    
    \item Setting \option{split.depth=TRUE} returns a matrix of size
    \code{n}x\code{p} where entry [i][j] is the minimal depth for variable
    [j] for case [i].  Used to select variables at the case-level.  See
    \code{max.subtree} for more details regarding minimal depth.
  }
}   
\value{
     An object of class \code{(rfsrc, grow)} with the following
     components:
     
     \item{call}{The original call to \command{rfsrc}.}
     \item{formula}{The formula used in the call.}
     \item{family}{The family used in the analysis.}
     \item{n}{Sample size of the data (depends upon \code{NA}'s, see \option{na.action}).}
     \item{ntree}{Number of trees grown.}
     \item{mtry}{Number of variables randomly selected for splitting at each node.}
     \item{nodesize}{Minimum size of terminal nodes.}
     \item{splitrule}{Splitting rule used.}
     \item{nsplit}{Number of randomly selected split points.}
     \item{yvar}{y-outcome values.}
     \item{yvar.names}{A character vector of the y-outcome names.}
     \item{xvar}{Data frame of x-variables.}
     \item{xvar.names}{A character vector of the x-variable names.}
     \item{xvar.wt}{Vector of non-negative weights specifying the
        probability used to select a variable for splitting a node.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         forest. Vector of length \option{ntree}.  A value of zero indicates
         a rejected tree (can occur when imputing missing data).
         Values of one indicate tree stumps.}
     \item{forest}{If \option{forest=TRUE}, the forest object is returned.
         This object is used for prediction with new test data
	 sets and is required for other R-wrappers.}
     \item{proximity}{If \option{proximity=TRUE}, a matrix of dimension
         \code{n}x\code{n} recording the frequency pairs of data points
         occur within the same terminal node.}
     \item{membership}{Matrix of dimension \code{n}x\code{ntree}
       recording terminal node membership where each column contains the
       node number that a case falls in for that tree.}
     \item{inbag}{Matrix of dimension \code{n}x\code{ntree}
       recording inbag membership where each column contains the number
     of times that a case appears in the bootstrap sample for that tree.}
     \item{var.used}{Count of the number of times a variable is used in
       growing the forest.}
     \item{imputed.indv}{Vector of indices for cases with missing
         values.}
     \item{imputed.data}{Data frame of the imputed data. The first
	 column(s) are reserved for the y-outcomes, after which the
	 x-variables are listed.}
     \item{split.depth}{Matrix of size \code{n}x\code{p} where entry
         [i][j] is the minimal depth for variable [j] for case [i].}
     \item{err.rate}{Tree cumulative OOB error rate.}
     \item{importance}{Variable importance (VIMP) for each x-variable.}
     \item{predicted}{In-bag predicted value.}
     \item{predicted.oob}{OOB predicted value.}\cr

\item{...... class}{for classification settings, additionally the following ......} \cr     
      \item{class}{In-bag predicted class labels.}
      \item{class.oob}{OOB predicted class labels.}\cr

\item{...... surv}{for survival settings, additionally the following ......} \cr
     
     \item{survival}{In-bag survival function.}
     \item{survival.oob}{OOB survival function.}
     \item{chf}{In-bag cumulative hazard function (CHF).}
     \item{chf.oob}{OOB CHF.}
     \item{time.interest}{Ordered unique death times.}
     \item{ndead}{Number of deaths.}\cr

\item{...... surv-CR}{for competing risks, additionally the following ......} \cr

     \item{chf}{In-bag cause-specific cumulative hazard function (CSCHF)
     for each event.}
     \item{chf.oob}{OOB CSCHF.}
     \item{cif}{In-bag cumulative incidence function (CIF) for each event.}
     \item{cif.oob}{OOB CIF.}
     \item{time.interest}{Ordered unique event times.}
     \item{ndead}{Number of events.}
}
\note{
  \enumerate{
    \item The returned object depends heavily on the family:
    \itemize{
      \item For regression, the predicted value is the predicted
      y-response.

      \item For classification, the predicted response is a matrix with
      columns containing the estimated class probability for each class.

      \item For survival, the predicted response is mortality (Ishwaran et
      al., 2008), a value representing estimated risk calibrated to the
      scale of the number of events.  As a specific example, if \emph{i} has
      a mortality value of 100, then if all individuals had the same
      x-values as \emph{i}, we expect an average of 100 events.  Also
      included in the grow object are matrices containing the CHF and
      survival function.  Each row corresponds to an individual's ensemble
      CHF or survival function evaluated at each time point in
      \code{time.interest}.

      \item For competing risks, the predicted response is a matrix with
      one column for each event recording the expected number of life
      years lost due to that event specific cause (Ishwaran et al.,
      2012).  The grow object also contains the cause-specific
      cumulative hazard function (CSCHF) and the cumulative incidence
      function (CIF) for each event type.  These are encoded as a 3-d
      array, with the third dimension used for the event type, each time
      point in \code{time.interest} making up the second dimension
      (columns), and the case (individual) being the first dimension
      (rows).
    }
  
    \item Different R-wrappers are provided to aid in parsing the grow
    object.
  }
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Breiman L., Friedman J.H., Olshen R.A. and Stone C.J. 
    \emph{Classification and Regression Trees}, Belmont, California, 1984.
    
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Cutler A. and Zhao G. (2001). Pert-Perfect random tree ensembles.
    \emph{Comp. Sci. Statist.}, 33: 490-497.

    Gray R.J. (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16: 1141-1154.

    Harrell et al. F.E. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247:2543-2546.

    Hothorn T. and Lausen B. (2003). On the exact distribution of maximally selected
    rank statistics, \emph{Comp. Statist. Data Anal.}, 43:121-137.    

    Ishwaran H. (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    Ishwaran H. and Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.
    
    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc.}, 105:205-217.

    Ishwaran H. (2012).  The effect of splitting on random forests.
    
    Ishwaran H., Gerds, T.A. Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2012). Random survival forests for competing risks.

    Lin Y. and Jeon Y. (2006).  Random forests and adaptive nearest
    neighbors, \emph{J. Amer. Statist. Assoc.}, 101:578-590.

    LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
    \emph{J. Amer. Statist. Assoc.}, 88:457-467.
    
    Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
    classification trees, \emph{Statist. Sinica}, 7:815-840.

    Mogensen, U.B, Ishwaran H. and Gerds T.A. (2012). Evaluating random
    forests for survival analysis using prediction error curves,
    \emph{J.  Statist. Software}, 50(11): 1-23.

    Segal M.R. (1988).  Regression trees for censored data, 
    \emph{Biometrics},  44:35-47.
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
  \command{\link{restore}},
  \command{\link{var.select}},
  \command{\link{vimp}}
}
\examples{
###------------------------------------------------------------
### Survival analysis

### Veteran data
### Randomized trial of two treatment regimens for lung cancer

data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100)

# print and plot the grow object
print(v.obj)
plot(v.obj)

# plot survival curves for first 10 individuals: direct way
matplot(v.obj$time.interest, 100 * t(v.obj$survival[1:10, ]),
    xlab = "Time", ylab = "Survival", type = "l", lty = 1)

# plot survival curves for first 10 individuals
# indirect way: using plot.survival (also generates hazard plots)
plot.survival(v.obj, subset = 1:10, haz.model = "ggamma")

\dontrun{
### Primary biliary cirrhosis (PBC) of the liver

data(pbc, package = "randomForestSRC") 
pbc.obj <- rfsrc(Surv(days, status) ~ ., pbc, nsplit = 10)
pbc.obj

# imputation
pbc.obj2 <- rfsrc(Surv(days, status)~., pbc, 
                nsplit=10, na.action="na.impute")
pbc.obj2

# here's a nice wrapper to combine original data + imputed data
combine.impute <- function(object) {
 impData <- cbind(object$yvar, object$xvar)
 if (!is.null(object$imputed.indv)) {
   impData[object$imputed.indv, ] <- object$imputed.data
 }
 impData
}

# combine original data + imputed data
pbc.imputed.data <- combine.impute(pbc.obj2)

# iterate the missing data algorithm
# compare to non-iterated algorithm
pbc.obj3 <- rfsrc(Surv(days, status) ~ ., pbc, nsplit=10,  
         na.action="na.impute", nimpute = 3)
pbc.iterate.imputed.data <- combine.impute(pbc.obj3)
tail(pbc.imputed.data)
tail(pbc.iterate.imputed.data)

# fast way to impute the data: no inference is done
# see impute.rfsc for more details
pbc.fast.impute.data <- impute.rfsrc(Surv(days, status) ~ ., pbc,
          nsplit = 10, nimpute = 3)
tail(pbc.fast.impute.data)


### Competing risks
### WIHS analysis
### Cumulative incidence function (CIF) for HAART and AIDS stratified by IDU

data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree = 100)
plot.competing.risk(wihs.obj)
cif <- wihs.obj$cif
Time <- wihs.obj$time.interest
idu <- wihs$idu
cif.haart <- cbind(apply(cif[,,1][idu == 0,], 2, mean), apply(cif[,,1][idu == 1,], 2, mean))
cif.aids  <- cbind(apply(cif[,,2][idu == 0,], 2, mean), apply(cif[,,2][idu == 1,], 2, mean))
matplot(Time, cbind(cif.haart, cif.aids), type = "l", 
        lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3,
        ylab = "Cumulative Incidence")  
legend("topleft",
       legend = c("HAART (Non-IDU)", "HAART (IDU)", "AIDS (Non-IDU)", "AIDS (IDU)"),
       lty = c(1,2,1,2), col = c(4, 4, 2, 2), lwd = 3, cex = 1.5)


### Compare RF-SRC to cox regression
### Assumes "pec" and "survival" libraries are loaded 

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE))
{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }

  ## data, formula specifications
  data(pbc, package = "randomForestSRC")
  pbc.na <- na.omit(pbc)  ##remove NA's
  surv.f <- as.formula(Surv(days, status) ~ .)
  pec.f <- as.formula(Hist(days,status) ~ 1)

  ## run cox/rfsrc models
  ## for illustration we use a small number of trees
  cox.obj <- coxph(surv.f, data=pbc.na)
  rfsrc.obj <- rfsrc(surv.f, pbc.na, nsplit = 10, ntree = 150)

  ## compute bootstrap cross-validation estimate of expected Brier score
  ## see Mogensen, Ishwaran and Gerds (2012) Journal of Statistical Software
  set.seed(17743)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = pbc.na, formula = pec.f,
                        splitMethod = "bootcv", B = 50)
  print(prederror.pbc)
  plot(prederror.pbc)

  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, pbc.na, nsplit = 10)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b\%\%10 == 0) cat("cox bootstrap:", b, "\n") 
    train <- sample(1:nrow(pbc.na), nrow(pbc.na), replace = TRUE) 
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    if (is.list(cox.obj)) {
      rcorr.cens(predict(cox.obj, pbc.na[-train, ]),
                 Surv(pbc.na$days[-train],
                 pbc.na$status[-train]))[1]
    } else NA 
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRandom Forests [S]RC: ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression      : ", mean(cox.err, na.rm = TRUE), "\n")
}

###------------------------------------------------------------
### Regression analysis

### New York air quality measurements

airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute")
airq.obj

# partial plot of variables (see plot.variable for more details)
plot.variable(airq.obj, partial = TRUE, smooth.lines = TRUE)

### motor trend cars

mtcars.obj <- rfsrc(mpg ~ ., data = mtcars)
mtcars.obj

# minimal depth variable selection via max.subtree
md.obj <- max.subtree(mtcars.obj)
cat("top variables:\n")
print(md.obj$topvars)

# equivalent way to select variables
# see var.select for more details
vs.obj <- var.select(object = mtcars.obj)


###------------------------------------------------------------
### Classification analysis

### Edgar Anderson's iris data

iris.obj <- rfsrc(Species ~., data = iris)
iris.obj

### Wisconsin prognostic breast cancer data

data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast, nsplit = 10)
breast.obj
plot(breast.obj)


}}
\keyword{forest}

