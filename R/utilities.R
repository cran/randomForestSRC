####**********************************************************************
####**********************************************************************
####
####  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
####  Version 1.4
####
####  Copyright 2012, University of Miami
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  ----------------------------------------------------------------
####  Project Partially Funded By: 
####  ----------------------------------------------------------------
####  Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
####  National Science Foundation and grant R01 CA163739 from the National
####  Cancer Institute.
####
####  Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
####  National Cancer Institute.
####  ----------------------------------------------------------------
####  Written by:
####  ----------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  commerce@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************


get.seed <- function (seed) {
  if ((is.null(seed)) || (abs(seed) < 1)) {
    seed <- runif(1,1,1e6)
  }
  seed <- -round(abs(seed))
  return (seed)
}
get.trace <- function (do.trace) {
  if (!is.logical(do.trace)) {
    if (do.trace >= 1) {
      do.trace <- round(do.trace)
    }
    else {
      do.trace <- 0
    }
  }
  else {
    do.trace <- 1 * do.trace
  }
  return (do.trace)
}
get.importance <-  function (importance) {
  if (!is.null(importance)) {
    if (importance == "none") {
      importance <- 0
    }
    else if (importance == "permute.ensemble") {
      importance <- 2^25 + 0
    }
    else if (importance == "random.ensemble") {
      importance <- 2^25 + 2^9
    }
    else if (importance == "permute.joint.ensemble") {
      importance <- 2^25 + 2^10 + 0
    }
    else if (importance == "random.joint.ensemble") {
      importance <- 2^25 + 2^10 + 2^9
    }
    else if (importance == "permute") {
      importance <- 2^25 + 2^24 + 0
    }
    else if (importance == "random") {
      importance <- 2^25 + 2^24 + 2^9
    }
    else if (importance == "permute.joint") {
      importance <- 2^25 + 2^24 + 2^10 + 0
    }
    else if (importance == "random.joint") {
      importance <- 2^25 + 2^24 + 2^10 + 2^9
    }
    else {
      stop("Invalid choice for 'importance' option:  ", importance)
    }
  }
  else {
    stop("Invalid choice for 'importance' option:  ", importance)
  }
  return (importance)
}
get.bootstrap <- function (bootstrap) {
  if (!is.null(bootstrap)) {
    if (bootstrap == "by.root") {
      bootstrap <- 0
    }
    else if (bootstrap == "by.node") {
      bootstrap <- 2^19
    }
    else if (bootstrap == "none") {
      bootstrap <- 2^20
    }
    else {
      stop("Invalid choice for 'bootstrap' option:  ", bootstrap)
    }
  }
  else {
    stop("Invalid choice for 'var.used' option:  ", bootstrap)
  }
  return (bootstrap)
}
get.forest <- function (forest) {
  if (!is.null(forest)) {
    if (forest == TRUE) {
      forest <- 2^5
    }
    else if (forest == FALSE) {
      forest <- 0
    }
    else {
      stop("Invalid choice for 'forest' option:  ", forest)
    }
  }
  else {
    stop("Invalid choice for 'forest' option:  ", forest)
  }
  return (forest)
}
get.proximity <- function (grow.equivalent, proximity) {
  if (!is.null(proximity)) {
    if (proximity == FALSE) {
      proximity <- 0
    }
    else if (grow.equivalent == TRUE) {
      if (proximity == TRUE) {
        proximity <- 2^28 
      }
      else if (proximity == "inbag") {
        proximity <- 2^28
      }
      else if (proximity == "oob") {
        proximity <- 2^29
      }
      else if (proximity == "all") {
        proximity <- 2^28 + 2^29
      }
      else {
        stop("Invalid choice for 'proximity' option:  ", proximity)
      }
    }
    else if (grow.equivalent == FALSE) {
      if (proximity == TRUE) {
        proximity <- 2^28 + 2^29
      }
      else if (proximity == "all") {
        proximity <- 2^28 + 2^29
      }
      else {
        stop("Invalid choice for 'proximity' option:  ", proximity)
      }
    }
    else {
      stop("Invalid choice for 'grow.equivalent' in proximity:  ", grow.equivalent)
    }
  }
  else {
    stop("Invalid choice for 'proximity' option:  ", proximity)
  }
  return (proximity)
}
get.split.null <- function (split.null) {
  if (!is.null(split.null)) {
    if (split.null == TRUE) {
      split.null <- 2^18
    }
    else if (split.null == FALSE) {
      split.null <- 0
    }
    else {
      stop("Invalid choice for 'split.null' option:  ", split.null)
    }
  }
  else {
    stop("Invalid choice for 'split.null' option:  ", split.null)
  }
  return (split.null)
}
get.outcome <- function (outcome) {
  if (is.null(outcome)) {
    outcome <- 0
  }
  else if (outcome == "train") {
    outcome <- 0
  }
  else if (outcome == "test") {
    outcome <- 2^17
  }
  else {
    stop("Invalid choice for 'outcome' option:  ", outcome)
  }
  return (outcome)
}
get.var.used <- function (var.used) {
  if (!is.null(var.used)) {
    if (var.used == "all.trees") {
      var.used <- 2^13 + 0
    }
    else if (var.used == "by.tree") {
      var.used <- 2^13 + 2^12
    }
    else if (var.used == FALSE) {
      var.used <- 0
    }
    else {
      stop("Invalid choice for 'var.used' option:  ", var.used)
    }
  }
  else {
    stop("Invalid choice for 'var.used' option:  ", var.used)
  }
  return (var.used)
}
get.split.depth <- function (split.depth) {
  if (!is.null(split.depth)) {
    if (split.depth == "all.trees") {
      split.depth <- 2^22
    }
    else if (split.depth == "by.tree") {
      split.depth <- 2^23
    }
    else if (split.depth == FALSE) {
      split.depth <- 0
    }
    else {
      stop("Invalid choice for 'split.depth' option:  ", split.depth)
    }
  }
  else {
    stop("Invalid choice for 'split.depth' option:  ", split.depth)
  }
  return (split.depth)
}
get.perf <-  function (perf, impute.only, family) {
  if (impute.only != TRUE) {
    if (!is.null(perf)) {
      if (perf == TRUE) {
        return (TRUE)
      }
      else if (perf == FALSE) {
        return (FALSE)
      }
      else {
        stop("Invalid choice for 'perf' option:  ", perf)
      }
    }
    else {
      return (TRUE)
    }
  }
  else {
    return (FALSE)
  }
}
get.cr.bits <- function (fmly) {
  if (fmly == "surv-CR") {
    return(2^21)
  } else {
    return(0)
  }
}
get.perf.bits <- function (perf) {
  if (perf) {
    return (2^2)
  }
  else {
    return (0)
  }
}
get.membership <- function (membership) {
  if (!is.null(membership)) {
    if (membership == TRUE) {
      membership <- 2^11
    }
    else if (membership == FALSE) {
      membership <- 0
    }
    else {
      stop("Invalid choice for 'membership' option:  ", membership)
    }
  }
  else {
    stop("Invalid choice for 'membership' option:  ", membership)
  }
  return (membership)
}
is.hidden.impute.only <-  function (user.option) {
  index = match("impute.only", names(user.option), 0)
  if(index == 0) {
    return (FALSE)
  }
  else {
    return (as.logical(as.character(user.option[index])))
  }
}
is.hidden.miss.tree.only <-  function (user.option) {
  index = match("miss.tree", names(user.option), 0)
  if(index == 0) {
    return (0)
  }
  else {
    miss.value <- as.character(user.option[index])
    if (is.na(as.logical(miss.value))) {
      return (as.numeric(miss.value))
    }
    else {
      if (as.logical(miss.value)) {
        return(1.0/3.0)
      }
      else {
        return (0)
      }
    } 
  }
}
get.impute.only <-  function (impute.only, nMiss) {
  if (impute.only) {
    if (nMiss > 0) {
      return (2^16)
    }
    else {
      stop("Data has no missing values, using 'impute' makes no sense.")
    }
  }
  else {
    return (0)
  }
}
get.restore.only <-  function (restore.only) {
  if (!is.null(restore.only)) {
    if (restore.only) {
      return (2^14)
    }
    else if (!restore.only) {
      return (0)
    }
    else {
      stop("Invalid choice for 'restore.only' option:  ", restore.only)    
    }
  }
  else {
    stop("Invalid choice for 'restore.only' option:  ", restore.only)
  }
}
get.vimp.only <-  function (vimp.only) {
  if (!is.null(vimp.only)) {
    if (vimp.only) {
      return (2^27)
    }
    else if (!vimp.only) {
      return (0)
    }
    else {
      stop("Invalid choice for 'vimp.only' option:  ", vimp.only)    
    }
  }
  else {
    stop("Invalid choice for 'vimp.only' option:  ", vimp.only)
  }
}
get.statistics <- function (statistics) {
  if (!is.null(statistics)) {
    if (statistics == TRUE) {
      statistics <- 2^27
    }
    else if (statistics == FALSE) {
      statistics <- 0
    }
    else {
      stop("Invalid choice for 'statistics' option:  ", statistics)
    }
  }
  else {
    stop("Invalid choice for 'statistics' option:  ", statistics)
  }
  return (statistics)
}
get.rf.cores <- function () {
  if (is.null(getOption("rf.cores"))) {
    if(!is.na(as.numeric(Sys.getenv("RF_CORES")))) {
      options(rf.cores = as.integer(Sys.getenv("RF_CORES")))
    }
  }
  return (getOption("rf.cores", -1L))
}
