////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.1.0
////
////  Copyright 2012, University of Miami
////
////  This program is free software; you can redistribute it and/or
////  modify it under the terms of the GNU General Public License
////  as published by the Free Software Foundation; either version 2
////  of the License, or (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public
////  License along with this program; if not, write to the Free
////  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
////  Boston, MA  02110-1301, USA.
////
////  ----------------------------------------------------------------
////  Project Partially Funded By: 
////  ----------------------------------------------------------------
////  Dr. Ishwaran's work was funded in part by DMS grant 1148991 from the
////  National Science Foundation and grant R01 CA163739 from the National
////  Cancer Institute.
////
////  Dr. Kogalur's work was funded in part by grant R01 CA163739 from the 
////  National Cancer Institute.
////  ----------------------------------------------------------------
////  Written by:
////  ----------------------------------------------------------------
////    Hemant Ishwaran, Ph.D.
////    Director of Statistical Methodology
////    Professor, Division of Biostatistics
////    Clinical Research Building, Room 1058
////    1120 NW 14th Street
////    University of Miami, Miami FL 33136
////
////    email:  hemant.ishwaran@gmail.com
////    URL:    http://web.ccs.miami.edu/~hishwaran
////    --------------------------------------------------------------
////    Udaya B. Kogalur, Ph.D.
////    Adjunct Staff
////    Dept of Quantitative Health Sciences
////    Cleveland Clinic Foundation
////    
////    Kogalur & Company, Inc.
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  ubk@kogalur.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include         "rfsrc.h"
#include         "entry.h"
SEXP rfsrcGrow(SEXP traceFlag,
               SEXP seedPtr,  
               SEXP opt,  
               SEXP splitRule,  
               SEXP splitRandomRule,  
               SEXP randomCovariateCount,  
               SEXP minimumNodeSize,
               SEXP maximumNodeDepth,
               SEXP crWeight,
               SEXP forestSize,  
               SEXP observationSize,
               SEXP rSize,
               SEXP rType,
               SEXP rLevels,
               SEXP rData,
               SEXP xSize,
               SEXP xType,
               SEXP xLevels,
               SEXP xWeight,
               SEXP xData,
               SEXP timeInterestSize,
               SEXP timeInterest,
               SEXP imputeSize,
               SEXP numThreads) {
  uint i;
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(opt)[0];
  RF_splitRule            = INTEGER(splitRule)[0];
  RF_splitRandomRule      = INTEGER(splitRandomRule)[0];
  RF_randomCovariateCount = INTEGER(randomCovariateCount)[0];
  RF_minimumNodeSize      = INTEGER(minimumNodeSize)[0];
  RF_maximumNodeDepth     = INTEGER(maximumNodeDepth)[0];
  RF_crWeight             = REAL(crWeight);  RF_crWeight--;
  RF_forestSize           = INTEGER(forestSize)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_rSize                = INTEGER(rSize)[0];
  RF_sexp_rType           = rType;
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_rData                = REAL(rData);
  RF_xSize                = INTEGER(xSize)[0];
  RF_sexp_xType           = xType;
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_randomCovariateWeight =REAL(xWeight);  RF_randomCovariateWeight--;
  RF_xData                = REAL(xData);
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest--;
  RF_imputeSize           = INTEGER(imputeSize)[0];
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_intrPredictorSize    = RF_xSize;
  RF_opt                  = RF_opt & (~OPT_VIMP_JOIN);
  RF_opt                  = RF_opt & (~OPT_OUTC_TYPE);
  RF_opt                  = RF_opt & (~OPT_COMP_RISK);
  if (RF_opt & OPT_IMPU_ONLY) {
    RF_opt                  = OPT_IMPU_ONLY | (RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE);
  }
  else {
    if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      RF_opt                  = RF_opt & (~OPT_OENS);  
      RF_opt                  = RF_opt & (~OPT_PERF);  
      RF_opt                  = RF_opt & (~OPT_PERF_CALB);  
      RF_opt                  = RF_opt & (~OPT_MORT_WGHT);  
      RF_opt                  = RF_opt & (~OPT_VIMP);  
    }
    else {
      RF_opt                  = RF_opt | OPT_OENS;  
      RF_opt                  = RF_opt | OPT_PERF;
    }
    RF_opt                  = RF_opt | OPT_FENS;  
    if (RF_opt & OPT_TREE) {
      RF_opt = RF_opt | OPT_SEED;
    }
    else {
      RF_opt = RF_opt & (~OPT_SEED);
    }
  }
  RF_opt                  = RF_opt | OPT_LEAF;  
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    RF_opt                  = RF_opt & (~OPT_MISS);
  }
  else {
    RF_opt                  = RF_opt | OPT_OMIS;  
  }
  RF_opt                  = RF_opt | OPT_MISS;
  if (seedValue >= 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Random seed must be less than zero.  \n");
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if ( RF_splitRule > RAND_SPLIT) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Invalid split rule:  %10d \n", RF_splitRule);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_splitRule == RAND_SPLIT) {
  }
  if ( ((RF_randomCovariateCount < 1) || (RF_randomCovariateCount > RF_xSize)) ) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of random covariate parameters must be greater");
    Rprintf("\nRF-SRC:  than zero and less than the total number of covariates:  %10d \n", RF_randomCovariateCount);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_minimumNodeSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Minimum node size must be greater than zero:  %10d \n", RF_minimumNodeSize);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  for (i = 1; i <= RF_xSize; i++) {
    if(RF_randomCovariateWeight[i] < 0) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Parameter verification failed.");
      Rprintf("\nRF-SRC:  Random covariate weight elements must be greater than or equal to zero:  %12.4f \n", RF_randomCovariateWeight[i]);
      Rprintf("\nRF-SRC:  The application will now exit.\n");
      return R_NilValue;
    }
  }
  return rfsrc(RF_GROW, seedValue, INTEGER(traceFlag)[0]);
}
SEXP rfsrcPredict(SEXP traceFlag,
                  SEXP seedPtr,  
                  SEXP opt,
                  SEXP forestSize, 
                  SEXP observationSize,
                  SEXP rSize,
                  SEXP rType,
                  SEXP rLevels,
                  SEXP rData,
                  SEXP xSize,
                  SEXP xType,
                  SEXP xLevels,
                  SEXP xData,
                  SEXP fobservationSize,
                  SEXP frSize,
                  SEXP frData,
                  SEXP fxData,
                  SEXP timeInterestSize,
                  SEXP timeInterest,
                  SEXP treeID,
                  SEXP nodeID,
                  SEXP parmID,
                  SEXP contPT,
                  SEXP mwcpSZ,
                  SEXP mwcpPT,
                  SEXP totalNodeCount,
                  SEXP minimumNodeSize,
                  SEXP seed,
                  SEXP intrPredictorSize,
                  SEXP intrPredictor,
                  SEXP numThreads) {
  int seedValue           = INTEGER(seedPtr)[0];
  RF_opt                  = INTEGER(opt)[0];
  RF_forestSize           = INTEGER(forestSize)[0];
  RF_observationSize      = INTEGER(observationSize)[0];
  RF_rSize                = INTEGER(rSize)[0];
  RF_sexp_rType           = rType;
  RF_rLevels              = INTEGER(rLevels); RF_rLevels--;
  RF_rData                = REAL(rData);
  RF_xSize                = INTEGER(xSize)[0];
  RF_sexp_xType           = xType;
  RF_xLevels              = INTEGER(xLevels); RF_xLevels--;
  RF_xData                = REAL(xData);
  RF_fobservationSize     = INTEGER(fobservationSize)[0];
  RF_frSize               = INTEGER(frSize)[0];
  RF_frData               = REAL(frData);
  RF_fxData               = REAL(fxData);
  RF_timeInterestSize     = INTEGER(timeInterestSize)[0];
  RF_timeInterest         = REAL(timeInterest);  RF_timeInterest --;
  RF_treeID_              = (uint*) INTEGER(treeID);  RF_treeID_ --;
  RF_nodeID_              = (uint*) INTEGER(nodeID);  RF_nodeID_ --;
  RF_parmID_              = (uint*) INTEGER(parmID);  RF_parmID_ --;
  RF_contPT_              = REAL(contPT);  RF_contPT_ --;
  RF_mwcpSZ_              = (uint*) INTEGER(mwcpSZ);  RF_mwcpSZ_ --;
  RF_mwcpPT_              = (uint*) INTEGER(mwcpPT);  RF_mwcpPT_ --;
  RF_totalNodeCount       = INTEGER(totalNodeCount)[0];
  RF_minimumNodeSize      = INTEGER(minimumNodeSize)[0];
  RF_seed_                = INTEGER(seed); RF_seed_ --;
  RF_intrPredictorSize    = INTEGER(intrPredictorSize)[0];
  RF_intrPredictor        = (uint*) INTEGER(intrPredictor);  RF_intrPredictor --;
  RF_numThreads           = INTEGER(numThreads)[0];
  RF_opt = RF_opt & (~OPT_TREE);  
  RF_opt = RF_opt & (~OPT_SEED);  
  RF_opt = RF_opt & (~OPT_IMPU_ONLY);
  RF_opt = RF_opt & (~OPT_SPLT_FAST);
  RF_imputeSize = 1;
  RF_opt                  = RF_opt | OPT_LEAF;  
  RF_opt                  = RF_opt | OPT_FENS;  
  RF_opt                  = RF_opt | OPT_MISS;  
  if (RF_opt & OPT_OUTC_TYPE) {
    RF_opt = RF_opt | OPT_REST;
    RF_opt = RF_opt & (~OPT_BOOT_NODE) & (~OPT_BOOT_NONE);
  }
  if (RF_opt & OPT_REST) {
    if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      RF_opt                  = RF_opt & (~OPT_OENS);  
      RF_opt                  = RF_opt & (~OPT_PERF);  
      RF_opt                  = RF_opt & (~OPT_PERF_CALB);  
      RF_opt                  = RF_opt & (~OPT_MORT_WGHT);  
      RF_opt                  = RF_opt & (~OPT_OMIS);  
      RF_opt                  = RF_opt & (~OPT_VIMP);  
    }
    else {
      RF_opt                  = RF_opt | OPT_OENS;  
      RF_opt                  = RF_opt | OPT_PERF;  
      RF_opt                  = RF_opt | OPT_OMIS;  
    }
  }
  else {
    RF_opt                  = RF_opt & (~OPT_OENS);  
    RF_opt                  = RF_opt & (~OPT_OMIS);  
    if (RF_fobservationSize < 1) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Parameter verification failed.");
      Rprintf("\nRF-SRC:  Number of individuals in prediction must be at least one:  %10d \n", RF_fobservationSize);
      Rprintf("\nRF-SRC:  The application will now exit.\n");
      return R_NilValue;
    }
  }
  if (seedValue >= 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  User random seed must be less than zero.  \n");
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_opt & OPT_REST) {
    return rfsrc(RF_REST, seedValue, INTEGER(traceFlag)[0]);
  }
  else {
    return rfsrc(RF_PRED, seedValue, INTEGER(traceFlag)[0]);
  }
}
