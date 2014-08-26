////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.5
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
////    email:  commerce@kogalur.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include     "factorOps.h"
#include     "splitUtil.h"
#include     "splitSurv.h"
#include     "splitRegr.h"
#include     "splitClas.h"
#include     "splitUspv.h"
#include    "regression.h"
#include         "split.h"
char getBestSplit(uint    treeID,
                  Node   *parent,
                  uint   *repMembrIndx,
                  uint    repMembrSize,
                  uint   *allMembrIndx,
                  uint    allMembrSize,
                  uint   *splitParameterMax,
                  double *splitValueMaxCont,
                  uint   *splitValueMaxFactSize,
                  uint  **splitValueMaxFactPtr,
                  double *splitStatistic,
                  char  **splitIndicator,
                  char    multImpFlag) {
  char  result;
  result = FALSE;  
  switch(RF_splitRule) {
  case SURV_LGRNK:
    result = logRankNCR(treeID,
                        parent,
                        repMembrIndx,
                        repMembrSize,
                        allMembrIndx,
                        allMembrSize,
                        splitParameterMax,
                        splitValueMaxCont,
                        splitValueMaxFactSize,
                        splitValueMaxFactPtr,
                        splitStatistic,
                        splitIndicator,
                        multImpFlag);
    break;
  case SURV_LRSCR:
    result = logRankNCR(treeID,
                        parent,
                        repMembrIndx,
                        repMembrSize,
                        allMembrIndx,
                        allMembrSize,
                        splitParameterMax,
                        splitValueMaxCont,
                        splitValueMaxFactSize,
                        splitValueMaxFactPtr,
                        splitStatistic,
                        splitIndicator,
                        multImpFlag);
    break;
  case SURV_CR_LAU:
    result = logRankCR(treeID,
                       parent,
                       repMembrIndx,
                       repMembrSize,
                       allMembrIndx,
                       allMembrSize,
                       splitParameterMax,
                       splitValueMaxCont,
                       splitValueMaxFactSize,
                       splitValueMaxFactPtr,
                       splitStatistic,
                       splitIndicator,
                       multImpFlag);
    break;
  case SURV_CR_LOG:
    result = logRankCR(treeID,
                       parent,
                       repMembrIndx,
                       repMembrSize,
                       allMembrIndx,
                       allMembrSize,
                       splitParameterMax,
                       splitValueMaxCont,
                       splitValueMaxFactSize,
                       splitValueMaxFactPtr,
                       splitStatistic,
                       splitIndicator,
                       multImpFlag);
    break;
  case RAND_SPLIT:
    result = randomSplit(treeID,
                         parent,
                         repMembrIndx,
                         repMembrSize,
                         allMembrIndx,
                         allMembrSize,
                         splitParameterMax,
                         splitValueMaxCont,
                         splitValueMaxFactSize,
                         splitValueMaxFactPtr,
                         splitStatistic,
                         splitIndicator,
                         multImpFlag);
    break;
  case REGR_WT_NRM:
    result = regressionXwghtSplit(treeID,
                                  parent,
                                  repMembrIndx,
                                  repMembrSize,
                                  allMembrIndx,
                                  allMembrSize,
                                  splitParameterMax,
                                  splitValueMaxCont,
                                  splitValueMaxFactSize,
                                  splitValueMaxFactPtr,
                                  splitStatistic,
                                  splitIndicator,
                                  multImpFlag);
       break;
  case REGR_WT_OFF:
       result = regressionXwghtSplit(treeID,
                                     parent,
                                     repMembrIndx,
                                     repMembrSize,
                                     allMembrIndx,
                                     allMembrSize,
                                     splitParameterMax,
                                     splitValueMaxCont,
                                     splitValueMaxFactSize,
                                     splitValueMaxFactPtr,
                                     splitStatistic,
                                     splitIndicator,
                                     multImpFlag);
       break;
  case REGR_WT_HVY:
       result = regressionXwghtSplit(treeID,
                                     parent,
                                     repMembrIndx,
                                     repMembrSize,
                                     allMembrIndx,
                                     allMembrSize,
                                     splitParameterMax,
                                     splitValueMaxCont,
                                     splitValueMaxFactSize,
                                     splitValueMaxFactPtr,
                                     splitStatistic,
                                     splitIndicator,
                                     multImpFlag);
       break;
  case CLAS_WT_NRM:
    result = classificationXwghtSplit(treeID,
                                      parent,
                                      repMembrIndx,
                                      repMembrSize,
                                      allMembrIndx,
                                      allMembrSize,
                                      splitParameterMax,
                                      splitValueMaxCont,
                                      splitValueMaxFactSize,
                                      splitValueMaxFactPtr,
                                      splitStatistic,
                                      splitIndicator,
                                      multImpFlag);
    break;
  case CLAS_WT_OFF:
    result = classificationXwghtSplit(treeID,
                                      parent,
                                      repMembrIndx,
                                      repMembrSize,
                                      allMembrIndx,
                                      allMembrSize,
                                      splitParameterMax,
                                      splitValueMaxCont,
                                      splitValueMaxFactSize,
                                      splitValueMaxFactPtr,
                                      splitStatistic,
                                      splitIndicator,
                                      multImpFlag);
    break;
  case CLAS_WT_HVY:
    result = classificationXwghtSplit(treeID,
                                      parent,
                                      repMembrIndx,
                                      repMembrSize,
                                      allMembrIndx,
                                      allMembrSize,
                                      splitParameterMax,
                                      splitValueMaxCont,
                                      splitValueMaxFactSize,
                                      splitValueMaxFactPtr,
                                      splitStatistic,
                                      splitIndicator,
                                      multImpFlag);
    break;
  case MVRG_SPLIT:
    result = multivariateSplit(treeID,
                               parent,
                               repMembrIndx,
                               repMembrSize,
                               allMembrIndx,
                               allMembrSize,
                               splitParameterMax,
                               splitValueMaxCont,
                               splitValueMaxFactSize,
                               splitValueMaxFactPtr,
                               splitStatistic,
                               splitIndicator,
                               multImpFlag);
    break;
  case MVCL_SPLIT:
    result = multivariateSplit(treeID,
                               parent,
                               repMembrIndx,
                               repMembrSize,
                               allMembrIndx,
                               allMembrSize,
                               splitParameterMax,
                               splitValueMaxCont,
                               splitValueMaxFactSize,
                               splitValueMaxFactPtr,
                               splitStatistic,
                               splitIndicator,
                               multImpFlag);
    break;
  case USPV_SPLIT:
    result = unsupervisedSplit(treeID,
                               parent,
                               repMembrIndx,
                               repMembrSize,
                               allMembrIndx,
                               allMembrSize,
                               splitParameterMax,
                               splitValueMaxCont,
                               splitValueMaxFactSize,
                               splitValueMaxFactPtr,
                               splitStatistic,
                               splitIndicator,
                               multImpFlag);
    break;
  default:
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid split rule:  %10d", RF_splitRule);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
    break;
  }
  return result;
}
char randomSplit(uint    treeID,
                 Node   *parent,
                 uint   *repMembrIndx,
                 uint    repMembrSize,
                 uint   *allMembrIndx,
                 uint    allMembrSize,
                 uint   *splitParameterMax,
                 double *splitValueMaxCont,
                 uint   *splitValueMaxFactSize,
                 uint  **splitValueMaxFactPtr,
                 double *splitStatistic,
                 char  **splitIndicator,
                 char    multImpFlag) {
  char    *randomCovariateFlag;
  uint    *randomCovariateIndex;
  uint    uniformCovariateIndex;
  uint    uniformCovariateSize;
  double *cdf;
  uint    cdfSize;
  uint   *cdfSort;
  uint   *density;
  uint    densitySize;
  uint  **densitySwap;
  uint     covariate;
  double  *permissibleSplit;
  uint     permissibleSplitSize;
  uint nonMissMembrSize, nonMissMembrSizeStatic;
  uint *nonMissMembrIndx, *nonMissMembrIndxStatic;
  uint   *indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint j;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             repMembrIndx,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if(result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    stackRandomCovariates(treeID,
                          parent,
                          repMembrSize,
                          multImpFlag,
                          & randomCovariateFlag,
                          & randomCovariateIndex,
                          & uniformCovariateSize,
                          & uniformCovariateIndex,
                          & cdf,
                          & cdfSize,
                          & cdfSort,
                          & density,
                          & densitySize,
                          & densitySwap);
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while ( ((*splitParameterMax) == 0) &&
            selectRandomCovariates(treeID,
                                   parent,
                                   repMembrIndx,
                                   repMembrSize,
                                   randomCovariateFlag,
                                   randomCovariateIndex,
                                   & uniformCovariateSize,
                                   & uniformCovariateIndex,
                                   cdf,
                                   & cdfSize,
                                   cdfSort,
                                   density,
                                   & densitySize,
                                   densitySwap,
                                   & covariate,
                                   & actualCovariateCount,
                                   & candidateCovariateCount,
                                   permissibleSplit,
                                   & permissibleSplitSize,
                                   & indxx,
                                   nonMissMembrSizeStatic,
                                   nonMissMembrIndxStatic,
                                   & nonMissMembrSize,
                                   & nonMissMembrIndx,
                                   multImpFlag)) {
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
      }
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 covariate,
                                                 permissibleSplit,
                                                 permissibleSplitSize,
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= nonMissMembrSize; j++) {
          localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
        }
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                              factorFlag,
                              mwcpSizeAbsolute,
                              covariate,
                              repMembrIndx,
                              repMembrSize,
                              nonMissMembrIndx,
                              nonMissMembrSize,
                              indxx,
                              permissibleSplitPtr,
                              j,
                              localSplitIndicator,
                              & leftSize,
                              priorMembrIter,
                              & currentMembrIter);
        rghtSize = nonMissMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          updateMaximumSplit(treeID,
                             0,  
                             covariate,
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             repMembrSize,
                             localSplitIndicator,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             permissibleSplitPtr,
                             splitIndicator);
          j = splitLength;
        }
      }  
      unstackSplitVector(treeID,
                         permissibleSplitSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
    }  
    unstackRandomCovariates(treeID,
                            randomCovariateFlag,
                            randomCovariateIndex,
                            uniformCovariateSize,
                            cdf,
                            cdfSize,
                            cdfSort,
                            density,
                            densitySize,
                            densitySwap,
                            repMembrSize,
                            nonMissMembrIndxStatic,
                            permissibleSplit);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
