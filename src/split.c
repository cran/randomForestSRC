////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.3
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
                  double *splitStatistic) {
  char  result;
  result = FALSE;  
  switch(RF_splitRule) {
  case SURV_LGRNK:
    result = logRank(treeID, 
                     parent, 
                     repMembrIndx, 
                     repMembrSize, 
                     allMembrIndx,
                     allMembrSize,
                     splitParameterMax, 
                     splitValueMaxCont, 
                     splitValueMaxFactSize, 
                     splitValueMaxFactPtr,
                     splitStatistic);
    break;
  case SURV_LRSCR:
    result = logRankScore(treeID, 
                          parent, 
                          repMembrIndx, 
                          repMembrSize, 
                          allMembrIndx,
                          allMembrSize,
                          splitParameterMax, 
                          splitValueMaxCont, 
                          splitValueMaxFactSize, 
                          splitValueMaxFactPtr,
                          splitStatistic);
    break;
  case SURV_CR_LAU:
    result = logRankLauCR(treeID, 
                          parent, 
                          repMembrIndx, 
                          repMembrSize, 
                          allMembrIndx,
                          allMembrSize,
                          splitParameterMax, 
                          splitValueMaxCont, 
                          splitValueMaxFactSize, 
                          splitValueMaxFactPtr,
                          splitStatistic);
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
                       splitStatistic);
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
                         splitStatistic);
    break;
  case REGR_WT_NRM:
       result = regressionSplit(treeID, 
                                parent, 
                                repMembrIndx, 
                                repMembrSize, 
                                allMembrIndx,
                                allMembrSize,
                                splitParameterMax, 
                                splitValueMaxCont, 
                                splitValueMaxFactSize, 
                                splitValueMaxFactPtr,
                                splitStatistic);
       break;
  case REGR_WT_OFF:
       result = regressionUwghtSplit(treeID, 
                                     parent, 
                                     repMembrIndx, 
                                     repMembrSize, 
                                     allMembrIndx,
                                     allMembrSize,
                                     splitParameterMax, 
                                     splitValueMaxCont, 
                                     splitValueMaxFactSize, 
                                     splitValueMaxFactPtr,
                                     splitStatistic);
       break;
  case REGR_WT_HVY:
       result = regressionHwghtSplit(treeID, 
                                      parent, 
                                      repMembrIndx, 
                                      repMembrSize, 
                                      allMembrIndx,
                                      allMembrSize,
                                      splitParameterMax, 
                                      splitValueMaxCont, 
                                      splitValueMaxFactSize, 
                                      splitValueMaxFactPtr,
                                      splitStatistic);
       break;
  case CLAS_WT_NRM:
    result = classificationSplit(treeID, 
                                 parent, 
                                 repMembrIndx, 
                                 repMembrSize, 
                                 allMembrIndx,
                                 allMembrSize,
                                 splitParameterMax, 
                                 splitValueMaxCont, 
                                 splitValueMaxFactSize, 
                                 splitValueMaxFactPtr,
                                 splitStatistic);
    break;
  case CLAS_WT_OFF:
    result = classificationUwghtSplit(treeID, 
                                       parent, 
                                       repMembrIndx, 
                                       repMembrSize, 
                                       allMembrIndx,
                                       allMembrSize,
                                       splitParameterMax, 
                                       splitValueMaxCont, 
                                       splitValueMaxFactSize, 
                                       splitValueMaxFactPtr,
                                       splitStatistic);
    break;
  case CLAS_WT_HVY:
    result = classificationHwghtSplit(treeID, 
                                       parent, 
                                       repMembrIndx, 
                                       repMembrSize, 
                                       allMembrIndx,
                                       allMembrSize,
                                       splitParameterMax, 
                                       splitValueMaxCont, 
                                       splitValueMaxFactSize, 
                                       splitValueMaxFactPtr,
                                       splitStatistic);
    break;
  case MVRG_SPLIT:
    result = mvRegressionSplit(treeID, 
                               parent, 
                               repMembrIndx, 
                               repMembrSize, 
                               allMembrIndx,
                               allMembrSize,
                               splitParameterMax, 
                               splitValueMaxCont, 
                               splitValueMaxFactSize, 
                               splitValueMaxFactPtr,
                               splitStatistic);
    break;
  case MVCL_SPLIT:
    result = mvClassificationSplit(treeID, 
                                   parent, 
                                   repMembrIndx, 
                                   repMembrSize, 
                                   allMembrIndx,
                                   allMembrSize,
                                   splitParameterMax, 
                                   splitValueMaxCont, 
                                   splitValueMaxFactSize, 
                                   splitValueMaxFactPtr,
                                   splitStatistic);
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
                 double *splitStatistic) {
  char result;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    result = randomSurvivalSplit(treeID, 
                                 parent, 
                                 repMembrIndx, 
                                 repMembrSize, 
                                 allMembrIndx,
                                 allMembrSize,
                                 splitParameterMax, 
                                 splitValueMaxCont, 
                                 splitValueMaxFactSize, 
                                 splitValueMaxFactPtr,
                                 splitStatistic);
  }
  else {
    result = randomNonSurvivalSplit(treeID, 
                                    parent, 
                                    repMembrIndx, 
                                    repMembrSize, 
                                    allMembrIndx,
                                    allMembrSize,
                                    splitParameterMax, 
                                    splitValueMaxCont, 
                                    splitValueMaxFactSize, 
                                    splitValueMaxFactPtr,
                                    splitStatistic);
  }
  return result;
}
char randomNonSurvivalSplit(uint    treeID, 
                            Node   *parent, 
                            uint   *repMembrIndx,
                            uint    repMembrSize,
                            uint   *allMembrIndx,
                            uint    allMembrSize,
                            uint   *splitParameterMax, 
                            double *splitValueMaxCont, 
                            uint   *splitValueMaxFactSize, 
                            uint  **splitValueMaxFactPtr,
                            double *splitStatistic) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    if (RF_maximumNodeDepth < 0) {
      result = TRUE;
    }
    else {
      if (parent -> depth < (uint) RF_maximumNodeDepth) {
        result = TRUE;
      }
      else {
        result = FALSE;
      }
    }
  }
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
  }
  if(result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    char *covariateStatus = NULL;  
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx,
                                                               repMembrSize,
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    if (actualCovariateCount > 0) {
      covariateStatus = cvector(1, actualCovariateCount);
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      covariateStatus[i] = TRUE;
    }
    i = getSelectableElement(treeID, actualCovariateCount, covariateStatus, NULL);
    while ((i != 0) && ((*splitParameterMax) == 0)) {
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      for (j = 1; j < splitLength; j++) {
        leftSize = virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           repMembrIndx,
                           repMembrSize,
                           permissibleSplitPtr,
                           j,
                           0,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           NULL,
                           localSplitIndicator);
        rghtSize = repMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          updateMaximumSplit(0,  
                             randomCovariateIndex[i],
                             j,
                             factorFlag,
                             mwcpSizeAbsolute,
                             & deltaMax,
                             splitParameterMax,
                             splitValueMaxCont,
                             splitValueMaxFactSize,
                             splitValueMaxFactPtr,
                             permissibleSplitPtr);
          j = splitLength;
        }  
      }  
      unstackSplitVector(treeID,
                         permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
      if(*splitParameterMax == 0) {
        covariateStatus[i] = FALSE;
        i = getSelectableElement(treeID, actualCovariateCount, covariateStatus, NULL);
      }
    }  
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    if (actualCovariateCount > 0) {
      free_cvector(covariateStatus, 1, actualCovariateCount);
    }
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
