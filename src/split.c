////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.2
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
////    email:  kogalurshear@gmail.com
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
#include         "split.h"
#include     "splitSurv.h"
char getBestSplit(uint    treeID, 
                  Node   *parent, 
                  uint   *repMembrIndx,
                  uint    repMembrSize,
                  uint   *allMembrIndx,
                  uint    allMembrSize,
                  uint   *splitParameterMax,
                  double *splitValueMaxCont,
                  uint   *splitValueMaxFactSize,
                  uint  **splitValueMaxFactPtr) {
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
                     splitValueMaxFactPtr);
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
                          splitValueMaxFactPtr);
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
                          splitValueMaxFactPtr);
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
                         splitValueMaxFactPtr);
    break;
  case REGR_SPLIT:
       result = regressionSplit(treeID, 
                             parent, 
                             repMembrIndx, 
                             repMembrSize, 
                             allMembrIndx,
                             allMembrSize,
                             splitParameterMax, 
                             splitValueMaxCont, 
                             splitValueMaxFactSize, 
                             splitValueMaxFactPtr);
    break;
  case CLAS_SPLIT:
    result = classificationSplit(treeID, 
                                 parent, 
                                 repMembrIndx, 
                                 repMembrSize, 
                                 allMembrIndx,
                                 allMembrSize,
                                 splitParameterMax, 
                                 splitValueMaxCont, 
                                 splitValueMaxFactSize, 
                                 splitValueMaxFactPtr);
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
                 uint  **splitValueMaxFactPtr) {
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
                                 splitValueMaxFactPtr);
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
                                    splitValueMaxFactPtr);
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
                            uint  **splitValueMaxFactPtr) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
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
  deltaMax               = - EPSILON;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    result = getStandardDeviation(repMembrSize, repMembrIndx, RF_response[treeID][1]);
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
                                                               & permissibleSplitSize);
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
                            permissibleSplit, 
                            permissibleSplitSize);
    if (actualCovariateCount > 0) {
      free_cvector(covariateStatus, 1, actualCovariateCount);
    }
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                deltaMax);
  return result;
}
char regressionSplit (uint    treeID, 
                      Node   *parent, 
                      uint   *repMembrIndx,
                      uint    repMembrSize,
                      uint   *allMembrIndx,
                      uint    allMembrSize,
                      uint   *splitParameterMax, 
                      double *splitValueMaxCont, 
                      uint   *splitValueMaxFactSize, 
                      uint  **splitValueMaxFactPtr) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumParent;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = - EPSILON;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    result = getStandardDeviation(repMembrSize, repMembrIndx, RF_response[treeID][1]);
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint   *nodeLeftIndex  = uivector(1, repMembrSize);
    uint   *nodeRghtIndex  = uivector(1, repMembrSize);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
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
          leftSize = rghtSize = 0;
          for (k=1; k <= repMembrSize; k++) {
            nodeLeftIndex[k] = nodeRghtIndex[k] = 0;
            if (localSplitIndicator[k] == LEFT) {
              nodeLeftIndex[++leftSize] = repMembrIndx[k];
            }
            else {
              nodeRghtIndex[++rghtSize] = repMembrIndx[k];
            }
          }
          sumLeft = 0;
          for (k=1; k <= leftSize; k++) {
            sumLeft += RF_response[treeID][1][nodeLeftIndex[k]];
          }
          sumParent = 0;
          for (k=1; k <= repMembrSize; k++) {
            sumParent += RF_response[treeID][1][repMembrIndx[k]];
          }
          sumRght = sumParent - sumLeft;
          sumLeft = pow (sumLeft, 2.0) / leftSize;
          sumRght = pow (sumRght, 2.0) / rghtSize;
          sumParent = pow (sumParent, 2.0) / repMembrSize;
          delta = sumLeft + sumRght - sumParent;
          updateMaximumSplit(delta,
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
        }  
      }  
      unstackSplitVector(treeID,
                         permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_uivector(nodeLeftIndex, 1, repMembrSize);
    free_uivector(nodeRghtIndex, 1, repMembrSize);
  }  
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                deltaMax);
  return result;
}
char weightedRegressionSplit (uint    treeID, 
                              Node   *parent, 
                              uint   *repMembrIndx,
                              uint    repMembrSize,
                              uint   *allMembrIndx,
                              uint    allMembrSize,
                              uint   *splitParameterMax, 
                              double *splitValueMaxCont, 
                              uint   *splitValueMaxFactSize, 
                              uint  **splitValueMaxFactPtr) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, normSize, meanLeft, meanRght;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = - DBL_MAX;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  if (result) {
    result = getStandardDeviation(repMembrSize, repMembrIndx, RF_response[treeID][1]);
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint   *nodeLeftIndex  = uivector(1, repMembrSize);
    uint   *nodeRghtIndex  = uivector(1, repMembrSize);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
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
          leftSize = rghtSize = 0;
          for (k=1; k <= repMembrSize; k++) {
            nodeLeftIndex[k] = nodeRghtIndex[k] = 0;
            if (localSplitIndicator[k] == LEFT) {
              nodeLeftIndex[++leftSize] = repMembrIndx[k];
            }
            else {
              nodeRghtIndex[++rghtSize] = repMembrIndx[k];
            }
          }
          meanLeft = 0;
          for (k=1; k <= leftSize; k++) {
            meanLeft += RF_response[treeID][1][nodeLeftIndex[k]];
          }
          meanLeft = meanLeft / leftSize;
          meanRght = 0;
          for (k=1; k <= rghtSize; k++) {
            meanRght += RF_response[treeID][1][nodeRghtIndex[k]];
          }
          meanRght = meanRght / rghtSize;
          sumLeft = 0;
          for (k=1; k <= leftSize; k++) {
            sumLeft += pow (RF_response[treeID][1][nodeLeftIndex[k]] - meanLeft, 2.0);
          }
          sumRght = 0;
          for (k=1; k <= rghtSize; k++) {
            sumRght += pow (RF_response[treeID][1][nodeRghtIndex[k]] - meanRght, 2.0);
          }
          normSize = pow (repMembrSize, 2.0);
          sumLeft = sumLeft * leftSize / normSize;
          sumRght = sumRght * rghtSize / normSize;
          delta = - (sumLeft + sumRght);
          updateMaximumSplit(delta,
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
        }  
      }  
      unstackSplitVector(treeID,
                         permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_uivector(nodeLeftIndex, 1, repMembrSize);
    free_uivector(nodeRghtIndex, 1, repMembrSize);
  }  
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                deltaMax);
  return result;
}
char classificationSplit (uint    treeID, 
                          Node   *parent, 
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *allMembrIndx,
                          uint    allMembrSize,
                          uint   *splitParameterMax, 
                          double *splitValueMaxCont, 
                          uint   *splitValueMaxFactSize, 
                          uint  **splitValueMaxFactPtr) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint responseClassCount;
  uint *parentClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRight;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k, p;
  mwcpSizeAbsolute = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = - EPSILON;
  if (repMembrSize >= (2 * RF_minimumNodeSize)) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  responseClassCount = RF_classLevelSize[1];
  parentClassProp = uivector(1, responseClassCount);
  if (result) {
    for (p=1; p <= responseClassCount; p++) {
      parentClassProp[p] = 0;
    }
    for (k=1; k <= repMembrSize; k++) {
      parentClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[k]]]] ++;
    }
    k = 0;
    for (p=1; p <= responseClassCount; p++) {
      if(parentClassProp[p] > 0) {
        k ++;
      }
    }
    if (k > 1) {
      result = TRUE;
    }
    else {
      result = FALSE;
    }
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint   *nodeLeftIndex  = uivector(1, repMembrSize);
    uint   *nodeRghtIndex  = uivector(1, repMembrSize);
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize);
    for (i = 1; i <= actualCovariateCount; i++) {
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
          leftSize = rghtSize = 0;
          for (k=1; k <= repMembrSize; k++) {
            nodeLeftIndex[k] = nodeRghtIndex[k] = 0;
            if (localSplitIndicator[k] == LEFT) {
              nodeLeftIndex[++leftSize] = repMembrIndx[k];
            }
            else {
              nodeRghtIndex[++rghtSize] = repMembrIndx[k];
            }
          }
          for (p=1; p <= responseClassCount; p++) {
            parentClassProp[p] = leftClassProp[p] = rghtClassProp[p] = 0;
          }
          for (k=1; k <= repMembrSize; k++) {
            parentClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[k]]]] ++;
          }
          for (k=1; k <= leftSize; k++) {
            leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][nodeLeftIndex[k]]]] ++;
          }
          sumLeft = 0;
          sumRight = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft   += (double) upower(leftClassProp[p], 2);
            sumRight += (double) upower(parentClassProp[p] - leftClassProp[p], 2);
          }
          sumLeft = sumLeft / leftSize;
          sumRight = sumRight / (repMembrSize - leftSize);
          delta = sumLeft + sumRight;
          updateMaximumSplit(delta,
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
        }  
      }  
      unstackSplitVector(treeID,
                         permissibleSplitSize[i],
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
    }  
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            permissibleSplit, 
                            permissibleSplitSize);
    free_uivector (leftClassProp, 1, responseClassCount);
    free_uivector (rghtClassProp, 1, responseClassCount);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_uivector(nodeLeftIndex, 1, repMembrSize);
    free_uivector(nodeRghtIndex, 1, repMembrSize);
  }  
  free_uivector (parentClassProp, 1, responseClassCount);
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 deltaMax);
  return result;
}
