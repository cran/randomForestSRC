////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.2
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
#include    "regression.h"
#include     "splitClas.h"
char classificationSplit (uint    treeID, 
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
  uint responseClassCount;
  uint *parentClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
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
  responseClassCount = RF_classLevelSize[1];
  parentClassProp = uivector(1, responseClassCount);
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
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
          sumLeft = sumRght = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(parentClassProp[p] - leftClassProp[p], 2);
          }
          sumLeft = sumLeft / leftSize;
          sumRght = sumRght / rghtSize;
          delta = (sumLeft + sumRght) / repMembrSize;
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
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char classificationUwghtSplit (uint    treeID, 
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
  uint responseClassCount;
  uint *parentClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
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
  responseClassCount = RF_classLevelSize[1];
  parentClassProp = uivector(1, responseClassCount);
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
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
          sumLeft = sumRght = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += pow((double) leftClassProp[p] / (double) leftSize, 2.0);
            sumRght += pow((double) (parentClassProp[p] - leftClassProp[p]) / (double) rghtSize, 2.0);
          }
          delta = sumLeft + sumRght;
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
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char classificationHwghtSplit (uint    treeID, 
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
  uint responseClassCount;
  uint *parentClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
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
  responseClassCount = RF_classLevelSize[1];
  parentClassProp = uivector(1, responseClassCount);
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
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
          sumLeft = sumRght = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(parentClassProp[p] - leftClassProp[p], 2);
          }
          delta = 
            (sumLeft / (double) (upower(repMembrSize, 2))) + 
            (sumRght / (double) (upower(repMembrSize, 2))) -
            pow((double) leftSize / repMembrSize, 2.0) -
            pow((double) rghtSize / repMembrSize, 2.0) + 2.0;
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
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char mvClassificationSplit (uint    treeID, 
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
  uint leftSize, rghtSize;
  char   *purity;
  double *mean;
  double *variance;
  uint   *parentClassProp;
  uint   *leftClassProp;
  uint   *rghtClassProp;
  char    puritySummary;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k, p, r;
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
  purity   = cvector(1, RF_rSize); 
  mean     = dvector(1, RF_rSize);
  variance = dvector(1, RF_rSize);
  if (result) {
    puritySummary = FALSE;
    for (r = 1; r <= RF_rSize; r++) {
      purity[r] = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][r], &mean[r], &variance[r]);
      puritySummary = puritySummary | purity[r];
    }
    result = puritySummary;
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
          delta = 0.0;
          for (r = 1; r <= RF_rSize; r++) {
            if (purity[r]) {
              parentClassProp = uivector(1, RF_classLevelSize[r]);
              leftClassProp   = uivector(1, RF_classLevelSize[r]);
              rghtClassProp   = uivector(1, RF_classLevelSize[r]);
              for (p=1; p <= RF_classLevelSize[r]; p++) {
                parentClassProp[p] = leftClassProp[p] = rghtClassProp[p] = 0;
              }
              for (k=1; k <= repMembrSize; k++) {
                parentClassProp[RF_classLevelIndex[r][(uint) RF_response[treeID][r][repMembrIndx[k]]]] ++;
              }
              for (k=1; k <= leftSize; k++) {
                leftClassProp[RF_classLevelIndex[r][(uint) RF_response[treeID][r][nodeLeftIndex[k]]]] ++;
              }
              sumLeft = sumRght = 0;
              for (p=1; p <= RF_classLevelSize[r]; p++) {
                sumLeft   += (double) upower(leftClassProp[p], 2);
                sumRght += (double) upower(parentClassProp[p] - leftClassProp[p], 2);
              }
              sumLeft = sumLeft / leftSize;
              sumRght = sumRght / (repMembrSize - leftSize);
              delta += sumLeft + sumRght;
              free_uivector(parentClassProp, 1, RF_classLevelSize[r]);
              free_uivector(leftClassProp,   1, RF_classLevelSize[r]);
              free_uivector(rghtClassProp,   1, RF_classLevelSize[r]);
            }  
          }  
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
  free_cvector(purity,   1, RF_rSize); 
  free_dvector(mean,     1, RF_rSize);
  free_dvector(variance, 1, RF_rSize);
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
