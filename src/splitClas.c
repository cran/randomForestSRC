////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.4
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
                          double *splitStatistic,
                          char  **splitIndicator,
                          char  **omitMembrFlag,
                          char    multImpFlag) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint responseClassCount;
  uint *parentClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k, p;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
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
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    for (p=1; p <= responseClassCount; p++) {
      parentClassProp[p] = 0;
    }
    for (j = 1; j <= repMembrSize; j++) {
      parentClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[j]]]] ++;
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = RIGHT;
        }
        for (p = 1; p <= responseClassCount; p++) {
          rghtClassProp[p] = parentClassProp[p];
          leftClassProp[p] = 0;
        }
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           repMembrIndx,
                           repMembrIndxx[i],
                           repMembrSize,
                           permissibleSplitPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        rghtSize = repMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (p=1; p <= responseClassCount; p++) {
              leftClassProp[p] = 0;
            }
            for (k=1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[k]]]] ++;
              }
            }
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p] - leftClassProp[p];
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] ++;
              rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] --;
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          sumLeft = sumRght = 0.0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(rghtClassProp[p], 2);
          }
          sumLeftSqr = sumLeft / leftSize;
          sumRghtSqr  = sumRght / rghtSize;
          delta = (sumLeftSqr + sumRghtSqr) / repMembrSize;
          updateMaximumSplit(treeID,
                             delta,
                             randomCovariateIndex[i],
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
        }  
        if (factorFlag == FALSE) {
          if (rghtSize  < RF_minimumNodeSize) {
            j = splitLength;
          }
          else {
            priorMembrIter = currentMembrIter - 1;
          }
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
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    free_uivector (leftClassProp, 1, responseClassCount);
    free_uivector (rghtClassProp, 1, responseClassCount);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
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
                               double *splitStatistic,
                               char  **splitIndicator,
                               char  **omitMembrFlag,
                               char    multImpFlag) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
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
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
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
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    for (p=1; p <= responseClassCount; p++) {
      parentClassProp[p] = 0;
    }
    for (j = 1; j <= repMembrSize; j++) {
      parentClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[j]]]] ++;
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = RIGHT;
        }
        for (p = 1; p <= responseClassCount; p++) {
          rghtClassProp[p] = parentClassProp[p];
          leftClassProp[p] = 0;
        }
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           repMembrIndx,
                           repMembrIndxx[i],
                           repMembrSize,
                           permissibleSplitPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        rghtSize = repMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (p=1; p <= responseClassCount; p++) {
              leftClassProp[p] = 0;
            }
            for (k=1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[k]]]] ++;
              }
            }
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p] - leftClassProp[p];
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] ++;
              rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] --;
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          sumLeft = sumRght = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += pow((double) leftClassProp[p] / (double) leftSize, 2.0);
            sumRght += pow((double) (rghtClassProp[p]) / (double) rghtSize, 2.0);
          }
          delta = sumLeft + sumRght;
          updateMaximumSplit(treeID,
                             delta,
                             randomCovariateIndex[i],
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
        }  
        if (factorFlag == FALSE) {
          if (rghtSize  < RF_minimumNodeSize) {
            j = splitLength;
          }
          else {
            priorMembrIter = currentMembrIter - 1;
          }
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
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    free_uivector (leftClassProp, 1, responseClassCount);
    free_uivector (rghtClassProp, 1, responseClassCount);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
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
                               double *splitStatistic,
                               char  **splitIndicator,
                               char  **omitMembrFlag,
                               char    multImpFlag) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
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
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
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
    uint *leftClassProp   = uivector(1, responseClassCount);
    uint *rghtClassProp   = uivector(1, responseClassCount);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    for (p=1; p <= responseClassCount; p++) {
      parentClassProp[p] = 0;
    }
    for (j = 1; j <= repMembrSize; j++) {
      parentClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[j]]]] ++;
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = RIGHT;
        }
        for (p = 1; p <= responseClassCount; p++) {
          rghtClassProp[p] = parentClassProp[p];
          leftClassProp[p] = 0;
        }
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           repMembrIndx,
                           repMembrIndxx[i],
                           repMembrSize,
                           permissibleSplitPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        rghtSize = repMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (p=1; p <= responseClassCount; p++) {
              leftClassProp[p] = 0;
            }
            for (k=1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[k]]]] ++;
              }
            }
            for (p=1; p <= responseClassCount; p++) {
              rghtClassProp[p] = parentClassProp[p] - leftClassProp[p];
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              leftClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] ++;
              rghtClassProp[RF_classLevelIndex[1][(uint) RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]]]] --;
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          sumLeft = sumRght = 0;
          for (p=1; p <= responseClassCount; p++) {
            sumLeft += (double) upower(leftClassProp[p], 2);
            sumRght += (double) upower(rghtClassProp[p], 2);
          }
          delta = 
            (sumLeft / (double) (upower(repMembrSize, 2))) + 
            (sumRght / (double) (upower(repMembrSize, 2))) -
            pow((double) leftSize / repMembrSize, 2.0) -
            pow((double) rghtSize / repMembrSize, 2.0) + 2.0;
          updateMaximumSplit(treeID,
                             delta,
                             randomCovariateIndex[i],
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
        }  
        if (factorFlag == FALSE) {
          if (rghtSize  < RF_minimumNodeSize) {
            j = splitLength;
          }
          else {
            priorMembrIter = currentMembrIter - 1;
          }
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
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    free_uivector (leftClassProp, 1, responseClassCount);
    free_uivector (rghtClassProp, 1, responseClassCount);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
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
                            double *splitStatistic,
                            char  **splitIndicator,
                            char  **omitMembrFlag,
                            char    multImpFlag) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  char   *purity;
  double *mean;
  double *variance;
  char    puritySummary;
  uint  **parentClassProp;
  uint  **leftClassProp;
  uint  **rghtClassProp;
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
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
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
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    parentClassProp = (uint**) vvector(1, RF_rSize);
    leftClassProp   = (uint**) vvector(1, RF_rSize);
    rghtClassProp   = (uint**) vvector(1, RF_rSize);
    for (r = 1; r <= RF_rSize; r++) {
      parentClassProp[r] = uivector(1, RF_classLevelSize[r]);
      leftClassProp[r]   = uivector(1, RF_classLevelSize[r]);
      rghtClassProp[r]   = uivector(1, RF_classLevelSize[r]);
      for (p=1; p <= RF_classLevelSize[r]; p++) {
        parentClassProp[r][p] = leftClassProp[r][p] = rghtClassProp[r][p] = 0;
      }
      for (j = 1; j <= repMembrSize; j++) {
        parentClassProp[r][RF_classLevelIndex[r][(uint) RF_response[treeID][r][repMembrIndx[j]]]] ++;
      }
    }
    for (i = 1; i <= actualCovariateCount; i++) {
      leftSize = 0;
      priorMembrIter = 0;
      splitLength = stackAndConstructSplitVector(treeID,
                                                 repMembrSize,
                                                 randomCovariateIndex[i], 
                                                 permissibleSplit[i], 
                                                 permissibleSplitSize[i],
                                                 & factorFlag,
                                                 & deterministicSplitFlag,
                                                 & mwcpSizeAbsolute,
                                                 & permissibleSplitPtr);
      if (factorFlag == FALSE) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = RIGHT;
        }
        for (r = 1; r <= RF_rSize; r++) {
          for (p=1; p <= RF_classLevelSize[r]; p++) {
            rghtClassProp[r][p] = parentClassProp[r][p];
            leftClassProp[r][p] = 0;
          }
        }
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNode(treeID,
                           factorFlag,
                           mwcpSizeAbsolute,
                           randomCovariateIndex[i],
                           repMembrIndx,
                           repMembrIndxx[i],
                           repMembrSize,
                           permissibleSplitPtr,
                           j,
                           localSplitIndicator,
                           & leftSize,
                           priorMembrIter,
                           & currentMembrIter);
        rghtSize = repMembrSize - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          delta = 0.0;
          for (r = 1; r <= RF_rSize; r++) {
            if (purity[r]) {
              if (factorFlag == TRUE) {
                for (p=1; p <= RF_classLevelSize[r]; p++) {
                  leftClassProp[r][p] = 0;
                }
                for (k=1; k <= repMembrSize; k++) {
                  if (localSplitIndicator[k] == LEFT) {
                    leftClassProp[r][RF_classLevelIndex[r][(uint) RF_response[treeID][r][repMembrIndx[k]]]] ++;
                  }
                }
                for (p=1; p <= RF_classLevelSize[r]; p++) {
                  rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                }
              }
              else {
                for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                  leftClassProp[r][RF_classLevelIndex[r][(uint) RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]]]] ++;
                  rghtClassProp[r][RF_classLevelIndex[r][(uint) RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]]]] --;
                }
              }
              sumLeft = sumRght = 0.0;
              for (p=1; p <= RF_classLevelSize[r]; p++) {
                sumLeft   += (double) upower(leftClassProp[r][p], 2);
                sumRght += (double) upower(rghtClassProp[r][p], 2);
              }
              sumLeft = sumLeft / leftSize;
              sumRght = sumRght / rghtSize;
              delta += sumLeft + sumRght; 
            }  
          }  
          if (factorFlag == FALSE) {
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          updateMaximumSplit(treeID,
                             delta,
                             randomCovariateIndex[i],
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
        }  
        if (factorFlag == FALSE) {
          if (rghtSize  < RF_minimumNodeSize) {
            j = splitLength;
          }
          else {
            priorMembrIter = currentMembrIter - 1;
          }
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
    for (r = 1; r <= RF_rSize; r++) {
      free_uivector(parentClassProp[r], 1, RF_classLevelSize[r]);
      free_uivector(leftClassProp[r], 1, RF_classLevelSize[r]);
      free_uivector(rghtClassProp[r], 1, RF_classLevelSize[r]);
    }
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_vvector(parentClassProp, 1, RF_rSize);
    free_vvector(leftClassProp,   1, RF_rSize);
    free_vvector(rghtClassProp,   1, RF_rSize);
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
