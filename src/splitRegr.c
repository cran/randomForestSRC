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
#include     "splitRegr.h"
char regressionSplit (uint    treeID, 
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
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
          sumRght = 0;
          for (k=1; k <= rghtSize; k++) {
            sumRght += RF_response[treeID][1][nodeRghtIndex[k]];
          }
          sumLeft = pow (sumLeft, 2.0) / leftSize;
          sumRght = pow (sumRght, 2.0) / rghtSize;
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
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_uivector(nodeLeftIndex, 1, repMembrSize);
    free_uivector(nodeRghtIndex, 1, repMembrSize);
  }  
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char regressionUwghtSplit (uint    treeID, 
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
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
          sumLeft = sumLeftSqr = 0.0;
          for (k=1; k <= leftSize; k++) {
            sumLeft += RF_response[treeID][1][nodeLeftIndex[k]];
            sumLeftSqr += pow(RF_response[treeID][1][nodeLeftIndex[k]], 2.0);
          }
          sumRght = sumRghtSqr = 0.0;
          for (k=1; k <= rghtSize; k++) {
            sumRght += RF_response[treeID][1][nodeRghtIndex[k]];
            sumRghtSqr += pow(RF_response[treeID][1][nodeRghtIndex[k]], 2.0);
          }
          sumLeft = pow(sumLeft, 2.0) / pow(leftSize, 2.0);
          sumRght = pow(sumRght, 2.0) / pow(rghtSize, 2.0);
          sumLeftSqr = sumLeftSqr / leftSize;
          sumRghtSqr = sumRghtSqr / rghtSize;
          delta = sumLeft + sumRght - sumLeftSqr - sumRghtSqr;
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
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char regressionHwghtSplit (uint    treeID, 
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumLeftSqr, sumRghtSqr;
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
          sumLeft = sumLeftSqr = 0;
          for (k=1; k <= leftSize; k++) {
            sumLeft += RF_response[treeID][1][nodeLeftIndex[k]];
            sumLeftSqr += pow(RF_response[treeID][1][nodeLeftIndex[k]], 2.0);
          }
          sumRght = sumRghtSqr = 0;
          for (k=1; k <= rghtSize; k++) {
            sumRght += RF_response[treeID][1][nodeRghtIndex[k]];
            sumRghtSqr += pow(RF_response[treeID][1][nodeRghtIndex[k]], 2.0);
          }
          sumLeft = pow(sumLeft, 2.0) / pow (repMembrSize, 2.0);
          sumRght = pow(sumRght, 2.0) / pow (repMembrSize, 2.0);
          sumLeftSqr = sumLeftSqr * leftSize / pow (repMembrSize, 2.0);
          sumRghtSqr = sumRghtSqr * rghtSize / pow (repMembrSize, 2.0);
          delta = sumLeft + sumRght - sumLeftSqr - sumRghtSqr;
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
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char mvRegressionSplit (uint    treeID, 
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
  char    puritySummary;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght;
  double sqrLeft, sqrRght;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k, r;
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
              sumLeft = 0;
              for (k=1; k <= leftSize; k++) {
                sumLeft += RF_response[treeID][r][nodeLeftIndex[k]] - mean[r];
              }
              sumRght = 0;
              for (k=1; k <= rghtSize; k++) {
                sumRght += RF_response[treeID][r][repMembrIndx[k]] - mean[r];
              }
              sqrLeft = pow (sumLeft, 2.0) / (leftSize * variance[r]);
              sqrRght = pow (sumRght, 2.0) / (rghtSize * variance[r]);
              delta += sqrLeft + sqrRght;
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
