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
                      double *splitStatistic,
                      char  **splitIndicator,
                      char  **omitMembrFlag,
                      char    multImpFlag) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  char **missMembrFlag;
  uint *nonMissMembrSize;
  uint **nonMissMembrIndx;
  uint   **indxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumRghtSave, sumLeftSqr, sumRghtSqr;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute       = 0;  
  sumLeft = sumRght      = 0;  
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
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint actualCovariateCount = stackAndSelectRandomCovariatesNew(treeID,
                                                                  parent, 
                                                                  repMembrIndx, 
                                                                  repMembrSize, 
                                                                  & randomCovariateIndex, 
                                                                  & permissibleSplit, 
                                                                  & permissibleSplitSize,
                                                                  & indxx,
                                                                  & missMembrFlag,
                                                                  & nonMissMembrSize,
                                                                  & nonMissMembrIndx,
                                                                  multImpFlag);
    for (i = 1; i <= actualCovariateCount; i++) {
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
      }
      sumRghtSave = 0.0;
      for (j = 1; j <= nonMissMembrSize[i]; j++) {
          sumRghtSave += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[i][indxx[i][j]]] ];
      }
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
        for (j = 1; j <= nonMissMembrSize[i]; j++) {
          localSplitIndicator[ nonMissMembrIndx[i][indxx[i][j]] ] = RIGHT;
        }
        sumRght      = sumRghtSave;
        sumLeft      = 0.0;
        leftSizeIter = 0;
        rghtSizeIter = nonMissMembrSize[i];
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        virtuallySplitNodeNew(treeID,
                              factorFlag,
                              mwcpSizeAbsolute,
                              randomCovariateIndex[i],
                              repMembrIndx,
                              repMembrSize,
                              nonMissMembrIndx[i],
                              nonMissMembrSize[i],
                              indxx[i],
                              permissibleSplitPtr,
                              j,
                              localSplitIndicator,
                              & leftSize,
                              missMembrFlag[i],                                 
                              priorMembrIter,
                              & currentMembrIter);
        rghtSize = nonMissMembrSize[i] - leftSize;
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            sumLeft = sumRght = 0.0;
            for (k = 1; k <= nonMissMembrSize[i]; k++) {
              if (localSplitIndicator[  nonMissMembrIndx[i][indxx[i][k]]  ] == LEFT) {
                sumLeft += RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[i][indxx[i][k]]]  ];
              }
              else {
                sumRght += RF_response[treeID][1][  repMembrIndx[nonMissMembrIndx[i][indxx[i][k]]]  ];
              }
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              sumLeft += RF_response[treeID][1][   repMembrIndx[nonMissMembrIndx[i][indxx[i][k]]]   ];
              sumRght -= RF_response[treeID][1][   repMembrIndx[nonMissMembrIndx[i][indxx[i][k]]]   ];
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          sumLeftSqr = pow (sumLeft, 2.0) / leftSize;
          sumRghtSqr = pow (sumRght, 2.0) / rghtSize;
          delta = sumLeftSqr + sumRghtSqr;
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
    unstackRandomCovariatesNew(treeID,
                               repMembrSize, 
                               randomCovariateIndex, 
                               actualCovariateCount,
                               permissibleSplit, 
                               permissibleSplitSize,
                               indxx,
                               missMembrFlag,
                               nonMissMembrSize,
                               nonMissMembrIndx);
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumRghtSave, sumRghtSqrSave, sumLeftSqr, sumRghtSqr;
  double leftTemp, rghtTemp, leftTempSqr, rghtTempSqr;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute       = 0;  
  sumLeft = sumRght      = 0;  
  sumLeftSqr             = 0;  
  sumRghtSqr             = 0;  
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
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
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
    sumRghtSave = sumRghtSqrSave = 0.0;
    for (j = 1; j <= repMembrSize; j++) {
      sumRghtSave += RF_response[treeID][1][repMembrIndx[j]];
      sumRghtSqrSave += pow(RF_response[treeID][1][repMembrIndx[j]], 2.0);
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
        sumRght      = sumRghtSave;
        sumRghtSqr   = sumRghtSqrSave; 
        sumLeft      = 0.0;
        sumLeftSqr   = 0.0;
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
            sumLeft = sumRght = 0.0;
            sumLeftSqr = sumRghtSqr = 0.0;
            for (k = 1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                sumLeft    += RF_response[treeID][1][repMembrIndx[k]];
                sumLeftSqr += pow(RF_response[treeID][1][repMembrIndx[k]], 2.0);
              }
              else {
                sumRght    += RF_response[treeID][1][repMembrIndx[k]];
                sumRghtSqr += pow(RF_response[treeID][1][repMembrIndx[k]], 2.0);
              }
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              sumLeft    += RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]];
              sumLeftSqr += pow(RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]], 2.0);
              sumRght    -= RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]];
              sumRghtSqr -= pow(RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]], 2.0);
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          leftTemp = pow(sumLeft, 2.0) / pow(leftSize, 2.0);
          rghtTemp = pow(sumRght, 2.0) / pow(rghtSize, 2.0);
          leftTempSqr = sumLeftSqr / leftSize;
          rghtTempSqr = sumRghtSqr / rghtSize;
          delta = leftTemp + rghtTemp - leftTempSqr - rghtTempSqr;
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double sumLeft, sumRght, sumRghtSave, sumRghtSqrSave, sumLeftSqr, sumRghtSqr;
  double leftTemp, rghtTemp, leftTempSqr, rghtTempSqr;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute = 0;        
  sumLeft = sumRght      = 0;  
  sumLeftSqr             = 0;  
  sumRghtSqr             = 0;  
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
  if (result) {
    result = getVariance(repMembrSize, repMembrIndx, RF_response[treeID][1], NULL, NULL);
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
    sumRghtSave = sumRghtSqrSave = 0.0;
    for (j = 1; j <= repMembrSize; j++) {
      sumRghtSave += RF_response[treeID][1][repMembrIndx[j]];
      sumRghtSqrSave += pow(RF_response[treeID][1][repMembrIndx[j]], 2.0);
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
        sumRght      = sumRghtSave;
        sumRghtSqr   = sumRghtSqrSave; 
        sumLeft      = 0.0;
        sumLeftSqr   = 0.0;
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
            sumLeft = sumRght = 0.0;
            sumLeftSqr = sumRghtSqr = 0.0;
            for (k = 1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                sumLeft    += RF_response[treeID][1][repMembrIndx[k]];
                sumLeftSqr += pow(RF_response[treeID][1][repMembrIndx[k]], 2.0);
              }
              else {
                sumRght    += RF_response[treeID][1][repMembrIndx[k]];
                sumRghtSqr += pow(RF_response[treeID][1][repMembrIndx[k]], 2.0);
              }
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              sumLeft    += RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]];
              sumLeftSqr += pow(RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]], 2.0);
              sumRght    -= RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]];
              sumRghtSqr -= pow(RF_response[treeID][1][repMembrIndx[repMembrIndxx[i][k]]], 2.0);
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          leftTemp = pow(sumLeft, 2.0) / pow (repMembrSize, 2.0);
          rghtTemp = pow(sumRght, 2.0) / pow (repMembrSize, 2.0);
          leftTempSqr = sumLeftSqr * leftSize / pow (repMembrSize, 2.0);
          rghtTempSqr = sumRghtSqr * rghtSize / pow (repMembrSize, 2.0);
          delta = leftTemp + rghtTemp - leftTempSqr - rghtTempSqr;
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
  char *localSplitIndicator;
  double delta, deltaMax;
  double *sumLeft, *sumRght, *sumRghtSave, sumLeftSqr, sumRghtSqr;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k, r;
  mwcpSizeAbsolute       = 0;  
  sumLeft = sumRght      = 0;  
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
      purity[r] = getVariance(repMembrSize, 
                              repMembrIndx, 
                              RF_response[treeID][r], 
                              &mean[r], 
                              &variance[r]);
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
    sumLeft      = dvector(1, RF_rSize);
    sumRght      = dvector(1, RF_rSize);
    sumRghtSave  = dvector(1, RF_rSize);
    for (r = 1; r <= RF_rSize; r++) {
      sumRghtSave[r] = 0.0;
      for (j = 1; j <= repMembrSize; j++) {
        sumRghtSave[r] += RF_response[treeID][r][repMembrIndx[j]] - mean[r];
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
          sumRght[r]      = sumRghtSave[r];
          sumLeft[r]      = 0.0;
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
                sumLeft[r] = sumRght[r] = 0.0;
                for (k = 1; k <= repMembrSize; k++) {
                  if (localSplitIndicator[k] == LEFT) {
                    sumLeft[r] += RF_response[treeID][r][repMembrIndx[k]] - mean[r];
                  }
                  else {
                    sumRght[r] += RF_response[treeID][r][repMembrIndx[k]] - mean[r];
                  }
                }
              }
              else {
                for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                  sumLeft[r] += RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                  sumRght[r] -= RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                }
              }
              sumLeftSqr = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
              sumRghtSqr = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
              delta += sumLeftSqr + sumRghtSqr;
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
    unstackRandomCovariates(treeID,
                            repMembrSize, 
                            randomCovariateIndex, 
                            actualCovariateCount,
                            permissibleSplit, 
                            permissibleSplitSize,
                            repMembrIndxx);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
    free_dvector(sumLeft, 1, RF_rSize);
    free_dvector(sumRght, 1, RF_rSize);
    free_dvector(sumRghtSave, 1, RF_rSize);
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
uint stackAndSelectRandomCovariatesNew(uint     treeID,
                                       Node     *parent,
                                       uint     *repMembrIndx,
                                       uint      repMembrSize,
                                       uint    **covariateIndex,
                                       double ***permissibleSplit,
                                       uint    **permissibleSplitSize,
                                       uint   ***indxx,
                                       char   ***missMembrFlag,
                                       uint    **nonMissMembrSize,
                                       uint   ***nonMissMembrIndx,
                                       char      multImpFlag) {
  uint i;
  uint actualCovariateCount;
  uint candidateCovariate;
  uint offset;
  uint indx;
  char mPredictorFlag;
  char splittable;
  *covariateIndex = uivector(1, RF_xSize);
  *permissibleSplit = dmatrix(1, RF_randomCovariateCount, 1, repMembrSize);
  *permissibleSplitSize = uivector(1, RF_randomCovariateCount);
  *indxx = (uint**) vvector(1, RF_randomCovariateCount);
  *missMembrFlag = (char **) vvector(1, RF_randomCovariateCount);
  *nonMissMembrSize = uivector(1, RF_randomCovariateCount);
  *nonMissMembrIndx = (uint**) vvector(1, RF_randomCovariateCount);
  char *randomSplitVector = cvector(1, RF_xSize);
  double *nonMissSplit = dvector(1, repMembrSize);
  if (repMembrSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid repMembrSize encountered in stackAndSelectRandomCovariates():  %10d", repMembrSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nrCopyVector(randomSplitVector, parent -> permissibleSplit, RF_xSize);
  for(i = 1; i <= RF_randomCovariateCount; i++) {
    (*covariateIndex)[i]    = 0;
    (*indxx)[i]             = NULL;
    (*missMembrFlag)[i]     = NULL;
    (*nonMissMembrSize)[i]  = 0;
    (*nonMissMembrIndx)[i]  = NULL;
  }
  actualCovariateCount =  1;
  candidateCovariate   = -1;
  while ((actualCovariateCount  <= RF_randomCovariateCount) && (candidateCovariate != 0)) {
    candidateCovariate = getSelectableElement(treeID, RF_xSize, randomSplitVector, RF_xWeight);
    if (candidateCovariate != 0) {
      (*indxx)[actualCovariateCount] = uivector(1, repMembrSize);
      (*missMembrFlag)[actualCovariateCount] = cvector(1, repMembrSize);
      (*nonMissMembrIndx)[actualCovariateCount] = uivector(1, repMembrSize);
      offset = RF_rSize + candidateCovariate;
      for (i = 1; i <= repMembrSize; i++) {
        mPredictorFlag = FALSE;
        if (RF_mRecordSize > 0) {
          if (RF_mRecordMap[repMembrIndx[i]] > 0) {
            if (RF_mpSign[offset][RF_mRecordMap[repMembrIndx[i]]] == 1) {
              if ((SPLIT_STATISTIC_OVERRIDE) && (!multImpFlag)) {
                mPredictorFlag = TRUE;
              }
            }
            if (RF_mpSign[1][RF_mRecordMap[repMembrIndx[i]]] == 1) {
              if ((SPLIT_STATISTIC_OVERRIDE) && (!multImpFlag)) {
                mPredictorFlag = TRUE;
              }
            }
          }
        }
        (*missMembrFlag)[actualCovariateCount][i] = mPredictorFlag;
        if (!mPredictorFlag) {
          (*nonMissMembrSize)[actualCovariateCount] ++;
          (*nonMissMembrIndx)[actualCovariateCount][(*nonMissMembrSize)[actualCovariateCount]] = i;
          nonMissSplit[(*nonMissMembrSize)[actualCovariateCount]] = RF_observation[treeID][candidateCovariate][repMembrIndx[i]];
        }
      }
      if ((*nonMissMembrSize)[actualCovariateCount] > 0) {
        splittable = TRUE;
      }
      else {
        splittable = FALSE;
      }
      if (splittable) {
        indexx((*nonMissMembrSize)[actualCovariateCount],
               nonMissSplit,
               (*indxx)[actualCovariateCount]);
        for (i = 1; i <= (*nonMissMembrSize)[actualCovariateCount]; i++) {
          indx = (*indxx)[actualCovariateCount][i];
          (*permissibleSplit)[actualCovariateCount][i] = nonMissSplit[indx];
        }
        (*permissibleSplitSize)[actualCovariateCount] = 1;
        for (i = 2; i <= (*nonMissMembrSize)[actualCovariateCount]; i++) {
          if ((*permissibleSplit)[actualCovariateCount][i] > (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]]) {
            (*permissibleSplitSize)[actualCovariateCount] ++;
            (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]] = (*permissibleSplit)[actualCovariateCount][i];
          }
        }
        if((*permissibleSplitSize)[actualCovariateCount] >= 2) {
          randomSplitVector[candidateCovariate] = ACTIVE;
          (*covariateIndex)[actualCovariateCount] = candidateCovariate;
          actualCovariateCount ++;
        }
        else {
          splittable = FALSE;
        }
      }  
      if (!splittable) {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        randomSplitVector[candidateCovariate] = FALSE;
        free_uivector((*indxx)[actualCovariateCount], 1, repMembrSize);
        (*indxx)[actualCovariateCount] = NULL;
        free_cvector((*missMembrFlag)[actualCovariateCount], 1, repMembrSize);
        (*missMembrFlag)[actualCovariateCount] = NULL;
        (*nonMissMembrSize)[actualCovariateCount] = 0;
        free_uivector((*nonMissMembrIndx)[actualCovariateCount], 1, repMembrSize);
        (*nonMissMembrIndx)[actualCovariateCount] = NULL;
      }
    }  
  }
  actualCovariateCount --;
  free_cvector(randomSplitVector, 1, RF_xSize);
  free_dvector(nonMissSplit, 1, repMembrSize);
  return actualCovariateCount;
}
void unstackRandomCovariatesNew(uint     treeID,
                                uint     repMembrSize, 
                                uint    *covariateIndex,
                                uint     covariateCount,
                                double **permissibleSplit,
                                uint    *permissibleSplitSize,
                                uint   **indxx,
                                char   **missMembrFlag,
                                uint    *nonMissMembrSize,
                                uint   **nonMissMembrIndx) {
  uint i;
  free_uivector(covariateIndex, 1, RF_xSize);
  free_dmatrix(permissibleSplit, 1, RF_randomCovariateCount, 1, repMembrSize);
  free_uivector(permissibleSplitSize, 1, RF_randomCovariateCount);
  for (i = 1; i <= covariateCount; i++) {
    if (indxx[i] != NULL) {
      free_uivector(indxx[i], 1, repMembrSize);
    }
    if (missMembrFlag[i] != NULL) {
      free_cvector(missMembrFlag[i], 1, repMembrSize);
    }
    if(nonMissMembrIndx[i] != NULL) {
      free_uivector(nonMissMembrIndx[i], 1, repMembrSize);
    }
  }
  free_vvector(indxx, 1, RF_randomCovariateCount);
  free_vvector(missMembrFlag, 1, RF_randomCovariateCount);
  free_uivector(nonMissMembrSize, 1, RF_randomCovariateCount);
  free_vvector(nonMissMembrIndx, 1, RF_randomCovariateCount);
}
uint virtuallySplitNodeNew(uint  treeID,
                           char  factorFlag,
                           uint  mwcpSizeAbsolute,
                           uint  randomCovariate,
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *nonMissMembrIndx,
                           uint  nonMissMembrSize,
                           uint *indxx,
                           void *permissibleSplitPtr,
                           uint  offset,
                           char *localSplitIndicator,
                           uint *leftSize,
                           char *missMembrFlag,
                           uint  priorMembrIter,
                           uint *currentMembrIter) {
  char daughterFlag;
  char iterFlag;
  iterFlag = TRUE;
  *currentMembrIter = priorMembrIter;
  while (iterFlag) {
    (*currentMembrIter) ++;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint)  RF_observation[treeID][randomCovariate][    repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]     ], 
                                   (uint*) permissibleSplitPtr + ((offset - 1) * mwcpSizeAbsolute));
      if ((*currentMembrIter) == nonMissMembrSize) {
        iterFlag = FALSE;
      }
    }
    else {
      if (RF_observation[treeID][randomCovariate][   repMembrIndx[nonMissMembrIndx[indxx[*currentMembrIter]]]    ] <= ((double*) permissibleSplitPtr)[offset]) {
        daughterFlag = LEFT;
      }
      else {
        daughterFlag = RIGHT;
        iterFlag = FALSE;
      }
    }
    localSplitIndicator[     nonMissMembrIndx[indxx[*currentMembrIter]]   ] = daughterFlag;
    if (daughterFlag == LEFT) {
      (*leftSize) ++;
    }  
    else {
    }
  }  
  return (*leftSize); 
 }
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexxx(unsigned int n, double *arr, char *omit, unsigned int *indx) {
  unsigned int i, j, k, l;
  unsigned int indxt, itemp, ir;
  unsigned int *istack, jstack;
  double a;
  char iFlag, jFlag;
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (j = 1; j <= n; j++) {
    indx[j] = j;
  }
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i = j-1; i >= l; i--) {
          if (omit[indx[i]]) {
            break;
          }
          else {
            if (omit[indxt]) {
            }
            else {
              if (arr[indx[i]] <= a) {
                break;
              }
              else {
              }
            }
          }
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1])
      if (omit[indx[ir]]) {
        if (!omit[indx[l]]) {
          SWAP(indx[l], indx[ir])
        }
        else {
        }
      }
      else {
        if (omit[indx[l]]) {
        }
        else {
          if (arr[indx[l]] > arr[indx[ir]]) {
            SWAP(indx[l], indx[ir])
          }
          else {
          }
        }
      }
      if (omit[indx[ir]]) {
        if (!omit[indx[l+1]]) {
          SWAP(indx[l+1], indx[ir])
        }
        else {
        }
      }
      else {
        if (omit[indx[l+1]]) {
        }
        else {
          if (arr[indx[l+1]] > arr[indx[ir]]) {
            SWAP(indx[l+1], indx[ir])
          }
          else {
          }
        }
      }
      if (omit[indx[l+1]]) {
        if (!omit[indx[l]]) {
          SWAP(indx[l], indx[l+1])
        }
        else {
        }
      }
      else {
        if (omit[indx[l]]) {
        }
        else {
          if (arr[indx[l]] > arr[indx[l+1]]) {
            SWAP(indx[l], indx[l+1])
          }
          else {
          }
        }
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        iFlag = TRUE;
        while(iFlag) {
          i++;
          if (omit[indx[i]]) {
            if (!omit[indxt]) {
              iFlag = TRUE;
            }
            else {
              iFlag = FALSE;
            }
          }
          else {
            if (omit[indxt]) {
              iFlag = FALSE;
            }
            else {
              if (arr[indx[i]] < a) {
                iFlag = TRUE;
              }
              else {
                iFlag = FALSE;
              }
            }
          }
        }
        jFlag = TRUE;
        while(jFlag) {
          j++;
          if (omit[indxt]) {
            if (!omit[indx[j]]) {
              jFlag = TRUE;
            }
            else {
              jFlag = FALSE;
            }
          }
          else {
            if (omit[indx[j]]) {
              jFlag = FALSE;
            }
            else {
              if (arr[indx[j]] > a) {
                jFlag = TRUE;
              }
              else {
                jFlag = FALSE;
              }
            }
          }
        }
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free_uivector(istack, 1, NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
