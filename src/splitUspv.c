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
char unsupervisedSplit (uint    treeID, 
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
  char   *purity;
  double *mean;
  double *variance;
  char    puritySummary;
  uint  **parentClassProp;
  uint  **leftClassProp;
  uint  **rghtClassProp;
  uint **pseudoResponse;
  uint *pseudoResponseClassSize;
  char *permissibleResponseFlag;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double *sumLeft, *sumRght, *sumRghtSave, sumLeftSqr, sumRghtSqr;
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
  sumLeft                = 0;  
  sumRght                = 0;  
  sumRghtSave            = 0;  
  parentClassProp        = NULL;   
  leftClassProp          = NULL;   
  rghtClassProp          = NULL;   
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
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
    if (actualCovariateCount > 0) {
      pseudoResponse = uimatrix(1, actualCovariateCount, 1, RF_randomResponseCount);
      permissibleResponseFlag = cvector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        permissibleResponseFlag[i] = TRUE;
      }
      for (i = 1; i <= actualCovariateCount; i++) {
        permissibleResponseFlag[randomCovariateIndex[i]] = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++) {
          pseudoResponse[i][r] = getSelectableElement(treeID, RF_xSize, permissibleResponseFlag, NULL);
        }
        permissibleResponseFlag[randomCovariateIndex[i]] = TRUE;
      }
      free_cvector(permissibleResponseFlag, 1, RF_xSize);
      purity   = cvector(1, RF_randomResponseCount); 
      mean     = dvector(1, RF_randomResponseCount);
      variance = dvector(1, RF_randomResponseCount);
      parentClassProp = (uint**) vvector(1, RF_randomResponseCount);
      leftClassProp   = (uint**) vvector(1, RF_randomResponseCount);
      rghtClassProp   = (uint**) vvector(1, RF_randomResponseCount);
      sumLeft      = dvector(1, RF_randomResponseCount);
      sumRght      = dvector(1, RF_randomResponseCount);
      sumRghtSave  = dvector(1, RF_randomResponseCount);
      pseudoResponseClassSize = uivector(1, RF_randomResponseCount);
      for (i = 1; i <= actualCovariateCount; i++) {
        puritySummary = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++) {
          purity[r] = getVariance(repMembrSize, 
                                  repMembrIndx, 
                                  RF_observation[treeID][pseudoResponse[i][r]],
                                  &mean[r], 
                                  &variance[r]);
          puritySummary = puritySummary | purity[r];
        }
        if (puritySummary) {
          for (r = 1; r <= RF_randomResponseCount; r++) {
            pseudoResponseClassSize[r] = 0;
            parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
            sumLeft[r] = sumRght[r] = sumRghtSave[r] = 0;
          }
          for (r = 1; r <= RF_randomResponseCount; r++) {
            if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
              pseudoResponseClassSize[r] = RF_xFactorSize[RF_xFactorMap[pseudoResponse[i][r]]];
              parentClassProp[r] = uivector(1, pseudoResponseClassSize[r]);
              leftClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              rghtClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                parentClassProp[r][p] = 0;
              }
              for (j = 1; j <= repMembrSize; j++) {
                parentClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[j]]] ++;
              }
            }
            else {
              sumRghtSave[r] = 0.0;
              for (j = 1; j <= repMembrSize; j++) {
                sumRghtSave[r] += RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[j]] - mean[r];
              }
            }
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
            for (j = 1; j <= repMembrSize; j++) {
              localSplitIndicator[j] = RIGHT;
            }
            for (r = 1; r <= RF_randomResponseCount; r++) {
              if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
                for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                  rghtClassProp[r][p] = parentClassProp[r][p];
                  leftClassProp[r][p] = 0;
                }
              }
              else {
                sumRght[r] = sumRghtSave[r];
                sumLeft[r] = 0.0;
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
              for (r = 1; r <= RF_randomResponseCount; r++) {
                if (purity[r]) {
                  if (factorFlag == TRUE) {
                    if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        leftClassProp[r][p] = 0;
                      }
                      for (k=1; k <= repMembrSize; k++) {
                        if (localSplitIndicator[k] == LEFT) {
                          leftClassProp[r][(uint) RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[k]]]  ++;
                        }
                      }
                      for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                        rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                      }
                    }
                    else {
                      sumLeft[r] = sumRght[r] = 0.0;
                      for (k = 1; k <= repMembrSize; k++) {
                        if (localSplitIndicator[k] == LEFT) {
                          sumLeft[r] += RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[k]] - mean[r];
                        }
                        else {
                          sumRght[r] += RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[k]] - mean[r];
                        }
                      }
                    }
                  }
                  else {
                    for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                      if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
                        leftClassProp[r][(uint) RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[repMembrIndxx[i][k]]]] ++;
                        rghtClassProp[r][(uint) RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[repMembrIndxx[i][k]]]] --;
                      }
                      else {
                        sumLeft[r] += RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                        sumRght[r] -= RF_observation[treeID][pseudoResponse[i][r]][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                      }
                    }
                  }  
                  if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
                    sumLeft[1] = sumRght[1] = 0;
                    for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                      sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                      sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                    }
                    sumLeftSqr = sumLeft[1] / leftSize;
                    sumRghtSqr  = sumRght[1] / rghtSize;
                  }
                  else {
                    sumLeftSqr = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
                    sumRghtSqr = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
                  }
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
          for (r = 1; r <= RF_randomResponseCount; r++) {
            if (strcmp(RF_xType[pseudoResponse[i][r]], "C") == 0) {
              free_uivector (parentClassProp[r], 1, pseudoResponseClassSize[r]);
              free_uivector (leftClassProp[r], 1, pseudoResponseClassSize[r]);
              free_uivector (rghtClassProp[r], 1, pseudoResponseClassSize[r]);
            }
            else {
            }
          }
        }  
      }  
      free_vvector(parentClassProp, 1, RF_randomResponseCount);
      free_vvector(leftClassProp,   1, RF_randomResponseCount);
      free_vvector(rghtClassProp,   1, RF_randomResponseCount);
      free_dvector(sumLeft,     1, RF_randomResponseCount);
      free_dvector(sumRght,     1, RF_randomResponseCount);
      free_dvector(sumRghtSave, 1, RF_randomResponseCount);
      free_uivector(pseudoResponseClassSize, 1, RF_randomResponseCount);
      free_cvector(purity,   1, RF_randomResponseCount); 
      free_dvector(mean,     1, RF_randomResponseCount);
      free_dvector(variance, 1, RF_randomResponseCount);
      free_uimatrix(pseudoResponse, 1, actualCovariateCount, 1, RF_randomResponseCount);
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
char multivariateSplit (uint    treeID, 
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
  char   *purity;
  double *mean;
  double *variance;
  char    puritySummary;
  uint  **parentClassProp;
  uint  **leftClassProp;
  uint  **rghtClassProp;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  double delta, deltaMax;
  double *sumLeft, *sumRght, *sumRghtSave, sumLeftSqr, sumRghtSqr;
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
  sumLeft                = 0;  
  sumRght                = 0;  
  sumRghtSave            = 0;  
  parentClassProp        = NULL;   
  leftClassProp          = NULL;   
  rghtClassProp          = NULL;   
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
    if (puritySummary) {
      stackSplitIndicator(repMembrSize, & localSplitIndicator);
      uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                                 parent, 
                                                                 repMembrIndx, 
                                                                 repMembrSize, 
                                                                 & randomCovariateIndex, 
                                                                 & permissibleSplit, 
                                                                 & permissibleSplitSize,
                                                                 & repMembrIndxx);
      if (TRUE) {
        parentClassProp = (uint**) vvector(1, RF_rSize);
        leftClassProp   = (uint**) vvector(1, RF_rSize);
        rghtClassProp   = (uint**) vvector(1, RF_rSize);
        sumLeft      = dvector(1, RF_rSize);
        sumRght      = dvector(1, RF_rSize);
        sumRghtSave  = dvector(1, RF_rSize);
        for (r = 1; r <= RF_rSize; r++) {
          parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
          sumLeft[r] = sumRght[r] = sumRghtSave[r] = 0;
        }
        for (r = 1; r <= RF_rSize; r++) {
          if (strcmp(RF_rType[r], "C") == 0) { 
            parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
            leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
            rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
            for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
              parentClassProp[r][p] = leftClassProp[r][p] = rghtClassProp[r][p] = 0;
            }
            for (j = 1; j <= repMembrSize; j++) {
              parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][repMembrIndx[j]]]] ++;
            }
          }
          else {
            sumRghtSave[r] = 0.0;
            for (j = 1; j <= repMembrSize; j++) {
              sumRghtSave[r] += RF_response[treeID][r][repMembrIndx[j]] - mean[r];
            }
          }
        }  
        for (i = 1; i <= actualCovariateCount; i++) {
          if (TRUE) {
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
                if (strcmp(RF_rType[r], "C") == 0) { 
                  for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                    rghtClassProp[r][p] = parentClassProp[r][p];
                    leftClassProp[r][p] = 0;
                  }
                }
                else {
                  sumRght[r] = sumRghtSave[r];
                  sumLeft[r] = 0.0;
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
                      if (strcmp(RF_rType[r], "C") == 0) { 
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                          leftClassProp[r][p] = 0;
                        }
                        for (k=1; k <= repMembrSize; k++) {
                          if (localSplitIndicator[k] == LEFT) {
                            leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][repMembrIndx[k]]]] ++;
                          }
                        }
                        for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                          rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                        }
                      }
                      else {
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
                    }
                    else {
                      for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                        if (strcmp(RF_rType[r], "C") == 0) { 
                          leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]]]] ++;
                          rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]]]] --;
                        }
                        else {
                          sumLeft[r] += RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                          sumRght[r] -= RF_response[treeID][r][repMembrIndx[repMembrIndxx[i][k]]] - mean[r];
                        }
                      }
                    }  
                    if (strcmp(RF_rType[r], "C") == 0) {
                      sumLeft[1] = sumRght[1] = 0;
                      for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                        sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                        sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                      }
                      sumLeftSqr = sumLeft[1] / leftSize;
                      sumRghtSqr  = sumRght[1] / rghtSize;
                    }
                    else {
                      sumLeftSqr = pow (sumLeft[r], 2.0) / (leftSize * variance[r]);
                      sumRghtSqr = pow (sumRght[r], 2.0) / (rghtSize * variance[r]);
                    }
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
        }  
        for (r = 1; r <= RF_rSize; r++) {
          if (strcmp(RF_rType[r], "C") == 0) {
            free_uivector (parentClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
            free_uivector (leftClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
            free_uivector (rghtClassProp[r], 1, RF_classLevelSize[RF_rFactorMap[r]]);
          }
          else {
          }
        }
        free_vvector(parentClassProp, 1, RF_rSize);
        free_vvector(leftClassProp,   1, RF_rSize);
        free_vvector(rghtClassProp,   1, RF_rSize);
        free_dvector(sumLeft, 1, RF_rSize);
        free_dvector(sumRght, 1, RF_rSize);
        free_dvector(sumRghtSave, 1, RF_rSize);
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
