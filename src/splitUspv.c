////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.6.0
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
#include   "splitCustom.h"
void getMemberCountOnly(uint treeID) {
  Terminal *parent;
  uint leaf, i;
  uint *membershipIndex;
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    membershipIndex = RF_identityMembershipIndex;
  }
  else {
    membershipIndex = RF_bootMembershipIndex[treeID];
  }
  for (leaf = 1; leaf <= RF_tLeafCount[treeID]; leaf++) {
    parent = RF_tTermList[treeID][leaf];
    parent -> membrCount = 0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_tTermMembership[treeID][membershipIndex[i]] == parent) {
        parent -> membrCount ++;
      }
    }
    if (parent -> membrCount == 0) {
      if (!(RF_opt & OPT_OUTC_TYPE)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Zero member count encountered in node during getMemberCountOnly():  %10d", leaf);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }  
}
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
                        char    multImpFlag) {
  char   *randomCovariateFlag;
  uint   *randomCovariateIndex;
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
  uint leftSizeIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  double delta, deltaMax;
  uint   deltaNorm;
  uint i, j, k, p, r;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             NULL,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if (result) {
    char   *impurity   = cvector(1, RF_randomResponseCount);
    double *mean       = dvector(1, RF_randomResponseCount);
    double *variance   = dvector(1, RF_randomResponseCount);
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
    uint  **parentClassProp = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    uint  **leftClassProp   = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    uint  **rghtClassProp   = (uint **) new_vvector(1, RF_randomResponseCount, NRUTIL_UPTR);
    double *sumLeft      = dvector(1, RF_randomResponseCount);
    double *sumRght      = dvector(1, RF_randomResponseCount);
    double *sumRghtSave  = dvector(1, RF_randomResponseCount);
    uint *pseudoResponseClassSize = uivector(1, RF_randomResponseCount);
    uint *pseudoResponse = uivector(1, RF_randomResponseCount);
    char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_randomResponseCount, NRUTIL_CPTR);
    uint  *secondNonMissMembrSize =              uivector(1, RF_randomResponseCount);
    uint  *secondNonMissMembrLeftSize =          uivector(1, RF_randomResponseCount);
    uint  *secondNonMissMembrRghtSize =          uivector(1, RF_randomResponseCount);
    char  *tempNonMissMembrFlag = 0;
    uint  *tempNonMissMembrIndx;
    char   mResponseFlag;
    uint   localIndex, localSize;
    char    nonMissImpuritySummary;
    double sumLeftSqr, sumRghtSqr;
    uint actualCovariateCount = 0;
    uint candidateCovariateCount = 0;
    while (selectRandomCovariates(treeID,
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
      uint *pseudoResponseIndex = uivector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        pseudoResponseIndex[i] = i;
      }
      pseudoResponseIndex[covariate] = pseudoResponseIndex[RF_xSize];
      localSize = RF_xSize - 1;
      for (r = 1; r <= RF_randomResponseCount; r++) {
        pseudoResponse[r] = sampleWithoutReplacement(treeID, pseudoResponseIndex, RF_xSize, & localIndex);
        pseudoResponseIndex[localIndex] = pseudoResponseIndex[localSize];
        localSize --;
      }
      free_uivector(pseudoResponseIndex, 1, RF_xSize);
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
        tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
        for (k = 1; k <= nonMissMembrSize; k++) {
          tempNonMissMembrFlag[k] = TRUE;
          tempNonMissMembrIndx[k] = k;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
          secondNonMissMembrSize[r] = nonMissMembrSize;
        }
        nonMissImpuritySummary = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    secondNonMissMembrSize[r],
                                    tempNonMissMembrIndx,
                                    RF_observation[treeID][pseudoResponse[r]],
                                    &mean[r],
                                    &variance[r]);
          nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
          secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
        }
        free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
      }
      else {
        tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
        nonMissImpuritySummary = FALSE;
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
          j = 0;
          for (k = 1; k <= nonMissMembrSize; k++) {
            mResponseFlag = FALSE;
            if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
              if (RF_mpSign[pseudoResponse[r]][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                mResponseFlag = TRUE;
              }
            }
            if (!mResponseFlag) {
              j ++;
              tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
              secondNonMissMembrFlag[r][k] = TRUE;
            }
            else {
              secondNonMissMembrFlag[r][k] = FALSE;
            }
          }  
          secondNonMissMembrSize[r] = j;
          impurity[r] = getVariance(repMembrSize,
                                    repMembrIndx,
                                    secondNonMissMembrSize[r],
                                    tempNonMissMembrIndx,
                                    RF_observation[treeID][pseudoResponse[r]],
                                    &mean[r],
                                    &variance[r]);
          nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
          secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
        }  
        free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
      }  
      if (nonMissImpuritySummary) {
        for (j = 1; j <= repMembrSize; j++) {
          localSplitIndicator[j] = NEITHER;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          pseudoResponseClassSize[r] = 0;
          parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
          sumLeft[r] = sumRght[r] = sumRghtSave[r] = 0.0;
        }
        for (r = 1; r <= RF_randomResponseCount; r++) {
          if (impurity[r]) {
            if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
              pseudoResponseClassSize[r] = RF_xFactorSize[RF_xFactorMap[pseudoResponse[r]]];
              parentClassProp[r] = uivector(1, pseudoResponseClassSize[r]);
              leftClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              rghtClassProp[r]   = uivector(1, pseudoResponseClassSize[r]);
              for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                parentClassProp[r][p] = 0;
              }
              for (j = 1; j <= nonMissMembrSize; j++) {
                if (secondNonMissMembrFlag[r][j] == TRUE) {
                  parentClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ]] ++;
                }
              }
            }
            else {
              sumRghtSave[r] = 0.0;
              for (j = 1; j <= nonMissMembrSize; j++) {
                if (secondNonMissMembrFlag[r][j] == TRUE) {
                  sumRghtSave[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] - mean[r];
                }
              }
            }
          }  
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
          for (r = 1; r <= RF_randomResponseCount; r++) {
            if (impurity[r]) {
              if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
                for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                  rghtClassProp[r][p] = parentClassProp[r][p];
                  leftClassProp[r][p] = 0;
                }
              }
              else {
                sumRght[r] = sumRghtSave[r];
                sumLeft[r] = 0.0;
              }
              secondNonMissMembrLeftSize[r] = 0;
              secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
            }
          }
          leftSizeIter = 0;
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
            for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrLeftSize[r] = 0;
              secondNonMissMembrRghtSize[r] = 0;
            }
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
            delta = 0.0;
            deltaNorm = 0;
            for (r = 1; r <= RF_randomResponseCount; r++) {
              if (impurity[r]) {
                if (factorFlag == TRUE) {
                  if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
                    for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                      leftClassProp[r][p] = 0;
                    }
                    for (k = 1; k <= nonMissMembrSize; k++) {
                      if (secondNonMissMembrFlag[r][k] == TRUE) {
                        if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                          leftClassProp[r][ (uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]  ++;
                          secondNonMissMembrLeftSize[r] ++;
                        }
                        else {
                        }
                      }
                    }
                    for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                      rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                    }
                  }
                  else {
                    sumLeft[r] = sumRght[r] = 0.0;
                    for (k = 1; k <= nonMissMembrSize; k++) {
                      if (secondNonMissMembrFlag[r][k] == TRUE) {
                        if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                          sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          secondNonMissMembrLeftSize[r] ++;
                        }
                        else {
                          sumRght[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        }
                      }
                    }
                  }
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                }
                else {
                  for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
                        leftClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] ++;
                        rghtClassProp[r][(uint) RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] --;
                      }
                      else {
                        sumLeft[r] += RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        sumRght[r] -= RF_observation[treeID][pseudoResponse[r]][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                      }
                      secondNonMissMembrLeftSize[r] ++;
                      secondNonMissMembrRghtSize[r] --;
                    }
                  }
                }  
                if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                  deltaNorm ++;
                  if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
                    sumLeft[1] = sumRght[1] = 0;
                    for (p = 1; p <= pseudoResponseClassSize[r]; p++) {
                      sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                      sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                    }
                    sumLeftSqr = sumLeft[1] / secondNonMissMembrLeftSize[r];
                    sumRghtSqr  = sumRght[1] / secondNonMissMembrRghtSize[r];
                  }
                  else {
                    sumLeftSqr = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                    sumRghtSqr = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                  }
                  delta += sumLeftSqr + sumRghtSqr;
                }
              }  
            }  
            if (deltaNorm > 0) {
              delta = delta / (double) deltaNorm;
            }
            else {
              delta = NA_REAL;
            }
            if (factorFlag == FALSE) {
              leftSizeIter = currentMembrIter - 1;
            }
            updateMaximumSplit(treeID,
                               delta,
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
                           permissibleSplitSize,
                           splitLength,
                           factorFlag,
                           deterministicSplitFlag,
                           mwcpSizeAbsolute,
                           permissibleSplitPtr);
        for (r = 1; r <= RF_randomResponseCount; r++) {
          if (impurity[r]) {
            if (strcmp(RF_xType[pseudoResponse[r]], "C") == 0) {
              free_uivector (parentClassProp[r], 1, pseudoResponseClassSize[r]);
              free_uivector (leftClassProp[r],   1, pseudoResponseClassSize[r]);
              free_uivector (rghtClassProp[r],   1, pseudoResponseClassSize[r]);
            }
            else {
            }
          }
        }
      }  
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
      }
      else {
        for (r = 1; r <= RF_randomResponseCount; r++)  {
          free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
        }
      }
    }  
    free_new_vvector(parentClassProp, 1, RF_randomResponseCount, NRUTIL_UPTR);
    free_new_vvector(leftClassProp,   1, RF_randomResponseCount, NRUTIL_UPTR);
    free_new_vvector(rghtClassProp,   1, RF_randomResponseCount, NRUTIL_UPTR);
    free_dvector(sumLeft,     1, RF_randomResponseCount);
    free_dvector(sumRght,     1, RF_randomResponseCount);
    free_dvector(sumRghtSave, 1, RF_randomResponseCount);
    free_uivector(pseudoResponseClassSize, 1, RF_randomResponseCount);
    free_uivector(pseudoResponse, 1, RF_randomResponseCount);
    free_new_vvector(secondNonMissMembrFlag,  1, RF_randomResponseCount, NRUTIL_CPTR);
    free_uivector(secondNonMissMembrSize,     1, RF_randomResponseCount);
    free_uivector(secondNonMissMembrLeftSize, 1, RF_randomResponseCount);
    free_uivector(secondNonMissMembrRghtSize, 1, RF_randomResponseCount);
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
    free_cvector(impurity,   1, RF_randomResponseCount);
    free_dvector(mean,     1, RF_randomResponseCount);
    free_dvector(variance, 1, RF_randomResponseCount);
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
                        char    multImpFlag) {
  char   *randomCovariateFlag;
  uint   *randomCovariateIndex;
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
  uint leftSizeIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  double delta, deltaMax;
  uint   deltaNorm;
  uint j, k, p, r;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             NULL,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if (result) {
    char   *impurity   = cvector(1, RF_rSize);
    double *mean       = dvector(1, RF_rSize);
    double *variance   = dvector(1, RF_rSize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_rSize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][r],
                                  &mean[r],
                                  &variance[r]);
        impuritySummary = impuritySummary | impurity[r];
      }
    }
    else {
      impuritySummary = TRUE;
    }
    if (impuritySummary) {
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
      uint **parentClassProp = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      uint **leftClassProp   = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      uint **rghtClassProp   = (uint **) new_vvector(1, RF_rSize, NRUTIL_UPTR);
      double *sumLeft      = dvector(1, RF_rSize);
      double *sumRght      = dvector(1, RF_rSize);
      double *sumRghtSave  = dvector(1, RF_rSize);
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_rSize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_rSize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_rSize);
      for (r = 1; r <= RF_rSize; r++) {
        parentClassProp[r] = leftClassProp[r] = rghtClassProp[r] = NULL;
        if (strcmp(RF_rType[r], "C") == 0) {
          parentClassProp[r] = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          leftClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
          rghtClassProp[r]   = uivector(1, RF_classLevelSize[RF_rFactorMap[r]]);
        }
        else {
        }
      }  
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      double sumLeftSqr, sumRghtSqr;
      uint actualCovariateCount = 0;
      uint candidateCovariateCount = 0;
      while (selectRandomCovariates(treeID,
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
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempNonMissMembrFlag[k] = TRUE;
          }
          for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
          }
          nonMissImpuritySummary = TRUE;
        }
        else {
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_rSize; r++)  {
            secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
            j = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              mResponseFlag = FALSE;
              if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                if (RF_mpSign[r][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                  mResponseFlag = TRUE;
                }
              }
              if (!mResponseFlag) {
                j ++;
                tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
                secondNonMissMembrFlag[r][k] = TRUE;
              }
              else {
                secondNonMissMembrFlag[r][k] = FALSE;
              }
            }  
            secondNonMissMembrSize[r] = j;
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_response[treeID][r],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }  
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }  
        if (nonMissImpuritySummary) {
          for (j = 1; j <= repMembrSize; j++) {
            localSplitIndicator[j] = NEITHER;
          }
          for (r = 1; r <= RF_rSize; r++) {
            if (impurity[r]) {
              if (strcmp(RF_rType[r], "C") == 0) {
                for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                  parentClassProp[r][p] = 0;
                }
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    parentClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ]]] ++;
                  }
                }
              }
              else {
                sumRghtSave[r] = 0.0;
                for (j = 1; j <= nonMissMembrSize; j++) {
                  if (secondNonMissMembrFlag[r][j] == TRUE) {
                    sumRghtSave[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] - mean[r];
                  }
                }
              }
            }  
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
            for (r = 1; r <= RF_rSize; r++) {
              if (impurity[r]) {
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
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r];
              }
            }
            leftSizeIter = 0;
          }
          for (j = 1; j < splitLength; j++) {
            if (factorFlag == TRUE) {
              priorMembrIter = 0;
              leftSize = 0;
              for (r = 1; r <= RF_rSize; r++) {
                secondNonMissMembrLeftSize[r] = 0;
                secondNonMissMembrRghtSize[r] = 0;
              }
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
              delta     = 0.0;
              deltaNorm = 0;
              for (r = 1; r <= RF_rSize; r++) {
                if (impurity[r]) {
                  if (factorFlag == TRUE) {
                    if (strcmp(RF_rType[r], "C") == 0) {
                      for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                        leftClassProp[r][p] = 0;
                      }
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        if (secondNonMissMembrFlag[r][k] == TRUE) {
                          if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                            leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                            secondNonMissMembrLeftSize[r] ++;
                          }
                          else {
                          }
                        }
                      }
                      for (p=1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                        rghtClassProp[r][p] = parentClassProp[r][p] - leftClassProp[r][p];
                      }
                    }
                    else {
                      sumLeft[r] = sumRght[r] = 0.0;
                      for (k = 1; k <= nonMissMembrSize; k++) {
                        if (secondNonMissMembrFlag[r][k] == TRUE) {
                          if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                            sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                            secondNonMissMembrLeftSize[r] ++;
                          }
                          else {
                            sumRght[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          }
                        }
                      }
                    }
                    secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                  }
                  else {
                    for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                      if (secondNonMissMembrFlag[r][k] == TRUE) {
                        if (strcmp(RF_rType[r], "C") == 0) {
                          leftClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] ++;
                          rghtClassProp[r][RF_classLevelIndex[RF_rFactorMap[r]][(uint) RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]] --;
                        }
                        else {
                          sumLeft[r] += RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                          sumRght[r] -= RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - mean[r];
                        }
                        secondNonMissMembrLeftSize[r] ++;
                        secondNonMissMembrRghtSize[r] --;
                      }
                    }
                  }  
                  if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                    deltaNorm ++;
                    if (strcmp(RF_rType[r], "C") == 0) {
                      sumLeft[1] = sumRght[1] = 0;
                      for (p = 1; p <= RF_classLevelSize[RF_rFactorMap[r]]; p++) {
                        sumLeft[1] += (double) upower(leftClassProp[r][p], 2);
                        sumRght[1] += (double) upower(rghtClassProp[r][p], 2);
                      }
                      sumLeftSqr = sumLeft[1] / secondNonMissMembrLeftSize[r];
                      sumRghtSqr = sumRght[1] / secondNonMissMembrRghtSize[r];
                    }
                    else {
                      sumLeftSqr = pow (sumLeft[r], 2.0) / (secondNonMissMembrLeftSize[r] * variance[r]);
                      sumRghtSqr = pow (sumRght[r], 2.0) / (secondNonMissMembrRghtSize[r] * variance[r]);
                    }
                    delta += sumLeftSqr + sumRghtSqr;
                  }
                }  
              }  
              if (deltaNorm > 0) {
                delta = delta / (double) deltaNorm;
              }
              else {
                delta = NA_REAL;
              }
              if (factorFlag == FALSE) {
                leftSizeIter = currentMembrIter - 1;
              }
              updateMaximumSplit(treeID,
                                 delta,
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
                             permissibleSplitSize,
                             splitLength,
                             factorFlag,
                             deterministicSplitFlag,
                             mwcpSizeAbsolute,
                             permissibleSplitPtr);
        }  
        unselectRandomCovariates(treeID,
                                 parent,
                                 repMembrSize,
                                 indxx,
                                 nonMissMembrSizeStatic,
                                 nonMissMembrIndx,
                                 multImpFlag);
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
        }
        else {
          for (r = 1; r <= RF_rSize; r++)  {
            free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
          }
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
      free_new_vvector(parentClassProp, 1, RF_rSize, NRUTIL_UPTR);
      free_new_vvector(leftClassProp,   1, RF_rSize, NRUTIL_UPTR);
      free_new_vvector(rghtClassProp,   1, RF_rSize, NRUTIL_UPTR);
      free_dvector(sumLeft,     1, RF_rSize);
      free_dvector(sumRght,     1, RF_rSize);
      free_dvector(sumRghtSave, 1, RF_rSize);
      free_new_vvector(secondNonMissMembrFlag,  1, RF_rSize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, RF_rSize);
      free_uivector(secondNonMissMembrLeftSize, 1, RF_rSize);
      free_uivector(secondNonMissMembrRghtSize, 1, RF_rSize);
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
    else {
      free_uivector(nonMissMembrIndxStatic, 1, repMembrSize);
      free_dvector(permissibleSplit, 1, repMembrSize);
    }
    free_cvector(impurity,   1, RF_rSize);
    free_dvector(mean,       1, RF_rSize);
    free_dvector(variance,   1, RF_rSize);
  }  
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char customSplitMultivariate (uint    treeID,
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
  char   *randomCovariateFlag;
  uint   *randomCovariateIndex;
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
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  double delta, deltaMax;
  uint   deltaNorm;
  uint j, k, m, r;
  mwcpSizeAbsolute       = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             NULL,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if (result) {
    char   *impurity   = cvector(1, RF_rSize);
    double *mean       = dvector(1, RF_rSize);
    double *variance   = dvector(1, RF_rSize);
    char impuritySummary;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      impuritySummary = FALSE;
      for (r = 1; r <= RF_rSize; r++)  {
        impurity[r] = getVariance(repMembrSize,
                                  repMembrIndx,
                                  0,
                                  NULL,
                                  RF_response[treeID][r],
                                  &mean[r],
                                  &variance[r]);
        impuritySummary = impuritySummary | impurity[r];
      }
    }
    else {
      impuritySummary = TRUE;
    }
    if (impuritySummary) {
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
      char **secondNonMissMembrFlag = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
      uint  *secondNonMissMembrSize =           uivector(1, RF_rSize);
      uint  *secondNonMissMembrLeftSize =       uivector(1, RF_rSize);
      uint  *secondNonMissMembrRghtSize =       uivector(1, RF_rSize);
      char  *tempNonMissMembrFlag = 0;
      uint  *tempNonMissMembrIndx;
      char   mResponseFlag;
      char   nonMissImpuritySummary;
      uint actualCovariateCount = 0;
      uint candidateCovariateCount = 0;
      while (selectRandomCovariates(treeID,
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
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          tempNonMissMembrFlag = cvector(1, nonMissMembrSize);
          for (k = 1; k <= nonMissMembrSize; k++) {
            tempNonMissMembrFlag[k] = TRUE;
          }
          for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrFlag[r] = tempNonMissMembrFlag;
              secondNonMissMembrSize[r] = nonMissMembrSize;
          }
          nonMissImpuritySummary = TRUE;
        }
        else {
          tempNonMissMembrIndx = uivector(1, nonMissMembrSize);
          nonMissImpuritySummary = FALSE;
          for (r = 1; r <= RF_rSize; r++)  {
            secondNonMissMembrFlag[r] = cvector(1, nonMissMembrSize);
            j = 0;
            for (k = 1; k <= nonMissMembrSize; k++) {
              mResponseFlag = FALSE;
              if (RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                if (RF_mpSign[r][RF_mRecordMap[ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]] == 1) {
                  mResponseFlag = TRUE;
                }
              }
              if (!mResponseFlag) {
                j ++;
                tempNonMissMembrIndx[j] = nonMissMembrIndx[indxx[k]];
                secondNonMissMembrFlag[r][k] = TRUE;
              }
              else {
                secondNonMissMembrFlag[r][k] = FALSE;
              }
            }  
            secondNonMissMembrSize[r] = j;
            impurity[r] = getVariance(repMembrSize,
                                      repMembrIndx,
                                      secondNonMissMembrSize[r],
                                      tempNonMissMembrIndx,
                                      RF_response[treeID][r],
                                      &mean[r],
                                      &variance[r]);
            nonMissImpuritySummary = nonMissImpuritySummary | impurity[r];
            secondNonMissMembrLeftSize[r] = secondNonMissMembrRghtSize[r] = 0;
          }  
          free_uivector(tempNonMissMembrIndx, 1, nonMissMembrSize);
        }  
        if (nonMissImpuritySummary) {
          for (j = 1; j <= repMembrSize; j++) {
            localSplitIndicator[j] = NEITHER;
          }
          splitLength = stackAndConstructSplitVector(treeID,
                                                     repMembrSize,
                                                     covariate,
                                                     permissibleSplit,
                                                     permissibleSplitSize,
                                                     & factorFlag,
                                                     & deterministicSplitFlag,
                                                     & mwcpSizeAbsolute,
                                                     & permissibleSplitPtr);
          double *userResponse = dvector(1, nonMissMembrSize);
          char   *userSplitIndicator = cvector(1, nonMissMembrSize);
          for (j = 1; j < splitLength; j++) {
            priorMembrIter = 0;
            leftSize = 0;
            for (r = 1; r <= RF_rSize; r++) {
              secondNonMissMembrLeftSize[r] = 0;
              secondNonMissMembrRghtSize[r] = 0;
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
            if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
              delta     = 0.0;
              deltaNorm = 0;
              for (r = 1; r <= RF_rSize; r++) {
                if (impurity[r]) {
                  m = 0;
                  for (k = 1; k <= nonMissMembrSize; k++) {
                    if (secondNonMissMembrFlag[r][k] == TRUE) {
                      if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                        secondNonMissMembrLeftSize[r] ++;
                      }
                      else {
                      }
                      userResponse[++m] = RF_response[treeID][r][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                      userSplitIndicator[m] = localSplitIndicator[ nonMissMembrIndx[indxx[k]] ];
                    }
                  }
                  secondNonMissMembrRghtSize[r] = secondNonMissMembrSize[r] - secondNonMissMembrLeftSize[r];
                  if ((secondNonMissMembrLeftSize[r] > 0) && (secondNonMissMembrRghtSize[r] > 0)) {
                    deltaNorm ++;
                    delta += getCustomSplitStatistic(m, userSplitIndicator, NULL, NULL, userResponse, mean[r], variance[r]);
                  }
                }  
              }  
              if (deltaNorm > 0) {
                delta = delta / (double) deltaNorm;
              }
              else {
                delta = NA_REAL;
              }
              updateMaximumSplit(treeID,
                                 delta,
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
            }  
          }  
          free_dvector (userResponse, 1, nonMissMembrSize);
          free_cvector (userSplitIndicator, 1, nonMissMembrSize);
          unstackSplitVector(treeID,
                             permissibleSplitSize,
                             splitLength,
                             factorFlag,
                             deterministicSplitFlag,
                             mwcpSizeAbsolute,
                             permissibleSplitPtr);
        }  
        unselectRandomCovariates(treeID,
                                 parent,
                                 repMembrSize,
                                 indxx,
                                 nonMissMembrSizeStatic,
                                 nonMissMembrIndx,
                                 multImpFlag);
        if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
          free_cvector(tempNonMissMembrFlag, 1, nonMissMembrSize);
        }
        else {
          for (r = 1; r <= RF_rSize; r++)  {
            free_cvector(secondNonMissMembrFlag[r], 1, nonMissMembrSize);
          }
        }
      }  
      free_new_vvector(secondNonMissMembrFlag,  1, RF_rSize, NRUTIL_CPTR);
      free_uivector(secondNonMissMembrSize,     1, RF_rSize);
      free_uivector(secondNonMissMembrLeftSize, 1, RF_rSize);
      free_uivector(secondNonMissMembrRghtSize, 1, RF_rSize);
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
    else {
      free_uivector(nonMissMembrIndxStatic, 1, repMembrSize);
      free_dvector(permissibleSplit, 1, repMembrSize);
    }
    free_cvector(impurity,   1, RF_rSize);
    free_dvector(mean,       1, RF_rSize);
    free_dvector(variance,   1, RF_rSize);
  }  
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
