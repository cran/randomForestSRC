////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.1
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
#include    "regression.h"
#include     "splitUtil.h"
void stackSplitIndicator(uint   nodeSize,
                         char **localSplitIndicator) {
  if (nodeSize > 0) {
    *localSplitIndicator = cvector(1, nodeSize);
  }
}
void unstackSplitIndicator(uint  nodeSize,
                           char *localSplitIndicator) {
  if (nodeSize > 0) {
    if (nodeSize > 0) {
      free_cvector(localSplitIndicator, 1, nodeSize);
    }
  }
}
void stackSplitEventTime(uint **localEventTimeCount,
                         uint **localEventTimeIndex) {
  *localEventTimeCount = uivector(1, RF_masterTimeSize);
  *localEventTimeIndex = uivector(1, RF_masterTimeSize);
}
void unstackSplitEventTime(uint *localEventTimeCount,
                           uint *localEventTimeIndex) {
  free_uivector(localEventTimeCount, 1, RF_masterTimeSize);
  free_uivector(localEventTimeIndex, 1, RF_masterTimeSize);
}
uint getSplitEventTime(uint   treeID,
                       uint   *repMembrIndx,
                       uint    repMembrSize,
                       uint   *nonMissMembrIndx,
                       uint    nonMissMembrSize,
                       uint   *localEventTimeCount,
                       uint   *localEventTimeIndex) {
  uint parentEventCount;
  uint i;
  uint eventTimeSize;
  parentEventCount = 0;
  eventTimeSize = 0;
  for (i=1; i <= RF_masterTimeSize; i++) {
    localEventTimeCount[i] = 0;
  }
  for (i = 1; i <= nonMissMembrSize; i++) {
    if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[i]] ] > 0) {
      localEventTimeCount[RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[i]] ]] ++;
      parentEventCount ++;
    }
  }
  for (i=1; i <= RF_masterTimeSize; i++) {
    if (localEventTimeCount[i] > 0) {
      localEventTimeIndex[++eventTimeSize] = i;
    }
  }
  return (eventTimeSize);
}
void stackSplitEventAndRisk(uint   eventTimeSize,
                            uint **nodeParentEvent,
                            uint **nodeParentAtRisk,
                            uint **nodeLeftEvent,
                            uint **nodeLeftAtRisk,
                            uint **nodeRightEvent,
                            uint **nodeRightAtRisk) {
  if (eventTimeSize > 0) {
    *nodeParentEvent  = uivector(1, eventTimeSize);
    *nodeParentAtRisk = uivector(1, eventTimeSize);
    *nodeLeftEvent  = uivector(1, eventTimeSize);
    *nodeLeftAtRisk = uivector(1, eventTimeSize);
    *nodeRightEvent  = uivector(1, eventTimeSize);
    *nodeRightAtRisk = uivector(1, eventTimeSize);
  }
  else {
    *nodeParentEvent = *nodeParentAtRisk = *nodeLeftEvent  = *nodeLeftAtRisk = *nodeRightEvent  = *nodeRightAtRisk = NULL;
  }
}
void unstackSplitEventAndRisk(uint  eventTimeSize,
                              uint *nodeParentEvent,
                              uint *nodeParentAtRisk,
                              uint *nodeLeftEvent,
                              uint *nodeLeftAtRisk,
                              uint *nodeRightEvent,
                              uint *nodeRightAtRisk) {
  if (eventTimeSize > 0) {
    free_uivector(nodeParentEvent, 1, eventTimeSize);
    free_uivector(nodeParentAtRisk, 1, eventTimeSize);
    free_uivector(nodeLeftEvent, 1, eventTimeSize);
    free_uivector(nodeLeftAtRisk, 1, eventTimeSize);
    free_uivector(nodeRightEvent, 1, eventTimeSize);
    free_uivector(nodeRightAtRisk, 1, eventTimeSize);
  }
}
void getSplitEventAndRisk(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *localEventTimeCount,
                          uint   *localEventTimeIndex,
                          uint    localEventTimeSize,
                          uint   *nodeParentEvent,
                          uint   *nodeParentAtRisk) {
  uint i, j;
  for (i=1; i <= localEventTimeSize; i++) {
    nodeParentAtRisk[i] = 0;
    nodeParentEvent[i] = localEventTimeCount[localEventTimeIndex[i]];
    for (j = 1; j <= nonMissMembrSize; j++) {
      if (localEventTimeIndex[i] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]) {
        nodeParentAtRisk[i] ++;
      }
    }
  }
}
void stackAndGetSplitSurv(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint  **localEventTimeCount,
                          uint  **localEventTimeIndex,
                          uint   *localEventTimeSize,
                          uint  **nodeParentEvent,
                          uint  **nodeParentAtRisk,
                          uint  **nodeLeftEvent,
                          uint  **nodeLeftAtRisk,
                          uint  **nodeRightEvent,
                          uint  **nodeRightAtRisk) {
  stackSplitEventTime(localEventTimeCount, localEventTimeIndex);
  *localEventTimeSize = getSplitEventTime( treeID,
                                           repMembrIndx,
                                           repMembrSize,
                                           nonMissMembrIndx,
                                           nonMissMembrSize,
                                          *localEventTimeCount,
                                          *localEventTimeIndex);
  stackSplitEventAndRisk(*localEventTimeSize,
                          nodeParentEvent,
                          nodeParentAtRisk,
                          nodeLeftEvent,
                          nodeLeftAtRisk,
                          nodeRightEvent,
                          nodeRightAtRisk);
  getSplitEventAndRisk( treeID,
                        repMembrIndx,
                        repMembrSize,
                        nonMissMembrIndx,
                        nonMissMembrSize,
                       *localEventTimeCount,
                       *localEventTimeIndex,
                       *localEventTimeSize,
                       *nodeParentEvent,
                       *nodeParentAtRisk);
}
void unstackSplitSurv(uint *localEventTimeCount,
                      uint *localEventTimeIndex,
                      uint  eventTimeSize,
                      uint *nodeParentEvent,
                      uint *nodeParentAtRisk,
                      uint *nodeLeftEvent,
                      uint *nodeLeftAtRisk,
                      uint *nodeRightEvent,
                      uint *nodeRightAtRisk) {
  unstackSplitEventTime(localEventTimeCount,
                        localEventTimeIndex);
  unstackSplitEventAndRisk(eventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
}
uint stackAndConstructSplitVector (uint     treeID,
                                   uint     repMembrSize,
                                   uint     randomCovariateIndex,
                                   double  *permissibleSplit,
                                   uint     permissibleSplitSize,
                                   char    *factorFlag,
                                   char    *deterministicSplitFlag,
                                   uint    *mwcpSizeAbsolute,
                                   void   **permissibleSplitPtr) {
  uint j, j2, k2;
  uint factorSizeAbsolute;
  uint offset;
  uint splitLength;
  uint relativePair;
  splitLength = 0;  
  (*permissibleSplitPtr) = NULL;  
  if (strcmp(RF_xType[randomCovariateIndex], "C") == 0) {
    *factorFlag = TRUE;
    if(RF_factorList[treeID][permissibleSplitSize] == NULL) {
      RF_factorList[treeID][permissibleSplitSize] = makeFactor(permissibleSplitSize, FALSE);
    }
    factorSizeAbsolute = RF_xFactorSize[RF_xFactorMap[randomCovariateIndex]];
    *mwcpSizeAbsolute = RF_factorList[treeID][factorSizeAbsolute] -> mwcpSize;
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomCount == 0) {
        *deterministicSplitFlag = TRUE;
        if ((RF_factorList[treeID][permissibleSplitSize] -> r) > MAX_EXACT_LEVEL) {
          *deterministicSplitFlag = FALSE;
        }
        else {
          if ( *((uint *) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) >= repMembrSize ) {
            *deterministicSplitFlag = FALSE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = repMembrSize + 1;
        }
        else {
          splitLength = *((uint*) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) + 1;
        }
      }
      else {
        *deterministicSplitFlag = FALSE;
        if ((RF_factorList[treeID][permissibleSplitSize] -> r) <= MAX_EXACT_LEVEL) {
          if (*((uint*) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) <= ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize)) {
            splitLength = *((uint*) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) + 1;
            *deterministicSplitFlag = TRUE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = 1 + ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize);
        }
      }  
    }  
    (*permissibleSplitPtr) = uivector(1, splitLength * (*mwcpSizeAbsolute));
    for (offset = 1; offset <= *mwcpSizeAbsolute; offset++) {
      ((uint*) (*permissibleSplitPtr) + ((splitLength - 1) * (*mwcpSizeAbsolute)))[offset] = 0;
    }
    if (*deterministicSplitFlag) {
      bookFactor(RF_factorList[treeID][permissibleSplitSize]);
      j2 = 0;
      for (j = 1; j <= RF_factorList[treeID][permissibleSplitSize] -> cardinalGroupCount; j++) {
        for (k2 = 1; k2 <= ((uint*) RF_factorList[treeID][permissibleSplitSize] -> cardinalGroupSize)[j]; k2++) {
          ++j2;
          relativePair = (RF_factorList[treeID][permissibleSplitSize] -> cardinalGroupBinary)[j][k2];
          convertRelToAbsBinaryPair(treeID,
                                    permissibleSplitSize,
                                    factorSizeAbsolute,
                                    relativePair,
                                    permissibleSplit,
                                    (uint*) (*permissibleSplitPtr) + ((j2 - 1) * (*mwcpSizeAbsolute)));
        }
      }
    }  
    else {
      for (j = 1; j < splitLength; j++) {
        getRandomPair(treeID, permissibleSplitSize, factorSizeAbsolute, permissibleSplit, (uint*) (*permissibleSplitPtr) + ((j - 1) * (*mwcpSizeAbsolute)));
      }
    }
  }  
  else {
    *factorFlag = FALSE;
    if (RF_splitRule == RAND_SPLIT) {
      splitLength = 1 + ((RF_splitRandomCount <= repMembrSize) ? RF_splitRandomCount : repMembrSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomCount == 0) {
        splitLength = permissibleSplitSize;
        (*permissibleSplitPtr) = permissibleSplit;
        *deterministicSplitFlag = TRUE;
      }
      else {
        if (permissibleSplitSize <= RF_splitRandomCount) {
          splitLength = permissibleSplitSize;
          (*permissibleSplitPtr) = permissibleSplit;
          *deterministicSplitFlag = TRUE;
        }
        else {
          splitLength = RF_splitRandomCount + 1;
          *deterministicSplitFlag = FALSE;
        }
      }  
    }  
    if (*deterministicSplitFlag == FALSE) {
      (*permissibleSplitPtr) = dvector(1, splitLength);
      ((double*) (*permissibleSplitPtr))[splitLength] = 0;
      for (j = 1; j < splitLength; j++) {
        ((double*) (*permissibleSplitPtr))[j]  = permissibleSplit[(uint) ceil(ran1B(treeID) * ((permissibleSplitSize - 1) * 1.0))];
      }
      hpsort(((double*) (*permissibleSplitPtr)), splitLength-1);
    }
  }  
  return splitLength;
}
void unstackSplitVector(uint   treeID,
                        uint   permissibleSplitSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        uint   mwcpSizeAbsolute,
                        void  *permissibleSplitPtr) {
  if (factorFlag == TRUE) {
    free_uivector(permissibleSplitPtr, 1, splitLength * mwcpSizeAbsolute);
    if (deterministicSplitFlag == FALSE) {
      if (permissibleSplitSize > SAFE_FACTOR_SIZE) {
        unbookFactor(RF_factorList[treeID][permissibleSplitSize]);
      }
    }
  }
  else {
    if (deterministicSplitFlag == FALSE) {
      free_dvector(permissibleSplitPtr, 1, splitLength);
    }
  }
}
void stackRandomCovariates(uint      treeID,
                           Node     *parent,
                           uint      repMembrSize,
                           char      multImpFlag,
                           char    **covariateFlag,
                           uint    **covariateIndex,
                           uint     *uniformCovariateSize,
                           uint     *uniformCovariateIndex,
                           double  **cdf,
                           uint     *cdfSize,
                           uint    **cdfSort,
                           uint    **density,
                           uint     *densitySize,
                           uint   ***densitySwap) {
  *covariateFlag = cvector(1, RF_xSize);
  *covariateIndex = uivector(1, RF_xSize);
  initializeCDF(treeID,
                parent,
                *covariateFlag,
                *covariateIndex,
                uniformCovariateSize,
                uniformCovariateIndex,
                cdf,
                cdfSize,
                cdfSort,
                density,
                densitySize,
                densitySwap);
}
void unstackRandomCovariates(uint     treeID,
                             char    *covariateFlag,
                             uint    *covariateIndex,
                             uint     uniformSize,
                             double  *cdf,
                             uint     cdfSize,
                             uint    *cdfSort,
                             uint    *density,
                             uint     densitySize,
                             uint   **densitySwap,
                             uint     repMembrSize,
                             uint    *nonMissMembrIndxStatic,
                             double  *permissibleSplit) {
  uint k;
  free_cvector(covariateFlag, 1, RF_xSize);
  free_uivector(covariateIndex, 1, RF_xSize);
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    free_uivector(density, 1, RF_xWeightDensitySize);
    for (k = 1; k <= RF_xSize; k++) {
      if (densitySwap[k] != NULL) {
        free_uivector(densitySwap[k], 1, (uint) RF_xWeight[k]);
        densitySwap[k] = NULL;
      }
    }
    free_vvector(densitySwap, 1, RF_xSize);
    break;
  case RF_WGHT_GENERIC:
    free_dvector(cdf, 1, RF_xSize);
    free_uivector(cdfSort, 1, RF_xSize);
    break;
  }
  if (TRUE) {
    free_uivector(nonMissMembrIndxStatic, 1, repMembrSize);
    free_dvector(permissibleSplit, 1, repMembrSize);
  }
}
void initializeCDF(uint     treeID,
                   Node    *parent,
                   char    *covariateFlag,
                   uint    *covariateIndex,
                   uint    *uniformCovariateSize,
                   uint     *uniformCovariateIndex,
                   double **cdf,
                   uint    *cdfSize,
                   uint   **cdfSort,
                   uint   **density,
                   uint    *densitySize,
                   uint  ***densitySwap) {
  uint i, j, k, kk;
  *uniformCovariateSize = 0;
  *uniformCovariateIndex = 0;
  *cdf = NULL;
  *cdfSize = 0;
  *cdfSort = NULL;
  *density = NULL;
  *densitySize = 0;
  *densitySwap = NULL;
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    *uniformCovariateSize = 0;
    for (k=1; k <= RF_xSize; k++) {
      covariateFlag[k] = (parent -> permissibleSplit)[k];
      if (covariateFlag[k]) {
        covariateIndex[++(*uniformCovariateSize)] = k;
      }
    }
    break;
  case RF_WGHT_INTEGER:
    (*density) = uivector(1, RF_xWeightDensitySize);
    (*densitySize) = 0;
    *densitySwap = (uint**) vvector(1, RF_xSize);
    for (k = RF_xSize; k >= 1; k--) {
      kk = RF_xWeightSorted[k];
      covariateFlag[kk] = (parent -> permissibleSplit)[kk];
      if (covariateFlag[kk]) {
        j = (uint) RF_xWeight[kk];
        if (j > 0) {
          (*densitySwap)[kk] = uivector(1, j);
          for (i = 1; i <= j; i++) {
            (*density)[++(*densitySize)] = kk;
            (*densitySwap)[kk][i] = (*densitySize);
          }
        }
        else {
          (*densitySwap)[kk] = NULL;
        }
      }
      else {
        (*densitySwap)[kk] = NULL;
      }
    }
    break;
  case RF_WGHT_GENERIC:
    i = 0;
    *cdf     = dvector(1, RF_xSize);
    *cdfSort = uivector(1, RF_xSize);
    *cdfSize = 0;
    for (k = 1; k <= RF_xSize; k++) {
      kk = RF_xWeightSorted[k];
      covariateFlag[kk] = (parent -> permissibleSplit)[kk];
      if (covariateFlag[kk]) {
        covariateIndex[kk] = ++ i;
        (*cdfSize) ++;
        (*cdfSort)[(*cdfSize)] = kk;
        (*cdf)[(*cdfSize)] = RF_xWeight[k];
      }
      else {
        covariateIndex[kk] = 0;
      }
    }
    for (k = 2; k <= (*cdfSize); k++) {
      (*cdf)[k] += (*cdf)[k-1];
    }
    break;
  }
}
void updateCDF(uint    treeID,
               char   *covariateFlag,
               uint   *covariateIndex,
               uint   *uniformCovariateSize,
               uint    uniformCovariateIndex,
               double *cdf,
               uint   *cdfSize,
               uint   *cdfSort,
               uint   *density,
               uint   *densitySize,
               uint  **densitySwap,
               uint    index) {
  double stepValue;
  uint sourcePt;
  uint stepIndex;
  uint currCov, nextCov;
  uint   i, j, k;
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    covariateFlag[covariateIndex[uniformCovariateIndex]] = FALSE;
    covariateIndex[uniformCovariateIndex] = covariateIndex[(*uniformCovariateSize)];
    (*uniformCovariateSize) --;
    break;
  case RF_WGHT_INTEGER:
    covariateFlag[index] = FALSE;
    currCov = nextCov = density[*densitySize];
    i = 0;
    j = (uint) RF_xWeight[currCov];
    k = (uint) RF_xWeight[index];
    while(i < k) {
      if (density[(*densitySize)] == index) {
        density[(*densitySize)] = 0;
        (*densitySize) --;
        densitySwap[index][k] = 0;
        k--;
        currCov = nextCov = density[(*densitySize)];
        j = (uint) RF_xWeight[currCov];
      }
      else {
        i++;
        sourcePt = densitySwap[index][i];
        density[sourcePt] = density[(*densitySize)];
        density[(*densitySize)] = 0;
        (*densitySize) --;
        densitySwap[currCov][j] = densitySwap[index][i];
        densitySwap[index][i] = 0;
        nextCov = density[(*densitySize)];
        if (nextCov == currCov) {
          j--;
        }
        else {
          hpsortui(densitySwap[currCov], (uint) RF_xWeight[currCov]);
          currCov = nextCov = density[(*densitySize)];
          j = (uint) RF_xWeight[currCov];
        }
      }
    }
    if (nextCov == currCov) {
      hpsortui(densitySwap[currCov], (uint) RF_xWeight[currCov]);
    }
    break;
  case RF_WGHT_GENERIC:
    covariateFlag[index] = FALSE;
    stepIndex = covariateIndex[index];
    stepValue = cdf[stepIndex];
    if (stepIndex > 1) {
      stepValue -= cdf[stepIndex-1];
    }
    for (k = stepIndex; k <= (*cdfSize); k++) {
      cdf[k] = cdf[k] - stepValue;
    }
    break;
  }
}
uint getSelectableElementNew (uint    treeID,
                              char   *covariateFlag,
                              uint   *covariateIndex,
                              uint    uniformCovariateSize,
                              uint   *uniformCovariateIndex,
                              double *cdf,
                              uint    cdfSize,
                              uint   *cdfSort,
                              uint   *density,
                              uint    densitySize) {
  double randomValue;
  uint low, mid, high, index, p;
  index = 0;  
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    if (uniformCovariateSize > 0) {
      (*uniformCovariateIndex) = (uint) ceil(ran1B(treeID) * (uniformCovariateSize * 1.0));
      index = covariateIndex[(*uniformCovariateIndex)];
    }
    else {
      index = 0;
    }
    break;
  case RF_WGHT_INTEGER:
    if (densitySize > 0) {
      p = (uint) ceil(ran1B(treeID) * (densitySize * 1.0));
      index = density[p];
    }
    else {
      index = 0;
    }
    break;
  case RF_WGHT_GENERIC:
    if (cdf[cdfSize] > 0) {
      randomValue = ran1B(treeID) * cdf[cdfSize];
      low  = 1;
      high = cdfSize;
      while (low < high) {
        mid  = (low + high) >> 1;
        if (randomValue > cdf[mid]) {
          if (low == mid) {
            low = high;
          }
          else {
            low = mid;
          }
        }
        else {
          if (low == mid) {
            low = high;
          }
          else {
            high = mid;
          }
        }
      }
      index = high;
    }
    else {
      index = 0;
    }
    break;
  }
  return index;
}
uint sampleWithoutReplacement (uint    treeID,
                               uint   *index,
                               uint    size,
                               uint   *sampleIndex) {
  uint absoluteIndex;
  if (size > 0) {
    (*sampleIndex) = (uint) ceil(ran1B(treeID) * (size * 1.0));
    absoluteIndex = index[*sampleIndex];
  }
  else {
    absoluteIndex = 0;
  }
  return absoluteIndex;
}
uint getSelectableElement (uint    treeID,
                           uint    length,
                           char   *permissible) {
  char   *localPermissible = NULL;  
  uint selectableCount;
  uint i, p, index;
  if (length > 0) {
    localPermissible = cvector(1, length);
  }
  selectableCount = 0;
  for (i=1; i <= length; i++) {
    if (permissible[i] == TRUE) {
      localPermissible[i] = TRUE;
      selectableCount ++;
    }
    else {
      localPermissible[i] = FALSE;
    }
  }
  if (selectableCount > 0) {
    p = (uint) ceil(ran1B(treeID) * (selectableCount * 1.0));
    index = 1;
    while (p > 0) {
      if (localPermissible[index] == TRUE) {
        p --;
      }
      index ++;
    }
    index --;
  }  
  else {
    index = 0;
  }
  if (length > 0) {
    free_cvector(localPermissible, 1, length);
  }
  return index;
}
char selectRandomCovariates(uint     treeID,
                            Node     *parent,
                            uint     *repMembrIndx,
                            uint      repMembrSize,
                            char     *covariateFlag,
                            uint     *covariateIndex,
                            uint     *uniformCovariateSize,
                            uint     *uniformCovariateIndex,
                            double   *cdf,
                            uint     *cdfSize,
                            uint     *cdfSort,
                            uint     *density,
                            uint     *densitySize,
                            uint    **densitySwap,
                            uint     *covariate,
                            uint     *actualCovariateCount,
                            uint     *candidateCovariateCount,
                            double   *permissibleSplit,
                            uint     *permissibleSplitSize,
                            uint    **indxx,
                            uint      nonMissMembrSizeStatic,
                            uint     *nonMissMembrIndxStatic,
                            uint     *nonMissMembrSize,
                            uint    **nonMissMembrIndx,
                            char      multImpFlag) {
  uint i;
  uint candidateCovariate;
  uint offset;
  uint indx;
  double *nonMissSplit;
  char mPredictorFlag;
  char splittable;
  if (nonMissMembrSizeStatic < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid nonMissMembrSizeStatic encountered in selectRandomCovariates():  %10d", nonMissMembrSizeStatic);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nonMissSplit = dvector(1, repMembrSize);
  (*covariate) = candidateCovariate = -1;
  splittable = FALSE;
  (*indxx) = uivector(1, repMembrSize);
  if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
    *nonMissMembrSize = nonMissMembrSizeStatic;
    *nonMissMembrIndx = nonMissMembrIndxStatic;
  }
  else {
    *nonMissMembrSize = 0;
    *nonMissMembrIndx = uivector(1, nonMissMembrSizeStatic);
  }
  while ( ((*candidateCovariateCount) < RF_randomCovariateCount) &&
          (candidateCovariate != 0) && (splittable == FALSE)) {
    candidateCovariate = getSelectableElementNew(treeID,
                                                 covariateFlag,
                                                 covariateIndex,
                                                 *uniformCovariateSize,
                                                 uniformCovariateIndex,
                                                 cdf,
                                                 *cdfSize,
                                                 cdfSort,
                                                 density,
                                                 *densitySize);
    if (candidateCovariate != 0) {
      updateCDF(treeID,
                covariateFlag,
                covariateIndex,
                uniformCovariateSize,
                *uniformCovariateIndex,
                cdf,
                cdfSize,
                cdfSort,
                density,
                densitySize,
                densitySwap,
                candidateCovariate);
      (*actualCovariateCount) ++;
      (*candidateCovariateCount) ++;
      splittable = TRUE;
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        for (i = 1; i <= repMembrSize; i++) {
          nonMissSplit[i] = RF_observation[treeID][candidateCovariate][repMembrIndx[i]];
        }
        *nonMissMembrSize = nonMissMembrSizeStatic;
        *nonMissMembrIndx = nonMissMembrIndxStatic;
      }
      else {
        offset = RF_rSize + candidateCovariate;
        (*nonMissMembrSize) = 0;
        for (i = 1; i <= nonMissMembrSizeStatic; i++) {
          mPredictorFlag = FALSE;
          if (RF_mRecordMap[repMembrIndx[nonMissMembrIndxStatic[i]]] > 0) {
            if (RF_mpSign[offset][RF_mRecordMap[repMembrIndx[nonMissMembrIndxStatic[i]]]] == 1) {
                mPredictorFlag = TRUE;
            }
          }
          if (!mPredictorFlag) {
            (*nonMissMembrSize) ++;
            (*nonMissMembrIndx)[*nonMissMembrSize] = nonMissMembrIndxStatic[i];
            nonMissSplit[*nonMissMembrSize] = RF_observation[treeID][candidateCovariate][repMembrIndx[(*nonMissMembrIndx)[*nonMissMembrSize]]];
          }
        }  
      }  
      if ((*nonMissMembrSize) == 0) {
        splittable = FALSE;
      }
      if (splittable) {
        indexx((*nonMissMembrSize),
               nonMissSplit,
               (*indxx));
        for (i = 1; i <= (*nonMissMembrSize); i++) {
          indx = (*indxx)[i];
          permissibleSplit[i] = nonMissSplit[indx];
        }
        (*permissibleSplitSize) = 1;
        for (i = 2; i <= (*nonMissMembrSize); i++) {
          if (permissibleSplit[i] > permissibleSplit[(*permissibleSplitSize)]) {
            (*permissibleSplitSize) ++;
            permissibleSplit[(*permissibleSplitSize)] = permissibleSplit[i];
          }
        }
        if((*permissibleSplitSize) >= 2) {
          (*covariate) = candidateCovariate;
        }
        else {
          splittable = FALSE;
        }
      }  
      if (!splittable) {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        (*actualCovariateCount) --;
      }
    }  
    else {
    }
  }  
  if (!splittable) {
    free_uivector(*indxx, 1, repMembrSize);
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      *nonMissMembrSize = 0;
      *nonMissMembrIndx = NULL;
    }
    else {
      *nonMissMembrSize = 0;
      free_uivector(*nonMissMembrIndx, 1, nonMissMembrSizeStatic);
    }
  }
  free_dvector(nonMissSplit, 1, repMembrSize);
  return splittable;
}
void unselectRandomCovariates(uint      treeID,
                              Node     *parent,
                              uint      repMembrSize,
                              uint     *indxx,
                              uint     nonMissMembrSizeStatic,
                              uint    *nonMissMembrIndx,
                              char      multImpFlag) {
  free_uivector((indxx), 1, repMembrSize);
  if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
    free_uivector(nonMissMembrIndx, 1, nonMissMembrSizeStatic);
  }
}
uint virtuallySplitNode(uint  treeID,
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
void getReweightedRandomPair (uint    treeID,
                              uint    relativeFactorSize,
                              uint    absoluteFactorSize,
                              double *absoluteLevel,
                              uint   *result) {
  uint randomGroupIndex;
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  randomGroupIndex = (uint) ceil(ran1B(treeID) * ((RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount) * 1.0));
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void getRandomPair (uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result) {
  uint randomGroupIndex;
  double randomValue;
  uint k;
  if(RF_factorList[treeID][relativeFactorSize] == NULL) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Factor not allocated for size:  %10d", relativeFactorSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  double *cdf = dvector(1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  if (relativeFactorSize <= MAX_EXACT_LEVEL) {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = (double) ((uint*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  else {
    for (k=1; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
      cdf[k] = ((double*) RF_factorList[treeID][relativeFactorSize] -> cardinalGroupSize)[k];
    }
  }
  for (k=2; k <= RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount; k++) {
    cdf[k] += cdf[k-1];
  }
  randomValue = ceil((ran1B(treeID) * cdf[RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount]));
  randomGroupIndex = 1;
  while (randomValue > cdf[randomGroupIndex]) {
    randomGroupIndex ++;
  }
  free_dvector(cdf, 1, RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount);
  createRandomBinaryPair(treeID, relativeFactorSize, absoluteFactorSize, randomGroupIndex, absoluteLevel, result);
}
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupIndex,
                            double *absoluteLevel,
                            uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint k, offset;
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  char *localPermissible = cvector(1, relativeFactorSize);
  uint *randomLevel = uivector(1, groupIndex);
  for (k = 1; k <= relativeFactorSize; k++) {
    localPermissible[k] = TRUE;
  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = getSelectableElement(treeID, relativeFactorSize, localPermissible);
    localPermissible[randomLevel[k]] = FALSE;
  }
  for (k = 1; k <= groupIndex; k++) {
    randomLevel[k] = (uint) absoluteLevel[randomLevel[k]];
  }
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= groupIndex; k++) {
    mwcpLevelIdentifier = (randomLevel[k] >> (3 + ulog2(SIZE_OF_INTEGER))) + ((randomLevel[k] & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
    pair[mwcpLevelIdentifier] += upower(2, randomLevel[k] - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
  }
  free_cvector(localPermissible, 1, relativeFactorSize);
  free_uivector(randomLevel, 1, groupIndex);
}
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair) {
  uint mwcpLevelIdentifier;
  uint mwcpSizeAbsolute;
  uint coercedAbsoluteLevel;
  uint k, offset;
  mwcpSizeAbsolute = RF_factorList[treeID][absoluteFactorSize] -> mwcpSize;
  for (offset = 1; offset <= mwcpSizeAbsolute; offset++) {
    pair[offset] = 0;
  }
  for (k = 1; k <= relativeFactorSize; k++) {
    if (relativePair & ((uint) 0x01)) {
      coercedAbsoluteLevel = (uint) absoluteLevel[k];
      mwcpLevelIdentifier = (coercedAbsoluteLevel >> (3 + ulog2(SIZE_OF_INTEGER))) + ((coercedAbsoluteLevel & (MAX_EXACT_LEVEL - 1)) ? 1 : 0);
      pair[mwcpLevelIdentifier] += upower(2, coercedAbsoluteLevel - ((mwcpLevelIdentifier - 1) * MAX_EXACT_LEVEL) - 1 );
    }
    relativePair = relativePair >> 1;
  }
}
char summarizeSplitResult(uint    splitParameterMax,
                          double  splitValueMaxCont,
                          uint    splitValueMaxFactSize,
                          uint   *splitValueMaxFactPtr,
                          double *splitStatistic,
                          double  deltaMax) {
  char result;
  if (splitParameterMax > 0) {
    *splitStatistic = deltaMax;
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
char getPreSplitResult (uint      treeID,
                        Node     *parent,
                        uint      repMembrSize,
                        uint     *repMembrIndx,
                        uint     *nonMissMembrSize,
                        uint    **nonMissMembrIndx,
                        double  **permissibleSplit,
                        char      multImpFlag) {
  uint i, r;
  char mPredictorFlag;
  char result;
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
    *nonMissMembrIndx = uivector(1, repMembrSize);
    *permissibleSplit = dvector(1, repMembrSize);
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)) || (repMembrIndx == NULL)) {
      (*nonMissMembrSize) = repMembrSize;
      for (i = 1; i <= repMembrSize; i++) {
        (*nonMissMembrIndx)[i] = i;
      }
    }
    else {
      (*nonMissMembrSize) = 0;
      for (i = 1; i <= repMembrSize; i++) {
        mPredictorFlag = FALSE;
        if (RF_mRecordMap[repMembrIndx[i]] > 0) {
          for (r = 1; r<= RF_rSize; r++) {
            if (RF_mpSign[r][RF_mRecordMap[repMembrIndx[i]]] == 1) {
              mPredictorFlag = TRUE;
            }
          }
        }
        if (!mPredictorFlag) {
          (*nonMissMembrSize) ++;
          (*nonMissMembrIndx)[(*nonMissMembrSize)] = i;
        }
      }  
    }  
    if (repMembrIndx != NULL) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        uint q,k,m;
        uint *evntProp = uivector(1, RF_eventTypeSize + 1);
        for (q=1; q <= RF_eventTypeSize + 1; q++) {
          evntProp[q] = 0;
        }
        for (i = 1; i <= (*nonMissMembrSize); i++) {
          m = (uint) RF_status[treeID][repMembrIndx[(*nonMissMembrIndx)[i]]];
          if (m > 0) {
            evntProp[RF_eventTypeIndex[m]] ++;
          }
          else {
            evntProp[RF_eventTypeSize + 1] ++;
          }
        }
        k = 0;
        for (q = 1; q <= RF_eventTypeSize + 1; q++) {
          if(evntProp[q] > 0) {
            k ++;
          }
        }
        if (k == 0) {
          result = FALSE;
        }
        else {
          if (k == 1) {
            if (evntProp[RF_eventTypeSize + 1] > 0) {
              result = FALSE;
            }
            else {
              result = getVariance(repMembrSize,
                                   repMembrIndx,
                                   *nonMissMembrSize,
                                   *nonMissMembrIndx,
                                   RF_time[treeID],
                                   NULL,
                                   NULL);
            }
          }
        }
        free_uivector(evntProp, 1, RF_eventTypeSize + 1);
      }
      else {
        if (RF_rSize == 1) {
          result = getVariance(repMembrSize,
                               repMembrIndx,
                               *nonMissMembrSize,
                               *nonMissMembrIndx,
                               RF_response[treeID][1],
                               NULL,
                               NULL);
        }
      }
    }
    if (!result) {
      (*nonMissMembrSize) = 0;
      free_uivector(*nonMissMembrIndx, 1, repMembrSize);
      free_dvector(*permissibleSplit, 1, repMembrSize);
    }
  }
  return result;
}
char updateMaximumSplit(uint    treeID,
                        double  delta,
                        uint    randomCovariate,
                        uint    index,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        uint    repMembrSize,
                        char   *localSplitIndicator,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        double *splitValueMaxCont,
                        uint   *splitValueMaxFactSize,
                        uint  **splitValueMaxFactPtr,
                        void   *permissibleSplitPtr,
                        char  **splitIndicator) {
  char flag;
  uint k;
  delta = delta * RF_splitWeight[randomCovariate];
  if(ISNA(*deltaMax)) {
    flag = TRUE;
  }
  else {
    if (delta > *deltaMax) {
      flag = TRUE;
    }
    else {
      flag = FALSE;
    }
  }
  if (flag) {
    *deltaMax = delta;
    *splitParameterMax = randomCovariate;
    if (factorFlag == TRUE) {
      if (*splitValueMaxFactSize > 0) {
        if (*splitValueMaxFactSize != mwcpSizeAbsolute) {
          free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
          *splitValueMaxFactSize = mwcpSizeAbsolute;
          *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
        }
      }
      else {
        *splitValueMaxFactSize = mwcpSizeAbsolute;
        *splitValueMaxFactPtr = uivector(1, *splitValueMaxFactSize);
      }
      *splitValueMaxCont = NA_REAL;
      for (k=1; k <= *splitValueMaxFactSize; k++) {
        (*splitValueMaxFactPtr)[k] =
          ((uint*) permissibleSplitPtr + ((index - 1) * (*splitValueMaxFactSize)))[k];
      }
    }
    else {
      if (*splitValueMaxFactSize > 0) {
        free_uivector(*splitValueMaxFactPtr, 1, *splitValueMaxFactSize);
        *splitValueMaxFactSize = 0;
        *splitValueMaxFactPtr = NULL;
      }
      else {
      }
      *splitValueMaxCont = ((double*) permissibleSplitPtr)[index];
    }
    if (*splitIndicator == NULL) {
     *splitIndicator = cvector(1, repMembrSize);
    }
   for (k=1; k <= repMembrSize; k++) {
     (*splitIndicator)[k] = localSplitIndicator[k];
   }
  }
  else {
  }
  return flag;
}
