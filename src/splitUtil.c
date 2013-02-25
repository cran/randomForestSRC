////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.1.0
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
////    email:  ubk@kogalur.com
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
void updateMaximumSplit(double  delta, 
                        uint    randomCovariate,
                        uint    index,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        double *splitValueMaxCont,
                        uint   *splitValueMaxFactSize,
                        uint  **splitValueMaxFactPtr,
                        void   *permissibleSplitPtr) {
  uint k;
  if (delta > *deltaMax) {
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
  }
  else {
  }
}
uint stackAndSelectRandomCovariates(uint     treeID,
                                    Node     *parent,
                                    uint     *repMembrIndx,
                                    uint      repMembrSize,
                                    uint    **covariateIndex,
                                    double ***permissibleSplit,
                                    uint    **permissibleSplitSize) {
  uint i;
  uint actualCovariateCount;
  uint candidateCovariate;
  *covariateIndex = uivector(1, RF_xSize);
  *permissibleSplit = dmatrix(1, RF_randomCovariateCount, 1, repMembrSize);
  *permissibleSplitSize = uivector(1, RF_randomCovariateCount);
  char *randomSplitVector = cvector(1, RF_xSize);
  if (repMembrSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid repMembrSize encountered in stackAndSelectRandomCovariates():  %10d", repMembrSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nrCopyVector(randomSplitVector, parent -> permissibleSplit, RF_xSize);
  for(i=1; i <= RF_randomCovariateCount; i++) {
    (*covariateIndex)[i] = 0;
  }
  actualCovariateCount =  1;
  candidateCovariate   = -1;
  while ((actualCovariateCount  <= RF_randomCovariateCount) && (candidateCovariate != 0)) {
    candidateCovariate = getSelectableElement(treeID, RF_xSize, randomSplitVector, RF_randomCovariateWeight);
    if (candidateCovariate != 0) {
      for (i=1; i <= repMembrSize; i++) {
        (*permissibleSplit)[actualCovariateCount][i] = RF_observation[treeID][candidateCovariate][repMembrIndx[i]];
      }
      if ( (strcmp(RF_xType[candidateCovariate], "C") != 0) && (RF_opt & OPT_SPLT_FAST) ) {
          (*permissibleSplitSize)[actualCovariateCount] = repMembrSize;
      }
      else {
        (*permissibleSplitSize)[actualCovariateCount] = 1;
        hpsort((*permissibleSplit)[actualCovariateCount], repMembrSize);
        for (i = 2; i <= repMembrSize; i++) {
          if ((*permissibleSplit)[actualCovariateCount][i] > (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]]) {
            (*permissibleSplitSize)[actualCovariateCount] ++;
            (*permissibleSplit)[actualCovariateCount][(*permissibleSplitSize)[actualCovariateCount]] = (*permissibleSplit)[actualCovariateCount][i];
          }
        }
      }
      if((*permissibleSplitSize)[actualCovariateCount] >= 2) {
        randomSplitVector[candidateCovariate] = ACTIVE;
        (*covariateIndex)[actualCovariateCount] = candidateCovariate;
        actualCovariateCount ++;
      }
      else {
        (parent -> permissibleSplit)[candidateCovariate] = FALSE;
        randomSplitVector[candidateCovariate] = FALSE;
      }
    }
  }
  actualCovariateCount --;
  free_cvector(randomSplitVector, 1, RF_xSize);
  return actualCovariateCount;
}
void unstackRandomCovariates(uint     treeID,
                             uint     repMembrSize, 
                             uint    *covariateIndex,
                             double **permissibleSplit,
                             uint    *permissibleSplitSize) {
  free_uivector(covariateIndex, 1, RF_xSize);
  free_dmatrix(permissibleSplit, 1, RF_randomCovariateCount, 1, repMembrSize);
  free_uivector(permissibleSplitSize, 1, RF_randomCovariateCount);
}
uint getSelectableElement (uint    treeID,
                           uint    length,
                           char   *permissible,
                           double *weight) {
  char   *localPermissible = NULL;  
  double *cdf = NULL;  
  uint selectableCount;
  uint covariateIndex;
  double randomValue;
  uint i, j, k, p, index;
  if (length > 0) {
    localPermissible = cvector(1, length);
    cdf = dvector(1, length);
  }
  selectableCount = 0;
  for (i=1; i <= length; i++) {
    if (permissible[i] == TRUE) {
      if (weight != NULL) {
        if (weight[i] > 0) {
          localPermissible[i] = TRUE;
          selectableCount ++;
        }
        else {
          localPermissible[i] = FALSE;
        }
      }
      else {
        localPermissible[i] = TRUE;
        selectableCount ++;
      }
    }
    else {
      localPermissible[i] = FALSE;
    }
  }
  if (selectableCount > 0) {
    if (weight != NULL) { 
      covariateIndex = 0;
      for (k=1; k <= RF_xSize; k++) {
        if (localPermissible[k] == TRUE) {
          cdf[++covariateIndex] = weight[k];
        }
      }
      for (k=2; k <= covariateIndex; k++) {
        cdf[k] += cdf[k-1];
      }
      randomValue = ran2(treeID) * cdf[covariateIndex];
      j=1;
      while (randomValue > cdf[j]) {
        j++;
      }
      for (index = 1; j > 0; index++) {
        if (localPermissible[index] == TRUE) {
          j--;
        }
      }
      index --;
    }
    else {
      p = (uint) ceil(ran2(treeID) * (selectableCount * 1.0));
      index = 1;
      while (p > 0) {
        if (permissible[index] == TRUE) {
          p --;
        }
        index ++;
      }
      index --;
    }
  }  
  else {
    index = 0;
  }
  if (length > 0) {
    free_cvector(localPermissible, 1, length);
    free_dvector(cdf, 1, length);
  }
  return index;
}
uint getEventTimeSize(uint   treeID, 
                      Node   *parent, 
                      uint   *repMembrIndx,
                      uint    repMembrSize,
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
  for (i=1; i <= repMembrSize; i++) {
    if (RF_status[treeID][repMembrIndx[i]] > 0) {
      localEventTimeCount[RF_masterTimeIndex[treeID][repMembrIndx[i]]] ++;
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
void stackSplitTime(uint **localEventTimeCount, 
                    uint **localEventTimeIndex) {
  if ((localEventTimeCount != NULL) && (localEventTimeIndex != NULL)) {
    *localEventTimeCount = uivector(1, RF_masterTimeSize);
    *localEventTimeIndex = uivector(1, RF_masterTimeSize);
  }
}
void unstackSplitTime(uint *localEventTimeCount, 
                      uint *localEventTimeIndex) {
  if ((localEventTimeCount != NULL) && (localEventTimeIndex != NULL)) {
    free_uivector(localEventTimeCount, 1, RF_masterTimeSize);
    free_uivector(localEventTimeIndex, 1, RF_masterTimeSize);
  }
}
void stackSplitCompactEventAndRisk(uint   eventTimeSize,
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
void unstackSplitCompactEventAndRisk(uint  eventTimeSize,
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
void getEventAndRisk(uint    treeID,
                uint   *repMembrIndx,
                uint    repMembrSize,
                uint   *localEventTimeCount,
                uint   *localEventTimeIndex,
                uint    localEventTimeSize,
                uint   *nodeParentEvent,
                uint   *nodeParentAtRisk) {
  uint i, j;
  for (i=1; i <= localEventTimeSize; i++) {
    nodeParentAtRisk[i] = 0;
    nodeParentEvent[i] = localEventTimeCount[localEventTimeIndex[i]];
    for (j=1; j <= repMembrSize; j++) {
      if (localEventTimeIndex[i] <= RF_masterTimeIndex[treeID][repMembrIndx[j]]) {
        nodeParentAtRisk[i] ++;
      }
    }
  }
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
  uint j, k, j2, k2;
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
      splitLength = 1 + ((RF_splitRandomRule <= repMembrSize) ? RF_splitRandomRule : repMembrSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomRule == 0) {
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
          if (*((uint*) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) <= ((RF_splitRandomRule <= repMembrSize) ? RF_splitRandomRule : repMembrSize)) {
            splitLength = *((uint*) RF_factorList[treeID][permissibleSplitSize] -> complementaryPairCount) + 1;
            *deterministicSplitFlag = TRUE;
          }
        }
        if (*deterministicSplitFlag == FALSE) {
          splitLength = 1 + ((RF_splitRandomRule <= repMembrSize) ? RF_splitRandomRule : repMembrSize);
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
      splitLength = 1 + ((RF_splitRandomRule <= repMembrSize) ? RF_splitRandomRule : repMembrSize);
      *deterministicSplitFlag = FALSE;
    }
    else {
      if(RF_splitRandomRule == 0) {
        splitLength = permissibleSplitSize;
        (*permissibleSplitPtr) = permissibleSplit;
        *deterministicSplitFlag = TRUE;
      }
      else {
        if (permissibleSplitSize <= RF_splitRandomRule) {
          splitLength = permissibleSplitSize;
          (*permissibleSplitPtr) = permissibleSplit;
          *deterministicSplitFlag = TRUE;
        }
        else {
          splitLength = RF_splitRandomRule + 1;
          *deterministicSplitFlag = FALSE;
        }
      }  
    }  
    if (*deterministicSplitFlag == FALSE) {
      (*permissibleSplitPtr) = dvector(1, splitLength);
      ((double*) (*permissibleSplitPtr))[splitLength] = 0;
      for (j = 1; j < splitLength; j++) {
        k = permissibleSplitSize - 1;
        ((double*) (*permissibleSplitPtr))[j]  = permissibleSplit[(uint) ceil(ran2(treeID) * (k * 1.0))];
      }
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
uint virtuallySplitNode(uint  treeID,
                        char  factorFlag,
                        uint  mwcpSizeAbsolute,
                        uint  randomCovariate,
                        uint *repMembrIndx,
                        uint  repMembrSize,
                        void *permissibleSplitPtr,
                        uint  offset,
                        uint  localEventTimeSize,
                        uint *localEventTimeIndex,
                        uint *nodeParentAtRisk,
                        uint *nodeParentEvent,
                        uint *nodeLeftAtRisk,
                        uint *nodeLeftEvent,
                        uint *leftEventTimeSize,
                        uint *nodeRightAtRisk,
                        uint *nodeRightEvent,
                        uint *rightEventTimeSize,
                        char *localSplitIndicator) {
  char daughterFlag;
  uint leftSize;
  uint index, k, m;
  leftSize = 0;
  if (localEventTimeSize > 0) {
    *leftEventTimeSize = *rightEventTimeSize = 0;
    for (k=1; k <= localEventTimeSize; k++) {
      nodeLeftEvent[k] = nodeLeftAtRisk[k] = 0;
    }
  }
  for (k=1; k <= repMembrSize; k++) {
    daughterFlag = RIGHT;
    if (factorFlag == TRUE) {
      daughterFlag = splitOnFactor((uint) RF_observation[treeID][randomCovariate][repMembrIndx[k]], (uint*) permissibleSplitPtr + ((offset - 1) * mwcpSizeAbsolute));
    }
    else {
      if (RF_observation[treeID][randomCovariate][repMembrIndx[k]] <= ((double*) permissibleSplitPtr)[offset]) {
        daughterFlag = LEFT;
      }
    }
    if (localSplitIndicator != NULL) {
      localSplitIndicator[k] = daughterFlag;
    }
    if (daughterFlag == LEFT) {
      leftSize ++;
      if (localEventTimeSize > 0) {
        index = 0;  
        for (m = 1; m <= localEventTimeSize; m++) {
          if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
            nodeLeftAtRisk[m] ++;
            index = m;
          }
          else {
            m = localEventTimeSize;
          }
        }
        if (RF_status[treeID][repMembrIndx[k]] > 0) {
          nodeLeftEvent[index] ++;
        }
      }
    }  
    else {
    }
  }  
  if (localEventTimeSize > 0) {
    for (k=1; k <= localEventTimeSize; k++) {
      nodeRightEvent[k] = nodeParentEvent[k] - nodeLeftEvent[k];
      nodeRightAtRisk[k] = nodeParentAtRisk[k] - nodeLeftAtRisk[k];
      if (nodeLeftEvent[k] > 0) {
        (*leftEventTimeSize) ++;
      }
      if (nodeRightEvent[k] > 0) {
        (*rightEventTimeSize) ++;
      }
    }
  }
  return (leftSize);
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
  randomGroupIndex = (uint) ceil(ran2(treeID) * ((RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount) * 1.0));
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
  randomValue = ceil((ran2(treeID) * cdf[RF_factorList[treeID][relativeFactorSize] -> cardinalGroupCount]));
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
    randomLevel[k] = getSelectableElement(treeID, relativeFactorSize, localPermissible, NULL);
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
char summarizeSplitResult(uint   splitParameterMax, 
                          double splitValueMaxCont,
                          uint   splitValueMaxFactSize,
                          uint  *splitValueMaxFactPtr,
                          double deltaMax) {
  char result;
  if (splitParameterMax > 0) {
    result = TRUE;
  }
  else {
    result = FALSE;
  }
  return result;
}
char getStandardDeviation(uint repSize, uint *repIndx, double *target) {
  uint i;
  double meanResult, sdResult;
  char result;
  meanResult = 0.0;
  for (i=1; i <= repSize; i++) {
      meanResult += target[repIndx[i]];
    }
    meanResult = meanResult / (double) repSize;
    sdResult = 0.0;
    for (i=1; i <= repSize; i++) {
      sdResult += pow(meanResult - target[repIndx[i]], 2.0);
    }
    sdResult = sdResult / (double) repSize;
    result = ((sdResult <= EPSILON) ? FALSE : TRUE);
    return(result);
}
