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
#include    "regression.h"
#include     "splitSurv.h"
char randomSurvivalSplit(uint    treeID, 
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
  uint priorMembrIter, currentMembrIter;
  uint leftSize, rghtSize;
  uint censProp, evntProp;
  char *localSplitIndicator;
  double deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute       = 0;  
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
    censProp = evntProp = 0;
    for (k=1; k <= repMembrSize; k++) {
      if (RF_status[treeID][repMembrIndx[k]] > 0) {
        evntProp ++;
      }
      else {
        censProp ++;
      }
    }
    if (censProp == 0) {
      result = getVariance(repMembrSize, repMembrIndx, RF_time[treeID], NULL, NULL);
    }
    else {
      if (evntProp == 0) {
        result = FALSE;
      }
    }
  }
  if (result) {
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
      for (j = 1; j < splitLength; j++) {
        virtuallySplitNodeNew(treeID,
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
        if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
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
char logRank (uint    treeID, 
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
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  uint *localEventTimeCount;
  uint *localEventTimeIndex;
  uint localEventTimeSize;
  uint censProp, evntProp;
  char *localSplitIndicator;
  double delta, deltaNum, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  stackSplitTime(& localEventTimeCount, & localEventTimeIndex);
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
    censProp = evntProp = 0;
    for (k=1; k <= repMembrSize; k++) {
      if (RF_status[treeID][repMembrIndx[k]] > 0) {
        evntProp ++;
      }
      else {
        censProp ++;
      }
    }
    if (censProp == 0) {
      result = getVariance(repMembrSize, repMembrIndx, RF_time[treeID], NULL, NULL);
    }
    else {
      if (evntProp == 0) {
        result = FALSE;
      }
    }
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = getEventTimeSize(treeID, 
                                          parent, 
                                          repMembrIndx, 
                                          repMembrSize, 
                                          localEventTimeCount, 
                                          localEventTimeIndex);
    stackSplitCompactEventAndRisk(localEventTimeSize,
                                  & nodeParentEvent,
                                  & nodeParentAtRisk,
                                  & nodeLeftEvent,
                                  & nodeLeftAtRisk,
                                  & nodeRightEvent,
                                  & nodeRightAtRisk);
    getEventAndRisk(treeID,
                    repMembrIndx, 
                    repMembrSize, 
                    localEventTimeCount,
                    localEventTimeIndex,
                    localEventTimeSize,
                    nodeParentEvent,
                    nodeParentAtRisk);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
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
        for (k = 1; k <= repMembrSize; k++) {
          localSplitIndicator[k] = RIGHT;
        }
        for (m = 1; m <= localEventTimeSize; m++) {
          nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
        }
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        leftSize = virtuallySplitNodeNew(treeID,
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
        if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
            }
            for (k = 1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                index = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
                    index = m;
                    nodeLeftAtRisk[index] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][repMembrIndx[k]] > 0) {
                  nodeLeftEvent[index] ++;
                }
              }
              else {
              }
            } 
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              index = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[repMembrIndxx[i][k]]]) {
                  index = m;
                  nodeLeftAtRisk[index] ++;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              if (RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]] > 0) {
                nodeLeftEvent[index] ++;
              }
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          delta = deltaNum = deltaDen =  0.0;
          for (k=1; k <= localEventTimeSize; k++) {
            deltaNum = deltaNum + ((double) nodeLeftEvent[k] - ((double) ( nodeLeftAtRisk[k] * nodeParentEvent[k]) / nodeParentAtRisk[k]));
            if (nodeParentAtRisk[k] >= 2) {
              deltaDen = deltaDen + (
                                     ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k]) *
                                     (1.0 - ((double) nodeLeftAtRisk[k] / nodeParentAtRisk[k])) *
                                     ((double) (nodeParentAtRisk[k] - nodeParentEvent[k]) / (nodeParentAtRisk[k] - 1)) * nodeParentEvent[k]
                                     );
            }
          }
          rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
          leftSizeIter = currentMembrIter - 1; 
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNum / deltaDen;
            }
          }
          else {
            delta = deltaNum / deltaDen;
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
    unstackSplitCompactEventAndRisk(localEventTimeSize,
                                    nodeParentEvent,
                                    nodeParentAtRisk,
                                    nodeLeftEvent,
                                    nodeLeftAtRisk,
                                    nodeRightEvent,
                                    nodeRightAtRisk);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackSplitTime(localEventTimeCount, localEventTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char logRankScore(uint    treeID, 
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
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  uint censProp, evntProp;
  char *localSplitIndicator;
  double delta, deltaNum, deltaNumAdj, deltaDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint i, j, k;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  deltaNum               = 0;  
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
    censProp = evntProp = 0;
    for (k=1; k <= repMembrSize; k++) {
      if (RF_status[treeID][repMembrIndx[k]] > 0) {
        evntProp ++;
      }
      else {
        censProp ++;
      }
    }
    if (censProp == 0 ) {
      result = getVariance(repMembrSize, repMembrIndx, RF_time[treeID], NULL, NULL);
    }
    else {
      if (evntProp == 0 ) {
        result = FALSE;
      }
    }
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    double  meanSurvRank;
    double  varSurvRank;
    uint   *survivalTimeIndexRank = uivector(1, repMembrSize);
    double *survivalRank = dvector(1, repMembrSize);
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
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
      for (k = 1; k <= repMembrSize; k++) {
        survivalTimeIndexRank[k] = 0;
        for (j = 1; j <= repMembrSize; j++) {
          if ( RF_masterTimeIndex[treeID][repMembrIndx[j]]  <= RF_masterTimeIndex[treeID][repMembrIndx[repMembrIndxx[i][k]]] ) {
            survivalTimeIndexRank[k] ++;
          }
        }
      }
      meanSurvRank = varSurvRank = 0;
      for (k = 1; k <= repMembrSize; k++) {
        survivalRank[k] = 0;
        for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
          survivalRank[k] = survivalRank[k] + (RF_status[treeID][repMembrIndx[repMembrIndxx[i][j]]] / (repMembrSize - survivalTimeIndexRank[j] + 1) );
        }
        survivalRank[k] = RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]] - survivalRank[k];
        meanSurvRank = meanSurvRank + survivalRank[k];
        varSurvRank = varSurvRank +  pow(survivalRank[k], 2.0);
      }
      varSurvRank = ( varSurvRank - (pow(meanSurvRank, 2.0) / repMembrSize) ) / (repMembrSize - 1);
      meanSurvRank = meanSurvRank / repMembrSize;
      if (factorFlag == FALSE) {
        for (k = 1; k <= repMembrSize; k++) {
          localSplitIndicator[k] = RIGHT;
        }
        deltaNum = 0.0;
        leftSizeIter = 0;
        rghtSizeIter = repMembrSize;
      }
      for (j = 1; j < splitLength; j++) {
        if (factorFlag == TRUE) {
          priorMembrIter = 0;
          leftSize = 0;
        }
        leftSize = virtuallySplitNodeNew(treeID,
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
        if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            deltaNum = 0.0;
            for (k = 1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                deltaNum = deltaNum + survivalRank[k];
              }
            }  
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              deltaNum = deltaNum + survivalRank[repMembrIndxx[i][k]];
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          deltaNumAdj  = deltaNum - (leftSize * meanSurvRank);
          deltaDen     = leftSize * (1.0 - (leftSize / repMembrSize)) * varSurvRank;
          deltaNumAdj = fabs(deltaNumAdj);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNumAdj <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNumAdj / deltaDen;
            }
          }
          else {
            delta = deltaNumAdj / deltaDen;
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
        if (factorFlag == FALSE) {
          if (rghtSize < RF_minimumNodeSize) {
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
    free_uivector(survivalTimeIndexRank, 1, repMembrSize);
    free_dvector(survivalRank, 1, repMembrSize);
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
char logRankLauCR (uint    treeID,
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
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  uint *localEventTimeCount;
  uint *localEventTimeIndex;
  uint localEventTimeSize;
  uint *evntProp;
  char *localSplitIndicator;
  double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m, q, r, s;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  stackSplitTime(& localEventTimeCount, & localEventTimeIndex);
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
    evntProp = uivector(1, RF_eventTypeSize + 1);
    for (q=1; q <= RF_eventTypeSize + 1; q++) {
      evntProp[q] = 0;
    }
    for (k=1; k <= repMembrSize; k++) {
      m = (uint) RF_status[treeID][repMembrIndx[k]];
      if (m > 0) {
        evntProp[RF_eventTypeIndex[m]] ++;
      }
      else {
        evntProp[RF_eventTypeSize + 1] ++;
      }
    }
    k = 0;
    for (q=1; q <= RF_eventTypeSize + 1; q++) {
      if(evntProp[q] > 0) {
        k ++;
      }
    }
    if (k == 1) {
      if (evntProp[RF_eventTypeSize + 1] > 0) {
        result = FALSE;
      }
      else {
        result = getVariance(repMembrSize, repMembrIndx, RF_time[treeID], NULL, NULL);
      }
    }
    free_uivector(evntProp, 1, RF_eventTypeSize + 1);
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint *nodeParentEvent, *nodeLeftEvent, *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = getEventTimeSize(treeID, 
                                          parent, 
                                          repMembrIndx, 
                                          repMembrSize, 
                                          localEventTimeCount, 
                                          localEventTimeIndex);
    stackSplitCompactEventAndRisk(localEventTimeSize,
                                  & nodeParentEvent,
                                  & nodeParentAtRisk,
                                  & nodeLeftEvent,
                                  & nodeLeftAtRisk,
                                  & nodeRightEvent,
                                  & nodeRightAtRisk);
    getEventAndRisk(treeID,
                    repMembrIndx, 
                    repMembrSize, 
                    localEventTimeCount,
                    localEventTimeIndex,
                    localEventTimeSize,
                    nodeParentEvent,
                    nodeParentAtRisk);
    uint **nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    uint **nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    uint **nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    uint **nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    for (m=1; m <= localEventTimeSize; m++) {
      for (q = 1; q <= RF_eventTypeSize; q++) {
        nodeParentEventCR[q][m] = 0;
      }
    }
    for (k=1; k <= repMembrSize; k++) {
      if (RF_status[treeID][repMembrIndx[k]] > 0) {
        index = 0;  
        for (m = 1; m <= localEventTimeSize; m++) {
          if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
            index = m;
          }
          else {
            m = localEventTimeSize;
          }
        }
        nodeParentEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[k]]]][index] ++;
      }
    }
    for (m=1; m <= localEventTimeSize; m++) {
      for (q = 1; q <= RF_eventTypeSize; q++) {
        if (RF_crWeight[q] > 0) {
          nodeParentInclusiveAtRisk[q][m] = nodeParentAtRisk[m];
          for (s = 1; s < m; s++) {
            for (r = 1; r <= RF_eventTypeSize; r++) {
              if (q != r) {
                nodeParentInclusiveAtRisk[q][m]  += nodeParentEventCR[r][s];
              }
            }
          }
        }
      }
    }
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
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
        for (k = 1; k <= repMembrSize; k++) {
          localSplitIndicator[k] = RIGHT;
        }
        for (m = 1; m <= localEventTimeSize; m++) {
          nodeLeftAtRisk[m] = 0;
          for (q = 1; q <= RF_eventTypeSize; q++) {
            nodeLeftEventCR[q][m] = 0;
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
        leftSize = virtuallySplitNodeNew(treeID,
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
        if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftAtRisk[m] = 0;
              for (q = 1; q <= RF_eventTypeSize; q++) {
                nodeLeftEventCR[q][m] = 0;
              }
            }
            for (k=1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                index = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
                    index = m;
                    nodeLeftAtRisk[index] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][repMembrIndx[k]] > 0) {
                  nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[k]]]][index] ++;
                }
              }
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              index = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[repMembrIndxx[i][k]]]) {
                  index = m;
                  nodeLeftAtRisk[index] ++;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              if (RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]] > 0) {
                nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]]]][index] ++;
              }
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          for (m=1; m <= localEventTimeSize; m++) {
            for (q = 1; q <= RF_eventTypeSize; q++) {
              if (RF_crWeight[q] > 0) {
                nodeLeftInclusiveAtRisk[q][m]   = nodeLeftAtRisk[m];
                for (s = 1; s < m; s++) {
                  for (r = 1; r <= RF_eventTypeSize; r++) {
                    if (q != r) {
                      nodeLeftInclusiveAtRisk[q][m]    += nodeLeftEventCR[r][s];
                    }
                  }
                }
              }
            }
          }
          delta = deltaNum = deltaDen =  0.0;
          for (q = 1; q <= RF_eventTypeSize; q++) {
            if (RF_crWeight[q] > 0) {
              deltaSubNum = 0;
              for (m=1; m <= localEventTimeSize; m++) {
                deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])));
              }
              deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
              deltaSubDen = 0;
              for (m=1; m <= localEventTimeSize; m++) {
                if (nodeParentAtRisk[m] >= 2) {
                  deltaSubDen = deltaSubDen  + (
                                                (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                (1.0 - ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])) *
                                                ((double) (nodeParentInclusiveAtRisk[q][m] - nodeParentEventCR[q][m]) / (nodeParentInclusiveAtRisk[q][m] - 1))
                                                );
                }
              }
              deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
            }
          }
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNum / deltaDen;
            }
          }
          else {
            delta = deltaNum / deltaDen;
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
    unstackSplitCompactEventAndRisk(localEventTimeSize,
                                    nodeParentEvent,
                                    nodeParentAtRisk,
                                    nodeLeftEvent,
                                    nodeLeftAtRisk,
                                    nodeRightEvent,
                                    nodeRightAtRisk);
    free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
    free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
    free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
    free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackSplitTime(localEventTimeCount, localEventTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
char logRankCR (uint    treeID,
                Node   *parent,
                uint   *repMembrIndx,
                uint    repMembrSize,
                uint   *allMembrIndx,
                uint    allMembrSize,
                uint   *splitParameterMax,
                double *splitValueMaxCont,
                uint   *splitValueMaxFactSize,
                uint  **splitValueMaxFactPtr,
                double  *splitStatistic) {
  uint    *randomCovariateIndex;
  double **permissibleSplit;
  uint    *permissibleSplitSize;
  uint   **repMembrIndxx;
  uint priorMembrIter, currentMembrIter;
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  uint *localEventTimeCount;
  uint *localEventTimeIndex;
  uint localEventTimeSize;
  uint *evntProp;
  char *localSplitIndicator;
  double delta, deltaNum, deltaSubNum, deltaDen, deltaSubDen, deltaMax;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  uint index;
  uint i, j, k, m, q;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  stackSplitTime(& localEventTimeCount, & localEventTimeIndex);
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
    evntProp = uivector(1, RF_eventTypeSize + 1);
    for (q=1; q <= RF_eventTypeSize + 1; q++) {
      evntProp[q] = 0;
    }
    for (k=1; k <= repMembrSize; k++) {
      m = (uint) RF_status[treeID][repMembrIndx[k]];
      if (m > 0) {
        evntProp[RF_eventTypeIndex[m]] ++;
      }
      else {
        evntProp[RF_eventTypeSize + 1] ++;
      }
    }
    k = 0;
    for (q=1; q <= RF_eventTypeSize + 1; q++) {
      if(evntProp[q] > 0) {
        k ++;
      }
    }
    if (k == 1) {
      if (evntProp[RF_eventTypeSize + 1] > 0) {
        result = FALSE;
      }
      else {
        result = getVariance(repMembrSize, repMembrIndx, RF_time[treeID], NULL, NULL);
      }
    }
    free_uivector(evntProp, 1, RF_eventTypeSize + 1);
  }
  if (result) {
    stackSplitIndicator(repMembrSize, & localSplitIndicator);
    uint *nodeParentEvent, *nodeLeftEvent, *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    localEventTimeSize = getEventTimeSize(treeID, 
                                          parent, 
                                          repMembrIndx, 
                                          repMembrSize, 
                                          localEventTimeCount, 
                                          localEventTimeIndex);
    stackSplitCompactEventAndRisk(localEventTimeSize,
                                  & nodeParentEvent,
                                  & nodeParentAtRisk,
                                  & nodeLeftEvent,
                                  & nodeLeftAtRisk,
                                  & nodeRightEvent,
                                  & nodeRightAtRisk);
    getEventAndRisk(treeID,
                    repMembrIndx, 
                    repMembrSize, 
                    localEventTimeCount,
                    localEventTimeIndex,
                    localEventTimeSize,
                    nodeParentEvent,
                    nodeParentAtRisk);
    uint **nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    uint **nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
    for (m=1; m <= localEventTimeSize; m++) {
      for (q = 1; q <= RF_eventTypeSize; q++) {
        nodeParentEventCR[q][m] = 0;
      }
    }
    for (k=1; k <= repMembrSize; k++) {
      if (RF_status[treeID][repMembrIndx[k]] > 0) {
        index = 0;  
        for (m = 1; m <= localEventTimeSize; m++) {
          if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
            index = m;
          }
          else {
            m = localEventTimeSize;
          }
        }
        nodeParentEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[k]]]][index] ++;
      }
    }
    uint actualCovariateCount = stackAndSelectRandomCovariates(treeID,
                                                               parent, 
                                                               repMembrIndx, 
                                                               repMembrSize, 
                                                               & randomCovariateIndex, 
                                                               & permissibleSplit, 
                                                               & permissibleSplitSize,
                                                               & repMembrIndxx);
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
        for (k = 1; k <= repMembrSize; k++) {
          localSplitIndicator[k] = RIGHT;
        }
        for (m = 1; m <= localEventTimeSize; m++) {
          nodeLeftAtRisk[m] = 0;
          for (q = 1; q <= RF_eventTypeSize; q++) {
            nodeLeftEventCR[q][m] = 0;
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
        leftSize = virtuallySplitNodeNew(treeID,
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
        if ((leftSize >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftAtRisk[m] = 0;
              for (q = 1; q <= RF_eventTypeSize; q++) {
                nodeLeftEventCR[q][m] = 0;
              }
            }
            for (k=1; k <= repMembrSize; k++) {
              if (localSplitIndicator[k] == LEFT) {
                index = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[k]]) {
                    index = m;
                    nodeLeftAtRisk[index] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][repMembrIndx[k]] > 0) {
                  nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[k]]]][index] ++;
                }
              }
            }
          }
          else {
            for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
              index = 0;  
              for (m = 1; m <= localEventTimeSize; m++) {
                if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][repMembrIndx[repMembrIndxx[i][k]]]) {
                  index = m;
                  nodeLeftAtRisk[index] ++;
                }
                else {
                  m = localEventTimeSize;
                }
              }
              if (RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]] > 0) {
                nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][repMembrIndx[repMembrIndxx[i][k]]]]][index] ++;
              }
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1; 
          }
          delta = deltaNum = deltaDen =  0.0;
          for (q = 1; q <= RF_eventTypeSize; q++) {
            if (RF_crWeight[q] > 0) {
              deltaSubNum = 0;
              for (m=1; m <= localEventTimeSize; m++) {
                deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])));
              }
              deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
              deltaSubDen = 0;
              for (m=1; m <= localEventTimeSize; m++) {
                if (nodeParentAtRisk[m] >= 2) {
                  deltaSubDen = deltaSubDen  + (
                                                (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                (1.0 - ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])) *
                                                ((double) (nodeParentAtRisk[m] - nodeParentEventCR[q][m]) / (nodeParentAtRisk[m] - 1))
                                                );
                }
              }
              deltaDen = deltaDen + (RF_crWeight[q] * RF_crWeight[q] * deltaSubDen);
            }
          }
          deltaNum = fabs(deltaNum);
          deltaDen = sqrt(deltaDen);
          if (deltaDen <= EPSILON) {
            if (deltaNum <= EPSILON) {
              delta = 0.0;
            }
            else {
              delta = deltaNum / deltaDen;
            }
          }
          else {
            delta = deltaNum / deltaDen;
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
    unstackSplitCompactEventAndRisk(localEventTimeSize,
                                    nodeParentEvent,
                                    nodeParentAtRisk,
                                    nodeLeftEvent,
                                    nodeLeftAtRisk,
                                    nodeRightEvent,
                                    nodeRightAtRisk);
    free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
    free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
  }  
  unstackSplitTime(localEventTimeCount, localEventTimeIndex);
  result = summarizeSplitResult(*splitParameterMax, 
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                splitStatistic,
                                deltaMax);
  return result;
}
