////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.0
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
char logRankNCR (uint    treeID,
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
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  double delta, deltaMax;
  uint j, k, m;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             repMembrIndx,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if (result) {
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
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    uint   *survivalTimeIndexRank;
    double *survivalRank;
    double  meanSurvRank, varSurvRank;
    double deltaNum, deltaNumAdj, deltaDen;
    uint   tIndx;
    meanSurvRank = varSurvRank = 0;  
    survivalTimeIndexRank = NULL;  
    survivalRank = NULL;  
    localEventTimeSize = 0;  
    delta = deltaNum = 0;  
    switch(RF_splitRule) {
    case SURV_LGRNK:
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndxStatic,
                             nonMissMembrSizeStatic,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
      }
      break;
    case SURV_LRSCR:
      survivalTimeIndexRank = uivector(1, repMembrSize);
      survivalRank = dvector(1, repMembrSize);
      break;
    default:
      break;
    }
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
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
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
      switch(RF_splitRule) {
      case SURV_LGRNK:
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          stackAndGetSplitSurv(treeID,
                               repMembrIndx,
                               repMembrSize,
                               nonMissMembrIndx,
                               nonMissMembrSize,
                               & localEventTimeCount,
                               & localEventTimeIndex,
                               & localEventTimeSize,
                               & nodeParentEvent,
                               & nodeParentAtRisk,
                               & nodeLeftEvent,
                               & nodeLeftAtRisk,
                               & nodeRightEvent,
                               & nodeRightAtRisk);
        }
        break;
      case SURV_LRSCR:
        break;
      default:
        break;
      }
      if (localEventTimeSize > 0) {
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          leftSizeIter = 0;
          rghtSizeIter = nonMissMembrSize;
          switch(RF_splitRule) {
          case SURV_LGRNK:
            for (m = 1; m <= localEventTimeSize; m++) {
              nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
            }
            break;
          case SURV_LRSCR:
            deltaNum =  0.0;
            break;
          default:
            break;
          }
        }
        switch(RF_splitRule) {
        case SURV_LGRNK:
          break;
        case SURV_LRSCR:
          for (k = 1; k <= nonMissMembrSize; k++) {
            survivalTimeIndexRank[k] = 0;
            for (j = 1; j <= nonMissMembrSize; j++) {
              if ( RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[j]] ]  <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] ) {
                survivalTimeIndexRank[k] ++;
              }
            }
          }
          meanSurvRank = varSurvRank = 0;
          for (k = 1; k <= nonMissMembrSize; k++) {
            survivalRank[k] = 0;
            for (j = 1; j <= survivalTimeIndexRank[k]; j++) {
              survivalRank[k] = survivalRank[k] + (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ] / (nonMissMembrSize - survivalTimeIndexRank[j] + 1) );
            }
            survivalRank[k] = RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] - survivalRank[k];
            meanSurvRank = meanSurvRank + survivalRank[k];
            varSurvRank = varSurvRank +  pow(survivalRank[k], 2.0);
          }
          varSurvRank = ( varSurvRank - (pow(meanSurvRank, 2.0) / nonMissMembrSize) ) / (nonMissMembrSize - 1);
          meanSurvRank = meanSurvRank / nonMissMembrSize;
          break;
        default:
          break;
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
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
            if (factorFlag == TRUE) {
              switch(RF_splitRule) {
              case SURV_LGRNK:
                for (m = 1; m <= localEventTimeSize; m++) {
                  nodeLeftEvent[m] = nodeLeftAtRisk[m] = 0;
                }
                for (k = 1; k <= nonMissMembrSize; k++) {
                  if (localSplitIndicator[  nonMissMembrIndx[indxx[k]]  ] == LEFT) {
                    tIndx = 0;  
                    for (m = 1; m <= localEventTimeSize; m++) {
                      if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                        tIndx = m;
                        nodeLeftAtRisk[tIndx] ++;
                      }
                      else {
                        m = localEventTimeSize;
                      }
                    }
                    if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                      nodeLeftEvent[tIndx] ++;
                    }
                  }
                  else {
                  }
                } 
                break;
              case SURV_LRSCR:
                deltaNum = 0.0;
                for (k = 1; k <= nonMissMembrSize; k++) {
                  if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                    deltaNum = deltaNum + survivalRank[k];
                  }
                }
                break;
              default:
                break;
              }
            }
            else {
              switch(RF_splitRule) {
              case SURV_LGRNK:
                for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                  tIndx = 0;  
                  for (m = 1; m <= localEventTimeSize; m++) {
                    if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                      tIndx = m;
                      nodeLeftAtRisk[tIndx] ++;
                    }
                    else {
                      m = localEventTimeSize;
                    }
                  }
                  if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                    nodeLeftEvent[tIndx] ++;
                  }
                }
                break;
              case SURV_LRSCR:
                for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                  deltaNum = deltaNum + survivalRank[ nonMissMembrIndx[indxx[k]] ];
                }
                break;
              default:
                break;
              }
              rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
              leftSizeIter = currentMembrIter - 1;
            }
            switch(RF_splitRule) {
            case SURV_LGRNK:
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
              break;
            case SURV_LRSCR:
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
              break;
            default:
              break;
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
      }  
      else {
      }
      unstackSplitVector(treeID,
                         permissibleSplitSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      switch(RF_splitRule) {
      case SURV_LGRNK:
        if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
          unstackSplitSurv(localEventTimeCount,
                           localEventTimeIndex,
                           localEventTimeSize,
                           nodeParentEvent,
                           nodeParentAtRisk,
                           nodeLeftEvent,
                           nodeLeftAtRisk,
                           nodeRightEvent,
                           nodeRightAtRisk);
        }
        break;
      case SURV_LRSCR:
        break;
      default:
        break;
      }
    }  
    switch(RF_splitRule) {
    case SURV_LGRNK:
      if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      }
    break;
    case SURV_LRSCR:
      free_uivector(survivalTimeIndexRank, 1, repMembrSize);
      free_dvector(survivalRank, 1, repMembrSize);
      break;
    default:
      break;
    }
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
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
  }  
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
  uint leftSizeIter, rghtSizeIter;
  uint leftSize, rghtSize;
  char *localSplitIndicator;
  uint splitLength;
  void *permissibleSplitPtr;
  char factorFlag;
  uint mwcpSizeAbsolute;
  char deterministicSplitFlag;
  char result;
  double delta, deltaMax;
  uint j, k, m;
  mwcpSizeAbsolute       = 0;  
  leftSizeIter           = 0;  
  rghtSizeIter           = 0;  
  *splitParameterMax     = 0;
  *splitValueMaxFactSize = 0;
  *splitValueMaxFactPtr  = NULL;
  *splitValueMaxCont     = NA_REAL;
  deltaMax               = NA_REAL;
  result = getPreSplitResult(treeID,
                             parent,
                             repMembrSize,
                             repMembrIndx,
                             & nonMissMembrSizeStatic,
                             & nonMissMembrIndxStatic,
                             & permissibleSplit,
                             multImpFlag);
  if (result) {
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
    uint *localEventTimeCount, *localEventTimeIndex;
    uint  localEventTimeSize;
    uint *nodeParentEvent,  *nodeLeftEvent,  *nodeRightEvent;
    uint *nodeParentAtRisk, *nodeLeftAtRisk, *nodeRightAtRisk;
    uint **nodeParentEventCR, **nodeLeftEventCR;
    uint **nodeParentInclusiveAtRisk, **nodeLeftInclusiveAtRisk;
    nodeParentInclusiveAtRisk = nodeLeftInclusiveAtRisk = NULL;  
    nodeParentEventCR = nodeLeftEventCR = NULL;  
    double deltaNum, deltaSubNum, deltaDen, deltaSubDen;
    uint   tIndx;
    uint   q, s, r;
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
      stackAndGetSplitSurv(treeID,
                           repMembrIndx,
                           repMembrSize,
                           nonMissMembrIndxStatic,
                           nonMissMembrSizeStatic,
                           & localEventTimeCount,
                           & localEventTimeIndex,
                           & localEventTimeSize,
                           & nodeParentEvent,
                           & nodeParentAtRisk,
                           & nodeLeftEvent,
                           & nodeLeftAtRisk,
                           & nodeRightEvent,
                           & nodeRightAtRisk);
      nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
      switch(RF_splitRule) {
      case SURV_CR_LAU:
        nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        break;
      case SURV_CR_LOG:
        break;
      default:
        break;
      }
    }
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
      for (j = 1; j <= repMembrSize; j++) {
        localSplitIndicator[j] = NEITHER;
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
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        stackAndGetSplitSurv(treeID,
                             repMembrIndx,
                             repMembrSize,
                             nonMissMembrIndx,
                             nonMissMembrSize,
                             & localEventTimeCount,
                             & localEventTimeIndex,
                             & localEventTimeSize,
                             & nodeParentEvent,
                             & nodeParentAtRisk,
                             & nodeLeftEvent,
                             & nodeLeftAtRisk,
                             & nodeRightEvent,
                             & nodeRightAtRisk);
        if (localEventTimeSize > 0) {
          nodeParentEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
          nodeLeftEventCR = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          if (localEventTimeSize > 0) {
            nodeParentInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
            nodeLeftInclusiveAtRisk = uimatrix(1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
      }
      if (localEventTimeSize > 0) {
        for (m=1; m <= localEventTimeSize; m++) {
          for (q = 1; q <= RF_eventTypeSize; q++) {
            nodeParentEventCR[q][m] = 0;
          }
        }
        for (k = 1; k <= nonMissMembrSize; k++) {
          if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
            tIndx = 0;  
            for (m = 1; m <= localEventTimeSize; m++) {
              if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                tIndx = m;
              }
              else {
                m = localEventTimeSize;
              }
            }
            nodeParentEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
          }
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          for (m = 1; m <= localEventTimeSize; m++) {
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
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
        if (factorFlag == FALSE) {
          for (j = 1; j <= nonMissMembrSize; j++) {
            localSplitIndicator[ nonMissMembrIndx[indxx[j]] ] = RIGHT;
          }
          for (m = 1; m <= localEventTimeSize; m++) {
            nodeLeftAtRisk[m] = 0;
            for (q = 1; q <= RF_eventTypeSize; q++) {
              nodeLeftEventCR[q][m] = 0;
            }
          }
          leftSizeIter = 0;
          rghtSizeIter = nonMissMembrSize;
        }
        for (j = 1; j < splitLength; j++) {
          if (factorFlag == TRUE) {
            priorMembrIter = 0;
            leftSize = 0;
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
            if (factorFlag == TRUE) {
              for (m = 1; m <= localEventTimeSize; m++) {
                nodeLeftAtRisk[m] = 0;
                for (q = 1; q <= RF_eventTypeSize; q++) {
                  nodeLeftEventCR[q][m] = 0;
                }
              }
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[  nonMissMembrIndx[indxx[k]]  ] == LEFT) {
                  tIndx = 0;  
                  for (m = 1; m <= localEventTimeSize; m++) {
                    if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                      tIndx = m;
                      nodeLeftAtRisk[tIndx] ++;
                    }
                    else {
                      m = localEventTimeSize;
                    }
                  }
                  if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                    nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
                  }
                }
              }
            }
            else {
              for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                tIndx = 0;  
                for (m = 1; m <= localEventTimeSize; m++) {
                  if (localEventTimeIndex[m] <= RF_masterTimeIndex[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]) {
                    tIndx = m;
                    nodeLeftAtRisk[tIndx] ++;
                  }
                  else {
                    m = localEventTimeSize;
                  }
                }
                if (RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ] > 0) {
                  nodeLeftEventCR[RF_eventTypeIndex[(uint) RF_status[treeID][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ]]][tIndx] ++;
                }
              }
              rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
              leftSizeIter = currentMembrIter - 1;
            }
            switch(RF_splitRule) {
            case SURV_CR_LAU:
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
              break;
            case SURV_CR_LOG:
              break;
            default:
              break;
            }
            delta = deltaNum = deltaDen =  0.0;
            switch(RF_splitRule) {
            case SURV_CR_LAU:
              for (q = 1; q <= RF_eventTypeSize; q++) {
                if (RF_crWeight[q] > 0) {
                  deltaSubNum = 0;
                  for (m = 1; m <= localEventTimeSize; m++) {
                    deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftInclusiveAtRisk[q][m] / nodeParentInclusiveAtRisk[q][m])));
                  }
                  deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                  deltaSubDen = 0;
                  for (m = 1; m <= localEventTimeSize; m++) {
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
              break;
            case SURV_CR_LOG:
              for (q = 1; q <= RF_eventTypeSize; q++) {
                if (RF_crWeight[q] > 0) {
                  deltaSubNum = 0;
                  for (m=1; m <= localEventTimeSize; m++) {
                    deltaSubNum = deltaSubNum + (nodeLeftEventCR[q][m] - (nodeParentEventCR[q][m] * ((double) nodeLeftAtRisk[m] / nodeParentAtRisk[m])));
                  }
                  deltaNum = deltaNum + (RF_crWeight[q] * deltaSubNum);
                  deltaSubDen = 0;
                  for (m = 1; m <= localEventTimeSize; m++) {
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
              break;
            default:
              break;
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
      }  
      else {
      }
      unstackSplitVector(treeID,
                         permissibleSplitSize,
                         splitLength,
                         factorFlag,
                         deterministicSplitFlag,
                         mwcpSizeAbsolute,
                         permissibleSplitPtr);
      unselectRandomCovariates(treeID,
                               parent,
                               repMembrSize,
                               indxx,
                               nonMissMembrSizeStatic,
                               nonMissMembrIndx,
                               multImpFlag);
      if (!((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP)))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
        if (localEventTimeSize > 0) {
          free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
          free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
        }
        switch(RF_splitRule) {
        case SURV_CR_LAU:
          if (localEventTimeSize > 0) {
            free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
            free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
          }
          break;
        case SURV_CR_LOG:
          break;
        default:
          break;
        }
      }
    }  
    if ((RF_mRecordSize == 0) || (multImpFlag) || (!(RF_optHigh & OPT_MISS_SKIP))) {
        unstackSplitSurv(localEventTimeCount,
                         localEventTimeIndex,
                         localEventTimeSize,
                         nodeParentEvent,
                         nodeParentAtRisk,
                         nodeLeftEvent,
                         nodeLeftAtRisk,
                         nodeRightEvent,
                         nodeRightAtRisk);
      free_uimatrix(nodeParentEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      free_uimatrix(nodeLeftEventCR, 1, RF_eventTypeSize, 1, localEventTimeSize);
      switch(RF_splitRule) {
      case SURV_CR_LAU:
        free_uimatrix(nodeParentInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
        free_uimatrix(nodeLeftInclusiveAtRisk, 1, RF_eventTypeSize, 1, localEventTimeSize);
        break;
      case SURV_CR_LOG:
        break;
      default:
        break;
      }
    }
    unstackSplitIndicator(repMembrSize, localSplitIndicator);
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
  }  
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                splitStatistic,
                                deltaMax);
  return result;
}
