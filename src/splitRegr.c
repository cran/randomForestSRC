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
#include     "splitRegr.h"
char regressionXwghtSplit (uint    treeID,
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
  uint j, k;
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
    double sumLeft, sumRght, sumRghtSave, sumLeftSqr, sumRghtSqr, sumRghtSqrSave;
    double leftTemp, rghtTemp, leftTempSqr, rghtTempSqr;
    sumLeft = sumRght = sumLeftSqr = sumRghtSqr = 0;  
    leftTempSqr = rghtTempSqr = 0;  
    delta = 0;  
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
      sumRghtSave = sumRghtSqrSave = 0.0;
      switch(RF_splitRule) {
      case REGR_WT_NRM:
        sumRghtSave = 0.0;
        for (j = 1; j <= nonMissMembrSize; j++) {
          sumRghtSave += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ];
        }
        break;
      default:
        sumRghtSave = sumRghtSqrSave = 0.0;
        for (j = 1; j <= nonMissMembrSize; j++) {
          sumRghtSave += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ];
          sumRghtSqrSave += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[j]]] ], 2.0);
        }
        break;
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
        sumRght      = sumRghtSave;
        sumRghtSqr   = sumRghtSqrSave;
        sumLeft      = 0.0;
        sumLeftSqr   = 0.0;
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
        if ((leftSize  >= (RF_minimumNodeSize)) && (rghtSize  >= (RF_minimumNodeSize))) {
          if (factorFlag == TRUE) {
            switch(RF_splitRule) {
            case REGR_WT_NRM:
              sumLeft = sumRght = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                  sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                }
                else {
                  sumRght += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                }
              } 
              break;
            default:
              sumLeft = sumRght = 0.0;
              sumLeftSqr = sumRghtSqr = 0.0;
              for (k = 1; k <= nonMissMembrSize; k++) {
                if (localSplitIndicator[ nonMissMembrIndx[indxx[k]] ] == LEFT) {
                  sumLeft    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                  sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                }
                else {
                  sumRght    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                  sumRghtSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                }
              } 
              break;
            }
          }
          else {
            switch(RF_splitRule) {
            case REGR_WT_NRM:
              for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                sumLeft += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumRght -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
              }
              break;
            default:
              for (k = leftSizeIter + 1; k < currentMembrIter; k++) {
                sumLeft    += RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumLeftSqr += pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
                sumRght    -= RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ];
                sumRghtSqr -= pow(RF_response[treeID][1][ repMembrIndx[nonMissMembrIndx[indxx[k]]] ], 2.0);
              }
              break;
            }
            rghtSizeIter = rghtSizeIter - (currentMembrIter - (leftSizeIter + 1));
            leftSizeIter = currentMembrIter - 1;
          }
          leftTemp = pow(sumLeft, 2.0) / pow(leftSize, 2.0);
          rghtTemp = pow(sumRght, 2.0) / pow(rghtSize, 2.0);
          switch(RF_splitRule) {
          case REGR_WT_NRM:
            sumLeftSqr = pow(sumLeft, 2.0) / leftSize;
            sumRghtSqr = pow(sumRght, 2.0) / rghtSize;
            delta = sumLeftSqr + sumRghtSqr;
            break;
          case REGR_WT_OFF:
            leftTempSqr = sumLeftSqr / leftSize;
            rghtTempSqr = sumRghtSqr / rghtSize;
            delta = leftTemp + rghtTemp - leftTempSqr - rghtTempSqr;
            break;
          case REGR_WT_HVY:
            leftTempSqr = sumLeftSqr * leftSize / pow (repMembrSize, 2.0);
            rghtTempSqr = sumRghtSqr * rghtSize / pow (repMembrSize, 2.0);
            delta = leftTemp + rghtTemp - leftTempSqr - rghtTempSqr;
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
    }  
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
  result = summarizeSplitResult(*splitParameterMax,
                                *splitValueMaxCont,
                                *splitValueMaxFactSize,
                                *splitValueMaxFactPtr,
                                 splitStatistic,
                                 deltaMax);
  return result;
}
