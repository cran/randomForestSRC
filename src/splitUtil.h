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


#ifndef RSFSPLITUTIL_H
#define RSFSPLITUTIL_H
#include "node.h"
void stackSplitIndicator(uint   nodeSize,
                         char **localSplitIndicator);
void unstackSplitIndicator(uint  nodeSize,
                           char *localSplitIndicator);
void stackSplitEventTime(uint **localEventTimeCount,
                         uint **localEventTimeIndex);
void unstackSplitEventTime(uint *localEventTimeCount,
                           uint *localEventTimeIndex);
uint getSplitEventTime(uint   treeID,
                       uint   *repMembrIndx,
                       uint    repMembrSize,
                       uint   *nonMissMembrIndx,
                       uint    nonMissMembrSize,
                       uint   *localEventTimeCount,
                       uint   *localEventTimeIndex);
void stackSplitEventAndRisk(uint   eventTimeSize,
                            uint **nodeParentEvent,
                            uint **nodeParentAtRisk,
                            uint **nodeLeftEvent,
                            uint **nodeLeftAtRisk,
                            uint **nodeRightEvent,
                            uint **nodeRightAtRisk);
void unstackSplitEventAndRisk(uint  eventTimeSize,
                              uint *nodeParentEvent,
                              uint *nodeParentAtRisk,
                              uint *nodeLeftEvent,
                              uint *nodeLeftAtRisk,
                              uint *nodeRightEvent,
                              uint *nodeRightAtRisk);
void getSplitEventAndRisk(uint    treeID,
                          uint   *repMembrIndx,
                          uint    repMembrSize,
                          uint   *nonMissMembrIndx,
                          uint    nonMissMembrSize,
                          uint   *localEventTimeCount,
                          uint   *localEventTimeIndex,
                          uint    localEventTimeSize,
                          uint   *nodeParentEvent,
                          uint   *nodeParentAtRisk);
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
                          uint  **nodeRightAtRisk);
void unstackSplitSurv(uint *localEventTimeCount,
                      uint *localEventTimeIndex,
                      uint  eventTimeSize,
                      uint *nodeParentEvent,
                      uint *nodeParentAtRisk,
                      uint *nodeLeftEvent,
                      uint *nodeLeftAtRisk,
                      uint *nodeRightEvent,
                      uint *nodeRightAtRisk);
uint stackAndConstructSplitVector(uint     treeID,
                                  uint     localMembershipSize,
                                  uint     randomCovariateIndex,
                                  double  *permissibleSplit,
                                  uint     permissibleSplitSize,
                                  char    *factorFlag,
                                  char    *deterministicSplitFlag,
                                  uint    *mwcpSizeAbsolute,
                                  void   **permissibleSplitPtr);
void unstackSplitVector(uint   treeID,
                        uint   permissibleSplitSize,
                        uint   splitLength,
                        char   factorFlag,
                        char   deterministicSplitFlag,
                        uint   mwcpSizeAbsolute,
                        void  *permissibleSplitPtr);
void stackRandomCovariates(uint      treeID,
                           Node     *parent,
                           uint      repMembrSize,
                           char      multImpFlag,
                           char    **covariateFlag,
                           uint    **covariateIndex,
                           uint     *covariateSize,
                           uint     *uniformCovariateIndex,
                           double  **cdf,
                           uint     *cdfSize,
                           uint    **cdfSort,
                           uint    **density,
                           uint     *densitySize,
                           uint   ***densitySwap);
void unstackRandomCovariates(uint     treeID,
                             char    *covariateFlag,
                             uint    *covariateIndex,
                             uint     covariateSize,
                             double  *cdf,
                             uint     cdfSize,
                             uint    *cdfSort,
                             uint    *density,
                             uint     densitySize,
                             uint   **densitySwap,
                             uint     repMembrSize,
                             uint    *nonMissMembrIndx,
                             double  *permissibleSplit);
void initializeCDF(uint     treeID,
                   Node    *parent,
                   char    *covariateFlag,
                   uint    *covariateIndex,
                   uint    *uniformCovariateSize,
                   uint    *uniformCovariateIndex,
                   double **cdf,
                   uint    *cdfSize,
                   uint   **cdfSort,
                   uint   **density,
                   uint    *densitySize,
                   uint  ***densitySwap);
void updateCDF(uint    treeID,
               char   *covariateFlag,
               uint   *covariateIndex,
               uint   *covariateSize,
               uint    uniformCovariateIndex,
               double *cdf,
               uint   *cdfSize,
               uint   *cdfSort,
               uint   *density,
               uint   *densitySize,
               uint  **densitySwap,
               uint    index);
uint getSelectableElementNew (uint    treeID,
                              char   *covariateFlag,
                              uint   *covariateIndex,
                              uint    uniformCovariateSize,
                              uint   *uniformCovariateIndex,
                              double *cdf,
                              uint    cdfSize,
                              uint   *cdfSort,
                              uint   *density,
                              uint    densitySize);
uint sampleWithoutReplacement (uint    treeID,
                               uint   *index,
                               uint    size,
                               uint   *sampleIndex);
uint getSelectableElement (uint    treeID,
                           uint    length,
                           char   *permissible);
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
                            char      multImpFlag);
void unselectRandomCovariates(uint      treeID,
                              Node     *parent,
                              uint      repMembrSize,
                              uint     *indxx,
                              uint      nonMissMembrSizeStatic,
                              uint     *nonMissMembrIndx,
                              char      multImpFlag);
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
                           uint *currentMembrIter);
void getReweightedRandomPair(uint    treeID,
                             uint    relativefactorSize,
                             uint    absoluteFactorSize,
                             double *absoluteLevel,
                             uint   *result);
void getRandomPair(uint treeID, uint relativeFactorSize, uint absoluteFactorSize, double *absoluteLevel, uint *result);
void createRandomBinaryPair(uint    treeID,
                            uint    relativeFactorSize,
                            uint    absoluteFactorSize,
                            uint    groupSize,
                            double *absolutelevel,
                            uint   *pair);
void convertRelToAbsBinaryPair(uint    treeID,
                               uint    relativeFactorSize,
                               uint    absoluteFactorSize,
                               uint    relativePair,
                               double *absoluteLevel,
                               uint   *pair);
char summarizeSplitResult(uint    splitParameterMax,
                          double  splitValueMaxCont,
                          uint    splitValueMaxFactSize,
                          uint   *splitValueMaxFactPtr,
                          double *splitStatistic,
                          double  deltaMax);
char getPreSplitResult (uint      treeID,
                        Node     *parent,
                        uint      repMembrSize,
                        uint     *repMembrIndx,
                        uint     *nonMissMembrSize,
                        uint    **nonMissMembrIndx,
                        double  **permissibleSplit,
                        char      multImpFlag);
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
                        char  **splitIndicator);
#endif
