////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.0
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
////    email:  kogalurshear@gmail.com
////    URL:    http://www.kogalur-shear.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************


#ifndef RSFSPLITUTIL_H
#define RSFSPLITUTIL_H
#include "node.h"
void updateMaximumSplit(double  delta, 
                        uint    randomCovariate,
                        uint    jLong,
                        char    factorFlag,
                        uint    mwcpSizeAbsolute,
                        double *deltaMax,
                        uint   *splitParameterMax,
                        double *splitValueMaxCont,
                        uint   *splitValueMaxFactSize,
                        uint  **splitValueMaxFactPtr,
                        void   *permissibleSplitPtr);
uint stackAndSelectRandomCovariates(uint      treeID,
                                    Node     *parent,
                                    uint     *repMembrIndx,
                                    uint      repMembrSize,
                                    uint    **covariateIndex,
                                    double ***permissibleSplit,
                                    uint    **permissibleSplitSize);
void unstackRandomCovariates(uint     treeID,
                             uint     nodeSize, 
                             uint    *covariateIndex,
                             double **permissibleSplit,
                             uint    *permissibleSplitSize);
uint getSelectableElement(uint    treeID,
                          uint    length,
                          char   *permissible,
                          double *weight);
void stackSplitTime(uint **localEventTimeCount, 
                    uint **localEventTimeIndex);
void unstackSplitTime(uint *localEventTimeCount, 
                      uint *localEventTimeIndex);
uint getEventTimeSize(uint   treeID, 
                      Node   *parent, 
                      uint   *repMembrIndx,
                      uint    repMembrSize,
                      uint   *localEventTimeCount, 
                      uint   *localEventTimeIndex);
void stackSplitCompactEventAndRisk(uint   eventTimeSize,
                                   uint **nodeParentEvent,
                                   uint **nodeParentAtRisk,
                                   uint **nodeLeftEvent,
                                   uint **nodeLeftAtRisk,
                                   uint **nodeRightEvent,
                                   uint **nodeRightAtRisk);
void unstackSplitCompactEventAndRisk(uint  eventTimeSize,
                                     uint *nodeParentEvent,
                                     uint *nodeParentAtRisk,
                                     uint *nodeLeftEvent,
                                     uint *nodeLeftAtRisk,
                                     uint *nodeRightEvent,
                                     uint *nodeRightAtRisk);
void stackSplitIndicator(uint   nodeSize,
                         char **localSplitIndicator);
void unstackSplitIndicator(uint  nodeSize,
                           char *localSplitIndicator);
void getEventAndRisk(uint    treeID,
                     uint   *repMembrIndx,
                     uint    repMembrSize,
                     uint   *localEventTimeCount,
                     uint   *localEventTimeIndex,
                     uint    localEventTimeSize,
                     uint   *nodeParentEvent,
                     uint   *nodeParentAtRisk);
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
                        char *localSplitIndicator);
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
char summarizeSplitResult(uint   splitParameterMax, 
                          double splitValueMaxCont,
                          uint   splitValueMaxFactSize,
                          uint  *splitValueMaxFactPtr,
                          double deltaMax);
char getStandardDeviation(uint repSize, uint *repIndx, double *target);
#endif
