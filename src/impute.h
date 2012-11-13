////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.1
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


#ifndef RSFIMPUTE_H
#define RSFIMPUTE_H
#include "node.h"
void imputeInteraction (uint treeID, Node *parent, uint *repMembrIndx, uint repMembrSize);
char imputeNode (uint     type,
                 char     chainFlag,
                 uint     treeID, 
                 Node    *parent,
                 uint    *repAbsIdx,
                 uint     repNodeSize,
                 uint    *iAbsIdx,
                 uint     iNodeSize);
char restoreNodeMembership(uint  mode, 
                           char  rootFlag,
                           uint  treeID, 
                           Node *parent, 
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *allMembrIndx,
                           uint  allMembrSize,
                           uint *ngAllMembrIndx,
                           uint  ngAllMembrSize,
                           uint *bootMembrIndxIter);
void imputeUpdateShadow (uint      mode, 
                         char      selectionFlag,
                         double  **shadowResponse, 
                         double  **shadowPredictor);
void imputeUpdateSummary (uint     mode, 
                          double **responsePtr, 
                          double **predictorPtr, 
                          uint     treeID);
void imputeUpdateSummaryNew (uint     mode, 
                          double **responsePtr, 
                          double **predictorPtr, 
                             uint     treeID);
void imputeSummary(uint      mode,
                   char      selectionFlag);
void imputeResponse(uint      mode,
                       uint      treeID,
                       double  **tempResponse);
void imputeCommon(uint      mode,
                  uint      treeID,
                  char      selectionFlag,
                  char      predictorFlag);
void imputeMultipleTime (char selectionFlag);
double getMaximalValue(double *value, uint size, uint treeID);
double getMedianValue(double *value, uint size);
double getMeanValue(double *value, uint size);
double getSampleValue(double *value, uint size, char chainFlag, uint treeID);
uint getRecordMap(uint     *map, 
                  uint      size, 
                  double  **resp, 
                  double  **data);
void updateTimeIndexArray(uint    treeID,
                          Node   *parent, 
                          uint   *allMemberIndx,
                          uint    allMembrSize,
                          double *time, 
                          char    naflag,
                          char    idFlag,
                          uint   *masterTimeIndex);
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mvSign,
                            uint   *mRecordIndex,
                            uint   *meIndividualSize,
                            uint  **eIndividual);
void stackShadow (uint mode, uint treeID);
void unstackShadow (uint mode, uint treeID);
#endif
