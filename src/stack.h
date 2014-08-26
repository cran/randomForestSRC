////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.5
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


#ifndef RSFSTACK_H
#define RSFSTACK_H
#include "node.h"
void initializeTimeArrays(char mode);
void stackFactorArrays();
void stackFactorGeneric(uint    size, 
                        char  **type, 
                        uint  **p_factorMap,
                        uint   *factorCount,
                        uint  **p_factorIndex,
                        uint  **p_factorSize);
void unstackFactorArrays();
void initializeFactorArrays(char mode);
char stackMissingArrays(char mode);
void unstackMissingArrays(char mode);
void stackMissingSignatures(uint     obsSize, 
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize, 
                            uint   **p_recordIndex, 
                            uint    *p_vSize,
                            int   ***p_vSign, 
                            int    **p_vIndex,
                            uint    *pRF_mrFactorSize,
                            uint   **pRF_mrFactorIndex,
                            uint    *pRF_mxFactorSize,
                            uint   **pRF_mxFactorIndex,
                            char    *pRF_mTimeFlag,
                            char    *pRF_mStatusFlag,
                            char    *pRF_mResponseFlag,
                            char    *pRF_mPredictorFlag);
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize, 
                              uint     *recordIndex, 
                              uint      vSize,
                              int     **vSign, 
                              int      *vIndex,
                              uint      mrFactorSize,
                              uint     *mrFactorIndex,
                              uint      mxFactorSize,
                              uint     *mxFactorIndex);
char stackCompetingArrays(char mode);
void unstackCompetingArrays(char mode);
char stackClassificationArrays(char mode);
void unstackClassificationArrays(char mode);
void getEventTypeSize(uint     obsSize, 
                      double  *status, 
                      uint    *mRecordMap, 
                      int    **mpSign,  
                      char     overWriteFlag,
                      uint    *eventTypeSize,
                      uint    *msize,
                      uint    *eventType);
void getClassLevelSize(uint      obsSize, 
                       double  **response, 
                       uint     *mRecordMap, 
                       int     **mpSign,  
                       uint     *classLevelSize,
                       uint    **classLevel);
#endif
