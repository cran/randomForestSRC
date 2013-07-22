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


#ifndef RSFUTIL_H
#define RSFUTIL_H
#include "node.h"
void updateEnsembleCalculations (char      multipleImputeFlag,
                                 uint      mode,
                                 uint      b);
void copyDenominator(uint mode, uint *denominatorCopy);
void copyEnsemble(uint mode, double **ensembleCopy);
void getVariablesUsed(uint treeID, Node *rootPtr, uint *varUsedVector);
char stackAndImputePerfResponse(uint      mode, 
                                char      multipleImputeFlag, 
                                uint      treeID,
                                uint      serialID, 
                                uint      rSize, 
                                double ***responsePtr);
double **stackAndImputeGenericResponse(char flag, uint mode, uint rSize, uint obsSize, uint treeID, uint serialID, double **responsePtr);
void unstackImputeResponse(char flag, uint rSize, uint obsSize, double **mResponsePtr);
void getPerformance(uint      b,
                    uint      mode,
                    uint      obsSize,
                    double  **responsePtr,
                    double   *outcome,
                    double  **conditionalOutcome,
                    uint     *denominatorCopy,
                    double   *performancePtr);
double *stackCondPerformance();
void    unstackCondPerformance(double *cpv);
void finalizeEnsembleEstimates(uint mode, uint rejectedTreeCount);
#endif
