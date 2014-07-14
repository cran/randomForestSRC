////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.4
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


#ifndef RSFIMPORTANCE_H
#define RSFIMPORTANCE_H
#include "node.h"
Node *identifyPerturbedMembership(Node    *parent,
                                  double **shadowVIMP,
                                  uint     index);
Node *randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter,
                          uint     treeID);
void permute(uint ranGenID, uint p, uint n, uint *indx);
void getRandomMembership(uint      mode,
                         uint      treeID,
                         Node    **vimpMembership,
                         uint      p);
void getPermuteMembership(uint      mode,
                          uint      treeID,
                          Node    **vimpMembership,
                          uint      p);
void getVimpMembership(uint      mode, 
                       uint      treeID,
                       Node    **vimpMembership,
                       uint      p);
void updateGenericVimpEnsemble (uint       mode,
                                uint       treeID,
                                uint       targetIndex,
                                Node     **noiseMembership,
                                char       ensembleFlag,
                                double   **outcome,
                                double  ***sOutcome,
                                double  ***mcOutcome);
void updateTreeEnsemble (uint       mode,
                         uint       treeID,
                         double   **outcome,
                         double  ***sOutcome,
                         double  ***mcOutcome);
void updateVimpEnsemble (uint       mode,
                         uint       treeID,
                         Node     **vimpMembership,
                         uint       p);
void summarizeVimpPerformance(uint       mode,
                              uint       treeID,
                              uint       p);
void finalizeVimpPerformance(uint mode, uint rejectedTreeCount);
void  stackVimpMembership(uint mode, Node ***membership);
void  unstackVimpMembership(uint mode, Node **membership);
void stackTreeEnsemble(uint         mode,
                       uint         treeID,
                       uint       **denomTree,
                       double    ***treeOutcome,
                       double   ****sTreeOutcome,
                       double   ****mcTreeOutcome);
void unstackTreeEnsemble(uint       mode,
                         uint       treeID,
                         uint      *denomTree,
                         double    **treeOutcome,
                         double   ***sTreeOutcome,
                         double   ***mcTreeOutcome);
void updateVimpCalculations (uint mode, uint b, uint intrIndex, Node **vimpMembership);
void summarizeTreePerformance(uint mode, uint treeID);
uint getEnsembleDim ();
#endif
