////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.4
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


#ifndef RSFTREEUTIL_H
#define RSFTREEUTIL_H
#include "node.h"
Node *getTerminalNode(uint treeID, uint leaf);
void getRawNodeSize(uint  type,
                    uint  treeID, 
                    Node *parent, 
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize);
char forkAndUpdate(uint    treeID, 
                   Node   *parent, 
                   uint   *repMembrIndx,
                   uint    repMembrSize,
                   uint   *allMembrIndx,
                   uint    allMembrSize,
                   uint    splitParameterMax,
                   double  splitValueMaxCont,
                   uint    splitValueMaxFactSize,
                   uint   *splitValueMaxFactPtr,
                   double  splitStatistic,
                   char   *localSplitIndicator,
                   char   *localOmitMembrFlag,
                   char    multImpFlag,
                   char   *membershipIndicator,
                   uint   *leftDaughterSize,
                   uint   *rightDaughterSize);
char growTree(char     rootFlag,
              char     multImpFlag,
              uint     b,
              Node    *parent,
              uint    *repMembrIndx,
              uint     repMembrSize,
              uint    *allMembrIndx,
              uint     allMembrSize,
              uint     depth,
              uint    *maximumDepth,
              uint    *bootMembrIndxIter);
char restoreTree(uint    b,
                 Node   *parent,
                 uint   *offset,
                 uint   *treeID,
                 uint   *nodeID,
                 uint   *parmID,
                 double *contPT,
                 uint   *mwcpSZ,
                 uint  **mwcpPtr,
                 uint    depth,
                 uint   *maximumDepth);
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *contPT,
              uint   *mwcpSZ,
              uint  **mwcpPtr);
void freeTree(uint treeID, Node *parent, char rootFlag);
void getSplitDepth(Node *parent, uint *maximumDepth);
void freeSplitDepth(uint treeID);
void saveStatistics(char    mode,
                    uint    b,
                    Node   *parent,
                    uint   *offset,
                    double *spltST,
                    double *spltVR);
uint getMaximumDepth(Node *parent);
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount);
void getTreeInfo(uint treeID, Node *parent);
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset);
#endif
