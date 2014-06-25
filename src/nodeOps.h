////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.3
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


#ifndef NODEOPS_H
#define NODEOPS_H
#include "node.h"
#include "terminal.h"
Terminal *makeTerminal();
void freeTerminal(Terminal *parent);
Node *makeNode(unsigned int xSize);
void freeNode(Node *parent);
void freeTerminalNodeLocalSurvivalStructures(Node *terminalNode);
void freeTerminalNodeSurvivalStructures(Node *terminalNode);
void getNodeInfo(Node *leaf);
void setParent(
  Node *daughter,
  Node *parent
);
void setLeftDaughter(
   Node *daughter,
   Node *parent
);
void setRightDaughter(
  Node *daughter,
  Node *parent
);
char forkNode(Node         *parent,
              unsigned int  splitParameter,
              double        splitValueMaxCont,
              unsigned int  splitValueMaxFactSize,
              unsigned int *splitValueMaxFactPtr);
void stackAtRiskAndEventCounts(Node *tNode, unsigned int eTypeSize, unsigned int mTimeSize);
void stackEventTimeIndex(Node *tNode, unsigned int mTimeSize);
void unstackAtRiskAndEventCounts(Node *tNode);
void unstackEventTimeIndex(Node *tNode);
void unstackAtRisk(Node *tNode);
void stackLocalRatio(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatio(Node *tNode);
void stackLocalSurvival(Node *tNode, unsigned int eTimeSize);
void unstackLocalSurvival(Node *tNode);
void stackLocalNelsonAalen(Node *tNode, unsigned int eTimeSize);
void unstackLocalNelsonAalen(Node *tNode);
void stackLocalCSH(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCSH(Node *tNode);
void stackLocalCIF(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCIF(Node *tNode);
void stackNelsonAalen(Node *tNode, unsigned int sTimeSize);
void unstackNelsonAalen(Node *tNode);
void stackSurvival(Node *tNode, unsigned int sTimeSize);
void unstackSurvival(Node *tNode);
void stackCSH(Node *tNode, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCSH(Node *tNode);
void stackCIF(Node *tNode, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCIF(Node *tNode);
void stackMortality(Node *tNode, unsigned int eTypeSize);
void unstackMortality(Node *tNode);
void stackMPSign(Node *node, unsigned int mpIndexSize);
void unstackMPSign(Node *node);
void stackFMPSign(Node *node, unsigned int fmpIndexSize);
void unstackFMPSign(Node *node);
void stackNodeLMPIndex(Node *node, unsigned int size);
void unstackNodeLMPIndex(Node *node);
void stackNodeFLMPIndex(Node *node, unsigned int size);
void unstackNodeFLMPIndex(Node *node);
void stackNodeLMRIndex(Node *node, unsigned int size);
void unstackNodeLMRIndex(Node *node);
void stackNodeFLMRIndex(Node *node, unsigned int size);
void unstackNodeFLMRIndex(Node *node);
void stackTermLMIIndex(Terminal *tNode, unsigned int size);
void unstackTermLMIIndex(Terminal *tNode);
void stackMultiClassProb(Node *tNode, unsigned int rfCount, unsigned int *rfSize);
void unstackMultiClassProb(Node *tNode);
void stackSplitDepth(Node *tNode, unsigned int depth);
void unstackSplitDepth(Node *tNode);
#endif
