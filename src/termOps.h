////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.6.0
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


#ifndef TERMOPS_H
#define TERMOPS_H
#include "node.h"
#include "terminal.h"
Terminal *makeTerminal();
void freeTerminal(Terminal *parent);
void stackTermLMIIndex(Terminal *tTerm, unsigned int size);
void unstackTermLMIIndex(Terminal *tTerm);
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresNonVimp(Terminal *tTerm);
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm);
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm);
void stackAtRiskAndEventCounts(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize);
void stackEventTimeIndex(Terminal *tTerm, unsigned int mTimeSize);
void unstackAtRiskAndEventCounts(Terminal *tTerm);
void unstackEventTimeIndex(Terminal *tTerm);
void unstackAtRisk(Terminal *tTerm);
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalRatio(Terminal *tTerm);
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalSurvival(Terminal *tTerm);
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize);
void unstackLocalNelsonAalen(Terminal *tTerm);
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCSH(Terminal *tTerm);
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize);
void unstackLocalCIF(Terminal *tTerm);
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize);
void unstackNelsonAalen(Terminal *tTerm);
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize);
void unstackSurvival(Terminal *tTerm);
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCSH(Terminal *tTerm);
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize);
void unstackCIF(Terminal *tTerm);
void stackMortality(Terminal *tTerm, unsigned int eTypeSize);
void unstackMortality(Terminal *tTerm);
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize);
void unstackMultiClassProb(Terminal *tTerm);
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount);
void unstackMeanResponse(Terminal *tTerm);
void getTerminalInfo(Terminal *termPtr);
#endif
