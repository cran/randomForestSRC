////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.2
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


#include        <stdlib.h>
#include        "nrutil.h"
#include       "nodeOps.h"
#include   <R_ext/Print.h>
#include   <Rdefines.h>
extern unsigned int getTraceFlag();
extern unsigned int getNodeDefTraceFlag();
extern unsigned int getForkDefTraceFlag();
extern unsigned int getTurnOffTraceFlag();
extern unsigned int getTurnOnTraceFlag();
#include <R_ext/Arith.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
Terminal *makeTerminal() {
  Terminal *parent = (Terminal*) gblock((size_t) sizeof(Terminal));
  parent -> lmiIndex      = NULL;
  parent -> lmiSize       = 0;
  parent -> lmiValue  = NULL;
  parent -> nodeID     = 0;
  return parent;
}
void freeTerminal(Terminal        *parent) {
  unstackTermLMIIndex(parent);
  free_gblock(parent, sizeof(Terminal));
}
Node *makeNode(unsigned int xSize) {
  unsigned int i;
  Node *parent = (Node*) gblock((size_t) sizeof(Node));
  parent -> xSize = xSize;
  parent -> permissibleSplit = cvector(1, xSize);
  for (i = 1; i <= xSize; i++) {
    (parent -> permissibleSplit)[i] = TRUE;
  }
  parent -> mate               = NULL;
  parent -> left               = NULL;
  parent -> right              = NULL;
  parent -> splitFlag            = TRUE;
  parent -> predictedOutcome     = NA_REAL;
  parent -> splitParameter       = 0;
  parent -> splitValueCont       = NA_REAL;
  parent -> splitValueFactSize   = 0;
  parent -> splitValueFactPtr    = NULL;
  parent -> splitStatistic       = NA_REAL;
  parent -> variance             = NA_REAL;
  parent -> nodeID               = 0;
  parent -> orderedNodeID        = 0;
  parent -> depth                = 0;
  parent -> splitDepth           = NULL;
  parent -> pseudoTerminal       = FALSE;
  parent -> eTypeSize            = 0;
  parent -> mTimeSize            = 0;
  parent -> eTimeSize            = 0;
  parent -> sTimeSize            = 0;
  parent -> atRiskCount          = NULL;
  parent -> eventCount           = NULL;
  parent -> eventTimeIndex       = NULL;
  parent -> localRatio           = NULL;
  parent -> localCSH             = NULL;
  parent -> localCIF             = NULL;
  parent -> localSurvival        = NULL;
  parent -> localNelsonAalen     = NULL;
  parent -> CSH                  = NULL;
  parent -> CIF                  = NULL;
  parent -> survival             = NULL;
  parent -> nelsonAalen          = NULL;
  parent -> rfCount              = 0;
  parent -> rfSize               = NULL;
  parent -> multiClassProb       = NULL;
  parent -> weight               = 0.0;
  parent -> membrCount           = 0;
  parent -> mpIndexSize          = 0;
  parent -> fmpIndexSize         = 0;
  parent -> mpSign               = NULL;
  parent -> fmpSign              = NULL;
  parent -> imputed              = FALSE;
  parent -> lmpIndex             = NULL;
  parent -> flmpIndex            = NULL;
  parent -> lmpIndexAllocSize    = 0;
  parent -> flmpIndexAllocSize   = 0;
  parent -> lmpIndexActualSize   = 0;
  parent -> flmpIndexActualSize  = 0;
  parent -> lmrIndex             = NULL;
  parent -> flmrIndex            = NULL;
  parent -> lmrIndexAllocSize    = 0;
  parent -> flmrIndexAllocSize   = 0;
  parent -> lmrIndexActualSize   = 0;
  parent -> flmrIndexActualSize  = 0;
  return parent;
}
void freeNode(Node         *parent) {
  if (parent -> permissibleSplit != NULL) {
    free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
    parent -> permissibleSplit = NULL;
  }
  if ((parent -> splitValueFactSize) > 0) {
    if (parent -> splitValueFactPtr != NULL) {
      free_uivector(parent -> splitValueFactPtr, 1, parent -> splitValueFactSize);
      parent -> splitValueFactPtr = NULL;
    }
  }
  if ((parent -> splitParameter) == 0) {
    if ((parent -> depth) > 0) {
      if (parent -> splitDepth != NULL) {
        free_uivector(parent -> splitDepth, 1, parent -> depth);
        parent -> splitDepth = NULL;
      }
    }
  }
  unstackMPSign(parent);
  unstackFMPSign(parent);
  unstackNodeLMPIndex(parent);
  unstackNodeFLMPIndex(parent);
  unstackNodeLMRIndex(parent);
  unstackNodeFLMRIndex(parent);
  if ((parent -> splitParameter) == 0) {
    freeTerminalNodeSurvivalStructures(parent);
  }
  free_gblock(parent, sizeof(Node));
}
void freeTerminalNodeLocalSurvivalStructures(Node *terminalNode) {
  unstackLocalRatio(terminalNode);
  unstackLocalSurvival(terminalNode);
  unstackLocalNelsonAalen(terminalNode);
  if (terminalNode -> eTypeSize > 1) {
    unstackLocalCSH(terminalNode);
    unstackLocalCIF(terminalNode);
  }
}
void freeTerminalNodeSurvivalStructures(Node *terminalNode) {
  if (terminalNode -> eTypeSize > 1) {
  }
  unstackEventTimeIndex(terminalNode);
  unstackNelsonAalen(terminalNode);
  unstackSurvival(terminalNode);
  if (terminalNode -> eTypeSize > 1) {
    unstackCSH(terminalNode);
    unstackCIF(terminalNode);
  }
}
void getNodeInfo(Node *nodePtr) {
  unsigned int i;
  Rprintf("\nNodeInfo:  %20x", nodePtr);
  Rprintf("\n   LeafCnt   SpltParm  ");
  Rprintf("\n%10d %10d \n", nodePtr -> nodeID, nodePtr -> splitParameter);
  if (nodePtr -> splitValueFactSize > 0) {
    Rprintf("FactorInfo %20x \n", nodePtr -> splitValueFactPtr);
    Rprintf("0x ");
    for (i = nodePtr -> splitValueFactSize; i >= 1; i--) {
      Rprintf("%8x ", (nodePtr -> splitValueFactPtr)[i]);
    }
  }
  else {
    Rprintf(" %12.4f \n", nodePtr -> splitValueCont);
  }
  Rprintf("\nSplit Statistic \n");
  Rprintf(" %12.4f \n", nodePtr -> splitStatistic);
  Rprintf("\nNode Variance \n");
  Rprintf(" %12.4f \n", nodePtr -> variance);
  if (nodePtr -> permissibleSplit != NULL) {
    Rprintf("\nPermissible Splits \n");
    for (i=1; i <= nodePtr -> xSize; i++) {
      Rprintf(" %10d", i);
    }
    Rprintf("\n");
    for (i=1; i <= nodePtr -> xSize; i++) {
      Rprintf(" %10d", (nodePtr -> permissibleSplit)[i]);
    }
    Rprintf("\n");
  }
  Rprintf("\n mpIndexSize   = %20d", nodePtr -> mpIndexSize);
  Rprintf("\n fmpIndexSize  = %20d", nodePtr -> fmpIndexSize);
  Rprintf("\n");
  Rprintf("\n mpSign       = %20x", nodePtr -> mpSign);
  Rprintf("\n fmpSign      = %20x", nodePtr -> fmpSign);
  Rprintf("\n");
  Rprintf("\n lmpIndexActualSize        = %20d", nodePtr -> lmpIndexActualSize);
  Rprintf("\n flmpIndexActualSize       = %20d", nodePtr -> flmpIndexActualSize);
  Rprintf("\n lmpIndexAllocSize         = %20d", nodePtr -> lmpIndexAllocSize);
  Rprintf("\n flmpIndexAllocSize        = %20d", nodePtr -> flmpIndexAllocSize);
  Rprintf("\n");
  Rprintf("\n lmpIndex            = %20x", nodePtr -> lmpIndex);
  Rprintf("\n flmpIndex           = %20x", nodePtr -> flmpIndex);
  Rprintf("\n");
  Rprintf("\n lmrIndexActualSize        = %20d", nodePtr -> lmrIndexActualSize);
  Rprintf("\n flmrIndexActualSize       = %20d", nodePtr -> flmrIndexActualSize);
  Rprintf("\n lmrIndexAllocSize         = %20d", nodePtr -> lmrIndexAllocSize);
  Rprintf("\n flmrIndexAllocSize        = %20d", nodePtr -> flmrIndexAllocSize);
  Rprintf("\n");
  Rprintf("\n lmrIndex            = %20x", nodePtr -> lmrIndex);
  Rprintf("\n flmrIndex           = %20x", nodePtr -> flmrIndex);
}
void setParent(Node *daughter, Node *parent) {
  daughter -> parent = parent;
}
void setLeftDaughter(Node *daughter, Node *parent) {
  parent -> left = daughter;
}
void setRightDaughter(Node *daughter, Node *parent) {
  parent -> right = daughter;
}
char forkNode(Node         *parent,
              unsigned int  splitParameter,
              double        splitValueMaxCont,
              unsigned int  splitValueMaxFactSize,
              unsigned int *splitValueMaxFactPtr) {
  if (parent == NULL) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    Rprintf("\nRF-SRC:  The parent node is NULL.");
    return FALSE;
  }
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    Rprintf("\nRF-SRC:  The daughter nodes are NON-NULL.");
    return FALSE;
  }
  if (parent -> splitFlag == FALSE) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    Rprintf("\nRF-SRC:  The split flag is FALSE.");
    return FALSE;
  }
  if (parent -> xSize < splitParameter) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    Rprintf("\nRF-SRC:  The split parameter index is out of range [1, xSize].");
    return FALSE;
  }
  if ((parent -> permissibleSplit)[splitParameter] == FALSE) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to forkNode().  ");
    Rprintf("\nRF-SRC:  The split parameter is marked unsplittable.");
    return FALSE;
  }
  Node *left  = makeNode(parent -> xSize);
  Node *right = makeNode(parent -> xSize);
  parent -> splitParameter = splitParameter;
  parent -> splitValueCont = splitValueMaxCont;
  parent -> splitValueFactSize = splitValueMaxFactSize;
  parent -> splitValueFactPtr = splitValueMaxFactPtr;
  setParent(left, parent);
  setParent(right, parent);
  setLeftDaughter(left, parent);
  setRightDaughter(right, parent);
  nrCopyVector(left  -> permissibleSplit, parent -> permissibleSplit, left -> xSize);
  nrCopyVector(right -> permissibleSplit, parent -> permissibleSplit, right -> xSize);
  free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
  parent -> permissibleSplit = NULL;
  parent -> splitFlag = FALSE;
  return TRUE;
}
void stackAtRiskAndEventCounts(Node *tNode, unsigned int eTypeSize, unsigned int mTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> mTimeSize > 0) {
    if (tNode -> mTimeSize != mTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  mTimeSize has been previously defined:  %10d vs %10d", tNode -> mTimeSize, mTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> mTimeSize = mTimeSize;
  }
  tNode -> atRiskCount     = uivector(1, mTimeSize);
  tNode -> eventCount      = uimatrix(1, eTypeSize, 1, mTimeSize);
}
void stackEventTimeIndex(Node *tNode, unsigned int eTimeSize) {
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> eventTimeIndex  = uivector(1, eTimeSize + 1);
}
void unstackAtRiskAndEventCounts(Node *tNode) {
  if (tNode -> atRiskCount != NULL) {
    free_uivector(tNode -> atRiskCount, 1, tNode -> mTimeSize);
    tNode -> atRiskCount = NULL;
  }
  if (tNode -> eventCount != NULL) {
    free_uimatrix(tNode -> eventCount, 1, tNode -> eTypeSize, 1, tNode -> mTimeSize);
    tNode -> eventCount = NULL;
  }
}
void unstackEventTimeIndex(Node *tNode) {
  if (tNode -> eventTimeIndex != NULL) {
    free_uivector(tNode -> eventTimeIndex, 1, tNode -> eTimeSize + 1);
    tNode -> eventTimeIndex = NULL;
  }
}
void stackLocalRatio(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> localRatio = dmatrix(1, eTypeSize, 1, tNode -> eTimeSize);
}
void unstackLocalRatio(Node *tNode) {
  if(tNode -> eTimeSize > 0) {
    if (tNode -> localRatio != NULL) {
      free_dmatrix(tNode -> localRatio, 1, tNode -> eTypeSize, 1, tNode -> eTimeSize);
      tNode -> localRatio = NULL;
    }
  }
}
void stackLocalSurvival(Node *tNode, unsigned int eTimeSize) {
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> localSurvival = dvector(1, tNode -> eTimeSize);
}
void unstackLocalSurvival(Node *tNode) {
  if(tNode -> eTimeSize > 0) {
    if (tNode -> localSurvival != NULL) {
      free_dvector(tNode -> localSurvival, 1, tNode -> eTimeSize);
      tNode -> localSurvival = NULL;
    }
  }
}
void stackLocalNelsonAalen(Node *tNode, unsigned int eTimeSize) {
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> localNelsonAalen = dvector(1, tNode -> eTimeSize);
}
void unstackLocalNelsonAalen(Node *tNode) {
  if(tNode -> eTimeSize > 0) {
    if (tNode -> localNelsonAalen != NULL) {
      free_dvector(tNode -> localNelsonAalen, 1, tNode -> eTimeSize);
      tNode -> localNelsonAalen = NULL;
    }
  }
}
void stackLocalCSH(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> localCSH = dmatrix(1, eTypeSize, 1, tNode -> eTimeSize);
}
void unstackLocalCSH(Node *tNode) {
  if(tNode -> eTimeSize > 0) {
    if (tNode -> localCSH != NULL) {
      free_dmatrix(tNode -> localCSH, 1, tNode -> eTypeSize, 1, tNode -> eTimeSize);
      tNode -> localCSH = NULL;
    }
  }
}
void stackLocalCIF(Node *tNode, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> eTimeSize > 0) {
    if (tNode -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tNode -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTimeSize = eTimeSize;
  }
  tNode -> localCIF = dmatrix(1, eTypeSize, 1, tNode -> eTimeSize);
}
void unstackLocalCIF(Node *tNode) {
  if(tNode -> eTimeSize > 0) {
    if (tNode -> localCIF != NULL) {
      free_dmatrix(tNode -> localCIF, 1, tNode -> eTypeSize, 1, tNode -> eTimeSize);
      tNode -> localCIF = NULL;
    }
  }
}
void stackNelsonAalen(Node *tNode, unsigned int sTimeSize) {
  if (tNode -> sTimeSize > 0) {
    if (tNode -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tNode -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> sTimeSize = sTimeSize;
  }
  tNode -> nelsonAalen = dvector(1, tNode -> sTimeSize);
}
void unstackNelsonAalen(Node *tNode) {
  if(tNode -> sTimeSize > 0) {
    if (tNode -> nelsonAalen != NULL) {
      free_dvector(tNode -> nelsonAalen, 1, tNode -> sTimeSize);
      tNode -> nelsonAalen = NULL;
    }
  }
}
void stackSurvival(Node *tNode, unsigned int sTimeSize) {
  if (tNode -> sTimeSize > 0) {
    if (tNode -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tNode -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> sTimeSize = sTimeSize;
  }
  tNode -> survival = dvector(1, tNode -> sTimeSize);
}
void unstackSurvival(Node *tNode) {
  if(tNode -> sTimeSize > 0) {
    if (tNode -> survival != NULL) {
      free_dvector(tNode -> survival, 1, tNode -> sTimeSize);
      tNode -> survival = NULL;
    }
  }
}
void stackCSH(Node *tNode, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> sTimeSize > 0) {
    if (tNode -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tNode -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> sTimeSize = sTimeSize;
  }
  tNode -> CSH = dmatrix(1, eTypeSize, 1, tNode -> sTimeSize);
}
void unstackCSH(Node *tNode) {
  if(tNode -> sTimeSize > 0) {
    if (tNode -> CSH != NULL) {
      free_dmatrix(tNode -> CSH, 1, tNode -> eTypeSize, 1, tNode -> sTimeSize);
      tNode -> CSH = NULL;
    }
  }
}
void stackCIF(Node *tNode, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  if (tNode -> sTimeSize > 0) {
    if (tNode -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tNode -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> sTimeSize = sTimeSize;
  }
  tNode -> CIF = dmatrix(1, eTypeSize, 1, tNode -> sTimeSize);
}
void unstackCIF(Node *tNode) {
  if(tNode -> sTimeSize > 0) {
    if (tNode -> CIF != NULL) {
      free_dmatrix(tNode -> CIF, 1, tNode -> eTypeSize, 1, tNode -> sTimeSize);
      tNode -> CIF = NULL;
    }
  }
}
void stackMortality(Node *tNode, unsigned int eTypeSize) {
  if (tNode -> eTypeSize > 0) {
    if (tNode -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tNode -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> eTypeSize = eTypeSize;
  }
  tNode -> mortality = dvector(1, eTypeSize);
}
void unstackMortality(Node *tNode) {
  if(tNode -> eTypeSize > 0) {
    if (tNode -> mortality != NULL) {
      free_dvector(tNode -> mortality, 1, tNode -> eTypeSize);
      tNode -> mortality = NULL;
    }
  }
}
void stackMPSign(Node *tNode, unsigned int mpIndexSize) {
  if (tNode -> mpIndexSize > 0) {
    if (tNode -> mpIndexSize != mpIndexSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  mpIndexSize has been previously defined:  %10d vs %10d", tNode -> mpIndexSize, mpIndexSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> mpIndexSize = mpIndexSize;
  }
  tNode -> mpSign = ivector(1, tNode -> mpIndexSize);
}
void unstackMPSign(Node *tNode) {
  if(tNode -> mpIndexSize > 0) {
    if (tNode -> mpSign != NULL) {
      free_ivector(tNode -> mpSign, 1, tNode -> mpIndexSize);
      tNode -> mpSign = NULL;
    }
  }
}
void stackFMPSign(Node *tNode, unsigned int fmpIndexSize) {
  if (tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpIndexSize != fmpIndexSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  fmpIndexSize has been previously defined:  %10d vs %10d", tNode -> fmpIndexSize, fmpIndexSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> fmpIndexSize = fmpIndexSize;
  }
  tNode -> fmpSign = ivector(1, tNode -> fmpIndexSize);
}
void unstackFMPSign(Node *tNode) {
  if(tNode -> fmpIndexSize > 0) {
    if (tNode -> fmpSign != NULL) {
      free_ivector(tNode -> fmpSign, 1, tNode -> fmpIndexSize);
      tNode -> fmpSign = NULL;
    }
  }
}
void stackNodeLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> lmpIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmpIndex has been previously defined:  %10d vs %10d", tNode -> lmpIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmpIndexAllocSize = size;
  }
  tNode -> lmpIndex = uivector(1, tNode -> lmpIndexAllocSize);
}
void unstackNodeLMPIndex(Node *tNode) {
  if(tNode -> lmpIndexAllocSize > 0) {
    if (tNode -> lmpIndex != NULL) {
      free_uivector(tNode -> lmpIndex, 1, tNode -> lmpIndexAllocSize);
      tNode -> lmpIndex = NULL;
      tNode -> lmpIndexAllocSize = 0;
    }
  }
}
void stackNodeFLMPIndex(Node *tNode, unsigned int size) {
  if (tNode -> flmpIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmpIndex has been previously defined:  %10d vs %10d", tNode -> flmpIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmpIndexAllocSize = size;
  }
  tNode -> flmpIndex = uivector(1, tNode -> flmpIndexAllocSize);
}
void unstackNodeFLMPIndex(Node *tNode) {
  if(tNode -> flmpIndexAllocSize > 0) {
    if (tNode -> flmpIndex != NULL) {
      free_uivector(tNode -> flmpIndex, 1, tNode -> flmpIndexAllocSize);
      tNode -> flmpIndex = NULL;
      tNode -> flmpIndexAllocSize = 0;
    }
  }
}
void stackNodeLMRIndex(Node *tNode, unsigned int size) {
  if (tNode -> lmrIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmrIndex has been previously defined:  %10d vs %10d", tNode -> lmrIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmrIndexAllocSize = size;
  }
  tNode -> lmrIndex = uivector(1, tNode -> lmrIndexAllocSize);
}
void unstackNodeLMRIndex(Node *tNode) {
  if(tNode -> lmrIndexAllocSize > 0) {
    if (tNode -> lmrIndex != NULL) {
      free_uivector(tNode -> lmrIndex, 1, tNode -> lmrIndexAllocSize);
      tNode -> lmrIndex = NULL;
      tNode -> lmrIndexAllocSize = 0;
    }
  }
}
void stackNodeFLMRIndex(Node *tNode, unsigned int size) {
  if (tNode -> flmrIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmrIndex has been previously defined:  %10d vs %10d", tNode -> flmrIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmrIndexAllocSize = size;
  }
  tNode -> flmrIndex = uivector(1, tNode -> flmrIndexAllocSize);
}
void unstackNodeFLMRIndex(Node *tNode) {
  if(tNode -> flmrIndexAllocSize > 0) {
    if (tNode -> flmrIndex != NULL) {
      free_uivector(tNode -> flmrIndex, 1, tNode -> flmrIndexAllocSize);
      tNode -> flmrIndex = NULL;
      tNode -> flmrIndexAllocSize = 0;
    }
  }
}
void stackTermLMIIndex(Terminal *tNode, unsigned int size) {
  if (tNode -> lmiSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmiIndex has been previously defined:  %10d vs %10d", tNode -> lmiSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmiSize = size;
  }
  tNode -> lmiIndex = uivector(1, tNode -> lmiSize);
  tNode -> lmiValue = dvector(1, tNode -> lmiSize);
}
void unstackTermLMIIndex(Terminal *tNode) {
  if(tNode -> lmiSize > 0) {
    if (tNode -> lmiIndex != NULL) {
      free_uivector(tNode -> lmiIndex, 1, tNode -> lmiSize);
      tNode -> lmiIndex = NULL;
    }
    if (tNode -> lmiValue != NULL) {
      free_dvector(tNode -> lmiValue, 1, tNode -> lmiSize);
      tNode -> lmiValue = NULL;
    }
  }
}
void stackMultiClassProb(Node *tNode, unsigned int rfCount, unsigned int *rfSize) {
  unsigned int j;
  if (tNode -> rfCount > 0) {
    if (tNode -> rfCount != rfCount) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tNode -> rfCount, rfCount);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> rfCount = rfCount;
  }
  tNode -> rfSize = uivector(1, tNode -> rfCount);
  tNode -> multiClassProb = (unsigned int **) new_vvector(1, tNode -> rfCount, NRUTIL_UPTR);
  for (j = 1; j <= tNode -> rfCount; j++) {
    (tNode -> rfSize)[j] = rfSize[j];
    (tNode -> multiClassProb)[j] = uivector(1, (tNode -> rfSize)[j]);
  }
}
void unstackMultiClassProb(Node *tNode) {
  unsigned int j;
  if (tNode -> rfCount > 0) {
    if (tNode -> rfSize != NULL) {
      if (tNode -> multiClassProb != NULL) {
        for (j = 1; j <= tNode -> rfCount; j++) {
          if (tNode -> multiClassProb[j] != NULL) {
            free_uivector(tNode -> multiClassProb[j], 1, tNode -> rfSize[j]);
            tNode -> multiClassProb[j] = NULL;
          }
        }
        free_new_vvector(tNode -> multiClassProb, 1, tNode -> rfCount, NRUTIL_UPTR);
        tNode -> multiClassProb = NULL;
      }
    }
    free_uivector(tNode -> rfSize, 1, tNode -> rfCount);
    tNode -> rfSize = NULL;
  }
}
void stackSplitDepth(Node *tNode, unsigned int depth) {
  if (tNode -> depth > 0) {
    if (tNode -> depth != depth) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  depth has been previously defined:  %10d vs %10d", tNode -> depth, depth);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> depth = depth;
  }
  tNode -> splitDepth = uivector(1, tNode -> depth);
}
void unstackSplitDepth(Node *tNode) {
  if (tNode -> splitDepth != NULL) {
    free_uivector(tNode -> splitDepth, 1, tNode -> depth);
    tNode -> splitDepth = NULL;
  }
}
void checkEventType(Node *tNode) {
  if (tNode -> eTypeSize <= 1) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize must be greater than one (1):  %10d ", tNode -> eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
  }
}
