////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.2
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
extern unsigned int getTraceFlag(unsigned int tree);
unsigned int getForkDefTraceFlag();
unsigned int getTurnOffTraceFlag();
unsigned int getTurnOnTraceFlag();
#include <R_ext/Arith.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
Terminal *makeTerminal() {
  Terminal *parent = (Terminal*) gblock((size_t) sizeof(Terminal));
  parent -> lmvIndex      = NULL;
  parent -> lmvIndexSize  = 0;
  parent -> flmvIndex     = NULL;
  parent -> flmvIndexSize = 0;;
  parent -> lmrIndex      = NULL;
  parent -> lmrIndexSize  = 0;
  parent -> flmrIndex     = NULL;
  parent -> flmrIndexSize = 0;;
  parent -> lmiIndex      = NULL;
  parent -> lmiValue      = NULL;
  parent -> lmiSizePtr    = NULL;
  parent -> lmiSize       = 0;
  parent -> flmiIndex     = NULL;
  parent -> flmiValue     = NULL;
  parent -> flmiSizePtr   = NULL;
  parent -> flmiSize      = 0;
  parent -> leafCount     = 0;
  parent -> mate          = NULL;
  parent -> dominant      = 0;
  parent -> fdominant     = 0;
  return parent;
}
void freeTerminal(Terminal        *parent) {
  unstackTermLMVIndex(parent);
  unstackTermFLMVIndex(parent);
  unstackTermLMRIndex(parent);
  unstackTermFLMRIndex(parent);
  unstackTermFLMISizePtr(parent);
  unstackTermLMISizePtr(parent);
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
  parent -> leafCount            = 0;
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
  parent -> membrCount           = 0;
  parent -> mvSignSize               = 0;
  parent -> fmvSignSize              = 0;
  parent -> mvSign               = NULL;
  parent -> fmvSign              = NULL;
  parent -> lmvIndex             = NULL;
  parent -> flmvIndex            = NULL;
  parent -> lmvIndexAllocSize    = 0;
  parent -> flmvIndexAllocSize   = 0;
  parent -> lmvIndexActualSize   = 0;
  parent -> flmvIndexActualSize  = 0;
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
  unstackMVSign(parent);
  unstackFMVSign(parent);
  unstackNodeLMVIndex(parent);
  unstackNodeFLMVIndex(parent);
  unstackNodeLMRIndex(parent);
  unstackNodeFLMRIndex(parent);
  if ((parent -> splitParameter) == 0) {
    freeTerminalNodeSurvivalStructures(parent);
  }
  free_gblock(parent, sizeof(Node));
}
void freeTerminalNodeSurvivalStructures(Node *terminalNode) {
  unstackAtRisk(terminalNode);
  unstackLocalRatio(terminalNode);
  unstackLocalSurvival(terminalNode);
  unstackLocalNelsonAalen(terminalNode);
  if (terminalNode -> eTypeSize > 1) {
    unstackLocalCSH(terminalNode);
    unstackLocalCIF(terminalNode);
  }
  unstackNelsonAalen(terminalNode);
  unstackSurvival(terminalNode);
  if (terminalNode -> eTypeSize > 1) {
    unstackCSH(terminalNode);
    unstackCIF(terminalNode);
  }
}
void getNodeInfo(Node *leaf) {
  unsigned int i;
  Rprintf("\nNodeInfo:  %20x", leaf);
  Rprintf("\n   LeafCnt   SpltParm  ");
  Rprintf("\n%10d %10d \n", leaf -> leafCount, leaf -> splitParameter);
  if (leaf -> splitValueFactSize > 0) {
    Rprintf("FactorInfo %20x \n", leaf -> splitValueFactPtr);
    Rprintf("0x ");
    for (i = leaf -> splitValueFactSize; i >= 1; i--) {
      Rprintf("%8x ", (leaf -> splitValueFactPtr)[i]);
    }
  }
  else {
    Rprintf(" %12.4f \n", leaf -> splitValueCont);
  }
  Rprintf("\nSplit Statistic \n");
  Rprintf(" %12.4f \n", leaf -> splitStatistic);
  Rprintf("\nNode Variance \n");
  Rprintf(" %12.4f \n", leaf -> variance);
  if (leaf -> permissibleSplit != NULL) {
    Rprintf("\nPermissible Splits \n");
    for (i=1; i <= leaf -> xSize; i++) {
      Rprintf(" %10d", i);
    }
    Rprintf("\n");
    for (i=1; i <= leaf -> xSize; i++) {
      Rprintf(" %10d", (leaf -> permissibleSplit)[i]);
    }
    Rprintf("\n");
  }
  Rprintf("\n mvSignSize   = %20d", leaf -> mvSignSize);
  Rprintf("\n fmvSignSize  = %20d", leaf -> fmvSignSize);
  Rprintf("\n");
  Rprintf("\n mvSign       = %20x", leaf -> mvSign);
  Rprintf("\n fmvSign      = %20x", leaf -> fmvSign);
  Rprintf("\n");
  Rprintf("\n lmvIndexActualSize        = %20d", leaf -> lmvIndexActualSize);
  Rprintf("\n flmvIndexActualSize       = %20d", leaf -> flmvIndexActualSize);
  Rprintf("\n lmvIndexAllocSize         = %20d", leaf -> lmvIndexAllocSize);
  Rprintf("\n flmvIndexAllocSize        = %20d", leaf -> flmvIndexAllocSize);
  Rprintf("\n");
  Rprintf("\n lmvIndex            = %20x", leaf -> lmvIndex);
  Rprintf("\n flmvIndex           = %20x", leaf -> flmvIndex);
  Rprintf("\n");
  Rprintf("\n lmrIndexActualSize        = %20d", leaf -> lmrIndexActualSize);
  Rprintf("\n flmrIndexActualSize       = %20d", leaf -> flmrIndexActualSize);
  Rprintf("\n lmrIndexAllocSize         = %20d", leaf -> lmrIndexAllocSize);
  Rprintf("\n flmrIndexAllocSize        = %20d", leaf -> flmrIndexAllocSize);
  Rprintf("\n");
  Rprintf("\n lmrIndex            = %20x", leaf -> lmrIndex);
  Rprintf("\n flmrIndex           = %20x", leaf -> flmrIndex);
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
void stackAtRisk(Node *tNode, unsigned int eTypeSize, unsigned int mTimeSize) {
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
  tNode -> eventTimeIndex  = uivector(1, mTimeSize);
}
void unstackAtRisk(Node *tNode) {
  if (tNode -> atRiskCount != NULL) {
    free_uivector(tNode -> atRiskCount, 1, tNode -> mTimeSize);
    tNode -> atRiskCount = NULL;
  }
  if (tNode -> eventCount != NULL) {
    free_uimatrix(tNode -> eventCount, 1, tNode -> eTypeSize, 1, tNode -> mTimeSize);
    tNode -> eventCount = NULL;
  }
  if (tNode -> eventTimeIndex != NULL) {
    free_uivector(tNode -> eventTimeIndex, 1, tNode -> mTimeSize);
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
void stackMVSign(Node *tNode, unsigned int mvSignSize) {
  if (tNode -> mvSignSize > 0) {
    if (tNode -> mvSignSize != mvSignSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  mvSignSize has been previously defined:  %10d vs %10d", tNode -> mvSignSize, mvSignSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> mvSignSize = mvSignSize;
  }
  tNode -> mvSign = ivector(1, tNode -> mvSignSize);
}
void unstackMVSign(Node *tNode) {
  if(tNode -> mvSignSize > 0) {
    if (tNode -> mvSign != NULL) {
      free_ivector(tNode -> mvSign, 1, tNode -> mvSignSize);
      tNode -> mvSign = NULL;
    }
  }
}
void stackFMVSign(Node *tNode, unsigned int fmvSignSize) {
  if (tNode -> fmvSignSize > 0) {
    if (tNode -> fmvSignSize != fmvSignSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  fmvSignSize has been previously defined:  %10d vs %10d", tNode -> fmvSignSize, fmvSignSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tNode -> fmvSignSize = fmvSignSize;
  }
  tNode -> fmvSign = ivector(1, tNode -> fmvSignSize);
}
void unstackFMVSign(Node *tNode) {
  if(tNode -> fmvSignSize > 0) {
    if (tNode -> fmvSign != NULL) {
      free_ivector(tNode -> fmvSign, 1, tNode -> fmvSignSize);
      tNode -> fmvSign = NULL;
    }
  }
}
void stackNodeLMVIndex(Node *tNode, unsigned int size) {
  if (tNode -> lmvIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmvIndex has been previously defined:  %10d vs %10d", tNode -> lmvIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmvIndexAllocSize = size;
  }
  tNode -> lmvIndex = uivector(1, tNode -> lmvIndexAllocSize);
}
void unstackNodeLMVIndex(Node *tNode) {
  if(tNode -> lmvIndexAllocSize > 0) {
    if (tNode -> lmvIndex != NULL) {
      free_uivector(tNode -> lmvIndex, 1, tNode -> lmvIndexAllocSize);
      tNode -> lmvIndex = NULL;
      tNode -> lmvIndexAllocSize = 0;
    }
  }
}
void stackNodeFLMVIndex(Node *tNode, unsigned int size) {
  if (tNode -> flmvIndexAllocSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmvIndex has been previously defined:  %10d vs %10d", tNode -> flmvIndexAllocSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmvIndexAllocSize = size;
  }
  tNode -> flmvIndex = uivector(1, tNode -> flmvIndexAllocSize);
}
void unstackNodeFLMVIndex(Node *tNode) {
  if(tNode -> flmvIndexAllocSize > 0) {
    if (tNode -> flmvIndex != NULL) {
      free_uivector(tNode -> flmvIndex, 1, tNode -> flmvIndexAllocSize);
      tNode -> flmvIndex = NULL;
      tNode -> flmvIndexAllocSize = 0;
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
void stackTermLMVIndex(Terminal *tNode, unsigned int size) {
  if (tNode -> lmvIndexSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmvIndex has been previously defined:  %10d vs %10d", tNode -> lmvIndexSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmvIndexSize = size;
  }
  tNode -> lmvIndex = uivector(1, tNode -> lmvIndexSize);
}
void unstackTermLMVIndex(Terminal *tNode) {
  if(tNode -> lmvIndexSize > 0) {
    if (tNode -> lmvIndex != NULL) {
      free_uivector(tNode -> lmvIndex, 1, tNode -> lmvIndexSize);
      tNode -> lmvIndex = NULL;
    }
  }
}
void stackTermFLMVIndex(Terminal *tNode, unsigned int size) {
  if (tNode -> flmvIndexSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmvIndex has been previously defined:  %10d vs %10d", tNode -> flmvIndexSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmvIndexSize = size;
  }
  tNode -> flmvIndex = uivector(1, tNode -> flmvIndexSize);
}
void unstackTermFLMVIndex(Terminal *tNode) {
  if(tNode -> flmvIndexSize > 0) {
    if (tNode -> flmvIndex != NULL) {
      free_uivector(tNode -> flmvIndex, 1, tNode -> flmvIndexSize);
      tNode -> flmvIndex = NULL;
    }
  }
}
void stackTermLMRIndex(Terminal *tNode, unsigned int size) {
  if (tNode -> lmrIndexSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmrIndex has been previously defined:  %10d vs %10d", tNode -> lmrIndexSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmrIndexSize = size;
  }
  tNode -> lmrIndex = uivector(1, tNode -> lmrIndexSize);
}
void unstackTermLMRIndex(Terminal *tNode) {
  if(tNode -> lmrIndexSize > 0) {
    if (tNode -> lmrIndex != NULL) {
      free_uivector(tNode -> lmrIndex, 1, tNode -> lmrIndexSize);
      tNode -> lmrIndex = NULL;
    }
  }
}
void stackTermFLMRIndex(Terminal *tNode, unsigned int size) {
  if (tNode -> flmrIndexSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmrIndex has been previously defined:  %10d vs %10d", tNode -> flmrIndexSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmrIndexSize = size;
  }
  tNode -> flmrIndex = uivector(1, tNode -> flmrIndexSize);
}
void unstackTermFLMRIndex(Terminal *tNode) {
  if(tNode -> flmrIndexSize > 0) {
    if (tNode -> flmrIndex != NULL) {
      free_uivector(tNode -> flmrIndex, 1, tNode -> flmrIndexSize);
      tNode -> flmrIndex = NULL;
    }
  }
}
void stackTermFLMISizePtr(Terminal *tNode, unsigned int size) {
  if (tNode -> flmiSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  flmiSize has been previously defined:  %10d vs %10d", tNode -> flmiSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> flmiSize = size;
  }
  tNode -> flmiSizePtr = uivector(1, tNode -> flmiSize);
}
void unstackTermFLMISizePtr(Terminal *tNode) {
  if(tNode -> flmiSize > 0) {
    if (tNode -> flmiSizePtr != NULL) {
      free_uivector(tNode -> flmiSizePtr, 1, tNode -> flmiSize);
      tNode -> flmiSizePtr = NULL;
    }
  }
}
void stackTermLMISizePtr(Terminal *tNode, unsigned int size) {
  if (tNode -> lmiSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmiSize has been previously defined:  %10d vs %10d", tNode -> lmiSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tNode -> lmiSize = size;
  }
  tNode -> lmiSizePtr = uivector(1, tNode -> lmiSize);
}
void unstackTermLMISizePtr(Terminal *tNode) {
  if(tNode -> lmiSize > 0) {
    if (tNode -> lmiSizePtr != NULL) {
      free_uivector(tNode -> lmiSizePtr, 1, tNode -> lmiSize);
      tNode -> lmiSizePtr = NULL;
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
  tNode -> multiClassProb = (unsigned int **) vvector(1, tNode -> rfCount);
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
        free_vvector(tNode -> multiClassProb, 1, tNode -> rfCount);
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
