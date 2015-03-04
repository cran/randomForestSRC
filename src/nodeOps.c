////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.6.1
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
#include       "termOps.h"
#include   <R_ext/Print.h>
#include   <Rdefines.h>
extern unsigned int getTraceFlag();
extern unsigned int getNodeDefTraceFlag();
extern unsigned int getForkDefTraceFlag();
extern unsigned int getTurnOffTraceFlag();
#include <R_ext/Arith.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
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
  free_gblock(parent, sizeof(Node));
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
