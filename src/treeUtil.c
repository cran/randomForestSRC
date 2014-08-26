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


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include       "nodeOps.h"
#include     "factorOps.h"
#include     "bootstrap.h"
#include        "impute.h"
#include         "split.h"
#include      "treeUtil.h"
Node *getTerminalNode(uint treeID, uint leaf) {
  uint i, j;
  Node *parent;
  parent = NULL;
  for (j = 1; j <= RF_observationSize; j++) {
    if ((RF_tNodeMembership[treeID][j] -> nodeID) == leaf) {
      parent = RF_tNodeMembership[treeID][j];
      j = RF_observationSize;
    }
  }
  if (parent == NULL) {
    Rprintf("\nDiagnostic Trace of (individual, boot, node, leaf) vectors in data set:  ");
    Rprintf("\n        index         boot         node         leaf \n");
    for (i = 1; i <= RF_observationSize; i++) {
      Rprintf(" %12d %12d %12x %12d \n", i,
              RF_bootMembershipFlag[treeID][i], RF_tNodeMembership[treeID][i],
              RF_tNodeMembership[treeID][i] -> nodeID);
    }
    Rprintf("\nDiagnostic State of TRAIN (SHADOW) data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, RF_status[treeID][j], RF_time[treeID][j]);
      for (i=1; i <= RF_xSize; i++) {
        Rprintf(" %12.4f", (RF_observation[treeID][i][j]));
      }
      Rprintf("\n");
    }
    Rprintf("\nDiagnostic State of TRAIN (INCOMING) data:  ");
    Rprintf("\n       index       status         time   observations -> \n");
    Rprintf("\n                                      ");
    for (i=1; i <= RF_xSize; i++) {
      Rprintf(" %12d", i);
    }
    Rprintf("\n");
    for (j = 1; j <= RF_observationSize; j++) {
      Rprintf("%12d %12.4f %12.4f", j, RF_responseIn[RF_statusIndex][j], RF_responseIn[RF_timeIndex][j]);
      for (i=1; i <= RF_xSize; i++) {
        Rprintf(" %12.4f", (RF_observationIn[i][j]));
      }
      Rprintf("\n");
    }
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Proxy member for (tree, node) = (%12d, %12d) not found.", treeID, leaf);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  return parent;
}
void getRawNodeSize(uint  type,
                    uint  treeID,
                    Node *parent,
                    uint *repMembrIndx,
                    uint *repMembrSize,
                    uint *allMembrIndx,
                    uint *allMembrSize) {
  uint      obsSize;
  Node   ***nodeMembershipPtr;
  uint i;
  obsSize           = 0;     
  nodeMembershipPtr = NULL;  
  switch (type) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    nodeMembershipPtr = RF_ftNodeMembership;
    break;
  default:
    obsSize = RF_observationSize;
    nodeMembershipPtr = RF_tNodeMembership;
    break;
  }
  *repMembrSize = 0;
  for (i=1; i <= RF_observationSize; i++) {
    if (RF_tNodeMembership[treeID][RF_bootMembershipIndex[treeID][i]] == parent) {
      repMembrIndx[++(*repMembrSize)] = RF_bootMembershipIndex[treeID][i];
    }
  }
  *allMembrSize = 0;
  for (i=1; i <= obsSize; i++) {
    if (nodeMembershipPtr[treeID][i] == parent) {
      allMembrIndx[++(*allMembrSize)] = i;
    }
  }
}
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
                   char    multImpFlag,
                   char   *membershipIndicator,
                   uint   *leftDaughterSize,
                   uint   *rghtDaughterSize) {
  char factorFlag;
  char daughterFlag;
  char *randomMembrFlag;
  char result;
  uint nonMissAllMembrSize;
  double leftProbability;
  char mPredictorFlag;
  uint offset;
  uint i;
  factorFlag = FALSE; 
  randomMembrFlag = cvector(1, allMembrSize);
  result = forkNode(parent,
                    splitParameterMax,
                    splitValueMaxCont,
                    splitValueMaxFactSize,
                    splitValueMaxFactPtr);
  if (result == TRUE) {
    if (RF_opt & OPT_TREE) {
      RF_nodeCount[treeID] += 2;
    }
    parent -> splitStatistic = splitStatistic;
    RF_tLeafCount[treeID]++;
    ((parent -> left) -> nodeID) = (parent -> nodeID);
    ((parent -> right) -> nodeID) = RF_tLeafCount[treeID];
    factorFlag = FALSE;
    if (strcmp(RF_xType[splitParameterMax], "C") == 0) {
      factorFlag = TRUE;
      if (RF_opt & OPT_TREE) {
        RF_mwcpCount[treeID] += parent -> splitValueFactSize;
      }
    }
    *leftDaughterSize = *rghtDaughterSize = 0;
    for (i = 1; i <= allMembrSize; i++) {
      membershipIndicator[allMembrIndx[i]] = NEITHER;
    }
    for (i = 1; i <= repMembrSize; i++) {
      membershipIndicator[repMembrIndx[i]] = localSplitIndicator[i];
    }
    offset = RF_rSize + splitParameterMax;
    for (i = 1; i <= allMembrSize; i++) {
      mPredictorFlag = FALSE;
      if (RF_mRecordSize > 0) {
        if (RF_mRecordMap[allMembrIndx[i]] > 0) {
          if (RF_mpSign[offset][RF_mRecordMap[allMembrIndx[i]]] == 1) {
            if ((RF_optHigh & OPT_MISS_RAND) && (!multImpFlag)) {
              mPredictorFlag = TRUE;
            }
          }
        }
      }
      randomMembrFlag[i] = mPredictorFlag;
    }
    for (i = 1; i <= allMembrSize; i++) {
      if (randomMembrFlag[i] == FALSE) {
        if(membershipIndicator[allMembrIndx[i]] == NEITHER) {
          daughterFlag = RIGHT;
          if (factorFlag == TRUE) {
            if (RF_observation[treeID][splitParameterMax][allMembrIndx[i]] != 0) {
              daughterFlag = splitOnFactor((uint) RF_observation[treeID][splitParameterMax][allMembrIndx[i]], splitValueMaxFactPtr);
            }
            else {
              Rprintf("\nRF-SRC:  *** ERROR *** ");
              Rprintf("\nRF-SRC:  Attempt to fork on NA factor value on (index, parameter):  (%10d, %10d)", allMembrIndx[i], splitParameterMax);
              Rprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
          else {
            if (!ISNA(RF_observation[treeID][splitParameterMax][allMembrIndx[i]])) {
              if (RF_observation[treeID][splitParameterMax][allMembrIndx[i]] <= splitValueMaxCont) {
                daughterFlag = LEFT;
              }
            }
            else {
              Rprintf("\nRF-SRC:  *** ERROR *** ");
              Rprintf("\nRF-SRC:  Attempt to fork on NA real value on (index, parameter):  (%10d, %10d)", allMembrIndx[i], splitParameterMax);
              Rprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
        }
        else {
          daughterFlag = membershipIndicator[allMembrIndx[i]];
        }
        membershipIndicator[allMembrIndx[i]] = daughterFlag;
        if (daughterFlag == LEFT) {
          (*leftDaughterSize) ++;
          RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> left;
        }
        else {
          (*rghtDaughterSize) ++;
          RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> right;
        }
      }
      else {
      }  
    }  
  }
  else {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  forkNode() failed.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nonMissAllMembrSize = (*leftDaughterSize) + (*rghtDaughterSize);
  leftProbability = (double) *leftDaughterSize / (double) nonMissAllMembrSize;
  for (i = 1; i <= allMembrSize; i++) {
    if (randomMembrFlag[i] == TRUE) {
      if (ran1A(treeID) <= leftProbability) {
        daughterFlag = LEFT;
        membershipIndicator[allMembrIndx[i]] = LEFT;
        (*leftDaughterSize) ++;
        RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> left;
      }
      else {
        daughterFlag = RIGHT;
        membershipIndicator[allMembrIndx[i]] = RIGHT;
        (*rghtDaughterSize) ++;
        RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> right;
      }
    }
  }
  if (localSplitIndicator != NULL) {
    free_cvector(localSplitIndicator, 1, repMembrSize);
  }
  else {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  NULL Local Split Indicator encountered in forkAndUpdate().");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_cvector(randomMembrFlag, 1, allMembrSize);
  return result;
}
char growTree (char     rootFlag,
               char     multImpFlag,
               uint     treeID,
               Node    *parent,
               uint    *repMembrIndx,
               uint     repMembrSize,
               uint    *allMembrIndx,
               uint     allMembrSize,
               uint     depth,
               uint    *maximumDepth,
               uint    *bootMembrIndxIter) {
  char  bootResult;
  char  splitResult;
  char  forkResult;
  char leftResult, rghtResult;
  char tnUpdateFlag;
  char bsUpdateFlag;
  uint *bootMembrIndx;
  uint *leftRepMembrIndx;
  uint *rghtRepMembrIndx;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint bootMembrSize;
  uint leftAllMembrSize;
  uint rghtAllMembrSize;
  uint leftRepMembrSize, jLeft;
  uint rghtRepMembrSize, jRght;
  uint     splitParameterMax;
  double   splitValueMaxCont;
  uint     splitValueMaxFactSize;
  uint    *splitValueMaxFactPtr;
  double   splitStatistic;
  char    *splitIndicator;
  uint i, p;
  parent -> depth = depth;
  bootResult = TRUE;
  tnUpdateFlag = TRUE;
  bsUpdateFlag = FALSE;
  splitIndicator = NULL;
  if (rootFlag | (RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    bootMembrIndx  = uivector(1, allMembrSize);
    bootMembrSize = allMembrSize;
    bootResult = bootstrap (RF_GROW,
                            treeID,
                            parent,
                            allMembrIndx,
                            allMembrSize,
                            bootMembrIndx);
    if (rootFlag & bootResult) {
      if (RF_opt & OPT_MEMB) {
        for (i=1; i <=  allMembrSize; i++) {
          RF_bootstrapMembershipPtr[treeID][bootMembrIndx[i]] ++;
        }
      }
      if (!(RF_opt & (OPT_BOOT_NODE | OPT_BOOT_NONE))) {
        bsUpdateFlag = TRUE;
      }
      if (RF_mRecordSize > 0) {
        for (p = 1; p <= RF_mpIndexSize; p++) {
          if (RF_mpIndex[p] > 0) {
            if (parent -> mpSign[p] == -1) {
              (parent -> permissibleSplit)[RF_mpIndex[p]] = FALSE;
            }
          }
        }
      }
    }
  }
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mpSign = (parent -> parent) -> mpSign;
  }
  if (bootResult) {
    if (!(RF_optHigh & OPT_MISS_RAND)) {
    if (multImpFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        imputeNode(RF_GROW,
                   FALSE,
                   TRUE,
                   treeID,
                   parent,
                   bootMembrIndx,
                   bootMembrSize,
                   allMembrIndx,
                   allMembrSize);
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(treeID,
                                 allMembrIndx,
                                 allMembrSize,
                                 RF_time[treeID],
                                 (RF_optHigh & OPT_MISS_RAND) ? TRUE : FALSE,
                                 FALSE,
                                 RF_masterTimeIndex[treeID]);
          }
        }
      }
    }  
    }
  }  
  if (bootResult) {
    if (rootFlag) {
      if (RF_opt & OPT_TREE) {
        RF_nodeCount[treeID] = 1;
      }
      RF_tLeafCount[treeID] = 1;
    }
    splitResult = getBestSplit(treeID,
                               parent,
                               bootMembrIndx,
                               bootMembrSize,
                               allMembrIndx,
                               allMembrSize,
                               & splitParameterMax,
                               & splitValueMaxCont,
                               & splitValueMaxFactSize,
                               & splitValueMaxFactPtr,
                               & splitStatistic,
                               & splitIndicator,
                               multImpFlag);
    if (splitResult == TRUE) {
      tnUpdateFlag = FALSE;
      char *membershipIndicator = cvector(1, RF_observationSize);
      forkResult = forkAndUpdate(treeID,
                                 parent,
                                 bootMembrIndx,
                                 bootMembrSize,
                                 allMembrIndx,
                                 allMembrSize,
                                 splitParameterMax,
                                 splitValueMaxCont,
                                 splitValueMaxFactSize,
                                 splitValueMaxFactPtr,
                                 splitStatistic,
                                 splitIndicator,
                                 multImpFlag,
                                 membershipIndicator,
                                 &leftAllMembrSize,
                                 &rghtAllMembrSize);
      if (forkResult == TRUE) {
        leftAllMembrIndx  = uivector(1, leftAllMembrSize);
        rghtAllMembrIndx  = uivector(1, rghtAllMembrSize);
        jLeft = jRght = 0;
        for (i = 1; i <= allMembrSize; i++) {
          if (membershipIndicator[allMembrIndx[i]] == LEFT) {
            leftAllMembrIndx[++jLeft] = allMembrIndx[i];
          }
          else {
            rghtAllMembrIndx[++jRght] = allMembrIndx[i];
          }
        }
        if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
          leftRepMembrIndx = leftAllMembrIndx;
          leftRepMembrSize = leftAllMembrSize;
          rghtRepMembrIndx = rghtAllMembrIndx;
          rghtRepMembrSize = rghtAllMembrSize;
        }
        else {
          leftRepMembrIndx  = uivector(1, bootMembrSize);
          rghtRepMembrIndx  = uivector(1, bootMembrSize);
          leftRepMembrSize = rghtRepMembrSize = 0;
          for (i = 1; i <= bootMembrSize; i++) {
            if (membershipIndicator[bootMembrIndx[i]] == LEFT) {
              leftRepMembrIndx[++leftRepMembrSize] = bootMembrIndx[i];
            }
            else {
              rghtRepMembrIndx[++rghtRepMembrSize] = bootMembrIndx[i];
            }
          }
        }
        free_cvector(membershipIndicator, 1, RF_observationSize);
        leftResult = growTree (FALSE,
                               multImpFlag,
                               treeID,
                               parent -> left,
                               leftRepMembrIndx,
                               leftRepMembrSize,
                               leftAllMembrIndx,
                               leftAllMembrSize,
                               (parent -> depth) + 1,
                               maximumDepth,
                               bootMembrIndxIter);
        if(!leftResult) {
        }
        rghtResult = growTree (FALSE,
                               multImpFlag,
                               treeID,
                               parent -> right,
                               rghtRepMembrIndx,
                               rghtRepMembrSize,
                               rghtAllMembrIndx,
                               rghtAllMembrSize,
                               (parent -> depth) + 1,
                               maximumDepth,
                               bootMembrIndxIter);
        if(!rghtResult) {
        }
        free_uivector(leftAllMembrIndx, 1, leftAllMembrSize);
        free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize);
        if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
        }
        else {
          free_uivector(leftRepMembrIndx, 1, bootMembrSize);
          free_uivector(rghtRepMembrIndx, 1, bootMembrSize);
        }
      }
      else {
        free_cvector(membershipIndicator, 1, RF_observationSize);
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  forkAndUpdate(%10d) failed.", treeID);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }  
    else {
      parent -> splitFlag = FALSE;
      free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
      parent -> permissibleSplit = NULL;
    }
  }  
  else {
    if (rootFlag) {
      if (!bootResult) {
        tnUpdateFlag = FALSE;
      }
    }
  }  
  if (tnUpdateFlag) {
    parent -> pseudoTerminal = TRUE;
    if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      bsUpdateFlag = TRUE;
    }
    if (!(RF_opt & OPT_IMPU_ONLY)) {
      if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
        getSplitDepth(parent, maximumDepth);
      }
    }
    parent -> orderedNodeID = ++RF_orderedLeafCount[treeID];
  }  
  if (bsUpdateFlag) {
    for (i = 1; i <= bootMembrSize; i++) {
      RF_bootMembershipIndex[treeID][++(*bootMembrIndxIter)] = bootMembrIndx[i];
      RF_bootMembershipFlag[treeID][bootMembrIndx[i]] = TRUE;
      RF_oobMembershipFlag[treeID][bootMembrIndx[i]]  = FALSE;
      RF_bootMembershipCount[treeID][bootMembrIndx[i]] ++;
    }
  }
  if (rootFlag | (RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    free_uivector(bootMembrIndx, 1, allMembrSize);
  }
  return bootResult;
}
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
                 uint   *maximumDepth) {
  char notTerminal;
  uint i;
  if (b != treeID[*offset]) {
    Rprintf("\nDiagnostic Trace of Tree Record:  \n");
    Rprintf("\n    treeID     nodeID     parmID       spltPT     mwcpSZ \n");
    Rprintf("%10d %10d %10d %12.4f %10d \n", treeID[*offset], nodeID[*offset], parmID[*offset], contPT[*offset], mwcpSZ[*offset]);
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid forest input record in tree:  %10d", b);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  parent -> depth = depth;
  parent -> left  = NULL;
  parent -> right = NULL;
  free_cvector(parent -> permissibleSplit, 1, parent -> xSize);
  parent -> permissibleSplit = NULL;
  parent -> splitFlag = FALSE;
  parent -> predictedOutcome = NA_REAL;
  parent -> nodeID = nodeID[*offset];
  parent -> splitParameter = parmID[*offset];
  if ((parent -> splitParameter) != 0) {
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      parent -> splitValueFactSize = mwcpSZ[*offset];
      parent -> splitValueFactPtr = uivector(1, mwcpSZ[*offset]);
      for (i = 1; i <= parent -> splitValueFactSize; i++) {
        (*mwcpPtr) ++;
        (parent -> splitValueFactPtr)[i] = **mwcpPtr;
      }
      parent -> splitValueCont = NA_REAL;
    }
    else {
      parent -> splitValueCont = contPT[*offset];
      parent -> splitValueFactSize = 0;
      parent -> splitValueFactPtr = NULL;
    }
  }
  else {
    parent -> splitValueCont     = NA_REAL;
    parent -> splitValueFactSize = 0;
    parent -> splitValueFactPtr  = NULL;
  }
  (*offset) ++;
  if ((parent -> splitParameter) != 0) {
    notTerminal = TRUE;
    parent -> left  = makeNode(parent -> xSize);
    setParent(parent ->  left, parent);
    restoreTree(b,
                parent -> left,
                offset,
                treeID,
                nodeID,
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                parent -> depth + 1,
                maximumDepth);
    parent -> right = makeNode(parent -> xSize);
    setParent(parent -> right, parent);
    restoreTree(b,
                parent -> right,
                offset,
                treeID,
                nodeID,
                parmID,
                contPT,
                mwcpSZ,
                mwcpPtr,
                parent -> depth + 1,
                maximumDepth);
  }
  else {
    notTerminal = FALSE;
  }
  if (!notTerminal) {
    parent -> pseudoTerminal = TRUE;
    RF_tNodeList[b][parent -> nodeID] = parent;
    if (!(RF_opt & OPT_IMPU_ONLY)) {
      if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
        getSplitDepth(parent, maximumDepth);
      }
    }
    parent -> orderedNodeID = ++RF_orderedLeafCount[b];
  }
  return notTerminal;
}
void saveTree(uint    b,
              Node   *parent,
              uint   *offset,
              uint   *treeID,
              uint   *nodeID,
              uint   *parmID,
              double *contPT,
              uint   *mwcpSZ,
              uint  **mwcpPtr) {
  uint i;
  treeID[*offset] = b;
  nodeID[*offset] = parent -> nodeID;
  parmID[*offset] = parent -> splitParameter;
  if ((parent -> splitParameter) != 0) {
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      mwcpSZ[*offset] = parent -> splitValueFactSize;
      for (i = 1; i <= mwcpSZ[*offset]; i++) {
        (*mwcpPtr) ++;
        **mwcpPtr = (parent -> splitValueFactPtr)[i];
      }
      contPT[*offset] = NA_REAL;
    }
    else {
      contPT[*offset] = parent -> splitValueCont;
      mwcpSZ[*offset] = 0;
    }
  }
  else {
    contPT[*offset] = NA_REAL;
    mwcpSZ[*offset] = 0;
  }
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveTree(b, parent ->  left, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr);
    saveTree(b, parent -> right, offset, treeID, nodeID, parmID, contPT, mwcpSZ, mwcpPtr);
  }
}
void freeTree(uint treeID, Node *parent, char rootFlag) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    freeTree(treeID, parent -> left, FALSE);
    freeTree(treeID, parent -> right, FALSE);
  }
  freeNode(parent);
}
void getSplitDepth(Node *parent, uint *maximumDepth) {
  Node *reversePtr;
  uint i;
  if (!(RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T))) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Call to calculate split depth without the option being active.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (parent -> depth > 0) {
    *maximumDepth = ((parent -> depth > *maximumDepth) ? parent -> depth : *maximumDepth);
    stackSplitDepth(parent, parent -> depth);
    reversePtr = parent;
    for (i = 1; i <= parent -> depth; i++) {
      if ((reversePtr -> parent) == NULL) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Reverse parsing of tree failed in restoreTree().");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      (parent -> splitDepth)[(parent -> depth) - i + 1] = (reversePtr -> parent) -> splitParameter;
      reversePtr = reversePtr -> parent;
    }
  }
}
void freeSplitDepth(uint treeID) {
  uint j;
  for (j = 1; j <= RF_tLeafCount[treeID]; j++) {
    unstackSplitDepth(RF_tNodeList[treeID][j]);
  }
}
void saveStatistics(char    mode,
                    uint    b,
                    Node   *parent,
                    uint   *offset,
                    double *spltST,
                    double *spltVR) {
  if (!(RF_opt & OPT_NODE_STAT)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Inconsistent call to saveStatistics().  The option is NOT active.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (mode == RF_GROW) {
    spltST[*offset] = parent -> splitStatistic;
  }
  else {
    if (RF_ptnCount == 0) {
      spltST[*offset] = parent -> variance;
    }
    else {
      spltST[*offset] = parent -> pseudoTerminal;
    }
  }
  (*offset) ++;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    saveStatistics(mode, b, parent ->  left, offset, spltST, spltVR);
    saveStatistics(mode, b, parent -> right, offset, spltST, spltVR);
  }
}
uint getMaximumDepth(Node *parent) {
  uint result, rLeft, rRight;
  result = parent -> depth;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    rLeft = getMaximumDepth(parent ->  left);
    rRight = getMaximumDepth(parent -> right);
    result = (rLeft > rRight) ? rLeft : rRight;
  }
  return result;
}
void getNodesAtDepth(Node *parent, uint tagDepth, Node **nodesAtDepth, uint *nadCount) {
  char recurseFlag;
  recurseFlag = TRUE;
  if (tagDepth == parent -> depth) {
    if ((parent -> splitParameter) != 0) {
      (*nadCount) ++;
      nodesAtDepth[*nadCount] = parent;
    }
    recurseFlag = FALSE;
  }
  else {
    if (((parent -> left) == NULL) && ((parent -> right) == NULL)) {
      recurseFlag = FALSE;
    }
  }
  if (recurseFlag) {
    getNodesAtDepth(parent ->  left, tagDepth, nodesAtDepth, nadCount);
    getNodesAtDepth(parent -> right, tagDepth, nodesAtDepth, nadCount);
  }
}
void getTreeInfo(uint treeID, Node *parent) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    getTreeInfo(treeID, parent ->  left);
    getTreeInfo(treeID, parent -> right);
  }
}
void getPTNodeList(Node    *parent,
                   Node   **list,
                   uint    *offset) {
  if (!(parent -> pseudoTerminal)) {
    getPTNodeList(parent ->  left, list, offset);
    getPTNodeList(parent -> right, list, offset);
  }
  else {
    (*offset) ++;
    list[*offset] = parent;
  }
}
