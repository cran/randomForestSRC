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


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include        "random.h"
#include        "impute.h"
#include     "bootstrap.h"
#include       "nodeOps.h"
#include      "treeUtil.h"
#include    "importance.h"
#include     "rfsrcUtil.h"
#include          "tree.h"
void acquireTree(uint mode, uint r, uint b) {
  Node  *rootPtr;
  uint **mwcpPtrPtr;
  uint  *mwcpPtr;
  uint   nodeOffset;
  char  termImputeFlag;
  char multipleImputeFlag;
  uint *allMembrIndx;
  uint *fallMembrIndx;
  uint *repMembrIndxImputed;
  uint *allMembrIndxImputed;
  uint *ngAllMembrIndxImputed;
  uint  repMembrSizeImputed;
  uint  allMembrSizeImputed;
  uint  ngAllMembrSizeImputed;
  Node *parent;
  uint bootMembrIndxIter;
  uint mRecordSize;
  uint *mRecordIndex;
  Node ***nodeMembershipPtr;
  char  result;
  Node  ***gNodeMembership;
  uint     obsSize;
  uint i, j;
  obsSize    = 0;  
  gNodeMembership = NULL;  
  multipleImputeFlag = FALSE;
  if (mode == RF_GROW) {
    if (r > 1) {
      multipleImputeFlag = TRUE;
    } 
  }
#ifdef SUPPORT_OPENMP
#endif
  rootPtr = makeNode(RF_xSize);  
  RF_tNodeMembership[b] = (Node **) vvector(1, RF_observationSize);
  RF_bootMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_bootMembershipFlag[b] = uivector(1, RF_observationSize);
  RF_oobMembershipFlag[b] = uivector(1, RF_observationSize);
  allMembrIndx = uivector(1, RF_observationSize);
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    obsSize = RF_observationSize; 
    gNodeMembership = RF_tNodeMembership;
  }
  else {
    if (mode == RF_PRED) {
      obsSize = RF_fobservationSize; 
      gNodeMembership = RF_ftNodeMembership;
    }
    else {
    }
  }
  if (mode == RF_PRED) {
    RF_ftNodeMembership[b] = (Node **) vvector(1, RF_fobservationSize);
  }
  if (mode == RF_PRED) {
    fallMembrIndx = uivector(1, RF_fobservationSize);
  }
  else {
    fallMembrIndx = NULL;
  }
  if (RF_ptnCount > 0) {
    RF_pNodeMembership[b] = (Node **) vvector(1, obsSize);
  }
  stackShadow(mode, b);
  if (mode == RF_GROW) {
    if (RF_nImpute > 1) {
      if (r > 1) {
        if (RF_mRecordSize > 0) {
          imputeUpdateShadow(RF_GROW, 
                             RF_response[b], 
                             RF_observation[b]);
        }
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(0, 
                                 NULL,
                                 RF_observationSize,
                                 RF_time[b],
                                 FALSE,
                                 TRUE,
                                 RF_masterTimeIndex[b]);
          }
        }
      }
    }
  }
  rootPtr -> parent = NULL;
  rootPtr -> nodeID = 1;
  RF_root[b] = rootPtr;
  RF_maxDepth[b] = 0;
  bootMembrIndxIter = 0;
  for (i=1; i <= RF_observationSize; i++) {
    allMembrIndx[i] = i;
    RF_tNodeMembership[b][i] = RF_root[b];
    RF_bootMembershipFlag[b][i] = FALSE;
    RF_oobMembershipFlag[b][i]  = TRUE;
  }
  if (RF_ptnCount > 0) {
    for (i = 1; i <= obsSize; i++) {
      RF_pNodeMembership[b][i] = gNodeMembership[b][i]; 
    }
  }
  RF_orderedLeafCount[b] = 0;
  if (mode == RF_GROW) {
    RF_tLeafCount[b] = 0;
    result = growTree (TRUE, 
                       multipleImputeFlag, 
                       b, 
                       rootPtr, 
                       NULL,
                       0, 
                       allMembrIndx, 
                       RF_observationSize, 
                       0, 
                       RF_maxDepth + b,
                       & bootMembrIndxIter);
    stackNodeList(b);
    initNodeList(b);
  }  
  else {
    if (mode != RF_REST) {
      for (i=1; i <= RF_fobservationSize; i++) {
        fallMembrIndx[i] = i;
        RF_ftNodeMembership[b][i] = RF_root[b];
      }
    }
    nodeOffset = 1;
    mwcpPtr = RF_mwcpPT_;
    mwcpPtrPtr = & mwcpPtr;
    for (j = 1; j < b; j++) {
      mwcpPtr += RF_mwcpCount[j];
      nodeOffset += RF_nodeCount[j];
    }
    stackNodeList(b);
    restoreTree(b,
                rootPtr,
                & nodeOffset,
                RF_treeID_,
                RF_nodeID_,
                RF_parmID_,
                RF_contPT_,
                RF_mwcpSZ_,
                mwcpPtrPtr,
                0,
                RF_maxDepth + b);
    result = restoreNodeMembership(mode, 
                                   TRUE,
                                   b, 
                                   rootPtr, 
                                   NULL,
                                   0,
                                   allMembrIndx, 
                                   RF_observationSize, 
                                   fallMembrIndx,
                                   RF_fobservationSize,
                                   & bootMembrIndxIter);
  }  
  if (result) {
    if (RF_opt & OPT_MISS) {
      switch (mode) {
      case RF_PRED:
        mRecordSize = RF_fmRecordSize;
        mRecordIndex = RF_fmRecordIndex;
        nodeMembershipPtr = RF_ftNodeMembership;
        break;
      default:
        mRecordSize = RF_mRecordSize;
        mRecordIndex = RF_mRecordIndex;
        nodeMembershipPtr = RF_tNodeMembership;
        break;
      } 
      RF_mTermList[b] = (Terminal **) vvector(1, RF_tLeafCount[b]);
      RF_mTermMembership[b] = (Terminal **) vvector(1, mRecordSize);
      for (j = 1; j <= RF_tLeafCount[b]; j++) {
        RF_mTermList[b][j] = makeTerminal();
        RF_mTermList[b][j] -> nodeID = RF_tNodeList[b][j] -> nodeID;
        RF_tNodeList[b][j] -> mate = RF_mTermList[b][j];
      }
      for (j = 1; j <= mRecordSize; j++) {
        RF_mTermMembership[b][j] = nodeMembershipPtr[b][mRecordIndex[j]] -> mate;
      }
    }  
    if (r == 1) {
      if (RF_mRecordSize > 0) {
        repMembrIndxImputed = uivector(1, RF_observationSize);
        allMembrIndxImputed = uivector(1, RF_observationSize);
        for (j = 1; j <= RF_tLeafCount[b]; j++) {
          unstackNodeLMPIndex(RF_tNodeList[b][j]);
          unstackNodeLMRIndex(RF_tNodeList[b][j]);
          parent = RF_tNodeList[b][j];
          getRawNodeSize(RF_GROW, b, parent, repMembrIndxImputed, & repMembrSizeImputed, allMembrIndxImputed,  & allMembrSizeImputed);
          imputeNode(RF_GROW,
                     TRUE,
                     TRUE,
                     b,
                     parent,
                     repMembrIndxImputed, 
                     repMembrSizeImputed, 
                     allMembrIndxImputed,
                     allMembrSizeImputed);
        }  
        free_uivector(repMembrIndxImputed, 1, RF_observationSize);
        free_uivector(allMembrIndxImputed, 1, RF_observationSize);
        if (mode != RF_PRED) {
          imputeUpdateSummary(mode, b);
        }
      }  
      if (mode == RF_PRED) {
        if (RF_fmRecordSize > 0) {
          repMembrIndxImputed = uivector(1, RF_observationSize);
          ngAllMembrIndxImputed = uivector(1, RF_fobservationSize);
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            unstackNodeFLMPIndex(RF_tNodeList[b][j]);
            unstackNodeFLMRIndex(RF_tNodeList[b][j]);
            parent = RF_tNodeList[b][j];
            getRawNodeSize(RF_PRED, b, parent, repMembrIndxImputed, & repMembrSizeImputed, ngAllMembrIndxImputed,  & ngAllMembrSizeImputed);
            imputeNode(RF_PRED,
                       TRUE,
                       FALSE,
                       b,
                       parent,
                       repMembrIndxImputed, 
                       repMembrSizeImputed, 
                       ngAllMembrIndxImputed,
                       ngAllMembrSizeImputed);
          }  
          free_uivector(repMembrIndxImputed, 1, RF_observationSize);
          free_uivector(ngAllMembrIndxImputed, 1, RF_fobservationSize);
          imputeUpdateSummary(RF_PRED, b);
        }  
      }  
    }  
    if (mode == RF_GROW) {
      termImputeFlag = FALSE;
      if (RF_opt & OPT_IMPU_ONLY) {
        if ((RF_nImpute > 1) && (r > 1) && (r <= RF_nImpute) ) {
          termImputeFlag = TRUE;
        }
      }
      else {
        if ((RF_nImpute > 1) && (r > 1) && (r < RF_nImpute) ) {
          termImputeFlag = TRUE;
        }
      }
      if (termImputeFlag) {
        if (RF_mRecordSize > 0) {
          repMembrIndxImputed = uivector(1, RF_observationSize);
          allMembrIndxImputed = uivector(1, RF_observationSize);
          for (j = 1; j <= RF_tLeafCount[b]; j++) {
            parent = RF_tNodeList[b][j];
            getRawNodeSize(RF_GROW, b, parent, repMembrIndxImputed, & repMembrSizeImputed, allMembrIndxImputed,  & allMembrSizeImputed);
            imputeNode(RF_GROW,
                       TRUE,
                       FALSE,
                       b,
                       parent,
                       repMembrIndxImputed, 
                       repMembrSizeImputed, 
                       allMembrIndxImputed,
                       allMembrSizeImputed);
          }  
          imputeUpdateSummary(RF_GROW, b);
          free_uivector(repMembrIndxImputed, 1, RF_observationSize);
          free_uivector(allMembrIndxImputed, 1, RF_observationSize);
        }  
      }  
    }  
  }  
  if (result) {
    if (RF_ptnCount > 0) {
      for (i = 1; i <= obsSize; i++) {
        RF_pNodeMembership[b][i] = gNodeMembership[b][i]; 
      }
      RF_pLeafCount[b] = pruneTree(mode, b, RF_ptnCount);
      RF_pNodeList[b] = (Node **) vvector(1, RF_pLeafCount[b] + 1);
      i = 0;
      getPTNodeList(RF_root[b], RF_pNodeList[b], &i);
      free_vvector(RF_pNodeList[b], 1, RF_pLeafCount[b] + 1);
    }
    RF_oobSize[b] = 0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_bootMembershipFlag[b][i] == FALSE) {
        RF_oobSize[b] ++;
      }
    }
    if ((mode == RF_GROW) || (mode == RF_REST)) {
      if (RF_mRecordSize > 0) {
        for (i = 1; i <= RF_mRecordSize; i++) {
          if (RF_bootMembershipFlag[b][RF_mRecordIndex[i]] == TRUE) {
            RF_dmRecordBootFlag[b][i] = TRUE;
          }
          else {
            RF_dmRecordBootFlag[b][i] = FALSE;
          }
        }
      }  
    }  
  }  
  if (result) {
    if (RF_opt & OPT_MEMB) {
      for (i=1; i <= obsSize; i++) {
        RF_tNodeMembershipIndexPtr[b][i] = gNodeMembership[b][i] -> nodeID;
        if (RF_ptnCount > 0) {
          RF_pNodeMembershipIndexPtr[b][i] = RF_pNodeMembership[b][i] -> nodeID;
        }
      }
    }
    if (r == RF_nImpute) {
      if (!(RF_opt & OPT_IMPU_ONLY)) {
        if ((RF_opt & OPT_PERF) | 
            (RF_opt & OPT_PERF_CALB) |
            (RF_opt & OPT_OENS) |
            (RF_opt & OPT_FENS)) {
          updateTerminalNodeOutcomes(b);
        }
      }
      if (RF_opt & OPT_VIMP) {
        uint vimpCount;
        if (RF_opt & OPT_VIMP_JOIN) {
          vimpCount = 1;
        }
        else {
          vimpCount = RF_intrPredictorSize;
        }
        for (uint intrIndex = 1; intrIndex <= vimpCount; intrIndex++) {
          uint pp;
          if (!(RF_opt & OPT_VIMP_JOIN)) {
            pp = RF_intrPredictor[intrIndex];
          }
          else {
            pp = 0;
          }
          stackVimpMembership(mode, & RF_vimpMembership[intrIndex][b]);
          getVimpMembership(mode, b, RF_vimpMembership[intrIndex][b], pp);
        }
      }
    }  
  }  
  unstackShadow(mode, b, FALSE, TRUE);
  free_uivector(allMembrIndx, 1, RF_observationSize);
  if (mode == RF_PRED) {
    free_uivector(fallMembrIndx, 1, RF_fobservationSize);
  }
}
void updateProximity(uint mode, uint treeID) {
  Node ***nodeMembership;
  uint    obsSize;
  char flag, selectionFlag;
  uint i, j, k;
  if (RF_tLeafCount[treeID] > 0) {
    if ((mode == RF_GROW) || (mode == RF_REST)) {
      nodeMembership = RF_tNodeMembership;
      obsSize = RF_observationSize;
    }
    else {
      nodeMembership = RF_ftNodeMembership;
      obsSize = RF_fobservationSize;
    }
    if (RF_opt & (OPT_PROX | OPT_PROX_TYPE)) {
      flag = NEITHER;
      if(!(RF_opt & OPT_PROX) && (RF_opt & OPT_PROX_TYPE)) {
        flag = TRUE;
      }
      else if((RF_opt & OPT_PROX) && !(RF_opt & OPT_PROX_TYPE)) {
        flag = FALSE;
      }
      else if((RF_opt & OPT_PROX) && (RF_opt & OPT_PROX_TYPE)) {
        flag = ACTIVE;
      }
      k = 0;
      for (i = 1; i <= obsSize; i++) {
        k += i - 1;
        for (j = 1; j <= i; j++) {
          selectionFlag = FALSE;
          if (flag == TRUE) {
            if ((RF_bootMembershipFlag[treeID][i] == TRUE) && (RF_bootMembershipFlag[treeID][j] == TRUE)) {
              selectionFlag = TRUE;
            }
          }
          else if (flag == FALSE) {
            if ((RF_bootMembershipFlag[treeID][i] == FALSE) && (RF_bootMembershipFlag[treeID][j] == FALSE)) {
              selectionFlag = TRUE;
            }
          }
          else if (flag == ACTIVE) {
            selectionFlag = TRUE;
          }
          if (selectionFlag) {
            if ( (nodeMembership[treeID][i] -> nodeID) == (nodeMembership[treeID][j] -> nodeID) ) {
              RF_proximity_[k + j] ++;
            }
          }
        }
      }
    }
  }
}
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth) {
  Node  *parent;
  double *localSplitDepth;
  uint index;
  uint i, j, k;
  if (RF_tLeafCount[treeID] > 0) {
    index = 0;  
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      if (RF_opt & OPT_SPLDPTH_F) {
        index = 1;
      }
      else {
        index = treeID;
      }
    }
    else {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Illegal updateSplitDepth() call.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
    localSplitDepth = dvector(1, RF_xSize);
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= RF_xSize; j++) {
        localSplitDepth[j] = NA_REAL;
      }
      parent = RF_tNodeMembership[treeID][i];
      for (k = 1; k <= parent -> depth; k++) {
        if (ISNA(localSplitDepth[(parent -> splitDepth)[k]])) {
          localSplitDepth[(parent -> splitDepth)[k]] = (double) k;
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        if (ISNA(localSplitDepth[j])) {
          localSplitDepth[j] = (double) maxDepth + 1;
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        RF_splitDepthPtr[index][j][i] += localSplitDepth[j];
      }
    }
    free_dvector(localSplitDepth, 1, RF_xSize);
    freeSplitDepth(treeID);
  }
}
char pruneBranch(uint mode, uint treeID, Node **nodesAtDepth, uint nadCount, uint ptnTarget, uint ptnCurrent) {
  char pruneFlag;
  uint obsSize;
  uint i, j;
  pruneFlag = TRUE;
  double *varianceAtDepth =  dvector(1, nadCount);
  uint   *vadSortedIndex  = uivector(1, nadCount);
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    obsSize = RF_observationSize; 
  }
  else {
    obsSize = RF_fobservationSize; 
  }
  for (i = 1; i <= nadCount; i++) {
    varianceAtDepth[i] = nodesAtDepth[i] -> variance;
  }
  indexx(nadCount, varianceAtDepth, vadSortedIndex);
  j = nadCount;
  while ((j >= 1) && pruneFlag) {
    nodesAtDepth[vadSortedIndex[j]] -> pseudoTerminal = TRUE;
    (nodesAtDepth[vadSortedIndex[j]] -> left)  -> pseudoTerminal = FALSE;
    (nodesAtDepth[vadSortedIndex[j]] -> right) -> pseudoTerminal = FALSE;
    for (i = 1; i <= obsSize; i++) {
      if ( (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> left) ||
           (RF_pNodeMembership[treeID][i] == nodesAtDepth[vadSortedIndex[j]] -> right)) {
        RF_pNodeMembership[treeID][i] = nodesAtDepth[vadSortedIndex[j]];
      }
    }
    j --;
    ptnCurrent --;
    if (ptnCurrent <= ptnTarget) {
      pruneFlag = FALSE;
    }
  }
  free_dvector(varianceAtDepth, 1, nadCount);
  free_uivector(vadSortedIndex, 1, nadCount);
  return pruneFlag;
}
uint pruneTree(uint mode, uint treeID, uint ptnTarget) {
  Node **nodesAtDepth;
  uint   ptnCurrent;
  uint   nadCount;
  uint   tagDepth;
  char   pruneFlag;
  uint   i;
  if (ptnTarget < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Illegal target PTN count in pruneTree():  %10d", ptnTarget);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (RF_tLeafCount[treeID] == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Illegal call to pruneTree() on a rejected tree:  %10d", treeID);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  nodesAtDepth = (Node **) vvector(1, RF_tLeafCount[treeID]);
  ptnCurrent = RF_tLeafCount[treeID];
  tagDepth = getMaximumDepth(RF_root[treeID]) - 1;
  pruneFlag = (ptnCurrent > ptnTarget) && (tagDepth > 0);
  while (pruneFlag) {
    for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
      nodesAtDepth[i] = NULL;
    }
    nadCount = 0;
    getNodesAtDepth(RF_root[treeID], tagDepth, nodesAtDepth, &nadCount);
    pruneFlag = pruneBranch(mode, treeID, nodesAtDepth, nadCount, ptnTarget, ptnCurrent);
    if(pruneFlag) {
      ptnCurrent -= nadCount;
      tagDepth --;
    }
    else {
      ptnCurrent = ptnTarget;
    }
  }
  free_vvector(nodesAtDepth, 1, RF_tLeafCount[treeID]);
  return ptnCurrent;
}
void unstackAuxiliary(uint mode, uint b) {
  uint obsSize;
  obsSize = 0;  
  free_vvector(RF_tNodeMembership[b], 1, RF_observationSize); 
  free_uivector(RF_bootMembershipIndex[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipFlag[b], 1, RF_observationSize);
  free_uivector(RF_oobMembershipFlag[b], 1, RF_observationSize);
  free_vvector(RF_tNodeList[b], 1, RF_tLeafCount[b] + 1);
  if (mode == RF_PRED) {
    free_vvector((Node **) RF_ftNodeMembership[b],  1, RF_fobservationSize);  
  }
  if (RF_ptnCount > 0) {
    if ((mode == RF_GROW) || (mode == RF_REST)) {
      obsSize = RF_observationSize; 
    }
    else {
      if (mode == RF_PRED) {
        obsSize = RF_fobservationSize; 
      }
      else {
      }
    }
    free_vvector(RF_pNodeMembership[b], 1, obsSize);
  }
}
void stackNodeList(uint treeID) {
  RF_tNodeList[treeID] = (Node **) vvector(1, RF_tLeafCount[treeID] + 1);
}
void initNodeList(uint treeID) {
  uint j;
  for (j = 1; j <= RF_tLeafCount[treeID]; j++) {
    RF_tNodeList[treeID][j] = getTerminalNode(treeID, j);
  }
}
