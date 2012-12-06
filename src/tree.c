////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.2
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
////    email:  kogalurshear@gmail.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************


#ifdef SUPPORT_OPENMP
#include           <omp.h>
#endif
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
  uint *maxDepthPtr;
  uint  maxDepth;
  uint **mwcpPtrPtr;
  uint  *mwcpPtr;
  uint   nodeOffset;
  double **responsePtr;
  double **predictorPtr;
  char  updateFlag;  
  char multipleImputeFlag;
  uint *allMembrIndx;
  uint *fallMembrIndx;
  uint *repMembrIndxImputed;
  uint *allMembrIndxImputed;
  uint  repMembrSizeImputed;
  uint  allMembrSizeImputed;
  Node *parent;
  uint bootMembrIndxIter;
  char  result;
  uint thisSerialTreeCount;
  Node  ***nodeMembershipPtr;
  uint     obsSize;
  uint i, j;
  obsSize    = 0;  
  responsePtr = NULL;  
  predictorPtr = NULL;  
  nodeMembershipPtr = NULL;  
  updateFlag  = FALSE;  
  thisSerialTreeCount = 0;  
  multipleImputeFlag = FALSE;
  if (mode == RF_GROW) {
    if (r > 1) {
      multipleImputeFlag = TRUE;
    } 
  }
#ifdef SUPPORT_OPENMP
#endif
  rootPtr = makeNode(RF_xSize);  
  RF_nodeMembership[b] = (Node **) vvector(1, RF_observationSize);
  RF_bootMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_bootMembershipFlag[b] = uivector(1, RF_observationSize);
  RF_oobMembershipFlag[b] = uivector(1, RF_observationSize);
  allMembrIndx = uivector(1, RF_observationSize);
  if (mode == RF_PRED) {
    RF_fnodeMembership[b] = (Node **) vvector(1, RF_fobservationSize);
    fallMembrIndx = uivector(1, RF_fobservationSize);
  }
  else {
    fallMembrIndx = NULL;
  }
  stackShadow(mode, b);
  if (mode == RF_GROW) {
    if (RF_imputeSize > 1) {
      if (r > 1) {
        if (r == 2) {
          if (RF_mRecordSize > 0) {
            imputeUpdateShadow(RF_GROW, 
                               FALSE, 
                               RF_response[b], 
                               RF_observation[b]);
          }
        }  
        if (r > 2) {
          if (RF_mRecordSize > 0) {
            imputeUpdateShadow(RF_GROW, 
                               ACTIVE, 
                               RF_response[b], 
                               RF_observation[b]);
          }
        }
        if (RF_timeIndex > 0) {
          if (RF_mTimeFlag == TRUE) {
            updateTimeIndexArray(0, 
                                 NULL,
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
  rootPtr -> leafCount = 1;
  RF_root[b] = rootPtr;
  maxDepth = 0;
  maxDepthPtr = &maxDepth;
  bootMembrIndxIter = 0;
  for (i=1; i <= RF_observationSize; i++) {
    allMembrIndx[i] = i;
    RF_nodeMembership[b][i] = RF_root[b];
    RF_bootMembershipFlag[b][i] = FALSE;
    RF_oobMembershipFlag[b][i]  = TRUE;
  }
  if (mode == RF_GROW) {
    RF_leafCount[b] = 0;
    result = growTree (TRUE, 
                       multipleImputeFlag, 
                       b, 
                       rootPtr, 
                       NULL,
                       0, 
                       allMembrIndx, 
                       RF_observationSize, 
                       0, 
                       maxDepthPtr,
                       & bootMembrIndxIter);
    RF_terminalNode[b] = (Node **) vvector(1, RF_leafCount[b] + 1);
    if ((RF_opt & OPT_MISS) | (RF_opt & OPT_OMIS)) {
      RF_mTerminalInfo[b] = (Terminal **) vvector(1, RF_leafCount[b] + 1);
      for (j = 1; j <= RF_leafCount[b]; j++) {
        RF_mTerminalInfo[b][j] = makeTerminal();
      }
    }
    for (j = 1; j <= RF_leafCount[b]; j++) {
      RF_terminalNode[b][j] = getTerminalNode(b, j);
    }
    if (result) {
      if ((RF_imputeSize > 1) && (r > 1) && (r < RF_imputeSize) ) {
        if (RF_mRecordSize > 0) {
          repMembrIndxImputed = uivector(1, RF_observationSize);
          allMembrIndxImputed = uivector(1, RF_observationSize);
          for (j = 1; j <= RF_leafCount[b]; j++) {
            parent = RF_terminalNode[b][j];
            getRawNodeSize(RF_GROW, b, parent, repMembrIndxImputed, & repMembrSizeImputed, allMembrIndxImputed,  & allMembrSizeImputed);
            imputeNode(RF_GROW,
                       FALSE,
                       b,
                       parent,
                       repMembrIndxImputed, 
                       repMembrSizeImputed, 
                       allMembrIndxImputed,
                       allMembrSizeImputed);
          }  
          free_uivector(repMembrIndxImputed, 1, RF_observationSize);
          free_uivector(allMembrIndxImputed, 1, RF_observationSize);
        }  
      }  
    }
  }  
  else {
    if (mode != RF_REST) {
      for (i=1; i <= RF_fobservationSize; i++) {
        fallMembrIndx[i] = i;
        RF_fnodeMembership[b][i] = RF_root[b];
      }
    }
    nodeOffset = 1;
    mwcpPtr = RF_mwcpPT_;
    mwcpPtrPtr = & mwcpPtr;
    for (j = 1; j < b; j++) {
      mwcpPtr += RF_mwcpCount[j];
      nodeOffset += RF_nodeCount[j];
    }
    RF_terminalNode[b] = (Node **) vvector(1, RF_leafCount[b] + 1);
    if ((RF_opt & OPT_MISS) | (RF_opt & OPT_OMIS)) {
      RF_mTerminalInfo[b] = (Terminal **) vvector(1, RF_leafCount[b] + 1);
    }
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
                maxDepthPtr);
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
  if (r == RF_imputeSize) {
    if (RF_opt & OPT_PROX) {
#ifdef SUPPORT_OPENMP
#pragma omp critical (RF_update_proximity_)
#endif
      {
        if (result) {
          updateProximity(mode, b);
        }
      }
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
#ifdef SUPPORT_OPENMP
#pragma omp critical (RF_update_split_depth_)
#endif
      {
        if (result) {
          updateSplitDepth(b, rootPtr, maxDepth);
        }
      }
      freeSplitDepth(b);
    }
  }  
#ifdef SUPPORT_OPENMP
#pragma omp critical (RF_update_serial_tree_count_)
#endif
  {    
    RF_serialTreeIndex[++RF_serialTreeCount] = b;
    thisSerialTreeCount = RF_serialTreeCount;
    if (result) {
      updateFlag = FALSE;
      switch (mode) {
      case RF_PRED:
        if (RF_fmRecordSize > 0) {
          updateFlag = TRUE;
          if(RF_frSize > 0) {
            responsePtr = RF_fresponse[b];
          }
          else {
            responsePtr = NULL;
          }
          predictorPtr = RF_fobservation[b];
        }
        break;
      default:
        if ((r == 1) || (r < RF_imputeSize)) {
          if (RF_mRecordSize > 0) {
            updateFlag = TRUE;
            responsePtr = RF_response[b];
            predictorPtr = RF_observation[b];
          }
        }
        break;
      }
      if (updateFlag == TRUE) {
        imputeUpdateSummary(mode, 
                            responsePtr, 
                            predictorPtr,
                            b);
      }
    }
  }
  if (result) {
    if (r == RF_imputeSize) {
      updateFlag = TRUE;
      if (mode == RF_GROW) {
        if (RF_opt & OPT_IMPU_ONLY) {
          updateFlag = FALSE;
        }
      }
      if (updateFlag) {
        updateEnsembleCalculations(multipleImputeFlag,
                                   mode,
                                   rootPtr,
                                   b,
                                   thisSerialTreeCount);
      }
      if (RF_opt & OPT_VUSE) {
        getVariablesUsed(rootPtr, RF_varUsedPtr[b]);
      }
    }
  }  
  if (RF_opt & OPT_MEMB) {
    if ((mode == RF_GROW) || (mode == RF_REST)) {
      obsSize = RF_observationSize; 
      nodeMembershipPtr = RF_nodeMembership;
    }
    else {
      if (mode == RF_PRED) {
        obsSize = RF_fobservationSize; 
        nodeMembershipPtr = RF_fnodeMembership;
      }
      else {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Illegal mode for SEXP membership output:  %10d", mode);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    for (i=1; i <= obsSize; i++) {
      RF_terminalNodeMembershipPtr[b][i] = nodeMembershipPtr[b][i] -> leafCount;
    }
  }
  if (!(RF_opt & OPT_TREE) || !(r == RF_imputeSize)) {
    freeTree(b, rootPtr, TRUE);
  }
  free_vvector(RF_nodeMembership[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipIndex[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipFlag[b], 1, RF_observationSize);
  free_uivector(RF_oobMembershipFlag[b], 1, RF_observationSize);
  free_uivector(allMembrIndx, 1, RF_observationSize);
  free_vvector(RF_terminalNode[b], 1, RF_leafCount[b] + 1);
  if (mode == RF_PRED) {
    free_vvector((Node **) RF_fnodeMembership[b],  1, RF_fobservationSize);
    free_uivector(fallMembrIndx, 1, RF_fobservationSize);
  }
  unstackShadow(mode, b);
}
void updateProximity(uint mode, uint treeID) {
  Node ***nodeMembership;
  uint    obsSize;
  uint i, j, k;
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    nodeMembership = RF_nodeMembership;
    obsSize = RF_observationSize;
  }
  else {
    nodeMembership = RF_fnodeMembership;
    obsSize = RF_fobservationSize;
  }
  k = 0;
  for (i = 1; i <= obsSize; i++) {
    k += i - 1;
    for (j = 1; j <= i; j++) {
      if ( (nodeMembership[treeID][i] -> leafCount) == (nodeMembership[treeID][j] -> leafCount) ) {
        RF_proximity_[k + j] ++;
      }
    }
  }
}
void updateSplitDepth(uint treeID, Node *rootPtr, uint maxDepth) {
  Node  *parent;
  double *localSplitDepth;
  uint index;
  uint i, j, k;
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
    parent = RF_nodeMembership[treeID][i];
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
}
