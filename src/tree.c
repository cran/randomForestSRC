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


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include        "random.h"
#include        "impute.h"
#include     "bootstrap.h"
#include       "nodeOps.h"
#include       "termOps.h"
#include      "treeUtil.h"
#include     "splitUspv.h"
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
  char  result;
  Node     ***gNodeMembership;
  Terminal ***gTermMembership;
  uint     obsSize;
  uint i, j;
  obsSize    = 0;  
  gNodeMembership = NULL;  
  gTermMembership = NULL;  
  multipleImputeFlag = FALSE;
  if (mode == RF_GROW) {
    if (r > 1) {
      multipleImputeFlag = TRUE;
    }
  }
#ifdef SUPPORT_OPENMP
#endif
  rootPtr = makeNode(RF_xSize);
  RF_tNodeMembership[b] = (Node **) new_vvector(1, RF_observationSize, NRUTIL_NPTR);
  RF_bootMembershipIndex[b] = uivector(1, RF_observationSize);
  RF_bootMembershipFlag[b] = cvector(1, RF_observationSize);
  RF_bootMembershipCount[b] = uivector(1, RF_observationSize);
  RF_oobMembershipFlag[b] = cvector(1, RF_observationSize);
  allMembrIndx = uivector(1, RF_observationSize);
  if (mode == RF_PRED) {
    RF_ftNodeMembership[b] = (Node **) new_vvector(1, RF_fobservationSize, NRUTIL_NPTR);
  }
  if (mode != RF_PRED) {
    obsSize = RF_observationSize;
    gNodeMembership = RF_tNodeMembership;
    gTermMembership = RF_tTermMembership;
  }
  else {
    obsSize = RF_fobservationSize;
    gNodeMembership = RF_ftNodeMembership;
    gTermMembership = RF_ftTermMembership;
  }
  if (mode == RF_PRED) {
    fallMembrIndx = uivector(1, RF_fobservationSize);
  }
  else {
    fallMembrIndx = NULL;
  }
  if (RF_ptnCount > 0) {
    RF_pNodeMembership[b] = (Node **) new_vvector(1, obsSize, NRUTIL_NPTR);
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
                                 FALSE,
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
    RF_bootMembershipFlag[b][i]  = FALSE;
    RF_bootMembershipCount[b][i] = 0;
    RF_oobMembershipFlag[b][i]   = TRUE;
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
#ifdef SUPPORT_OPENMP
#else
    if (getTraceFlag(b) & SUMM_USR_TRACE) {
        Rprintf("\nTree Complete:  (MII:  %10d, ID:  %10d)", r, b);
    }
#endif
    if (result) {
      stackNodeList(b);
      initNodeList(b);
      stackTermList(b);
      initTermList(b);
    }
  }  
  else {
    if (mode == RF_PRED) {
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
    if (RF_tLeafCount[b] > 0) {
      stackNodeList(b);
      stackTermList(b);
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
    if (result) {
      if ((RF_optHigh & OPT_TERM) && !(RF_opt & (OPT_BOOT_NODE | OPT_BOOT_NONE)) && RF_fmRecordSize == 0) {
        restoreNodeMembershipGrow(b);
      }
    }
  }  
  if (result) {
    stackAndInitTermMembership(mode, b);
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
      RF_pLeafCount[b] = pruneTree(obsSize, b, RF_ptnCount);
      RF_pNodeList[b] = (Node **) new_vvector(1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
      i = 0;
      getPTNodeList(RF_root[b], RF_pNodeList[b], &i);
      free_new_vvector(RF_pNodeList[b], 1, RF_pLeafCount[b] + 1, NRUTIL_NPTR);
    }
    RF_oobSize[b] = 0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_bootMembershipFlag[b][i] == FALSE) {
        RF_oobSize[b] ++;
      }
    }
    if (mode != RF_PRED) {
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
    if (mode == RF_REST) {
      if(RF_sobservationSize > 0) {
        RF_soobSize[b] = 0;
        for (i = 1; i <= RF_sobservationSize; i++) {
          if (RF_bootMembershipFlag[b][RF_sobservationIndv[i]] == FALSE) {
            RF_soobSize[b] ++;
          }
        }
      }
    }
  }  
  if (result) {
    if (RF_opt & OPT_MEMB) {
      for (i=1; i <= obsSize; i++) {
        RF_tTermMembershipIndexPtr[b][i] = gTermMembership[b][i] -> nodeID;
        if (RF_ptnCount > 0) {
          RF_pNodeMembershipIndexPtr[b][i] = RF_pNodeMembership[b][i] -> nodeID;
        }
      }
    }
    if (r == RF_nImpute) {
      if ((RF_opt & OPT_PERF) |
          (RF_opt & OPT_PERF_CALB) |
          (RF_opt & OPT_OENS) |
          (RF_opt & OPT_FENS)) {
        char multipleImputeFlag;
        multipleImputeFlag = FALSE;
        if (mode == RF_GROW) {
          if (r > 1) {
            multipleImputeFlag = TRUE;
          }
        }
        updateEnsembleCalculations(multipleImputeFlag, mode, b);
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
          updateVimpCalculations(mode, b, intrIndex, RF_vimpMembership[intrIndex][b]);
          unstackVimpMembership(mode, RF_vimpMembership[intrIndex][b]);
        }
      }
    }  
  }  
  else {
    if (RF_opt & OPT_PROX) {
      anticipateProximity(mode, b);
    }
  }
  unstackShadow(mode, b, FALSE, TRUE);
  free_uivector(allMembrIndx, 1, RF_observationSize);
  if (mode == RF_PRED) {
    free_uivector(fallMembrIndx, 1, RF_fobservationSize);
  }
}
void getWeight(uint mode) {
  Terminal ***gTermMembership;
  uint    obsSize;
  char    flag;
  uint   *weightDenom;
  uint i,j;
  int b;
  if (mode != RF_PRED) {
    gTermMembership = RF_tTermMembership;
    obsSize         = RF_observationSize;
  }
  else {
    gTermMembership = RF_ftTermMembership;
    obsSize         = RF_fobservationSize;
  }
  if(RF_optHigh & OPT_WGHT_TYP2) {
    flag = ACTIVE;
  }
  else {
    if(!(RF_optHigh & OPT_WGHT_TYP1)  && !(RF_optHigh & OPT_WGHT_TYP2)) {
      flag = TRUE;
    }
    else if(RF_optHigh & OPT_WGHT_TYP1) {
      flag = FALSE;
    }
    else {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Illegal getWeight() call.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  if (!((RF_opt & OPT_OENS) | (RF_opt & OPT_FENS))) {
    if (RF_numThreads > 0) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (b = 1; b <= RF_forestSize; b++) {
        getMemberCountOnly(b);
      }
    }
    else {
      for (b = 1; b <= RF_forestSize; b++) {
        getMemberCountOnly(b);
      }
    }
  }
  if (RF_numThreads > 0) {
    for (b = 1; b <= RF_forestSize; b++) {
      updateWeight(b, flag, obsSize, gTermMembership);
    }
  }
  else {
    for (b = 1; b <= RF_forestSize; b++) {
      updateWeight(b, flag, obsSize, gTermMembership);
    }
  }
  if (flag != ACTIVE) {
    if (flag == TRUE) {
      for (i = 1; i <= RF_observationSize; i++) {
        for (j = 1; j <= obsSize; j++) {
          RF_weightPtr[j][i] = RF_weightPtr[j][i] / RF_validTreeCount;
        }
      }
    }
    else {
      weightDenom = uivector(1, RF_observationSize);
      for (i = 1; i <= RF_observationSize; i++) {
        weightDenom[i] = 0;
      }
      for (b = 1; b <= RF_forestSize; b++) {
        for (i = 1; i <= RF_observationSize; i++) {
          if(RF_bootMembershipCount[b][i] == 0) {
            weightDenom[i] ++;
          }
        }
      }
      for (i = 1; i <= RF_observationSize; i++) {
        for (j = 1; j <= RF_observationSize; j++) {
          if (weightDenom[j] > 0) {
            RF_weightPtr[j][i] = RF_weightPtr[j][i] / (double) weightDenom[j];
          }
          else {
            RF_weightPtr[j][i] = NA_REAL;
          }
        }
      }
      free_uivector (weightDenom, 1, RF_observationSize);
    }
  }
  else {
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= obsSize; j++) {
        RF_weightPtr[j][i] = RF_weightPtr[j][i] / RF_validTreeCount;
      }
    }
  }
}
void updateWeight(uint b, char flag, uint obsSize, Terminal ***gTermMembership) {
  Terminal *parent;
  uint leaf;
  uint i, j;
  if (RF_tLeafCount[b] > 0) {
    if (flag != ACTIVE) {
      if (flag == TRUE) {
        for (i = 1; i <= RF_observationSize; i++) {
          for (j = 1; j <= obsSize; j++) {
            if ( RF_tTermMembership[b][i] == gTermMembership[b][j] ) {
                RF_weightPtr[j][i] +=  (double) RF_bootMembershipCount[b][i] / (double) (gTermMembership[b][j] -> membrCount);
            }
          }
        }
      }
      else {
        for (i = 1; i <= RF_observationSize; i++) {
          for (j = 1; j <= RF_observationSize; j++) {
            if (RF_bootMembershipCount[b][j] == 0) {
              if ( RF_tTermMembership[b][i] == RF_tTermMembership[b][j] ) {
                RF_weightPtr[j][i] +=  (double) RF_bootMembershipCount[b][i] / (double) (RF_tTermMembership[b][j] -> membrCount);
              }
            }
          }
        }
      }
    }
    else {
      for (leaf = 1; leaf <= RF_tLeafCount[b]; leaf++) {
        parent = RF_tTermList[b][leaf];
        parent -> weight = 0.0;
        for (i = 1; i <= RF_observationSize; i++) {
          if ((RF_tTermMembership[b][i] -> nodeID)  == leaf) {
            (parent -> weight) ++;
          }
        }
      }
      for (i = 1; i <= RF_observationSize; i++) {
        for (j = 1; j <= obsSize; j++) {
          if ( RF_tTermMembership[b][i] == gTermMembership[b][j] ) {
            RF_weightPtr[j][i] +=  1.0 / (double) (gTermMembership[b][j] -> weight);
          }
        }
      }
    }
  }  
}
void getProximity(uint mode, double *proximityPtr) {
  uint  obsSize;
  char  flag;
  uint  i, j;
  int   b;
  if (mode != RF_PRED) {
    obsSize = RF_observationSize;
  }
  else {
    obsSize = RF_fobservationSize;
  }
  uint *offset  = uivector(1, obsSize);
  offset[1] = 0;
  for (i = 2; i <= obsSize; i++) {
    offset[i] = offset[i-1] + i - 1;
  }
  if(RF_opt & OPT_PROX_TYP2) {
    flag = ACTIVE;
  }
  else {
    if(!(RF_opt & OPT_PROX_TYP1)  && !(RF_opt & OPT_PROX_TYP2)) {
      flag = TRUE;
    }
    else if(RF_opt & OPT_PROX_TYP1) {
      flag = FALSE;
    }
    else {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Illegal getProximity() call.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  if (RF_numThreads > 0) {
    for (b = 1; b <= RF_forestSize; b++) {
      updateProximity(b, offset, obsSize, flag, proximityPtr);
    }
  }
  else {
    for (b = 1; b <= RF_forestSize; b++) {
      updateProximity(b, offset, obsSize, flag, proximityPtr);
    }
  }
  if (flag != ACTIVE) {
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= i; j++) {
        if (RF_proximityDen[offset[i] + j] > 0) {
          proximityPtr[offset[i] + j] = proximityPtr[offset[i] + j] /  RF_proximityDen[offset[i] + j];
        }
      }
    }
  }
  else {
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= i; j++) {
        proximityPtr[offset[i] + j] = proximityPtr[offset[i] + j] /  RF_validTreeCount;
      }
    }
  }
  free_uivector(offset, 1, obsSize);
}
void updateProximity(uint b, uint *offset, uint obsSize, char flag, double *proximityPtr) {
  uint i, j;
  if (flag != ACTIVE) {
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= i; j++) {
        if ((RF_bootMembershipFlag[b][i] == flag) && (RF_bootMembershipFlag[b][j] == flag)) {
          RF_proximityDen[offset[i] + j] ++;
          if ( RF_tTermMembershipIndexPtr[b][i] == RF_tTermMembershipIndexPtr[b][j] ) {
            proximityPtr[offset[i] + j] ++;
          }
        }
      }
    }
  }
  else {
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= i; j++) {
        if ( RF_tTermMembershipIndexPtr[b][i] == RF_tTermMembershipIndexPtr[b][j] ) {
          proximityPtr[offset[i] + j] ++;
        }
      }
    }
  }
}
void anticipateProximity(char mode, uint b) {
  uint obsSize;
  uint i;
  if (mode != RF_PRED) {
    obsSize = RF_observationSize;
    for (i = 1; i <= RF_observationSize; i++) {
      RF_bootMembershipFlag[b][i] = ACTIVE;
    }
  }
  else {
    obsSize = RF_fobservationSize;
  }
  for (i=1; i <= obsSize; i++) {
    RF_tTermMembershipIndexPtr[b][i] = i;
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
char pruneBranch(uint obsSize, uint treeID, Node **nodesAtDepth, uint nadCount, uint ptnTarget, uint ptnCurrent) {
  char pruneFlag;
  uint i, j;
  pruneFlag = TRUE;
  double *varianceAtDepth =  dvector(1, nadCount);
  uint   *vadSortedIndex  = uivector(1, nadCount);
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
uint pruneTree(uint obsSize, uint treeID, uint ptnTarget) {
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
  nodesAtDepth = (Node **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  ptnCurrent = RF_tLeafCount[treeID];
  tagDepth = getMaximumDepth(RF_root[treeID]) - 1;
  pruneFlag = (ptnCurrent > ptnTarget) && (tagDepth > 0);
  while (pruneFlag) {
    for (i = 1; i <= RF_tLeafCount[treeID]; i++) {
      nodesAtDepth[i] = NULL;
    }
    nadCount = 0;
    getNodesAtDepth(RF_root[treeID], tagDepth, nodesAtDepth, &nadCount);
    pruneFlag = pruneBranch(obsSize, treeID, nodesAtDepth, nadCount, ptnTarget, ptnCurrent);
    if(pruneFlag) {
      ptnCurrent -= nadCount;
      tagDepth --;
    }
    else {
      ptnCurrent = ptnTarget;
    }
  }
  free_new_vvector(nodesAtDepth, 1, RF_tLeafCount[treeID], NRUTIL_NPTR);
  return ptnCurrent;
}
void unstackAuxiliary(uint mode, uint b) {
  uint obsSize;
  obsSize = 0;  
  free_new_vvector(RF_tNodeMembership[b], 1, RF_observationSize, NRUTIL_NPTR);
  free_uivector(RF_bootMembershipIndex[b], 1, RF_observationSize);
  free_cvector(RF_bootMembershipFlag[b], 1, RF_observationSize);
  free_uivector(RF_bootMembershipCount[b], 1, RF_observationSize);
  free_cvector(RF_oobMembershipFlag[b], 1, RF_observationSize);
  if (mode == RF_PRED) {
    free_new_vvector(RF_ftNodeMembership[b],  1, RF_fobservationSize, NRUTIL_NPTR);
  }
  if (RF_ptnCount > 0) {
    if (mode != RF_PRED) {
      obsSize = RF_observationSize;
    }
    else {
      obsSize = RF_fobservationSize;
    }
    free_new_vvector(RF_pNodeMembership[b], 1, obsSize, NRUTIL_NPTR);
  }
}
void stackNodeList(uint treeID) {
  RF_tNodeList[treeID] = (Node **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_NPTR);
}
void unstackNodeList(uint treeID) {
  free_new_vvector(RF_tNodeList[treeID], 1, RF_tLeafCount[treeID], NRUTIL_NPTR);
}
void initNodeList(uint treeID) {
  uint j;
  for (j = 1; j <= RF_tLeafCount[treeID]; j++) {
    RF_tNodeList[treeID][j] = getTerminalNode(treeID, j);
  }
}
void stackTermList(uint treeID) {
  RF_tTermList[treeID] = (Terminal **) new_vvector(1, RF_tLeafCount[treeID], NRUTIL_TPTR);
}
void unstackTermList(uint treeID) {
  free_new_vvector(RF_tTermList[treeID], 1, RF_tLeafCount[treeID], NRUTIL_NPTR);
}
void initTermList(uint treeID) {
  uint j;
  for (j = 1; j <= RF_tLeafCount[treeID]; j++) {
    RF_tTermList[treeID][j] = makeTerminal();
    RF_tTermList[treeID][j] -> nodeID = j;
  }
}
void stackAndInitTermMembership(uint mode, uint treeID) {
  uint j;
  RF_tTermMembership[treeID] = (Terminal **) new_vvector(1, RF_observationSize, NRUTIL_TPTR);
  for (j = 1; j <= RF_observationSize; j++) {
    RF_tTermMembership[treeID][j] = RF_tTermList[treeID][RF_tNodeMembership[treeID][j] -> nodeID];
  }
  if (mode == RF_PRED) {
    RF_ftTermMembership[treeID] = (Terminal **) new_vvector(1, RF_fobservationSize, NRUTIL_TPTR);
    for (j = 1; j <= RF_fobservationSize; j++) {
      RF_ftTermMembership[treeID][j] = RF_tTermList[treeID][RF_ftNodeMembership[treeID][j] -> nodeID];
    }
  }
}
void unstackTermMembership(uint mode, uint treeID) {
  free_new_vvector(RF_tTermMembership[treeID], 1, RF_observationSize, NRUTIL_TPTR);
  if (mode == RF_PRED) {
    free_new_vvector(RF_ftTermMembership[treeID], 1, RF_fobservationSize, NRUTIL_TPTR);
  }
}
