////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.3
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
#include       "nodeOps.h"
#include     "bootstrap.h"
#include    "splitUtil.h"
#include      "treeUtil.h"
#include     "factorOps.h"
#include    "regression.h"
#include        "impute.h"
char imputeNode (uint     type,
                 char     lmvFlag,
                 char     chainFlag,
                 uint     treeID, 
                 Node    *nodePtr,
                 uint    *repMembrIndx,
                 uint     repMembrSize,
                 uint    *allMembrIndx,
                 uint     allMembrSize) {
  double  **response;
  double  **predictor;
  uint    *mRecordMap;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  int     *mvNSptr;
  uint     mRecordSize;
  double *valuePtr, *imputePtr;
  char mPredictorFlag;
  int  signedSignatureIndex;
  uint unsignedIndexSource;
  uint unsignedIndexTarget;
  char result;
  uint  *glmpIndexPtr;
  uint  *glmpIndexSize;
  uint  *glmpIndexParentPtr;
  uint   glmpIndexParentSize;
  uint  *glmrIndexPtr;
  uint  *glmrIndexSize;
  uint  *glmrIndexParentPtr;
  uint   glmrIndexParentSize;
  char mvFlag;
  uint i,p;
  uint localDistributionSize;
  mvNSptr = NULL;  
  mpIndex = NULL;  
  mpSign  = NULL;  
  mpIndexSize  = 0;  
  mRecordMap = NULL;  
  mRecordSize = 0;    
  predictor  = NULL;  
  response   = NULL;  
  glmpIndexPtr = NULL;
  glmpIndexSize = NULL;
  glmpIndexParentPtr = NULL;
  glmpIndexParentSize = 0;
  glmrIndexPtr = NULL;
  glmrIndexSize = NULL;
  glmrIndexParentPtr = NULL;
  glmrIndexParentSize = 0;
  result = FALSE;
  switch (type) {
  case RF_PRED:
    mRecordSize = RF_fmRecordSize;
    if (mRecordSize > 0) {
      if (RF_frSize > 0) {
        response = RF_fresponse[treeID];
      } 
      else {
        response = NULL;
      }
      predictor = RF_fobservation[treeID];
      mRecordMap = RF_fmRecordMap;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      mvNSptr = nodePtr -> fmpSign;
      if((nodePtr -> parent) == NULL) {
        glmpIndexParentPtr = uivector(1, mpIndexSize);
        glmpIndexParentSize = mpIndexSize;
        for (p = 1; p <= glmpIndexParentSize; p++) {
          glmpIndexParentPtr[p] = p;
        }
        glmrIndexParentPtr = uivector(1, mRecordSize);
        glmrIndexParentSize = mRecordSize;
        for (i = 1; i <= glmrIndexParentSize; i++) {
          glmrIndexParentPtr[i] = i;
        }
        stackNodeFLMPIndex(nodePtr, glmpIndexParentSize);
        glmpIndexPtr  = nodePtr -> flmpIndex;
        glmpIndexSize = & (nodePtr -> flmpIndexActualSize);
        *glmpIndexSize = 0;
        stackNodeFLMRIndex(nodePtr, glmrIndexParentSize);
        glmrIndexPtr  = nodePtr -> flmrIndex;
        glmrIndexSize = & (nodePtr -> flmrIndexActualSize);
        *glmrIndexSize = 0;
      }
      else {
        if((nodePtr -> parent) -> flmpIndexActualSize > 0) {
          glmpIndexParentPtr = (nodePtr -> parent) -> flmpIndex;
          glmpIndexParentSize = (nodePtr -> parent) -> flmpIndexActualSize;
          stackNodeFLMPIndex(nodePtr, glmpIndexParentSize);
          glmpIndexPtr  = nodePtr -> flmpIndex;
          glmpIndexSize = & (nodePtr -> flmpIndexActualSize);
          *glmpIndexSize = 0;
        }
        else {
          glmpIndexParentPtr  = NULL;
          glmpIndexParentSize = 0;
          glmpIndexPtr = glmpIndexSize = NULL;
        }
        if((nodePtr -> parent) -> flmrIndexActualSize > 0) {
          glmrIndexParentPtr = (nodePtr -> parent) -> flmrIndex;
          glmrIndexParentSize = (nodePtr -> parent) -> flmrIndexActualSize;
          stackNodeFLMRIndex(nodePtr, glmrIndexParentSize);
          glmrIndexPtr  = nodePtr -> flmrIndex;
          glmrIndexSize = & (nodePtr -> flmrIndexActualSize);
          *glmrIndexSize = 0;
        }
        else {
          glmrIndexParentPtr = NULL;
          glmrIndexParentSize = 0;
          glmrIndexPtr = glmrIndexSize = NULL;
        }
      }
      result = TRUE;
    }
    break;
  default:
    mRecordSize = RF_mRecordSize;
    if (mRecordSize > 0) {
      response = RF_response[treeID];
      predictor = RF_observation[treeID];
      mRecordMap = RF_mRecordMap;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      mvNSptr = nodePtr -> mpSign;
      if (lmvFlag) {
        if((nodePtr -> parent) == NULL) {
          glmpIndexParentPtr = uivector(1, mpIndexSize);
          glmpIndexParentSize = mpIndexSize;
          for (p = 1; p <= glmpIndexParentSize; p++) {
            glmpIndexParentPtr[p] = p;
          }
          glmrIndexParentPtr = uivector(1, mRecordSize);
          glmrIndexParentSize = mRecordSize;
          for (i = 1; i <= glmrIndexParentSize; i++) {
            glmrIndexParentPtr[i] = i;
          }
          stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
          glmpIndexPtr  = nodePtr -> lmpIndex;
          glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
          *glmpIndexSize = 0;
          stackNodeLMRIndex(nodePtr, glmrIndexParentSize);
          glmrIndexPtr  = nodePtr -> lmrIndex;
          glmrIndexSize = & (nodePtr -> lmrIndexActualSize);
          *glmrIndexSize = 0;
        }
        else {
          if((nodePtr -> parent) -> lmpIndexActualSize > 0) {
            glmpIndexParentPtr = (nodePtr -> parent) -> lmpIndex;
            glmpIndexParentSize = (nodePtr -> parent) -> lmpIndexActualSize;
            stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
            glmpIndexPtr  = nodePtr -> lmpIndex;
            glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
            *glmpIndexSize = 0;
          }
          else {
            glmpIndexParentPtr = NULL;
            glmpIndexParentSize = 0;
            glmpIndexPtr = glmpIndexSize = NULL;
          }
          if((nodePtr -> parent) -> lmrIndexActualSize > 0) {
            glmrIndexParentPtr = (nodePtr -> parent) -> lmrIndex;
            glmrIndexParentSize = (nodePtr -> parent) -> lmrIndexActualSize;
            stackNodeLMRIndex(nodePtr, glmrIndexParentSize);
            glmrIndexPtr  = nodePtr -> lmrIndex;
            glmrIndexSize = & (nodePtr -> lmrIndexActualSize);
            *glmrIndexSize = 0;
          }
          else {
            glmrIndexParentPtr = NULL;
            glmrIndexParentSize = 0;
            glmrIndexPtr = glmrIndexSize = NULL;
          }
        }
      }  
      else {
        glmpIndexParentPtr = uivector(1, mpIndexSize);
        glmpIndexParentSize = mpIndexSize;
        for (p = 1; p <= glmpIndexParentSize; p++) {
          glmpIndexParentPtr[p] = p;
        }
        glmrIndexParentPtr = uivector(1, mRecordSize);
        glmrIndexParentSize = mRecordSize;
        for (i = 1; i <= glmrIndexParentSize; i++) {
          glmrIndexParentPtr[i] = i;
        }
        stackNodeLMPIndex(nodePtr, glmpIndexParentSize);
        glmpIndexPtr  = nodePtr -> lmpIndex;
        glmpIndexSize = & (nodePtr -> lmpIndexActualSize);
        *glmpIndexSize = 0;
        stackNodeLMRIndex(nodePtr, glmrIndexParentSize);
        glmrIndexPtr  = nodePtr -> lmrIndex;
        glmrIndexSize = & (nodePtr -> lmrIndexActualSize);
        *glmrIndexSize = 0;
      }
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to impute node with no missingness in type:  %10d", type);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  double *localDistribution = dvector(1, repMembrSize + 1);
  for (p = 1; p <= glmpIndexParentSize; p++) {
    if (mvNSptr[glmpIndexParentPtr[p]] != -1) {
      signedSignatureIndex = mpIndex[glmpIndexParentPtr[p]];
      if (signedSignatureIndex < 0) {
        unsignedIndexSource = unsignedIndexTarget = (uint) abs(signedSignatureIndex);
        valuePtr = RF_response[treeID][(uint) abs(signedSignatureIndex)];
        imputePtr = response[(uint) abs(signedSignatureIndex)];
      }
      else { 
        unsignedIndexSource = RF_rSize + (uint) signedSignatureIndex;
        if (type == RF_PRED) {
          if (RF_frSize > 0) {
            unsignedIndexTarget = RF_rSize + (uint) signedSignatureIndex;
          }
          else {
            unsignedIndexTarget = (uint) signedSignatureIndex;
          }
        }
        else {
          unsignedIndexTarget = RF_rSize + (uint) signedSignatureIndex;
        }
        valuePtr = RF_observation[treeID][(uint) signedSignatureIndex];
        imputePtr = predictor[(uint) signedSignatureIndex];
      }
      localDistributionSize = 0;
      for (i = 1; i <= repMembrSize; i++) {
        mPredictorFlag = TRUE;
        if (RF_mRecordMap[repMembrIndx[i]] == 0) {
          mPredictorFlag = FALSE;
        }
        else if (RF_mpSign[unsignedIndexSource][RF_mRecordMap[repMembrIndx[i]]] == 0) {
          mPredictorFlag = FALSE;
        }
        if (mPredictorFlag == FALSE) {
          localDistributionSize ++;
          localDistribution[localDistributionSize] = valuePtr[repMembrIndx[i]];
        }
      }  
      mvFlag = FALSE;
      for (i = 1; i <= allMembrSize; i++) {
        if (mRecordMap[allMembrIndx[i]] > 0) {
          if(mpSign[unsignedIndexTarget][mRecordMap[allMembrIndx[i]]] == 1) {
            mvFlag = TRUE;
            if (localDistributionSize > 0) {
              imputePtr[allMembrIndx[i]] = getSampleValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }
          }  
        }
      }  
      if (mvFlag) {
        glmpIndexPtr[++(*glmpIndexSize)] = glmpIndexParentPtr[p]; 
      }
      if (localDistributionSize == 0) {
      }
    }  
  }  
  for (i = 1; i <= allMembrSize; i++) {
    if (mRecordMap[allMembrIndx[i]] > 0) {
      glmrIndexPtr[++(*glmrIndexSize)] = mRecordMap[allMembrIndx[i]];
    }
  }
  free_dvector(localDistribution, 1, repMembrSize + 1);
  if (lmvFlag) {
    if((nodePtr -> parent) == NULL) {
      free_uivector(glmpIndexParentPtr, 1, mpIndexSize);
      free_uivector(glmrIndexParentPtr, 1, mRecordSize);
    }
  }
  else {
    free_uivector(glmpIndexParentPtr, 1, mpIndexSize);
    free_uivector(glmrIndexParentPtr, 1, mRecordSize);
  }
  return TRUE;
}
char restoreNodeMembership(uint  mode, 
                           char  rootFlag,
                           uint  treeID, 
                           Node *parent, 
                           uint *repMembrIndx,
                           uint  repMembrSize,
                           uint *allMembrIndx,
                           uint  allMembrSize,
                           uint *ngAllMembrIndx,
                           uint  ngAllMembrSize,
                           uint *bootMembrIndxIter) {
  char  bootResult;
  char leftResult, rghtResult;
  char tnUpdateFlag;
  char bsUpdateFlag;
  uint *bootMembrIndx;  
  uint *leftRepMembrIndx; 
  uint *rghtRepMembrIndx;
  uint *leftAllMembrIndx;
  uint *rghtAllMembrIndx;
  uint *ngLeftAllMembrIndx;  
  uint *ngRghtAllMembrIndx;  
  uint bootMembrSize; 
  uint leftRepMembrSize, rghtRepMembrSize;
  uint leftAllMembrSize, ngLeftAllMembrSize;
  uint rghtAllMembrSize, ngRghtAllMembrSize;
  uint jLeft;
  uint jRght;
  char daughterFlag;
  uint i;
  bootResult = TRUE;
  tnUpdateFlag = TRUE;
  bsUpdateFlag = FALSE;
  if (rootFlag | (RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    bootMembrIndx  = uivector(1, allMembrSize);
    bootMembrSize = allMembrSize;
    bootResult = bootstrap (mode,
                            treeID, 
                            parent,
                            allMembrIndx,
                            allMembrSize,
                            bootMembrIndx);
    if (rootFlag & bootResult) {
      if (!(RF_opt & (OPT_BOOT_NODE | OPT_BOOT_NONE))) {
        if (RF_opt & OPT_MEMB) {
          if (mode == RF_REST) {
            for (i=1; i <=  allMembrSize; i++) {
              RF_bootstrapMembershipPtr[treeID][bootMembrIndx[i]] ++;
            }
          }
          else {
            for (i=1; i <= ngAllMembrSize; i++) {
              RF_bootstrapMembershipPtr[treeID][i] ++;
            }
          }
        }
        bsUpdateFlag = TRUE;
      }
      repMembrIndx = bootMembrIndx;
      repMembrSize = bootMembrSize;
    }
  }
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mpSign = (parent -> parent) -> mpSign;
    parent -> fmpSign = (parent -> parent) -> fmpSign;
  }
  if (bootResult) {
    if (RF_mRecordSize > 0) {
      imputeNode(RF_GROW,
                 TRUE,
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
                               parent, 
                               allMembrIndx, 
                               allMembrSize, 
                               RF_time[treeID], 
                               FALSE, 
                               FALSE, 
                               RF_masterTimeIndex[treeID]);
        }
      }
    }
    switch (mode) {
    case RF_PRED:
      if (RF_fmRecordSize > 0) {
        imputeNode(RF_PRED,
                   TRUE,
                   FALSE,
                   treeID,
                   parent,
                   bootMembrIndx, 
                   bootMembrSize, 
                   ngAllMembrIndx,
                   ngAllMembrSize);
      }
      break;
    case RF_REST:
      break;
    default:
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Unused GROW case in switch encountered. ");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
      break;
    }
  }  
  if (bootResult) {
    if ((RF_opt & OPT_NODE_STAT) || (RF_ptnCount > 0)) {
      getVariance(repMembrSize, repMembrIndx, RF_response[treeID][RF_rTarget], NULL, & (parent -> variance));  
    }
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      tnUpdateFlag = FALSE;
      uint *membershipIndicator = uivector(1, RF_observationSize);
      leftAllMembrSize = rghtAllMembrSize = 0;
      for (i=1; i <= allMembrSize; i++) {
        daughterFlag = RIGHT;
        if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
          daughterFlag = splitOnFactor((uint) RF_observation[treeID][parent -> splitParameter][allMembrIndx[i]], parent -> splitValueFactPtr);
        }
        else {
          if ( RF_observation[treeID][parent -> splitParameter][allMembrIndx[i]] <= (parent -> splitValueCont) ) {
            daughterFlag = LEFT;
          }
        }
        if (daughterFlag == LEFT) {
          membershipIndicator[allMembrIndx[i]] = LEFT;
          leftAllMembrSize ++;
          RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> left;
        }
        else {
          membershipIndicator[allMembrIndx[i]] = RIGHT;
          rghtAllMembrSize ++;
          RF_tNodeMembership[treeID][allMembrIndx[i]] = parent -> right;
        }
      }
      leftAllMembrIndx  = uivector(1, leftAllMembrSize + 1);
      rghtAllMembrIndx  = uivector(1, rghtAllMembrSize + 1);
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
        leftRepMembrIndx  = uivector(1, bootMembrSize + 1);
        rghtRepMembrIndx  = uivector(1, bootMembrSize + 1);
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
      ngLeftAllMembrIndx = ngRghtAllMembrIndx = NULL;
      ngLeftAllMembrSize = ngRghtAllMembrSize = 0;
      if (mode == RF_PRED) {
        uint *ngMembershipIndicator = uivector(1, RF_fobservationSize);
        for (i=1; i <= ngAllMembrSize; i++) {
          daughterFlag = RIGHT;
          if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
            daughterFlag = splitOnFactor((uint) RF_fobservation[treeID][parent -> splitParameter][ngAllMembrIndx[i]], parent -> splitValueFactPtr);
          }
          else {
            if ( RF_fobservation[treeID][parent -> splitParameter][ngAllMembrIndx[i]] <= (parent -> splitValueCont) ) {
              daughterFlag = LEFT;
            }
          }
          if (daughterFlag == LEFT) {
            ngMembershipIndicator[ngAllMembrIndx[i]] = LEFT;
            ngLeftAllMembrSize ++;
            RF_ftNodeMembership[treeID][ngAllMembrIndx[i]] = parent -> left;
          }
          else {
            ngMembershipIndicator[ngAllMembrIndx[i]] = RIGHT;
            ngRghtAllMembrSize ++;
            RF_ftNodeMembership[treeID][ngAllMembrIndx[i]] = parent -> right;
          }
        }
        ngLeftAllMembrIndx  = uivector(1, ngLeftAllMembrSize + 1);
        ngRghtAllMembrIndx  = uivector(1, ngRghtAllMembrSize + 1);
        jLeft = jRght = 0;
        for (i = 1; i <= ngAllMembrSize; i++) {
          if (ngMembershipIndicator[ngAllMembrIndx[i]] == LEFT) {
            ngLeftAllMembrIndx[++jLeft] = ngAllMembrIndx[i];
          }
          else {
            ngRghtAllMembrIndx[++jRght] = ngAllMembrIndx[i];
          }
        }
        free_uivector(ngMembershipIndicator, 1, RF_fobservationSize);
      }  
      leftResult = restoreNodeMembership(mode, 
                                         FALSE,
                                         treeID, 
                                         parent -> left, 
                                         leftRepMembrIndx,
                                         leftRepMembrSize,
                                         leftAllMembrIndx,
                                         leftAllMembrSize,
                                         ngLeftAllMembrIndx,
                                         ngLeftAllMembrSize,
                                         bootMembrIndxIter);
      if(!leftResult) {
      }
      rghtResult = restoreNodeMembership(mode, 
                                         FALSE,
                                         treeID, 
                                         parent -> right, 
                                         rghtRepMembrIndx,
                                         rghtRepMembrSize,
                                         rghtAllMembrIndx,
                                         rghtAllMembrSize,
                                         ngRghtAllMembrIndx,
                                         ngRghtAllMembrSize,
                                         bootMembrIndxIter);
      if(!rghtResult) {
      }
      free_uivector(leftAllMembrIndx, 1, leftAllMembrSize + 1);
      free_uivector(rghtAllMembrIndx, 1, rghtAllMembrSize + 1);
      if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      }
      else {
        free_uivector(leftRepMembrIndx, 1, bootMembrSize + 1);
        free_uivector(rghtRepMembrIndx, 1, bootMembrSize + 1);
      }
      if (mode == RF_PRED) {
        free_uivector(ngLeftAllMembrIndx, 1, ngLeftAllMembrSize + 1);
        free_uivector(ngRghtAllMembrIndx, 1, ngRghtAllMembrSize + 1);
      }
      free_uivector(membershipIndicator, 1, RF_observationSize);
    }  
    else {
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
    if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      bsUpdateFlag = TRUE;
    }
  }  
  if (bsUpdateFlag) {
    for (i = 1; i <= bootMembrSize; i++) {
      RF_bootMembershipIndex[treeID][++(*bootMembrIndxIter)] = bootMembrIndx[i];
      RF_bootMembershipFlag[treeID][bootMembrIndx[i]] = TRUE;
      RF_oobMembershipFlag[treeID][bootMembrIndx[i]]  = FALSE;
    }
  }
  if (rootFlag | (RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    free_uivector(bootMembrIndx, 1, allMembrSize);
  }  
  return bootResult;
}
void imputeUpdateSummary (uint     mode, 
                          double **responsePtr, 
                          double **predictorPtr, 
                          uint     treeID) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  int     *mvNodeSign;
  double  *valuePtr;
  int      signedSignatureIndex;
  uint     unsignedIndex;
  char result;
  uint rspSize;
  uint i, p;
  rspSize = 0;  
  mpIndex = 0;  
  mpSign  = 0;  
  mpIndexSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  result = FALSE;
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      rspSize = RF_rSize;
      result = TRUE;
    }
  }
  else {
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      rspSize = RF_frSize;
      result = TRUE;
    }
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to update forest impute data with no missingness in mode:  %10d", mode);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (p = 1; p <= mpIndexSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      if ((mode == RF_GROW) || (mode == RF_REST)) {
          mvNodeSign = RF_tNodeMembership[treeID][mRecordIndex[i]] -> mpSign;
      }
      else {
          mvNodeSign = RF_ftNodeMembership[treeID][mRecordIndex[i]] -> fmpSign;
      }
      if (mvNodeSign[p] != -1) {
        signedSignatureIndex = mpIndex[p];
        if (signedSignatureIndex < 0) {
          unsignedIndex = (uint) abs(signedSignatureIndex);
          valuePtr = responsePtr[(uint) abs(signedSignatureIndex)];
        }
        else {
          unsignedIndex = (uint) signedSignatureIndex + rspSize;
          valuePtr = predictorPtr[(uint) signedSignatureIndex];
        }
        if (mpSign[unsignedIndex][i] == 1) {
          if (ISNA(valuePtr[mRecordIndex[i]])) {
            Rprintf("\nDiagnostic Trace of Shadowed Data:  ");
            Rprintf("\n       index   imputation -> \n");
            Rprintf(  "            ");
            for (p=1; p <= mpIndexSize; p++) {
              Rprintf(" %12d", mpIndex[p]);
            }
            Rprintf("\n");
            for (i = 1; i <= mRecordSize; i++) {
              Rprintf("%12d", mRecordIndex[i]);
              for (p = 1; p <= mpIndexSize; p++) {
                if (mpIndex[p] < 0) {
                  valuePtr = responsePtr[(uint) abs(mpIndex[p])];
                }
                else {
                  valuePtr = predictorPtr[(uint) mpIndex[p]];
                }
                Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
              }
              Rprintf("\n");
            }
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Attempt to update forest impute data with invalid shadowed value, NA. ");
            Rprintf("\nRF-SRC:  Invalid value for:  [indv][outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mpIndex[p]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }  
          else {
            RF_dmvImputation[treeID][i][p] = valuePtr[mRecordIndex[i]];
          }
        }  
      }  
    }  
  }  
}
void imputeUpdateShadow (uint      mode, 
                         char      selectionFlag,
                         double  **shadowResponse, 
                         double  **shadowPredictor) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  double **outResponse;
  double **outPredictor;
  double  *valuePtr;
  double  *outputPtr;
  uint unsignedIndex;
  char outcomeFlag;
  uint rspSize;
  uint i, p;
  mRecordSize  = 0;     
  mRecordIndex = NULL;  
  mpIndexSize       = 0;     
  mpSign       = NULL;  
  mpIndex      = NULL;  
  outResponse  = NULL;  
  outPredictor = NULL;  
  valuePtr     = NULL;  
  outputPtr    = NULL;  
  unsignedIndex = 0;    
  switch (mode) {
  case RF_PRED:
    mRecordSize = RF_fmRecordSize;
    mRecordIndex = RF_fmRecordIndex;
    mpIndexSize = RF_fmpIndexSize;
    mpSign = RF_fmpSign;
    mpIndex = RF_fmpIndex;
    outResponse  = RF_sImputeResponsePtr;
    outPredictor = RF_sImputePredictorPtr;
    rspSize = RF_frSize;
    break;
  default:
    mRecordSize = RF_mRecordSize;
    mRecordIndex = RF_mRecordIndex;
    mpIndexSize = RF_mpIndexSize;
    mpSign = RF_mpSign;
    mpIndex = RF_mpIndex;
    if ((selectionFlag == TRUE) || (selectionFlag == ACTIVE)) {
      outResponse  = RF_sImputeResponsePtr;
      outPredictor = RF_sImputePredictorPtr;
    }
    else {
      outResponse  = RF_sOOBImputeResponsePtr;
      outPredictor = RF_sOOBImputePredictorPtr;
    }
    rspSize = RF_rSize;
    break;
  }
  if (mRecordSize == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to update shadow data with no missingness in mode:  %10d", mode);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  outcomeFlag = TRUE;
  for (p = 1; p <= mpIndexSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      if (mpIndex[p] < 0) {
        unsignedIndex = (uint) abs(mpIndex[p]);
        valuePtr = shadowResponse[(uint) abs(mpIndex[p])];
        outputPtr = outResponse[(uint) abs(mpIndex[p])];
      }
      else {
        if (shadowPredictor != NULL) {
          unsignedIndex = (uint) mpIndex[p] + rspSize;
          valuePtr = shadowPredictor[(uint) mpIndex[p]];
          outputPtr = outPredictor[(uint) mpIndex[p]];
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || (shadowPredictor != NULL)) {
        if (mpSign[unsignedIndex][i] == 1) {
          if (ISNA(outputPtr[i])) {
          }
          valuePtr[mRecordIndex[i]] = outputPtr[i];
        }  
      }
    }  
  }  
}
void imputeSummary(uint      mode,
                   char      selectionFlag) {
  imputeCommonNew(mode,
               0,
               selectionFlag, 
               TRUE);
}
void imputeResponse(uint      mode,
                    uint      serialTreeID,
                    double  **tempResponse) {
  switch(mode) {
  case RF_PRED:
    imputeCommonNew(mode, serialTreeID, ACTIVE, FALSE);
    imputeUpdateShadow(mode, ACTIVE, tempResponse, NULL);
    break;
  default:
    imputeCommonNew(mode, serialTreeID, FALSE, FALSE);
    imputeUpdateShadow(mode, FALSE, tempResponse, NULL);
    break;
  }
}
void imputeCommon(uint      mode,
                  uint      serialTreeID,
                  char      selectionFlag,
                  char      predictorFlag) {
  uint *localSerialIndex;
  uint  localSerialCount;
  uint *serialPtr;
  char mFlag;
  char outcomeFlag;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  double **outResponse;
  double **outPredictor;
  double *valuePtr;
  double *naivePtr;
  uint    unsignedIndex;
  double imputedValue;
  uint localDistributionSize;
  uint maxDistributionSize;
  uint rspSize;
  char result;
  uint i,p,tree;
  valuePtr      = NULL;  
  naivePtr      = NULL;  
  unsignedIndex = 0;     
  maxDistributionSize = 0;  
  outResponse         = 0;  
  outPredictor        = 0;  
  rspSize = 0;  
  mpIndex = 0;  
  mpSign  = 0;  
  mpIndexSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  if ((selectionFlag != TRUE) && (selectionFlag != FALSE) && (selectionFlag != ACTIVE)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid selectionFlag in imputeCommon():  %10d", selectionFlag);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      outResponse  = RF_sImputeResponsePtr;
      outPredictor = RF_sImputePredictorPtr;
      rspSize = RF_frSize;
      result = TRUE;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      if ((selectionFlag == TRUE) || (selectionFlag == ACTIVE)) {
        outResponse  = RF_sImputeResponsePtr;
        outPredictor = RF_sImputePredictorPtr;
      }
      else {
        outResponse  = RF_sOOBImputeResponsePtr;
        outPredictor = RF_sOOBImputePredictorPtr;
      }
      rspSize = RF_rSize;
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to impute in imputeCommon() with no missingness in mode:  %10d", mode);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (serialTreeID == 0) {
    localSerialIndex = uivector(1, RF_forestSize);
    for (tree = 1; tree <= RF_forestSize; tree++) {
      localSerialIndex[tree] = tree;
    }
    serialPtr = localSerialIndex;
    localSerialCount = RF_forestSize;
  }
  else {
    serialPtr = RF_serialTreeIndex;
    localSerialCount = serialTreeID;
  }
  imputedValue = 0.0;  
  double *localDistribution = dvector(1, maxDistributionSize);
  char  *naiveFlag = cvector(1, mpIndexSize);
  char **naiveSign = cmatrix(1, mRecordSize, 1, mpIndexSize);
  for (p = 1; p <= mpIndexSize; p++) {
    naiveFlag[p] = FALSE;
  }
  for (i = 1; i <= mRecordSize; i++) {
    outcomeFlag = TRUE;
    for (p = 1; p <= mpIndexSize; p++) {
      naiveSign[i][p] = FALSE;  
      if (mpIndex[p] < 0) {
        unsignedIndex = (uint) abs(mpIndex[p]);
      }
      else {
        if (predictorFlag == TRUE) {
          unsignedIndex = (uint) mpIndex[p] + rspSize;
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || predictorFlag) {
        if (mpSign[unsignedIndex][i] == 1) {
          localDistributionSize = 0;
          for (tree = 1; tree <= localSerialCount; tree++) {
            if (RF_tLeafCount[serialPtr[tree]] > 0) {
              if ((RF_dmRecordBootFlag[serialPtr[tree]][i] == selectionFlag) || (selectionFlag == ACTIVE)) {
                if (!ISNA(RF_dmvImputation[serialPtr[tree]][i][p])) {
                  localDistribution[++localDistributionSize] = RF_dmvImputation[serialPtr[tree]][i][p];
                }
                else {
                }  
              }  
            }  
          }  
          if (localDistributionSize > 0) {
            if (mpIndex[p] < 0) {
              if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "T") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "S") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "I") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "C") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              outResponse[(uint) abs(mpIndex[p])][i] = imputedValue;
            }  
            else {
              if (strcmp(RF_xType[(uint) mpIndex[p]], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              outPredictor[(uint) mpIndex[p]][i] = imputedValue;
            }
          }  
          else {
            naiveFlag[p] = TRUE;
            naiveSign[i][p] = TRUE;
          }
        }  
      }  
      else {
        p = mpIndexSize;
      }
    }  
  }  
  outcomeFlag = TRUE;
  for (p = 1; p <= mpIndexSize; p++) {
    if (mpIndex[p] < 0) {
      unsignedIndex = (uint) abs(mpIndex[p]);
      valuePtr = RF_responseIn[(uint) abs(mpIndex[p])];
      naivePtr = outResponse[(uint) abs(mpIndex[p])];
    }
    else {
      if (predictorFlag == TRUE) {
        unsignedIndex = (uint) mpIndex[p] + rspSize;
        valuePtr = RF_observationIn[(uint) mpIndex[p]];
        naivePtr = outPredictor[(uint) mpIndex[p]];
      }
      outcomeFlag = FALSE;
    }
    if (outcomeFlag || predictorFlag) {
      if (naiveFlag[p] == TRUE) {
        localDistributionSize = 0;
        for (i=1; i <= RF_observationSize; i++) {
          mFlag = TRUE;
          if (RF_mRecordMap[i] == 0) {
            mFlag = FALSE;
          }
          else if (RF_mpSign[unsignedIndex][RF_mRecordMap[i]] == 0) {
            mFlag = FALSE;
          }
          if (mFlag == FALSE) {
            localDistribution[++localDistributionSize] = valuePtr[i];
          }
        }  
        if (localDistributionSize > 0) {
          for (i=1; i <= mRecordSize; i++) {
            if (naiveSign[i][p] == TRUE) {
              naivePtr[i] = getSampleValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
            }
          }
        }  
        else {
          if (mpIndex[p] < 0) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Naive imputation failed for [indv, outcome] = [%10d, %10d] \n", mRecordIndex[i], mpIndex[p]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
          else {
          }
        }
      }  
    }  
    else {
      p = mpIndexSize;
    }
  }  
  if (serialTreeID == 0) {
    free_uivector(localSerialIndex, 1, RF_forestSize);
  }
  free_dvector(localDistribution, 1, maxDistributionSize);
  free_cvector(naiveFlag, 1, mpIndexSize);
  free_cmatrix(naiveSign, 1, mRecordSize, 1, mpIndexSize);
}
void imputeMultipleTime (char selectionFlag) {
  double  *outTime;
  double meanValue;
  double leftDistance, rightDistance;
  uint minimumIndex;
  char result;
  uint i,j;
  result = FALSE;
    if (RF_timeIndex > 0) {
      if (RF_mRecordSize > 0) {
      if (RF_mTimeFlag == TRUE) {
        result = TRUE;
      }
      else {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Attempt to multiply impute time with no missingness in time vector.");
      }
    }
  }
  else {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to multiply impute time in a !SURV environment.");
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (selectionFlag == FALSE) {
    outTime   = RF_sOOBImputeResponsePtr[RF_timeIndex];
  }
  else {
    outTime  = RF_sImputeResponsePtr[RF_timeIndex];
  }    
  for (i=1; i <= RF_mRecordSize; i++) {
    if(RF_mpSign[RF_timeIndex][i] == 1) { 
      meanValue = outTime[i];
      if ((meanValue < RF_masterTime[1]) || (meanValue > RF_masterTime[RF_masterTimeSize])) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  The summary mean value for time is out of range:  indv %10d, value %12.4f", RF_mRecordIndex[i], meanValue);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      leftDistance = rightDistance = 0;
      minimumIndex = 0;
      for (j = 1; j <= RF_masterTimeSize; j++) {
        if (meanValue <= RF_masterTime[j]) {
          minimumIndex = j;
          j = RF_masterTimeSize;
        }
      }
      if (minimumIndex == 1) {
      }
      else {
        leftDistance = meanValue - RF_masterTime[minimumIndex-1];
        rightDistance = RF_masterTime[minimumIndex] - meanValue;
        if (leftDistance < rightDistance) {
          minimumIndex = minimumIndex - 1;            
        }
        else {
          if (fabs(leftDistance - rightDistance) < EPSILON) {
            if (ran2(1) <= 0.5) {
              minimumIndex = minimumIndex - 1;
            }
          }
        }
      }
      outTime[i] = RF_masterTime[minimumIndex];
    }
  }
}
double getMaximalValue(double *value, uint size, uint treeID) {
  double result;
  uint classCount, maximalClassSize, maximalClassCount;
  uint randomIndex;
  uint j;
  uint   *classSize  = uivector(1, size);
  for (j = 1; j <= size; j++) {
    classSize[j] = 0;
  }
  hpsort(value, size);
  classCount = 1;
  classSize[1] = 1;
  for (j = 2; j <= size; j++) {
    if (value[j] > value[classCount]) {
      classCount ++;
      value[classCount] = value[j];
    }
    classSize[classCount] ++;
  }
  maximalClassSize = maximalClassCount = 0;
  for (j=1; j <= classCount; j++) {
    if (classSize[j] > maximalClassSize) {
      maximalClassSize = classSize[j];
    }
  }
  for (j=1; j <= classCount; j++) {
    if (classSize[j] == maximalClassSize) {
      maximalClassCount ++;
    }
  }
  if (maximalClassCount > 1) {
    randomIndex = (uint) ceil(ran2(treeID)*((maximalClassCount)*1.0));
  }
  else {
    randomIndex = 1;
  }
  j = 0;
  while (randomIndex > 0) {
    j++;
    if (classSize[j] == maximalClassSize) {
      randomIndex --;
    }
  }
  result = value[j];
  free_uivector(classSize, 1, size);
  return result;
}
double getMedianValue(double *value, uint size) {
  double result;
  uint medianIndex;
  hpsort(value, size);
  if (size > 1) {
    medianIndex = (uint) ceil(size/2);
  }
  else {
    medianIndex = 1;
  }
  result = value[medianIndex];
  return result;
}
double getMeanValue(double *value, uint size) {
  double result;
  uint j;
  result = 0.0;
  for (j = 1; j <= size; j++) {
    result = result + value[j];
  }
  result = result / size;
  return result;
}
double getSampleValue(double *value, uint size, char chainFlag, uint treeID) {
  uint randomIndex;
  if(chainFlag) {
    randomIndex = (uint) ceil(ran1(treeID)*((size)*1.0));
  }
  else {
    randomIndex = (uint) ceil(ran2(treeID)*((size)*1.0));
  }
  return value[randomIndex];
}
uint getRecordMap(uint    *map, 
                  uint     obsSize, 
                  double **resp,
                  double **data) {
  uint i, p, r;
  uint mSize;
  char mFlag;
  mSize  = 0;
  for (i = 1; i <= obsSize; i++) {
    mFlag = FALSE;
    if (resp != NULL) {
      for (r = 1; r <= RF_rSize; r++) {
        if (ISNA(resp[r][i])) {
          mFlag = TRUE;
          r = RF_rSize;
        }
      }
    }
    if (mFlag == FALSE) {
      for (p = 1; p <= RF_xSize; p++) {
        if (ISNA(data[p][i])) {
          mFlag = TRUE;
          p = RF_xSize;
        }
      }
    }
    if (mFlag == TRUE) {
      mSize ++;
      map[i] = mSize;
    }
    else {
      map[i] = 0;
    }
  }
  return mSize;
}
void updateTimeIndexArray(uint    treeID,
                          Node   *parent, 
                          uint   *allMembrIndx,
                          uint    allMembrSize,
                          double *time, 
                          char    naAllowFlag,
                          char    noIdxAllowFlag,
                          uint   *masterTimeIndex) {
  uint *membrIndx;
  char idxFoundFlag;
  uint i,k;
  if (allMembrIndx == NULL) {
    membrIndx = uivector(1, allMembrSize);
    for (i = 1; i <= allMembrSize; i++) {
      membrIndx[i] = i;
    }
  }
  else {
    membrIndx = allMembrIndx;
  }
  for (i=1; i <= allMembrSize; i++) {
    idxFoundFlag = FALSE;
    if (!ISNA(time[membrIndx[i]])) {
      k = 1;
      while (k <= RF_masterTimeSize) {
        if (time[membrIndx[i]] == RF_masterTime[k]) {
          masterTimeIndex[membrIndx[i]] = k;
          idxFoundFlag = TRUE;
          k = RF_masterTimeSize;
        }
        k++;
      }
    }
    else {
      if (naAllowFlag == FALSE) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Missing event time encountered for individual:  %10d, %12.4f", i, time[membrIndx[i]]);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      else {
        masterTimeIndex[membrIndx[i]] = 0;
        idxFoundFlag = TRUE;
      }
    }
    if (idxFoundFlag == FALSE) {
      if (noIdxAllowFlag == FALSE) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Invalid event time encountered for individual:  %10d, %12.4f", i, time[membrIndx[i]]);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      else {
        masterTimeIndex[membrIndx[i]] = 0;
      }
    }
  }
  if (allMembrIndx == NULL) {
    free_uivector(membrIndx, 1, allMembrSize);
  }
}
void updateEventTypeSubsets(double *summaryStatus, 
                            uint    mRecordSize,
                            int   **mpSign,
                            uint   *mRecordIndex,
                            uint   *meIndividualSize,
                            uint  **eIndividual) {
  uint i, j;
  if (RF_eventTypeSize == 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to update event type subsets in a non-CR analysis.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (j = 1; j <= RF_eventTypeSize; j++) {
    for (i = 1; i <= RF_eIndividualSize[j]; i++) {
      eIndividual[j][i] = RF_eIndividualIn[j][i];
    }
  }
  if (RF_mStatusSize > 0) {
    uint *eventCounter = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eventCounter[j] = RF_eIndividualSize[j]; 
    }
    for (i = 1; i <= mRecordSize; i++) {
      if (mpSign[RF_statusIndex][i] == 1) {
        if ((uint) summaryStatus[mRecordIndex[i]] > 0) {
          j = RF_eventTypeIndex[(uint) summaryStatus[mRecordIndex[i]]];
          eventCounter[j] ++;
          eIndividual[j][eventCounter[j]] = mRecordIndex[i];
        }
        else {
          for (j=1; j <= RF_eventTypeSize; j++) {
            eventCounter[j] ++;
            eIndividual[j][eventCounter[j]] = mRecordIndex[i];
          }
        }
      }
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      meIndividualSize[j] = eventCounter[j];
    }
    free_uivector(eventCounter, 1, RF_eventTypeSize);
  }
}
void stackShadow (uint mode, uint treeID) {
  char vimpShadowFlag;
  uint i, j, p;
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    RF_factorList[treeID] = (Factor **) vvector(1, RF_maxFactorLevel);
    for (j = 1; j <= RF_maxFactorLevel; j++) {
      RF_factorList[treeID][j] = NULL;
    }
    for (j = 1; j <= RF_xFactorCount; j++) {
      if (RF_factorList[treeID][RF_xFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_xFactorSize[j]] = makeFactor(RF_xFactorSize[j], FALSE);
      }
    }
    for (j = 1; j <= RF_rFactorCount; j++) {
      if (RF_factorList[treeID][RF_rFactorSize[j]] == NULL) {
        RF_factorList[treeID][RF_rFactorSize[j]] = makeFactor(RF_rFactorSize[j], FALSE);
      }
    }
  }
  vimpShadowFlag = FALSE;
  if ((RF_opt & OPT_VIMP) && !(RF_opt & OPT_VIMP_TYPE)) {
    vimpShadowFlag = TRUE;
  }
  if (RF_mResponseFlag == TRUE) {
    RF_response[treeID] = (double **) vvector(1, RF_rSize);
    for (p = 1; p <= RF_rSize; p++) {
      RF_response[treeID][p] = RF_responseIn[p];
    }
    for (p = 1; p <= RF_mpIndexSize; p++) {
      if (RF_mpIndex[p] < 0) {
        RF_response[treeID][(uint) abs(RF_mpIndex[p])] = dvector(1, RF_observationSize);
        for (i = 1; i <= RF_observationSize; i++) {
          RF_response[treeID][(uint) abs(RF_mpIndex[p])][i] = RF_responseIn[(uint) abs(RF_mpIndex[p])][i];
        }
      }
      else {
        p = RF_mpIndexSize;
      }
    }
    if (RF_timeIndex > 0) {
      RF_time[treeID] = RF_response[treeID][RF_timeIndex];
      if (RF_mTimeFlag == TRUE) {
        RF_masterTimeIndex[treeID] = uivector(1, RF_observationSize);
        for (i = 1; i <= RF_observationSize; i++) {
          RF_masterTimeIndex[treeID][i] = RF_masterTimeIndexIn[i];
        }
      }
      else {
        RF_masterTimeIndex[treeID] = RF_masterTimeIndexIn;
      }
    }
    if (RF_statusIndex > 0) {
      RF_status[treeID] =  RF_response[treeID][RF_statusIndex];
    }
  }
  if(vimpShadowFlag == TRUE) {
    RF_observation[treeID] = dmatrix(1, RF_xSize, 1, RF_observationSize);
    for (p = 1; p <= RF_xSize; p++) {
      for (i = 1; i <= RF_observationSize; i++) {
        RF_observation[treeID][p][i] = RF_observationIn[p][i];
      }
    }
  }
  else {
    if(RF_mPredictorFlag == TRUE) {
      RF_observation[treeID] = (double **) vvector(1, RF_xSize);
      for (p = 1; p <= RF_xSize; p++) {
        RF_observation[treeID][p] = RF_observationIn[p];
      }
      for (p = 1; p <= RF_mpIndexSize; p++) {
        if (RF_mpIndex[p] > 0) {
          RF_observation[treeID][(uint) RF_mpIndex[p]] = dvector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            RF_observation[treeID][(uint) RF_mpIndex[p]][i] = RF_observationIn[(uint) RF_mpIndex[p]][i];
          }
        }
      }
    }
  }
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_fmResponseFlag == TRUE) {
        RF_fresponse[treeID] = (double **) vvector(1, RF_rSize);
        for (p = 1; p <= RF_frSize; p++) {
          RF_fresponse[treeID][p] = RF_fresponseIn[p];
        }
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] < 0) {
            RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])][i] = RF_fresponseIn[(uint) abs(RF_fmpIndex[p])][i];
            }
          }
          else {
            p = RF_fmpIndexSize;
          }
        }
        if (RF_timeIndex > 0) {
        }
        if (RF_statusIndex > 0) {
        }
      }
    }
    if(vimpShadowFlag == TRUE) {
      RF_fobservation[treeID] = dmatrix(1, RF_xSize, 1, RF_fobservationSize);
      for (p = 1; p <= RF_xSize; p++) {
        for (i = 1; i <= RF_fobservationSize; i++) {
          RF_fobservation[treeID][p][i] = RF_fobservationIn[p][i];
        }
      }
    }
    else {
      if(RF_fmPredictorFlag == TRUE) {
        RF_fobservation[treeID] = (double **) vvector(1, RF_xSize);
        for (p = 1; p <= RF_xSize; p++) {
          RF_fobservation[treeID][p] = RF_fobservationIn[p];
        }
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] > 0) {
            RF_fobservation[treeID][(uint) RF_fmpIndex[p]] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fobservation[treeID][(uint) RF_fmpIndex[p]][i] = RF_fobservationIn[(uint) RF_fmpIndex[p]][i];
            }
          }
        }
      }
    }
  }  
}
void unstackShadow (uint mode, uint treeID) {
  char vimpShadowFlag;
  uint k, p;
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    if (RF_factorList[treeID] != NULL) {
      for (k = 1; k <= RF_maxFactorLevel; k++) {
        if (RF_factorList[treeID][k] != NULL) {
          free_Factor(RF_factorList[treeID][k]);
        }
      }
      free_vvector(RF_factorList[treeID], 1, RF_maxFactorLevel);
      RF_factorList[treeID] = NULL;
    }
  }
  vimpShadowFlag = FALSE;
  if ((RF_opt & OPT_VIMP) && !(RF_opt & OPT_VIMP_TYPE)) {
    vimpShadowFlag = TRUE;
  }
  if (RF_mResponseFlag == TRUE) {
    for (p = 1; p <= RF_mpIndexSize; p++) {
      if (RF_mpIndex[p] < 0) {
        free_dvector(RF_response[treeID][(uint) abs(RF_mpIndex[p])], 1, RF_observationSize);
      }
      else {
        p = RF_mpIndexSize;
      }
    }
    free_vvector(RF_response[treeID], 1, RF_rSize);
    if (RF_timeIndex > 0) {
      if (RF_mTimeFlag == TRUE) {
        free_uivector(RF_masterTimeIndex[treeID], 1, RF_observationSize);
      }
    }
  }
  if(vimpShadowFlag == TRUE) {
    free_dmatrix(RF_observation[treeID], 1, RF_xSize, 1, RF_observationSize);
  }
  else {
    if(RF_mPredictorFlag == TRUE) {
      for (p = 1; p <= RF_mpIndexSize; p++) {
        if (RF_mpIndex[p] > 0) {
          free_vvector(RF_observation[treeID][(uint) RF_mpIndex[p]], 1, RF_observationSize);
        }
      }
      free_vvector(RF_observation[treeID], 1, RF_xSize);
    }
  }
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_fmResponseFlag == TRUE) {
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] < 0) {
            free_vvector(RF_fresponse[treeID][(uint) abs(RF_fmpIndex[p])], 1, RF_fobservationSize);
          }
          else {
            p = RF_fmpIndexSize;
          }
        }
        free_vvector(RF_fresponse[treeID], 1, RF_rSize);
      }
    }
    if(vimpShadowFlag == TRUE) {
      free_dmatrix(RF_fobservation[treeID], 1, RF_xSize, 1, RF_fobservationSize);
    }
    else {
      if(RF_fmPredictorFlag == TRUE) {
        for (p = 1; p <= RF_fmpIndexSize; p++) {
          if (RF_fmpIndex[p] > 0) {
            free_vvector(RF_fobservation[treeID][(uint) RF_fmpIndex[p]], 1, RF_fobservationSize);
          }
        }
        free_vvector(RF_fobservation[treeID], 1, RF_xSize);
      }
    }
  }
}
void imputeUpdateSummaryNew (uint     mode,
                             uint     treeID) {
  Node     *leafNodePtr;
  Terminal *infoNodePtr;
  double **response;
  double **predictor;
  double  *valuePtr;
  int    **mpSign;
  int     *mpIndex;
  uint    *mRecordIndex;
  uint *glmdIndex;
  uint  glmdSize;
  uint *lmiIndex;
  uint  lmiSize;
  uint *lmiRaggedSize;
  int   signedSignatureIndex;
  uint  unsignedSignatureIndex;
  uint  absoluteTargetIndex;
  uint  iterator;
  uint i, j, r, p, t;
  for (t = 1; t <= RF_tLeafCount[treeID]; t++) {
    leafNodePtr = RF_tNodeList[treeID][t];
    infoNodePtr = RF_mTermList[treeID][t];
    if (xferMissingness(mode, leafNodePtr, infoNodePtr)) {
      if ((mode == RF_GROW) || (mode == RF_REST)) {
        mRecordIndex = RF_mRecordIndex;
        mpSign  = RF_mpSign;
        mpIndex = RF_mpIndex;
        response = RF_response[treeID];
        predictor = RF_observation[treeID];
        if (infoNodePtr -> dominant == P_DOMINANT) {
          glmdIndex = leafNodePtr -> lmpIndex;
          glmdSize = leafNodePtr -> lmpIndexActualSize;
        }
        else {
          glmdIndex = leafNodePtr -> lmrIndex;
          glmdSize = leafNodePtr -> lmrIndexActualSize;
        }
      }
      else {
        mRecordIndex = RF_fmRecordIndex;
        mpSign  = RF_fmpSign;
        mpIndex = RF_fmpIndex;
        if (RF_frSize > 0) {
          response = RF_fresponse[treeID];
        } 
        else {
          response = NULL;
        }
        predictor = RF_fobservation[treeID];
        if (infoNodePtr -> dominant == P_DOMINANT) {
          glmdIndex = leafNodePtr -> flmpIndex;
          glmdSize = leafNodePtr -> flmpIndexActualSize;
        }
        else {
          glmdIndex = leafNodePtr -> flmrIndex;
          glmdSize = leafNodePtr -> flmrIndexActualSize;
        }
      }
      lmiSize = infoNodePtr -> lmiSize;
      lmiIndex = infoNodePtr -> lmiIndex;
      lmiRaggedSize = infoNodePtr -> lmiRaggedSize;
      if (infoNodePtr -> dominant == P_DOMINANT) {
        for (r = 1; r <= lmiSize; r++) {
          lmiRaggedSize[r] = 0;
          for (p = 1; p <= glmdSize; p++) {
            signedSignatureIndex = mpIndex[glmdIndex[p]];
            if (signedSignatureIndex < 0) {
              unsignedSignatureIndex = (uint) abs(signedSignatureIndex);
            }
            else { 
              if ((mode == RF_GROW) || (mode == RF_REST)) {
                unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
              }
              else {
                if (RF_frSize > 0) {
                  unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
                }
                else {
                  unsignedSignatureIndex = (uint) signedSignatureIndex;
                }
              }
            }
            if(mpSign[unsignedSignatureIndex][lmiIndex[r]] == 1) {
              lmiRaggedSize[r] ++;
            }
          }
        }
      }
      else {
        for (p = 1; p <= lmiSize; p++) {
          signedSignatureIndex = mpIndex[lmiIndex[p]];
          if (signedSignatureIndex < 0) {
            unsignedSignatureIndex = (uint) abs(signedSignatureIndex);
          }
          else { 
            if ((mode == RF_GROW) || (mode == RF_REST)) {
              unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
            }
            else {
              if (RF_frSize > 0) {
                unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
              }
              else {
                unsignedSignatureIndex = (uint) signedSignatureIndex;
              }
            }
          }
          lmiRaggedSize[p] = 0;
          for (r = 1; r <= glmdSize; r++) {
            if(mpSign[unsignedSignatureIndex][glmdIndex[r]] == 1) {
              lmiRaggedSize[p] ++;
            }
          }
        }
      }
      stackTermLMIRagged(infoNodePtr);
      for (i = 1; i <= infoNodePtr -> lmiSize; i++) {
        for (j = 1; j <= (infoNodePtr -> lmiRaggedSize)[i]; j++) {
          (infoNodePtr -> lmiRaggedValue)[i][j] = NA_REAL;
        }
      }
      if (infoNodePtr -> dominant == P_DOMINANT) {
        for (r = 1; r <= lmiSize; r++) {
          iterator = 0;
          for (p = 1; p <= glmdSize; p++) {
            signedSignatureIndex = mpIndex[glmdIndex[p]];
            if (signedSignatureIndex < 0) {
              unsignedSignatureIndex = (uint) abs(signedSignatureIndex);
              absoluteTargetIndex = (uint) abs(signedSignatureIndex);
              valuePtr = response[absoluteTargetIndex];
            }
            else { 
              if ((mode == RF_GROW) || (mode == RF_REST)) {
                unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
              }
              else {
                if (RF_frSize > 0) {
                  unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
                }
                else {
                  unsignedSignatureIndex = (uint) signedSignatureIndex;
                }
              }
              absoluteTargetIndex = (uint) signedSignatureIndex;
              valuePtr = predictor[absoluteTargetIndex];
            }
            if(mpSign[unsignedSignatureIndex][lmiIndex[r]] == 1) {
              (infoNodePtr -> lmiRaggedIndex)[r][++iterator] = glmdIndex[p];
              (infoNodePtr -> lmiRaggedValue)[r][iterator] = valuePtr[mRecordIndex[lmiIndex[r]]];
            }
          }
        }
      }
      else {
        for (p = 1; p <= lmiSize; p++) {
          iterator = 0;
          signedSignatureIndex = mpIndex[lmiIndex[p]];
          if (signedSignatureIndex < 0) {
            unsignedSignatureIndex = (uint) abs(signedSignatureIndex);
            absoluteTargetIndex = (uint) abs(signedSignatureIndex);
            valuePtr = response[absoluteTargetIndex];
          }
          else { 
            if ((mode == RF_GROW) || (mode == RF_REST)) {
              unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
            }
            else {
              if (RF_frSize > 0) {
                unsignedSignatureIndex = RF_rSize + (uint) signedSignatureIndex;
              }
              else {
                unsignedSignatureIndex = (uint) signedSignatureIndex;
              }
            }
            absoluteTargetIndex = (uint) signedSignatureIndex;
            valuePtr = predictor[absoluteTargetIndex];
          }
          for (r = 1; r <= glmdSize; r++) {
            if(mpSign[unsignedSignatureIndex][glmdIndex[r]] == 1) {
              (infoNodePtr -> lmiRaggedIndex)[p][++iterator] = glmdIndex[r];
              (infoNodePtr -> lmiRaggedValue)[p][iterator] = valuePtr[mRecordIndex[glmdIndex[r]]];
            }
          }
        }
      }
    }  
  }  
}
void imputeCommonNew(uint      mode,
                  uint      serialTreeID,
                  char      selectionFlag,
                  char      predictorFlag) {
  uint *localSerialIndex;
  uint  localSerialCount;
  uint *serialPtr;
  char mFlag;
  char outcomeFlag;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mpIndexSize;
  int    **mpSign;
  int     *mpIndex;
  double **outResponse;
  double **outPredictor;
  double *valuePtr;
  double *naivePtr;
  uint    unsignedSignatureIndex;
  Terminal *info;
  double imputedValue;
  uint localDistributionSize;
  uint maxDistributionSize;
  uint rspSize;
  char result;
  uint i, p, r, v, tree;
  valuePtr      = NULL;  
  naivePtr      = NULL;  
  unsignedSignatureIndex = 0;     
  maxDistributionSize = 0;  
  outResponse         = 0;  
  outPredictor        = 0;  
  rspSize = 0;  
  mpIndex = 0;  
  mpSign  = 0;  
  mpIndexSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  if ((selectionFlag != TRUE) && (selectionFlag != FALSE) && (selectionFlag != ACTIVE)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid selectionFlag in imputeCommon():  %10d", selectionFlag);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mpIndexSize = RF_fmpIndexSize;
      mpSign = RF_fmpSign;
      mpIndex = RF_fmpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      outResponse  = RF_sImputeResponsePtr;
      outPredictor = RF_sImputePredictorPtr;
      rspSize = RF_frSize;
      result = TRUE;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mpIndexSize = RF_mpIndexSize;
      mpSign = RF_mpSign;
      mpIndex = RF_mpIndex;
      maxDistributionSize = ((RF_observationSize) > (RF_forestSize)) ? (RF_observationSize) : (RF_forestSize);
      if ((selectionFlag == TRUE) || (selectionFlag == ACTIVE)) {
        outResponse  = RF_sImputeResponsePtr;
        outPredictor = RF_sImputePredictorPtr;
      }
      else {
        outResponse  = RF_sOOBImputeResponsePtr;
        outPredictor = RF_sOOBImputePredictorPtr;
      }
      rspSize = RF_rSize;
      result = TRUE;
    }
    break;
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to impute in imputeCommon() with no missingness in mode:  %10d", mode);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (serialTreeID == 0) {
    localSerialIndex = uivector(1, RF_forestSize);
    for (tree = 1; tree <= RF_forestSize; tree++) {
      localSerialIndex[tree] = tree;
    }
    serialPtr = localSerialIndex;
    localSerialCount = RF_forestSize;
  }
  else {
    serialPtr = RF_serialTreeIndex;
    localSerialCount = serialTreeID;
  }
  imputedValue = 0.0;  
  double *localDistribution = dvector(1, maxDistributionSize);
  char  *naiveFlag = cvector(1, mpIndexSize);
  char **naiveSign = cmatrix(1, mRecordSize, 1, mpIndexSize);
  for (p = 1; p <= mpIndexSize; p++) {
    naiveFlag[p] = FALSE;
  }
  for (i = 1; i <= mRecordSize; i++) {
    outcomeFlag = TRUE;
    for (p = 1; p <= mpIndexSize; p++) {
      naiveSign[i][p] = FALSE;  
      if (mpIndex[p] < 0) {
        unsignedSignatureIndex = (uint) abs(mpIndex[p]);
      }
      else {
        if (predictorFlag == TRUE) {
          unsignedSignatureIndex = (uint) mpIndex[p] + rspSize;
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || predictorFlag) {
        if (mpSign[unsignedSignatureIndex][i] == 1) {
          localDistributionSize = 0;
          for (tree = 1; tree <= localSerialCount; tree++) {
            if (RF_tLeafCount[serialPtr[tree]] > 0) {
              if ((RF_dmRecordBootFlag[serialPtr[tree]][i] == selectionFlag) || (selectionFlag == ACTIVE)) {
                info = RF_mTermMembership[serialPtr[tree]][i];
                if (info -> dominant == P_DOMINANT) {
                  for (r = 1; r <= info -> lmiSize; r++) {
                    if ((info -> lmiIndex)[r] == i) {
                      for (v = 1; v <= info -> lmiRaggedSize[r]; v++) {
                        if ((info -> lmiRaggedIndex)[r][v] == p) {
                          if (!ISNA((info -> lmiRaggedValue)[r][v])) {
                            localDistribution[++localDistributionSize] = (info -> lmiRaggedValue)[r][v];
                          }
                          else {
                          }  
                          v = (info -> lmiRaggedSize)[r];
                        }
                      }
                      r = info -> lmiSize;
                    }
                  } 
                }
                else {
                  for (v = 1; v <= info -> lmiSize; v++) {
                    if ((info -> lmiIndex)[v] == p) {
                      for (r = 1; r <= info -> lmiRaggedSize[v]; r++) {
                        if ((info -> lmiRaggedIndex)[v][r] == i) {
                          if (!ISNA((info -> lmiRaggedValue)[v][r])) {
                            localDistribution[++localDistributionSize] = (info -> lmiRaggedValue)[v][r];
                          }
                          else {
                          }  
                          r = (info -> lmiRaggedSize)[v];
                        }
                      }
                      v = info -> lmiSize;
                    }
                  } 
                }
              }  
            }  
            else {
            }
          }  
          if (localDistributionSize > 0) {
            if (mpIndex[p] < 0) {
              if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "T") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "S") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "I") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              else if (strcmp(RF_rType[(uint) abs(mpIndex[p])], "C") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              outResponse[(uint) abs(mpIndex[p])][i] = imputedValue;
            }  
            else {
              if (strcmp(RF_xType[(uint) mpIndex[p]], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, localSerialCount);
              }
              outPredictor[(uint) mpIndex[p]][i] = imputedValue;
            }
          }  
          else {
            naiveFlag[p] = TRUE;
            naiveSign[i][p] = TRUE;
          }
        }  
      }  
      else {
        p = mpIndexSize;
      }
    }  
  }  
  outcomeFlag = TRUE;
  for (p = 1; p <= mpIndexSize; p++) {
    if (mpIndex[p] < 0) {
      unsignedSignatureIndex = (uint) abs(mpIndex[p]);
      valuePtr = RF_responseIn[(uint) abs(mpIndex[p])];
      naivePtr = outResponse[(uint) abs(mpIndex[p])];
    }
    else {
      if (predictorFlag == TRUE) {
        unsignedSignatureIndex = (uint) mpIndex[p] + rspSize;
        valuePtr = RF_observationIn[(uint) mpIndex[p]];
        naivePtr = outPredictor[(uint) mpIndex[p]];
      }
      outcomeFlag = FALSE;
    }
    if (outcomeFlag || predictorFlag) {
      if (naiveFlag[p] == TRUE) {
        localDistributionSize = 0;
        for (i=1; i <= RF_observationSize; i++) {
          mFlag = TRUE;
          if (RF_mRecordMap[i] == 0) {
            mFlag = FALSE;
          }
          else if (RF_mpSign[unsignedSignatureIndex][RF_mRecordMap[i]] == 0) {
            mFlag = FALSE;
          }
          if (mFlag == FALSE) {
            localDistribution[++localDistributionSize] = valuePtr[i];
          }
        }  
        if (localDistributionSize > 0) {
          for (i=1; i <= mRecordSize; i++) {
            if (naiveSign[i][p] == TRUE) {
              naivePtr[i] = getSampleValue(localDistribution, localDistributionSize, FALSE, localSerialCount);
            }
          }
        }  
        else {
          if (mpIndex[p] < 0) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Naive imputation failed for [indv, outcome] = [%10d, %10d] \n", mRecordIndex[i], mpIndex[p]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
          else {
          }
        }
      }  
    }  
    else {
      p = mpIndexSize;
    }
  }  
  if (serialTreeID == 0) {
    free_uivector(localSerialIndex, 1, RF_forestSize);
  }
  free_dvector(localDistribution, 1, maxDistributionSize);
  free_cvector(naiveFlag, 1, mpIndexSize);
  free_cmatrix(naiveSign, 1, mRecordSize, 1, mpIndexSize);
}
char xferMissingness(uint mode, Node *source, Terminal *destination) {
  uint *sourcePtr;
  uint sourceLen;
  uint p;
  char result;
  char xferFlag;
  sourcePtr = NULL;  
  sourceLen = 0;     
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    if (RF_mRecordSize > 0) {
      result = TRUE;
        destination -> dominant = N_DOMINANT;
        sourcePtr = source -> lmpIndex;
        sourceLen = source -> lmpIndexActualSize;
    }
  }
  else {
    if (RF_fmRecordSize > 0) {
      result = TRUE;
        destination -> dominant = N_DOMINANT;
        sourcePtr = source -> flmpIndex;
        sourceLen = source -> flmpIndexActualSize;
    }
  }
  if (result == FALSE) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to update forest impute data with no missingness in mode:  %10d", mode);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (sourceLen > 0) {
    stackTermLMIIndex(destination, sourceLen);
    for (p = 1; p <= sourceLen; p++) {
      (destination -> lmiIndex)[p] = sourcePtr[p]; 
    } 
    xferFlag = TRUE;
  }
  else {
    xferFlag = FALSE;
  }
  return xferFlag;
}
