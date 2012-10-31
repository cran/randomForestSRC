////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.0
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
////    URL:    http://www.kogalur-shear.com
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
#include      "treeUtil.h"
#include     "factorOps.h"
#include        "impute.h"
char imputeNode (uint     type,
                 char     chainFlag,
                 uint     treeID, 
                 Node    *parent,
                 uint    *repMembrIndx,
                 uint     repMembrSize,
                 uint    *allMembrIndx,
                 uint     allMembrSize) {
  double  **response;
  double  **predictor;
  uint    *mRecordMap;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int     *mvNSptr;
  double *valuePtr, *imputePtr;
  char mPredictorFlag;
  uint unsignedIndexSource, unsignedIndexTarget;
  char result;
  uint i,p;
  uint localDistributionSize;
  mvNSptr = 0;  
  mvIndex = 0;  
  mvSign  = 0;  
  mvSize  = 0;  
  mRecordMap = NULL;  
  predictor  = NULL;  
  response   = NULL;  
  result = FALSE;
  switch (type) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      if (RF_frSize > 0) {
        response = RF_fresponse[treeID];
      } 
      else {
        response = NULL;
      }
      predictor = RF_fobservation[treeID];
      mRecordMap = RF_fmRecordMap;
      mvSize = RF_fmvSize;
      mvSign = RF_fmvSign;
      mvIndex = RF_fmvIndex;
      mvNSptr = parent -> fmvSign;
      result = TRUE;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      response = RF_response[treeID];
      predictor = RF_observation[treeID];
      mRecordMap = RF_mRecordMap;
      mvSize = RF_mvSize;
      mvSign = RF_mvSign;
      mvIndex = RF_mvIndex;
      mvNSptr = parent -> mvSign;
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
  for (p = 1; p <= mvSize; p++) {
    if (mvNSptr[p] != -1) {
      if (mvIndex[p] < 0) {
        unsignedIndexSource = unsignedIndexTarget = (uint) abs(mvIndex[p]);
        valuePtr = RF_response[treeID][(uint) abs(mvIndex[p])];
        imputePtr = response[(uint) abs(mvIndex[p])];
      }
      else { 
        unsignedIndexSource = RF_rSize + (uint) mvIndex[p];
        if (type == RF_PRED) {
          if (RF_frSize > 0) {
            unsignedIndexTarget = RF_rSize + (uint) mvIndex[p];
          }
          else {
            unsignedIndexTarget = (uint) mvIndex[p];
          }
        }
        else {
          unsignedIndexTarget = RF_rSize + (uint) mvIndex[p];
        }
        valuePtr = RF_observation[treeID][(uint) mvIndex[p]];
        imputePtr = predictor[(uint) mvIndex[p]];
      }
      localDistributionSize = 0;
      for (i = 1; i <= repMembrSize; i++) {
          mPredictorFlag = TRUE;
          if (RF_mRecordMap[repMembrIndx[i]] == 0) {
            mPredictorFlag = FALSE;
          }
          else if (RF_mvSign[unsignedIndexSource][RF_mRecordMap[repMembrIndx[i]]] == 0) {
            mPredictorFlag = FALSE;
          }
          if (mPredictorFlag == FALSE) {
            localDistributionSize ++;
            localDistribution[localDistributionSize] = valuePtr[repMembrIndx[i]];
          }
      }  
      if (localDistributionSize > 0) {
        for (i = 1; i <= allMembrSize; i++) {
          if (mRecordMap[allMembrIndx[i]] > 0) {
            if(mvSign[unsignedIndexTarget][mRecordMap[allMembrIndx[i]]] == 1) {
              imputePtr[allMembrIndx[i]] = getSampleValue(localDistribution, localDistributionSize, chainFlag, treeID);
            }  
          }
        }  
      }  
      else {
      }
    }  
  }  
  free_dvector(localDistribution, 1, repMembrSize + 1);
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
    }
  }
  else {
    bootMembrIndx = repMembrIndx;
    bootMembrSize = repMembrSize;
    parent -> mvSign = (parent -> parent) -> mvSign;
    parent -> fmvSign = (parent -> parent) -> fmvSign;
  }
  if (bootResult) {
    if (RF_mRecordSize > 0) {
      imputeNode(RF_GROW,
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
          RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> left;
        }
        else {
          membershipIndicator[allMembrIndx[i]] = RIGHT;
          rghtAllMembrSize ++;
          RF_nodeMembership[treeID][allMembrIndx[i]] = parent -> right;
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
            RF_fnodeMembership[treeID][ngAllMembrIndx[i]] = parent -> left;
          }
          else {
            ngMembershipIndicator[ngAllMembrIndx[i]] = RIGHT;
            ngRghtAllMembrSize ++;
            RF_fnodeMembership[treeID][ngAllMembrIndx[i]] = parent -> right;
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
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int     *mvNodeSign;
  double  *valuePtr;
  uint     unsignedIndex;
  char result;
  uint rspSize;
  uint i, p;
  rspSize = 0;  
  mvIndex = 0;  
  mvSign  = 0;  
  mvSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  result = FALSE;
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mvSize = RF_mvSize;
      mvSign = RF_mvSign;
      mvIndex = RF_mvIndex;
      rspSize = RF_rSize;
      result = TRUE;
    }
  }
  else {
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mvSize = RF_fmvSize;
      mvSign = RF_fmvSign;
      mvIndex = RF_fmvIndex;
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
  for (p = 1; p <= mvSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      if ((mode == RF_GROW) || (mode == RF_REST)) {
          mvNodeSign = RF_nodeMembership[treeID][mRecordIndex[i]] -> mvSign;
      }
      else {
          mvNodeSign = RF_fnodeMembership[treeID][mRecordIndex[i]] -> fmvSign;
      }
      if (mvNodeSign[p] != -1) {
        if (mvIndex[p] < 0) {
          unsignedIndex = (uint) abs(mvIndex[p]);
          valuePtr = responsePtr[(uint) abs(mvIndex[p])];
        }
        else {
          unsignedIndex = (uint) mvIndex[p] + rspSize;
          valuePtr = predictorPtr[(uint) mvIndex[p]];
        }
        if (mvSign[unsignedIndex][i] == 1) {
          if (ISNA(valuePtr[mRecordIndex[i]])) {
            Rprintf("\nDiagnostic Trace of Shadowed Data:  ");
            Rprintf("\n       index   imputation -> \n");
            Rprintf(  "            ");
            for (p=1; p <= mvSize; p++) {
              Rprintf(" %12d", mvIndex[p]);
            }
            Rprintf("\n");
            for (i = 1; i <= mRecordSize; i++) {
              Rprintf("%12d", mRecordIndex[i]);
              for (p = 1; p <= mvSize; p++) {
                if (mvIndex[p] < 0) {
                  valuePtr = responsePtr[(uint) abs(mvIndex[p])];
                }
                else {
                  valuePtr = predictorPtr[(uint) mvIndex[p]];
                }
                Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
              }
              Rprintf("\n");
            }
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Attempt to update forest impute data with invalid shadowed value, NA. ");
            Rprintf("\nRF-SRC:  Invalid value for:  [indv][outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mvIndex[p]);
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
void imputeUpdateSummaryNew (uint     mode, 
                          double **responsePtr, 
                          double **predictorPtr, 
                          uint     treeID) {
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
  int     *mvNodeSign;
  double  *valuePtr;
  uint     unsignedIndex;
  char result;
  uint rspSize;
  uint i, p;
  rspSize = 0;  
  mvIndex = 0;  
  mvSign  = 0;  
  mvSize  = 0;  
  mRecordIndex = 0;  
  mRecordSize  = 0;  
  result = FALSE;
  if ((mode == RF_GROW) || (mode == RF_REST)) {
    if (RF_mRecordSize > 0) {
      mRecordSize = RF_mRecordSize;
      mRecordIndex = RF_mRecordIndex;
      mvSize = RF_mvSize;
      mvSign = RF_mvSign;
      mvIndex = RF_mvIndex;
      rspSize = RF_rSize;
      result = TRUE;
    }
  }
  else {
    if (RF_fmRecordSize > 0) {
      mRecordSize = RF_fmRecordSize;
      mRecordIndex = RF_fmRecordIndex;
      mvSize = RF_fmvSize;
      mvSign = RF_fmvSign;
      mvIndex = RF_fmvIndex;
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
  for (p = 1; p <= mvSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      if ((mode == RF_GROW) || (mode == RF_REST)) {
          mvNodeSign = RF_nodeMembership[treeID][mRecordIndex[i]] -> mvSign;
      }
      else {
          mvNodeSign = RF_fnodeMembership[treeID][mRecordIndex[i]] -> fmvSign;
      }
      if (mvNodeSign[p] != -1) {
        if (mvIndex[p] < 0) {
          unsignedIndex = (uint) abs(mvIndex[p]);
          valuePtr = responsePtr[(uint) abs(mvIndex[p])];
        }
        else {
          unsignedIndex = (uint) mvIndex[p] + rspSize;
          valuePtr = predictorPtr[(uint) mvIndex[p]];
        }
        if (mvSign[unsignedIndex][i] == 1) {
          if (ISNA(valuePtr[mRecordIndex[i]])) {
            Rprintf("\nDiagnostic Trace of Shadowed Data:  ");
            Rprintf("\n       index   imputation -> \n");
            Rprintf(  "            ");
            for (p=1; p <= mvSize; p++) {
              Rprintf(" %12d", mvIndex[p]);
            }
            Rprintf("\n");
            for (i = 1; i <= mRecordSize; i++) {
              Rprintf("%12d", mRecordIndex[i]);
              for (p = 1; p <= mvSize; p++) {
                if (mvIndex[p] < 0) {
                  valuePtr = responsePtr[(uint) abs(mvIndex[p])];
                }
                else {
                  valuePtr = predictorPtr[(uint) mvIndex[p]];
                }
                Rprintf(" %12.4f", valuePtr[mRecordIndex[i]]);
              }
              Rprintf("\n");
            }
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Attempt to update forest impute data with invalid shadowed value, NA. ");
            Rprintf("\nRF-SRC:  Invalid value for:  [indv][outcome/predictor] = [%10d][%10d] ", mRecordIndex[i], mvIndex[p]);
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
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
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
  mvSize       = 0;     
  mvSign       = NULL;  
  mvIndex      = NULL;  
  outResponse  = NULL;  
  outPredictor = NULL;  
  valuePtr     = NULL;  
  outputPtr    = NULL;  
  unsignedIndex = 0;    
  switch (mode) {
  case RF_PRED:
    mRecordSize = RF_fmRecordSize;
    mRecordIndex = RF_fmRecordIndex;
    mvSize = RF_fmvSize;
    mvSign = RF_fmvSign;
    mvIndex = RF_fmvIndex;
    outResponse  = RF_sImputeResponsePtr;
    outPredictor = RF_sImputePredictorPtr;
    rspSize = RF_frSize;
    break;
  default:
    mRecordSize = RF_mRecordSize;
    mRecordIndex = RF_mRecordIndex;
    mvSize = RF_mvSize;
    mvSign = RF_mvSign;
    mvIndex = RF_mvIndex;
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
  for (p = 1; p <= mvSize; p++) {
    for (i = 1; i <= mRecordSize; i++) {
      if (mvIndex[p] < 0) {
        unsignedIndex = (uint) abs(mvIndex[p]);
        valuePtr = shadowResponse[(uint) abs(mvIndex[p])];
        outputPtr = outResponse[(uint) abs(mvIndex[p])];
      }
      else {
        if (shadowPredictor != NULL) {
          unsignedIndex = (uint) mvIndex[p] + rspSize;
          valuePtr = shadowPredictor[(uint) mvIndex[p]];
          outputPtr = outPredictor[(uint) mvIndex[p]];
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || (shadowPredictor != NULL)) {
        if (mvSign[unsignedIndex][i] == 1) {
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
  imputeCommon(mode,
               RF_forestSize,
               selectionFlag, 
               TRUE);
}
void imputeResponse(uint      mode,
                    uint      serialTreeID,
                    double  **tempResponse) {
  switch(mode) {
  case RF_PRED:
    imputeCommon(mode, serialTreeID, ACTIVE, FALSE);
    imputeUpdateShadow(mode, ACTIVE, tempResponse, NULL);
    break;
  default:
    imputeCommon(mode, serialTreeID, FALSE, FALSE);
    imputeUpdateShadow(mode, FALSE, tempResponse, NULL);
    break;
  }
}
void imputeCommon(uint      mode,
                  uint      serialTreeID,
                  char      selectionFlag,
                  char      predictorFlag) {
  uint i,p,tree;
  char mFlag;
  char outcomeFlag;
  uint     mRecordSize;
  uint    *mRecordIndex;
  uint     mvSize;
  int    **mvSign;
  int     *mvIndex;
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
  valuePtr      = NULL;  
  naivePtr      = NULL;  
  unsignedIndex = 0;     
  maxDistributionSize = 0;  
  outResponse         = 0;  
  outPredictor        = 0;  
  rspSize = 0;  
  mvIndex = 0;  
  mvSign  = 0;  
  mvSize  = 0;  
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
      mvSize = RF_fmvSize;
      mvSign = RF_fmvSign;
      mvIndex = RF_fmvIndex;
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
      mvSize = RF_mvSize;
      mvSign = RF_mvSign;
      mvIndex = RF_mvIndex;
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
  imputedValue = 0.0;  
  double *localDistribution = dvector(1, maxDistributionSize);
  char  *naiveMvFlag = cvector(1, mvSize);
  char **naiveSign = cmatrix(1, mRecordSize, 1, mvSize);
  for (p = 1; p <= mvSize; p++) {
    naiveMvFlag[p] = FALSE;
  }
  for (i = 1; i <= mRecordSize; i++) {
    outcomeFlag = TRUE;
    for (p = 1; p <= mvSize; p++) {
      naiveSign[i][p] = FALSE;  
      if (mvIndex[p] < 0) {
        unsignedIndex = (uint) abs(mvIndex[p]);
      }
      else {
        if (predictorFlag == TRUE) {
          unsignedIndex = (uint) mvIndex[p] + rspSize;
        }
        outcomeFlag = FALSE;
      }
      if (outcomeFlag || predictorFlag) {
        if (mvSign[unsignedIndex][i] == 1) {
          localDistributionSize = 0;
          for (tree = 1; tree <= serialTreeID; tree++) {
            if (RF_leafCount[RF_serialTreeIndex[tree]] > 0) {
              if ((RF_dmRecordBootFlag[RF_serialTreeIndex[tree]][i] == selectionFlag) || (selectionFlag == ACTIVE)) {
                if (!ISNA(RF_dmvImputation[RF_serialTreeIndex[tree]][i][p])) {
                  localDistribution[++localDistributionSize] = RF_dmvImputation[RF_serialTreeIndex[tree]][i][p];
                }
                else {
                }  
              }  
            }  
          }  
          if (localDistributionSize > 0) {
            if (mvIndex[p] < 0) {
              if (strcmp(RF_rType[(uint) abs(mvIndex[p])], "T") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mvIndex[p])], "S") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, serialTreeID);
              }
              else if (strcmp(RF_rType[(uint) abs(mvIndex[p])], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else if (strcmp(RF_rType[(uint) abs(mvIndex[p])], "I") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, serialTreeID);
              }
              else if (strcmp(RF_rType[(uint) abs(mvIndex[p])], "C") == 0) {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, serialTreeID);
              }
              outResponse[(uint) abs(mvIndex[p])][i] = imputedValue;
            }  
            else {
              if (strcmp(RF_xType[(uint) mvIndex[p]], "R") == 0) {
                imputedValue = getMeanValue(localDistribution, localDistributionSize);
              }
              else {
                imputedValue = getMaximalValue(localDistribution, localDistributionSize, serialTreeID);
              }
              outPredictor[(uint) mvIndex[p]][i] = imputedValue;
            }
          }  
          else {
            naiveMvFlag[p] = TRUE;
            naiveSign[i][p] = TRUE;
          }
        }  
      }  
      else {
        p = mvSize;
      }
    }  
  }  
  outcomeFlag = TRUE;
  for (p = 1; p <= mvSize; p++) {
    if (mvIndex[p] < 0) {
      unsignedIndex = (uint) abs(mvIndex[p]);
      valuePtr = RF_responseIn[(uint) abs(mvIndex[p])];
      naivePtr = outResponse[(uint) abs(mvIndex[p])];
    }
    else {
      if (predictorFlag == TRUE) {
        unsignedIndex = (uint) mvIndex[p] + rspSize;
        valuePtr = RF_observationIn[(uint) mvIndex[p]];
        naivePtr = outPredictor[(uint) mvIndex[p]];
      }
      outcomeFlag = FALSE;
    }
    if (outcomeFlag || predictorFlag) {
      if (naiveMvFlag[p] == TRUE) {
        localDistributionSize = 0;
        for (i=1; i <= RF_observationSize; i++) {
          mFlag = TRUE;
          if (RF_mRecordMap[i] == 0) {
            mFlag = FALSE;
          }
          else if (RF_mvSign[unsignedIndex][RF_mRecordMap[i]] == 0) {
            mFlag = FALSE;
          }
          if (mFlag == FALSE) {
            localDistribution[++localDistributionSize] = valuePtr[i];
          }
        }  
        if (localDistributionSize > 0) {
          for (i=1; i <= mRecordSize; i++) {
            if (naiveSign[i][p] == TRUE) {
              naivePtr[i] = getSampleValue(localDistribution, localDistributionSize, FALSE, serialTreeID);
            }
          }
        }  
        else {
          if (mvIndex[p] < 0) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Naive imputation failed for [indv, outcome] = [%10d, %10d] \n", mRecordIndex[i], mvIndex[p]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
          else {
          }
        }
      }  
    }  
    else {
      p = mvSize;
    }
  }  
  free_dvector(localDistribution, 1, maxDistributionSize);
  free_cvector(naiveMvFlag, 1, mvSize);
  free_cmatrix(naiveSign, 1, mRecordSize, 1, mvSize);
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
    if(RF_mvSign[RF_timeIndex][i] == 1) { 
      meanValue = outTime[i];
      if ((meanValue < RF_masterTime[1]) || (meanValue > RF_masterTime[RF_masterTimeSize])) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  The summary mean value for time is out of range:  indv %10d, value %12.4f", RF_mRecordIndex[i], meanValue);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      leftDistance = rightDistance = minimumIndex = 0;
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
          minimumIndex = j-1;            
        }
        else {
          if (fabs(leftDistance - rightDistance) < EPSILON) {
            if (ran2(1) <= 0.5) {
              minimumIndex = j-1;
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
                            int   **mvSign,
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
      if (mvSign[RF_statusIndex][i] == 1) {
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
    for (p = 1; p <= RF_mvSize; p++) {
      if (RF_mvIndex[p] < 0) {
        RF_response[treeID][(uint) abs(RF_mvIndex[p])] = dvector(1, RF_observationSize);
        for (i = 1; i <= RF_observationSize; i++) {
          RF_response[treeID][(uint) abs(RF_mvIndex[p])][i] = RF_responseIn[(uint) abs(RF_mvIndex[p])][i];
        }
      }
      else {
        p = RF_mvSize;
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
      for (p = 1; p <= RF_mvSize; p++) {
        if (RF_mvIndex[p] > 0) {
          RF_observation[treeID][(uint) RF_mvIndex[p]] = dvector(1, RF_observationSize);
          for (i = 1; i <= RF_observationSize; i++) {
            RF_observation[treeID][(uint) RF_mvIndex[p]][i] = RF_observationIn[(uint) RF_mvIndex[p]][i];
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
        for (p = 1; p <= RF_fmvSize; p++) {
          if (RF_fmvIndex[p] < 0) {
            RF_fresponse[treeID][(uint) abs(RF_fmvIndex[p])] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fresponse[treeID][(uint) abs(RF_fmvIndex[p])][i] = RF_fresponseIn[(uint) abs(RF_fmvIndex[p])][i];
            }
          }
          else {
            p = RF_fmvSize;
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
        for (p = 1; p <= RF_fmvSize; p++) {
          if (RF_fmvIndex[p] > 0) {
            RF_fobservation[treeID][(uint) RF_fmvIndex[p]] = dvector(1, RF_fobservationSize);
            for (i = 1; i <= RF_fobservationSize; i++) {
              RF_fobservation[treeID][(uint) RF_fmvIndex[p]][i] = RF_fobservationIn[(uint) RF_fmvIndex[p]][i];
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
    for (p = 1; p <= RF_mvSize; p++) {
      if (RF_mvIndex[p] < 0) {
        free_dvector(RF_response[treeID][(uint) abs(RF_mvIndex[p])], 1, RF_observationSize);
      }
      else {
        p = RF_mvSize;
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
      for (p = 1; p <= RF_mvSize; p++) {
        if (RF_mvIndex[p] > 0) {
          free_vvector(RF_observation[treeID][(uint) RF_mvIndex[p]], 1, RF_observationSize);
        }
      }
      free_vvector(RF_observation[treeID], 1, RF_xSize);
    }
  }
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_fmResponseFlag == TRUE) {
        for (p = 1; p <= RF_fmvSize; p++) {
          if (RF_fmvIndex[p] < 0) {
            free_vvector(RF_fresponse[treeID][(uint) abs(RF_fmvIndex[p])], 1, RF_fobservationSize);
          }
          else {
            p = RF_fmvSize;
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
        for (p = 1; p <= RF_fmvSize; p++) {
          if (RF_fmvIndex[p] > 0) {
            free_vvector(RF_fobservation[treeID][(uint) RF_fmvIndex[p]], 1, RF_fobservationSize);
          }
        }
        free_vvector(RF_fobservation[treeID], 1, RF_xSize);
      }
    }
  }
}
