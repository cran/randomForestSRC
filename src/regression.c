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


#include          "global.h"
#include          "extern.h"
#include           "trace.h"
#include          "nrutil.h"
#include  "regression.h"
void getMeanResponse(uint treeID) {
  Node *parent;
  uint leaf, i;
  double sumResponse;
  uint *membershipIndex;
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    membershipIndex = RF_identityMembershipIndex;
  }
  else {
    membershipIndex = RF_bootMembershipIndex[treeID];
  }
  for (leaf = 1; leaf <= RF_tLeafCount[treeID]; leaf++) {
    parent = RF_tNodeList[treeID][leaf];
    parent -> membrCount = 0;
    sumResponse = 0.0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_tNodeMembership[treeID][membershipIndex[i]] == parent) {
        sumResponse += RF_response[treeID][RF_rTarget][membershipIndex[i]];
        parent -> membrCount ++;
      }
    }
    if (parent -> membrCount > 0) {
      parent -> predictedOutcome =  sumResponse / (double) (parent -> membrCount);
    }
    else {
      parent -> predictedOutcome =  NA_REAL;
      if (!(RF_opt & OPT_OUTC_TYPE)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Zero regression count encountered in node during getMeanResponse():  %10d %10d", treeID, leaf);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }  
}
void updateEnsembleMean(uint     mode,
                        uint     treeID,
                        double  *ensembleOutcome) {
  uint  obsSize;
  unsigned char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  double ***ensemblePtr;
  Node   ***nodeMembershipPtr;
  uint     *ensembleDenPtr;
  Node* parent;
  uint i;
  ensemblePtr    = NULL;  
  ensembleDenPtr = NULL;  
  oobFlag = fullFlag = FALSE;
  outcomeFlag = TRUE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    nodeMembershipPtr = RF_ftNodeMembership;
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    nodeMembershipPtr = RF_tNodeMembership;
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensemblePtr = RF_oobEnsemblePtr;
      ensembleDenPtr = RF_oobEnsembleDen;
    }
    else {
      if (fullFlag == TRUE) {
        ensemblePtr = RF_fullEnsemblePtr;
        ensembleDenPtr = RF_fullEnsembleDen;
      }
    }
    for (i = 1; i <= obsSize; i++) {
      selectionFlag = TRUE;
      if (oobFlag == TRUE) {
        if (RF_bootMembershipFlag[treeID][i] == FALSE) {
          selectionFlag = TRUE;
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        parent = nodeMembershipPtr[treeID][i];
        if (RF_opt & OPT_OUTC_TYPE) {
          if (!ISNA(parent -> predictedOutcome)) {
          }
          else {
            selectionFlag = FALSE;
          }
        }
      }
      if (selectionFlag) {
        ensemblePtr[1][1][i] += parent -> predictedOutcome;
        ensembleDenPtr[i] ++;
      }
      if (outcomeFlag == TRUE) {
        if (ensembleDenPtr[i] != 0) {
          ensembleOutcome[i] = ensemblePtr[1][1][i] / ensembleDenPtr[i];
        }
        else {
          ensembleOutcome[i] = NA_REAL;
        }
      }  
    }  
    if (outcomeFlag == TRUE) {
      outcomeFlag = FALSE;
    }
    if (oobFlag == TRUE) {
      oobFlag = FALSE;
    }
    else {
      if (fullFlag == TRUE) {
        fullFlag = FALSE;
      }
    }
  }  
}
double getMeanSquareError(uint    size,
                          double *responsePtr,
                          double *predictedOutcome,
                          uint   *denomCount) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i=1; i <= size; i++) {
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      result += pow (responsePtr[i] - predictedOutcome[i], 2.0);
    }
  }  
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = result / (double) cumDenomCount;
  }
  return result;
}
char getVariance(uint    repMembrSize,
                 uint   *repMembrIndx,
                 uint    nonMissMembrSize,
                 uint   *nonMissMembrIndx,
                 double *targetResponse,
                 double *mean,
                 double *variance) {
  uint i;
  uint denom;
  double meanResult, varResult;
  char result;
  uint *genIndx;
  uint  genSize;
  if ((nonMissMembrSize == 0) || (nonMissMembrIndx == NULL)) {
    genIndx = RF_identityMembershipIndex;
    genSize = repMembrSize;
  }
  else {
    genIndx = nonMissMembrIndx;
    genSize = nonMissMembrSize;
  }
  denom      = 0;
  meanResult = 0.0;
  for (i = 1; i <= genSize; i++) {
    if(!ISNA(targetResponse[repMembrIndx[genIndx[i]]])) {
      denom ++;
      meanResult += targetResponse[repMembrIndx[genIndx[i]]];
    }
  }
  if (denom > 0) {
    meanResult = meanResult / (double) denom;
  }
  else {
    meanResult = NA_REAL;
  }
  if (mean != NULL) {
    *mean = meanResult;
  }
  varResult = 0.0;
  if(!ISNA(meanResult)) {
    for (i = 1; i <= genSize; i++) {
      if(!ISNA(targetResponse[repMembrIndx[genIndx[i]]])) {
        varResult += pow(meanResult - targetResponse[repMembrIndx[genIndx[i]]], 2.0);
      }
    }
    varResult = varResult / (double) denom;
    result = ((varResult <= EPSILON) ? FALSE : TRUE);
  }
  else {
    varResult = NA_REAL;
    result = FALSE;
  }
  if (variance != NULL) {
    *variance = varResult;
  }
  return(result);
}
