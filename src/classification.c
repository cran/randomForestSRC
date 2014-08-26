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


#include          "global.h"
#include          "extern.h"
#include           "trace.h"
#include          "nrutil.h"
#include         "nodeOps.h"
#include  "classification.h"
void getMultiClassProb (uint treeID) {
  Node *parent;
  double maxValue, maxClass;
  uint leaf, i, j, k;
  uint *membershipIndex;
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    membershipIndex = RF_identityMembershipIndex;
  }
  else {
    membershipIndex = RF_bootMembershipIndex[treeID];
  }
  for (leaf = 1; leaf <= RF_tLeafCount[treeID]; leaf++) {
    parent = RF_tNodeList[treeID][leaf];
    stackMultiClassProb(parent, RF_rFactorCount, RF_rFactorSize);
    for (j=1; j <= RF_rFactorCount; j++) {
      for (k=1; k <= RF_rFactorSize[j]; k++) {
        (parent -> multiClassProb)[j][k] = 0;
      }
    }
    parent -> membrCount = 0;
    for (i = 1; i <= RF_observationSize; i++) {
      if (RF_tNodeMembership[treeID][membershipIndex[i]] == parent) {
        for (j=1; j <= RF_rFactorCount; j++) {
          (parent -> multiClassProb)[j][(uint) RF_response[treeID][RF_rFactorIndex[j]][membershipIndex[i]]] ++;
        }
        parent -> membrCount ++;
      }
    }
    if ((parent -> membrCount) > 0) {
      maxValue = 0;
      maxClass = 0;
      for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; k++) {
        if (maxValue < (parent -> multiClassProb[RF_rFactorMap[RF_rTarget]][k])) {
          maxValue = parent -> multiClassProb[RF_rFactorMap[RF_rTarget]][k];
          maxClass = (double) k;
        }
      }
      parent -> predictedOutcome = maxClass;
    }
    else {
      parent -> predictedOutcome = NA_REAL;
      if (!(RF_opt & OPT_OUTC_TYPE)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Zero classification count encountered in node during getMultiClassProb():  %10d", leaf);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }  
}
void updateEnsembleMultiClass(uint     mode,
                              uint     treeID,
                              double  *ensembleOutcome) {
  uint obsSize;
  unsigned char oobFlag, fullFlag, selectionFlag, outcomeFlag;
  double ***ensemblePtr;
  Node   ***nodeMembershipPtr;
  uint     *ensembleDenPtr;
  uint i, k;
  Node *parent;
  double maxValue;
  double maxClass;
  ensemblePtr    = NULL;  
  ensembleDenPtr = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    outcomeFlag = TRUE;
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
    outcomeFlag = TRUE;
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
    for (i=1; i <= obsSize; i++) {
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
        for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; k++) {
          ensemblePtr[1][k][i] += (double) (parent -> multiClassProb)[RF_rFactorMap[RF_rTarget]][k] / (double) (parent -> membrCount);
        }
        ensembleDenPtr[i] ++;
      }
      if (outcomeFlag == TRUE) {
        if (ensembleDenPtr[i] != 0) {
          maxValue = 0;
          maxClass = 0;
          for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; k++) {
            if (maxValue < ensemblePtr[1][k][i]) {
              maxValue = ensemblePtr[1][k][i];
              maxClass = (double) k;
            }
          }
          ensembleOutcome[i] = maxClass;
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
double getBrierScore(uint     obsSize,
                     double  *responsePtr,
                     double **ensemblePtr,
                     uint    *denomCount,
                     double  *condPerformance) {
  uint k;
  uint against;
  uint *oaaResponse;
  uint cumDenomCount;
  double result;
  oaaResponse = uivector(1, obsSize);
  result = 0.0;
  cumDenomCount = 0;
  for (k = 1; k <= obsSize; k ++) {
    if (denomCount[k] != 0) {
      cumDenomCount += 1;
    }
  }
  for (against = 1; against < RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; against++) {
    for (k = 1; k <= obsSize; k ++) {
      if ((uint) responsePtr[k] == against) {
        oaaResponse[k] = 1;
      }
      else {
        oaaResponse[k] = 0;
      }
    }
    condPerformance[against] = 0.0;
    for (k = 1; k <= obsSize; k ++) {
      if (denomCount[k] != 0) {
        condPerformance[against] += pow(((double) oaaResponse[k] - (ensemblePtr[against][k] / (double) denomCount[k])), 2.0);
      }
    }
    if (cumDenomCount == 0) {
      condPerformance[against] = NA_REAL;
    }
    else {
      condPerformance[against] = condPerformance[against] / (double) cumDenomCount;
      result += condPerformance[against];
    }
  }
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = result / RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
  }
  free_uivector(oaaResponse, 1, obsSize);
  return result;
}
void getConditionalClassificationIndex(uint    size,
                                       double *responsePtr,
                                       double *predictedOutcome,
                                       uint   *denomCount,
                                       double *condPerformance) {
  uint i, j;
  uint cumDenomCount;
  uint *condClassificationCount;
  cumDenomCount = 0;
  condClassificationCount = uivector(1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]]);
  for (j=1; j <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
    condClassificationCount[j] = 0;
    condPerformance[j] = 0;
  }
  for (i=1; i <= size; i++) {
    condClassificationCount[(uint) responsePtr[i]] ++;
    if (denomCount[i] != 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == predictedOutcome[i]) {
        condPerformance[(uint) responsePtr[i]] += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    for (j=1; j <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
        condPerformance[j] = NA_REAL;
    }
  }
  else {
    for (j=1; j <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
      if (condClassificationCount[j] != 0) {
        condPerformance[j] = 1.0 - condPerformance[j] / (double) condClassificationCount[j];
      }
      else {
        condPerformance[j] = NA_REAL;
      }
    }
  }
  free_uivector(condClassificationCount, 1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]]);
  return;
}
double getClassificationIndex(uint    size,
                              double *responsePtr,
                              double *predictedOutcome,
                              uint   *denomCount) {
  uint i;
  uint cumDenomCount;
  double result;
  cumDenomCount = 0;
  result = 0.0;
  for (i=1; i <= size; i++) {
    if (denomCount[i] > 0) {
      cumDenomCount += 1;
      if (responsePtr[i] == predictedOutcome[i]) {
        result += 1.0;
      }
    }
  }  
  if (cumDenomCount == 0) {
    result = NA_REAL;
  }
  else {
    result = 1.0 - result / (double) cumDenomCount;
  }
  return result;
}
