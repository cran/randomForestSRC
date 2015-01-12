////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.6.0
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
#include        "random.h"
#include     "factorOps.h"
#include        "impute.h"
#include      "survival.h"
#include     "survivalE.h"
#include     "rfsrcUtil.h"
#include    "importance.h"
Node *identifyPerturbedMembership (Node    *parent,
                                   double **shadowVIMP,
                                   uint     index) {
  char daughterFlag;
  Node *result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    daughterFlag = RIGHT;
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      daughterFlag = splitOnFactor((uint) shadowVIMP[parent -> splitParameter][index], parent -> splitValueFactPtr);
    }
    else {
      if (shadowVIMP[parent -> splitParameter][index] <= (parent -> splitValueCont)) {
        daughterFlag = LEFT;
      }
    }
    if (daughterFlag == LEFT) {
      result = identifyPerturbedMembership(parent ->  left, shadowVIMP, index);
    }
    else {
      result = identifyPerturbedMembership(parent -> right, shadowVIMP, index);
    }
  }
  return result;
}
Node *randomizeMembership(Node    *parent,
                          double **predictor,
                          uint     individual,
                          uint     splitParameter,
                          uint     treeID) {
  char daughterFlag;
  char randomSplitFlag;
  Node *result;
  result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    randomSplitFlag = FALSE;
    if (splitParameter > 0) {
      if ((parent -> splitParameter) == splitParameter) {
        randomSplitFlag = TRUE;
      }
    }
    else {
      if(RF_importanceFlag[parent -> splitParameter] == TRUE) {
        randomSplitFlag = TRUE;
      }
    }
    if(randomSplitFlag == TRUE) {
      if (ran1C(treeID) <= 0.5) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter, treeID);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter, treeID);
      }
    }
    else {
      daughterFlag = RIGHT;
      if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
        daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][individual], parent -> splitValueFactPtr);
      }
      else {
        if (predictor[parent -> splitParameter][individual] <= (parent -> splitValueCont)) {
          daughterFlag = LEFT;
        }
      }
      if (daughterFlag == LEFT) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter, treeID);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter, treeID);
      }
    }
  }
  return result;
}
void permute(uint ranGenID, uint parallelID, uint n, uint *indx) {
  float (*ranX) (uint);
  uint i,j,k;
  ranX = NULL;  
  if ((ranGenID != 1) && (ranGenID != 2) && (ranGenID != 3)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Invalid random generator selected:  %10d", ranGenID);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  switch(ranGenID) {
  case 1:
    ranX = ran1A;
    break;
  case 2:
    ranX = ran1B;
    break;
  case 3:
    ranX = ran1C;
    break;
  }
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ranX(parallelID)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
}
void getRandomMembership (uint       mode,
                          uint       treeID,
                          Terminal **vimpMembership,
                          uint       p) {
  Node    *rootPtr;
  uint     obsSize;
  double **predictorPtr;
  char    *membershipFlag;
  char     selectionFlag;
  uint     i;
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    predictorPtr = RF_fobservation[treeID];
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    predictorPtr = RF_observation[treeID];
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    break;
  }
  if (RF_sobservationSize > 0) {
    for (i = 1; i <= obsSize; i++) {
      if (membershipFlag[i] == selectionFlag) {
        vimpMembership[i] = RF_tTermMembership[treeID][i];
      }
      else {
        vimpMembership[i] = NULL;
      }
    }
    for (i = 1; i <= RF_sobservationSize; i++) {
      if (membershipFlag[RF_sobservationIndv[i]] == selectionFlag) {
        vimpMembership[RF_sobservationIndv[i]] = RF_tTermList[treeID][ randomizeMembership(rootPtr, predictorPtr, RF_sobservationIndv[i], p, treeID) -> nodeID ];
      }
    }
  }
  else {
    for (i = 1; i <= obsSize; i++) {
      if (membershipFlag[i] == selectionFlag) {
        vimpMembership[i] = RF_tTermList[treeID][ randomizeMembership(rootPtr, predictorPtr, i, p, treeID) -> nodeID ];
      }
      else {
        vimpMembership[i] = NULL;
      }
    }
  }  
}
void getPermuteMembership (uint       mode,
                           uint       treeID,
                           Terminal **vimpMembership,
                           uint       p) {
  Node    *rootPtr;
  uint     obsSize;
  double **predictorPtr;
  char    *membershipFlag;
  char     selectionFlag;
  uint     permuteObsSize;
  uint    *indexVIMP;
  uint    *permuteVIMP;
  double **shadowVIMP;
  uint     pInnerCount, pIn;
  uint     i, j, k, targetCov;
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    predictorPtr = RF_fobservation[treeID];
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    permuteObsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    if(RF_sobservationSize > 0) {
      permuteObsSize = RF_soobSize[treeID];
    }
    else {
      permuteObsSize = RF_oobSize[treeID];
    }
    predictorPtr = RF_observation[treeID];
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    break;
  }
  indexVIMP = uivector(1, permuteObsSize + 1);
  permuteVIMP = uivector(1, permuteObsSize + 1);
  if (RF_sobservationSize > 0) {
    k = 0;
    for (i = 1; i <= RF_sobservationSize; i++) {
      if (membershipFlag[RF_sobservationIndv[i]] == selectionFlag) {
        k++;
        indexVIMP[k] = RF_sobservationIndv[i];
      }
    }
  }
  else {
    k = 0;
    for (i = 1; i <= obsSize; i++) {
      if (membershipFlag[i] == selectionFlag) {
        k++;
        indexVIMP[k] = i;
      }
    }
  }
  if (k != permuteObsSize) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  VIMP candidate selection failed.");
    Rprintf("\nRF-SRC:  %10d available, %10d selected.", permuteObsSize, k);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (p > 0) {
    pInnerCount = 1;
  }
  else {
    pInnerCount = RF_intrPredictorSize;
  }
  shadowVIMP = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
  for (j = 1; j <= RF_xSize; j++) {
    shadowVIMP[j] = predictorPtr[j];
  }
  for (pIn = 1; pIn <= pInnerCount; pIn++) {
    if (p > 0) {
      targetCov = p;
    }
    else {
      targetCov = RF_intrPredictor[pIn];
    }
    shadowVIMP[targetCov] = dvector(1, obsSize);
    for (i = 1; i <= obsSize; i++) {
      shadowVIMP[targetCov][i] = predictorPtr[targetCov][i];
    }
    permute(3, treeID, permuteObsSize, permuteVIMP);
    for (k = 1; k <= permuteObsSize; k++) {
      shadowVIMP[targetCov][indexVIMP[k]] = predictorPtr[targetCov][indexVIMP[permuteVIMP[k]]];
    }
  }
  if (RF_sobservationSize > 0) {
    for (i = 1; i <= obsSize; i++) {
      if (membershipFlag[i] == selectionFlag) {
        vimpMembership[i] = RF_tTermMembership[treeID][i];
      }
      else {
        vimpMembership[i] = NULL;
      }
    }
    for (i = 1; i <= RF_sobservationSize; i++) {
      if (membershipFlag[RF_sobservationIndv[i]] == selectionFlag) {
        vimpMembership[RF_sobservationIndv[i]] = RF_tTermList[treeID][ identifyPerturbedMembership(rootPtr, shadowVIMP, RF_sobservationIndv[i]) -> nodeID ];
      }
    }
  }
  else {
    for (i = 1; i <= obsSize; i++) {
      if (membershipFlag[i] == selectionFlag) {
        vimpMembership[i] = RF_tTermList[treeID][ identifyPerturbedMembership(rootPtr, shadowVIMP, i) -> nodeID ];
      }
      else {
        vimpMembership[i] = NULL;
      }
    }
  }  
  for (pIn = 1; pIn <= pInnerCount; pIn++) {
    if (p > 0) {
      targetCov = p;
    }
    else {
      targetCov = RF_intrPredictor[pIn];
    }
    free_dvector(shadowVIMP[targetCov], 1, obsSize);
  }
  free_new_vvector(shadowVIMP, 1, RF_xSize, NRUTIL_DPTR);
  free_uivector(indexVIMP, 1, permuteObsSize + 1);
  free_uivector(permuteVIMP, 1, permuteObsSize + 1);
}
void getVimpMembership (uint mode, uint treeID, Terminal **vimpMembership, uint p) {
  char result;
  if (!(RF_opt & OPT_VIMP)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to compute variable importance though not requested.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  result = FALSE;
  switch (mode) {
  case RF_PRED:
    result = TRUE;
    break;
  default:
    if (RF_oobSize[treeID] > 0) {
      result = TRUE;
    }
    break;
  }
  if (result == TRUE) {
    if (RF_opt & OPT_VIMP_TYPE) {
      getRandomMembership(mode, treeID, vimpMembership, p);
    }
    else {
      getPermuteMembership(mode, treeID, vimpMembership, p);
    }
  }
  else {
  }
}
void updateVimpEnsemble (uint       mode,
                         uint       treeID,
                         Terminal **vimpMembership,
                         uint       p) {
  char   ensembleFlag;
  if (RF_opt & OPT_VIMP_LEOB) {
    ensembleFlag = FALSE;
  }
  else {
    ensembleFlag = TRUE;
  }
  updateGenericVimpEnsemble(mode,
                            treeID,
                            p,
                            vimpMembership,
                            ensembleFlag,
                            RF_vimpOutcome,
                            RF_sVimpOutcome,    
                            RF_cVimpEnsemble);  
}
void updateTreeEnsemble (uint       mode,
                         uint       treeID,
                         double   **treeOutcome,
                         double  ***sTreeOutcome,
                         double  ***mcTreeEnsemble) {
  char ensembleFlag;
  Terminal **termMembership;
  ensembleFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    termMembership = RF_ftTermMembership[treeID];
    break;
  default:
    termMembership = RF_tTermMembership[treeID];
    break;
  }
  updateGenericVimpEnsemble(mode,
                            treeID,
                            1,
                            termMembership,
                            ensembleFlag,
                            treeOutcome,
                            sTreeOutcome,
                            mcTreeEnsemble);
}
void updateGenericVimpEnsemble (uint       mode,
                                uint       treeID,
                                uint       targetIndex,
                                Terminal **noiseMembership,
                                char       ensembleFlag,
                                double   **genOutcome,
                                double  ***sGenOutcome,
                                double  ***mcGenEnsemble) {
  Terminal *terminalNode;
  uint   ensembleDim;
  uint   obsSize;
  char  *membershipFlag;
  char   selectionFlag;
  uint   i, j;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    break;
  }
  ensembleDim = getEnsembleDim();
  for (i=1; i <= obsSize; i++) {
    if (membershipFlag[i] == selectionFlag) {
      terminalNode = noiseMembership[i];
      if (!ISNA(terminalNode -> predictedOutcome)) {
        if (!ensembleFlag) {
          genOutcome[targetIndex][i] = terminalNode -> predictedOutcome;
          RF_vimpEnsembleDen[targetIndex][i] = 1;
        }
        else {
            genOutcome[targetIndex][i] += terminalNode -> predictedOutcome;
            RF_vimpEnsembleDen[targetIndex][i] ++;
        }
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          for (j=1; j <= ensembleDim; j++) {
            if (!ensembleFlag) {
              sGenOutcome[targetIndex][j][i] = terminalNode -> mortality[j];
            }
            else {
                sGenOutcome[targetIndex][j][i] += terminalNode -> mortality[j];
            }
          }
        }
        else {
          if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
            for (j=1; j<= ensembleDim; j++) {
              if (!ensembleFlag) {
                mcGenEnsemble[targetIndex][j][i] = (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTarget]][j] / (double) (terminalNode -> membrCount);
              }
              else {
                  mcGenEnsemble[targetIndex][j][i] += (double) (terminalNode -> multiClassProb)[RF_rFactorMap[RF_rTarget]][j] / (double) (terminalNode -> membrCount);
              }
            }
          }
        }
      }
      else {
        if (!(RF_opt & OPT_OUTC_TYPE)) {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  NA encountered for VIMP outcome in terminal node:  %10d", terminalNode -> nodeID);
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    else {
      if (!ensembleFlag) {
        RF_vimpEnsembleDen[targetIndex][i] = 0;
      }
    }
  }  
}
void summarizeVimpPerformance(uint       mode,
                              uint       treeID,
                              uint       p) {
  uint      obsSize;
  double    **responsePtr;
  double    **importancePtr;
  uint      *vimpDenomPtr;
  double     *vimpOutcomePtr;
  double    **subVimpOutcomePtr;
  char        responseImputeFlag;
  double    maxValue, maxClass;
  uint      normalizationFlag;
  uint i, j;
  if (!(RF_opt & OPT_VIMP)) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  VIMP analysis requested while OPT bit not set.  \n");
    return;
  }
  vimpOutcomePtr     = NULL;
  subVimpOutcomePtr  = NULL;
  responseImputeFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    responsePtr = RF_fresponseIn;
    if (RF_fmRecordSize > 0) {
      if(RF_fmResponseFlag == TRUE) {
        responseImputeFlag = TRUE;
      }
    }
    break;
  default:
    obsSize = RF_observationSize;
    responsePtr = RF_responseIn;
    if (RF_mRecordSize > 0) {
      if(RF_mResponseFlag == TRUE) {
        responseImputeFlag = TRUE;
      }
    }
    break;
  }
  responsePtr = stackAndImputeGenericResponse(responseImputeFlag, mode, RF_rSize, obsSize, 0, RF_forestSize, responsePtr);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    vimpDenomPtr = RF_vimpEnsembleDen[p];
    if (!(RF_opt & OPT_COMP_RISK)) {
      getEnsembleMortality(mode, treeID, obsSize, RF_sVimpOutcome[p], vimpDenomPtr, RF_sVimpOutcome[p][1]);
    }
    else {
      getEnsembleMortalityCR(mode, treeID, obsSize, RF_sVimpOutcome[p], vimpDenomPtr, RF_sVimpOutcome[p]);
    }
  }  
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      vimpDenomPtr = RF_vimpEnsembleDen[p];
      for (i = 1; i <= obsSize; i++) {
        if(vimpDenomPtr[i] > 0) {
          maxValue = 0;
          maxClass = 0;
          for (j=1; j <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
            if (maxValue < RF_cVimpEnsemble[p][j][i]) {
              maxValue = RF_cVimpEnsemble[p][j][i];
              maxClass = (double) j;
            }
          }
          RF_vimpOutcome[p][i] = maxClass;
        }
        else {
          RF_vimpOutcome[p][i] = NA_REAL;
        }
      }  
    }  
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    normalizationFlag = FALSE;
  }
  else {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      normalizationFlag = FALSE;
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        normalizationFlag = FALSE;
      }
      else {
        normalizationFlag = TRUE;
      }
    }
  }
  vimpDenomPtr = RF_vimpEnsembleDen[p];
  if (normalizationFlag == TRUE) {
    for (i = 1; i <= obsSize; i++) {
      if (vimpDenomPtr[i] != 0) {
        RF_vimpOutcome[p][i] = RF_vimpOutcome[p][i] / vimpDenomPtr[i];
      }
    }
  }
  vimpOutcomePtr = RF_vimpOutcome[p];
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_COMP_RISK) {
      vimpOutcomePtr = NULL;
      subVimpOutcomePtr = RF_sVimpOutcome[p];
    }
    else {
      subVimpOutcomePtr = NULL;
      vimpOutcomePtr = RF_sVimpOutcome[p][1];
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      subVimpOutcomePtr = RF_cVimpEnsemble[p];
    }
    else {
      subVimpOutcomePtr = NULL;
    }
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    importancePtr = RF_vimpLeo[treeID];
  }
  else {
    importancePtr = RF_importancePtr;
  }
  getPerformance(treeID,
                 mode,
                 obsSize,
                 responsePtr,
                 vimpOutcomePtr,
                 subVimpOutcomePtr,
                 vimpDenomPtr,
                 importancePtr[p]);
  unstackImputeResponse(responseImputeFlag, RF_rSize, obsSize, responsePtr);
}
void finalizeVimpPerformance(uint       mode,
                             uint       rejectedTreeCount) {
  uint varCount, perfDimOne;
  double result;
  uint cumDenomCount;
  uint i,j,k;
  switch (mode) {
  case RF_PRED:
    if (RF_opt & OPT_VIMP_JOIN) {
      varCount = 1;
    }
    else {
      varCount = RF_intrPredictorSize;
    }
    break;
  default:
    if (RF_opt & OPT_VIMP_JOIN) {
      varCount = 1;
    }
    else {
      varCount = RF_intrPredictorSize;
    }
    break;
  }
  perfDimOne = getEnsembleDim();
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      perfDimOne ++;
    }
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    for(j = 1; j <= varCount; j++) {
      for (k = 1; k <= perfDimOne; k++) {
        cumDenomCount = 0;
        result = 0.0;
        for (i = 1; i <= RF_forestSize; i++) {
          if(!ISNA(RF_vimpLeo[i][j][k])) {
            if(!ISNA(RF_perfLeo[i][k])) {
              result += RF_vimpLeo[i][j][k] - RF_perfLeo[i][k];
              cumDenomCount ++;
            }
          }
        }
        if (cumDenomCount != 0) {
          RF_importancePtr[j][k] = result / (double) cumDenomCount;
        }
        else {
          RF_importancePtr[j][k] = NA_REAL;
        }
      }
    }
  }
  else {
    for(j = 1; j <= varCount; j++) {
      for (k = 1; k <= perfDimOne; k++) {
        if(!ISNA(RF_importancePtr[j][k])) {
          if(!ISNA(RF_performancePtr[RF_forestSize][k])) {
            RF_importancePtr[j][k] = RF_importancePtr[j][k] - RF_performancePtr[RF_forestSize][k];
          }
          else {
          }
        }
      }
    }
  }
}
void stackVimpMembership(uint mode, Terminal ***membership) {
  uint obsSize;
  (*membership) = NULL;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    *membership = (Terminal **) new_vvector(1, obsSize, NRUTIL_NPTR);
  }
}
void unstackVimpMembership(uint mode, Terminal **membership) {
  uint obsSize;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    free_new_vvector(membership, 1, obsSize, NRUTIL_NPTR);
  }
}
void stackTreeEnsemble(uint         mode,
                       uint         treeID,
                       uint       **denomTree,
                       double    ***treeOutcome,
                       double   ****sTreeOutcome,
                       double   ****mcTreeEnsemble) {
  uint  obsSize;
  char *denomPtr;
  uint  ensembleDim;
  uint i,j,k;
  *sTreeOutcome     = NULL;
  *mcTreeEnsemble   = NULL;
  ensembleDim = getEnsembleDim();
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    denomPtr = RF_testMembershipFlag;
    break;
  default:
    obsSize = RF_observationSize;
    denomPtr = RF_oobMembershipFlag[treeID];
    break;
  }
  *denomTree = uivector(1, obsSize);
  for (i = 1; i <= obsSize; i++) {
    (*denomTree)[i] = (uint) denomPtr[i];
  }
  *treeOutcome = dmatrix(1, 1, 1, obsSize);
  for (i = 1; i <= obsSize; i++) {
    (*treeOutcome)[1][i] = 0.0;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    *sTreeOutcome = dmatrix3(1, 1, 1, ensembleDim, 1, obsSize);
    for (j = 1; j <= ensembleDim; j++) {
      for (i = 1; i <= obsSize; i++) {
        (*sTreeOutcome)[1][j][i] = 0.0;
      }
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      *mcTreeEnsemble = dmatrix3(1, 1, 1, ensembleDim, 1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        for (k = 1; k <= ensembleDim; k++) {
          (*mcTreeEnsemble)[1][k][i] = 0.0;
        }
      }
    }
  }
}
void unstackTreeEnsemble(uint        mode,
                         uint        treeID,
                         uint       *denomTree,
                         double    **treeOutcome,
                         double   ***sTreeOutcome,
                         double   ***mcTreeEnsemble) {
  uint obsSize;
  uint ensembleDim;
  ensembleDim = getEnsembleDim();
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  free_uivector(denomTree, 1, obsSize);
  free_dmatrix(treeOutcome, 1, 1, 1, obsSize);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    free_dmatrix3(sTreeOutcome, 1, 1, 1, ensembleDim, 1, obsSize);
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      free_dmatrix3(mcTreeEnsemble, 1, 1, 1, ensembleDim, 1, obsSize);
    }
  }
}
void updateVimpCalculations (uint mode, uint b, uint intrIndex, Terminal **vimpMembership) {
  if (RF_tLeafCount[b] == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to compute importance on a rejected tree:  %10d", b);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (!(RF_opt & OPT_VIMP)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to compute importance when the option is not active:  %10d", b);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
#ifdef SUPPORT_OPENMP
#pragma omp critical (_update_gve)
#endif
  {  
    updateVimpEnsemble(mode, b, vimpMembership, intrIndex);
  }  
  if (RF_opt & OPT_VIMP_LEOB) {
    summarizeVimpPerformance(mode, b, intrIndex);
  }
}
void summarizeTreePerformance(uint mode, uint treeID) {
  uint i, j;
  uint  obsSize;
  uint *denomTree;
  double    **responsePtr;
  double    **treeOutcome;
  double   ***sTreeOutcome;
  double   ***mcTreeEnsemble;
  double     *treeOutcomePtr;
  double    **subTreeOutcomePtr;
  uint ensembleDim;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    responsePtr = RF_fresponse[treeID];
    break;
  default:
    obsSize = RF_observationSize;
    responsePtr = RF_response[treeID];
    break;
  }
  stackTreeEnsemble(mode, treeID, &denomTree, &treeOutcome, &sTreeOutcome, &mcTreeEnsemble);
  updateTreeEnsemble(mode, treeID, treeOutcome, sTreeOutcome, mcTreeEnsemble);
  ensembleDim = getEnsembleDim();
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    for (j = 1; j <= ensembleDim; j ++) {
      for (i = 1; i <= obsSize; i++) {
        if (denomTree[i] != 0) {
        }
        else {
          sTreeOutcome[1][j][i] = NA_REAL;
        }
      }
    }
  }
  treeOutcomePtr = treeOutcome[1];
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_COMP_RISK) {
      treeOutcomePtr = NULL;
      subTreeOutcomePtr = sTreeOutcome[1];
    }
    else {
      treeOutcomePtr = sTreeOutcome[1][1];
      subTreeOutcomePtr = NULL;
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      subTreeOutcomePtr = mcTreeEnsemble[1];
    }
    else {
      subTreeOutcomePtr = NULL;
    }
  }
  getPerformance(treeID,
                 mode,
                 obsSize,
                 responsePtr,
                 treeOutcomePtr,
                 subTreeOutcomePtr,
                 denomTree,
                 RF_perfLeo[treeID]);
  unstackTreeEnsemble(mode, treeID, denomTree, treeOutcome, sTreeOutcome, mcTreeEnsemble);
}
uint getEnsembleDim () {
  uint result = 0;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_COMP_RISK) {
      result = RF_eventTypeSize;
    }
    else {
      result = 1;
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      result = RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
    }
    else {
      result = 1;
    }
  }
  return result;
}
