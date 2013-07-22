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
#include       "nodeOps.h"
#include        "random.h"
#include     "factorOps.h"
#include        "impute.h"
#include      "survival.h"
#include     "survivalE.h"
#include     "rfsrcUtil.h"
#include    "importance.h"
Node *getPerturbedMembership (Node    *parent, 
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
      result = getPerturbedMembership(parent ->  left, shadowVIMP, index);
    }
    else {
      result = getPerturbedMembership(parent -> right, shadowVIMP, index);
    }
  }
  return result;
}
Node *randomizeMembership(Node    *parent, 
                          double **predictor, 
                          uint     individual, 
                          uint     splitParameter) {
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
      if (ran3(parent -> splitParameter) <= 0.5) {
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
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
        result = randomizeMembership(parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(parent -> right, predictor, individual, splitParameter);
      }
    }
  }
  return result;
}
void permute(uint p, uint n, uint *indx) {
  uint i,j,k;
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ran3(p)*(i*1.0));
    for (j = 1; k > 0; j++) {
      if (indx[j] == 0) {
        k--;
      }
    }
    indx[j-1] = i;
  }
}
void getRandomMembership (uint      mode,
                          uint      treeID,
                          Node    **vimpMembership,
                          uint      p) {
  Node    *rootPtr;
  uint     obsSize;
  double **predictorPtr;
  uint    *membershipFlag;
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
  for (i=1; i <= obsSize; i++) {
    if ((membershipFlag[i] == selectionFlag)) {
      vimpMembership[i] = randomizeMembership(rootPtr, predictorPtr, i, p);
    }
    else {
      vimpMembership[i] = NULL;
    }
  }
}
void getPermuteMembership (uint      mode,
                           uint      treeID,
                           Node    **vimpMembership,
                           uint      p) {
  Node    *rootPtr;
  uint     obsSize;
  double **predictorPtr;
  uint    *membershipFlag;
  char     selectionFlag;
  uint     permuteObsSize;
  uint    *indexVIMP;
  uint    *permuteVIMP;
  double **originalVIMP;
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
    predictorPtr = RF_observation[treeID];
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    permuteObsSize = RF_oobSize[treeID];
    break;
  }
  indexVIMP = uivector(1, permuteObsSize);
  permuteVIMP = uivector(1, permuteObsSize);
  k = 0;
  for (i=1; i <= obsSize; i++) {
    if ((membershipFlag[i] == selectionFlag)) {
      k++;
      indexVIMP[k] = i;
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
  originalVIMP = dmatrix(1, pInnerCount, 1, permuteObsSize);
  shadowVIMP = (double**) vvector(1, RF_xSize);
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
    for (k=1; k<= permuteObsSize; k++) {
      originalVIMP[pIn][k] = predictorPtr[targetCov][indexVIMP[k]];
    }
    permute(targetCov, permuteObsSize, permuteVIMP);
    for (k=1; k <= obsSize; k++) {
      shadowVIMP[targetCov][k] = predictorPtr[targetCov][k];
    }
    for (k=1; k <= permuteObsSize; k++) {
      shadowVIMP[targetCov][indexVIMP[k]] = predictorPtr[targetCov][permuteVIMP[k]];
    }
  }
  for (i=1; i <= obsSize; i++) {
    if ((membershipFlag[i] == selectionFlag)) {
      vimpMembership[i] = getPerturbedMembership(rootPtr, shadowVIMP, i);
    }
    else {
      vimpMembership[i] = NULL;
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
  free_dmatrix(originalVIMP, 1, pInnerCount, 1, permuteObsSize);
  free_vvector(shadowVIMP, 1, RF_xSize);
  free_uivector(indexVIMP, 1, permuteObsSize);
  free_uivector(permuteVIMP, 1, permuteObsSize);
} 
 void getVIMPmembership (uint mode, uint treeID, Node **vimpMembership, uint p) {
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
                         Node     **vimpMembership,
                         uint       p) {
  uint   obsSize;
  uint  *membershipFlag;
  char   selectionFlag;
  char   ensembleFlag;
  selectionFlag = ACTIVE;  
  if (RF_opt & OPT_VIMP_LEOB) {
    ensembleFlag = FALSE;
  }
  else {
    ensembleFlag = TRUE;
  }
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
  updateGenericVimpEnsemble(treeID,
                            p,
                            obsSize,
                            selectionFlag,
                            membershipFlag,
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
  uint   obsSize;
  char   selectionFlag;
  char   ensembleFlag;
  uint      *membershipFlag;
  Node     **nodeMembership;
  selectionFlag = ACTIVE;  
  ensembleFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    membershipFlag = RF_testMembershipFlag;
    nodeMembership = RF_ftNodeMembership[treeID];
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    membershipFlag = RF_bootMembershipFlag[treeID];
    nodeMembership = RF_tNodeMembership[treeID];
    selectionFlag = FALSE;
    break;
  }
  updateGenericVimpEnsemble(treeID,
                            1,
                            obsSize,
                            selectionFlag,
                            membershipFlag,
                            nodeMembership,
                            ensembleFlag,
                            treeOutcome,
                            sTreeOutcome,
                            mcTreeEnsemble);
}
void updateGenericVimpEnsemble (uint       treeID,
                                uint       targetIndex,
                                uint       obsSize,
                                uint       selectionFlag,
                                uint      *membershipFlag,
                                Node     **noiseMembership,
                                char       ensembleFlag,
                                double   **genOutcome,
                                double  ***sGenOutcome,
                                double  ***mcGenEnsemble) {
  Node  *terminalNode;
  uint   ensembleDim;
  uint   i, j;
  ensembleDim = getEnsembleDim();
    for (i=1; i <= obsSize; i++) {
      if ((membershipFlag[i] == selectionFlag)) {
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
            if (RF_rFactorCount > 0) {
              for (j=1; j<= ensembleDim; j++) {
                if (!ensembleFlag) {
                  mcGenEnsemble[targetIndex][j][i] = (double) (terminalNode -> multiClassProb)[1][j] / (double) (terminalNode -> membrCount);
                }
                else {
                  mcGenEnsemble[targetIndex][j][i] += (double) (terminalNode -> multiClassProb)[1][j] / (double) (terminalNode -> membrCount);
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
    if (RF_rFactorCount > 0) {
      vimpDenomPtr = RF_vimpEnsembleDen[p];
      for (i = 1; i <= obsSize; i++) {
        if(vimpDenomPtr[i] > 0) {
          maxValue = 0;
          maxClass = 0;
          for (j=1; j <= RF_rFactorSize[1]; j++) {
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
      if (RF_rFactorCount > 0) {
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
    if (RF_rFactorCount > 0) {
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
  if (RF_rFactorCount > 0) {
    perfDimOne ++;
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
Node **stackVimpMembership(uint mode) {
  Node **membership;
  uint obsSize;
  membership = NULL;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      break;
    default:
      obsSize = RF_observationSize;
      break;
    }
    membership = (Node **) vvector(1, obsSize);
  }
  return membership;
}
void unstackVimpMembership(uint mode, Node **membership) {
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
    free_vvector(membership, 1, obsSize);
  }
}
void stackTreeEnsemble(uint         mode,
                       double    ***treeOutcome,
                       double   ****sTreeOutcome,
                       double   ****mcTreeEnsemble) {
  uint obsSize;
  uint ensembleDim;
  uint i,j,k;
  *sTreeOutcome     = NULL;
  *mcTreeEnsemble   = NULL;
  ensembleDim = getEnsembleDim();
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
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
    if (RF_rFactorCount > 0) {
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
  free_dmatrix(treeOutcome, 1, 1, 1, obsSize);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    free_dmatrix3(sTreeOutcome, 1, 1, 1, ensembleDim, 1, obsSize);
  }
  else {
    if (RF_rFactorCount > 0) {
      free_dmatrix3(mcTreeEnsemble, 1, 1, 1, ensembleDim, 1, obsSize);
    }
  }
}
void updateVimpCalculations (uint mode, uint b, uint intrIndex) {
  Node   **vimpMembership;
  uint     p;
  vimpMembership = NULL;  
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
  if (!(RF_opt & OPT_VIMP_JOIN)) {
    p = RF_intrPredictor[intrIndex];
  }
  else {
    p = 0;
  }
  vimpMembership = stackVimpMembership(mode);
  getVIMPmembership(mode, b, vimpMembership, p);
  updateVimpEnsemble(mode, b, vimpMembership, intrIndex);
  if (RF_opt & OPT_VIMP_LEOB) {
    summarizeVimpPerformance(mode, b, intrIndex);
  }
  unstackVimpMembership(mode, vimpMembership);
}
void summarizeTreePerformance(uint mode, uint treeID) {
  uint i, j;
  uint obsSize;
  uint *denomPtr;
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
    denomPtr = RF_testMembershipFlag;
    break;     
  default:
    obsSize = RF_observationSize;
    responsePtr = RF_response[treeID];
    denomPtr = RF_oobMembershipFlag[treeID];
    break;
  }
  stackTreeEnsemble(mode, &treeOutcome, &sTreeOutcome, &mcTreeEnsemble);
  updateTreeEnsemble(mode, treeID, treeOutcome, sTreeOutcome, mcTreeEnsemble);
  ensembleDim = getEnsembleDim();
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    for (j = 1; j <= ensembleDim; j ++) {
      for (i = 1; i <= obsSize; i++) {
        if (denomPtr[i] != 0) {
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
    if (RF_rFactorCount > 0) {
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
                 denomPtr,
                 RF_perfLeo[treeID]);
  unstackTreeEnsemble(mode, treeOutcome, sTreeOutcome, mcTreeEnsemble);
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
    if (RF_rFactorCount > 0) {
      result = RF_rFactorSize[1];
    }
    else {
      result = 1;
    }
  }
  return result;
}
