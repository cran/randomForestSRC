////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.0.1
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
#include       "nodeOps.h"
#include        "random.h"
#include     "factorOps.h"
#include        "impute.h"
#include      "survival.h"
#include     "survivalE.h"
#include     "rfsrcUtil.h"
#include    "importance.h"
Node *getProxyMember (Node    *parent, 
                      double **predictor, 
                      uint     index) {
  char daughterFlag;
  Node *result = parent;
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    daughterFlag = RIGHT;
    if (strcmp(RF_xType[parent -> splitParameter], "C") == 0) {
      daughterFlag = splitOnFactor((uint) predictor[parent -> splitParameter][index], parent -> splitValueFactPtr);
    }
    else {
      if (predictor[parent -> splitParameter][index] <= (parent -> splitValueCont)) {
        daughterFlag = LEFT;
      }
    }
    if (daughterFlag == LEFT) {
      result = getProxyMember(parent ->  left, predictor, index);
    }
    else {
      result = getProxyMember(parent -> right, predictor, index);
    }
  }
  return result;
}
Node *randomizeMembership(uint     treeID,
                          Node    *parent, 
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
      if (ran2(treeID) <= 0.5) {
        result = randomizeMembership(treeID, parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(treeID, parent -> right, predictor, individual, splitParameter);
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
        result = randomizeMembership(treeID, parent ->  left, predictor, individual, splitParameter);
      }
      else {
        result = randomizeMembership(treeID, parent -> right, predictor, individual, splitParameter);
      }
    }
  }
  return result;
}
void permute(uint treeID, uint n, uint *indx) {
  uint i,j,k;
  for (i=1; i<= n; i++) {
    indx[i] = 0;
  }
  for (i=n; i > 0; i--) {
    k = (uint) ceil(ran2(treeID)*(i*1.0));
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
                          Node   ***vimpMembership) {
  Node    *rootPtr;
  uint     obsSize;
  uint     varCount;
  double **predictorPtr;
  uint    *membershipFlag;
  char     selectionFlag;
  uint     pOuterCount;
  uint     i, p;
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    varCount = RF_intrPredictorSize;
    predictorPtr = RF_fobservation[treeID];
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    varCount = RF_intrPredictorSize;
    predictorPtr = RF_observation[treeID];
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    break;
  }
  if (!(RF_opt & OPT_VIMP_JOIN)) {
    pOuterCount = varCount;
  }
  else {
    pOuterCount = 1;
  }
  for (p=1; p <= pOuterCount; p++) {
    for (i=1; i <= obsSize; i++) {
      if ((membershipFlag[i] == selectionFlag)) {
        if (!(RF_opt & OPT_VIMP_JOIN)) {
          vimpMembership[p][i] = randomizeMembership(treeID, rootPtr, predictorPtr, i, RF_intrPredictor[p]);
        }
        else {
          vimpMembership[p][i] = randomizeMembership(treeID, rootPtr, predictorPtr, i, 0);
        }
      }
      else {
        vimpMembership[p][i] = NULL;
      }
    }
  }
}
void getPermuteMembership (uint      mode,
                           uint      treeID,
                           Node   ***vimpMembership) {
  Node    *rootPtr;
  uint     obsSize;
  uint     varCount;
  double **predictorPtr;
  uint    *membershipFlag;
  char     selectionFlag;
  uint     permuteObsSize;
  uint    *indexVIMP;
  uint    *permuteVIMP;
  double **originalVIMP;
  uint     pInnerCount, pIn;
  uint     pOuterCount, pOut;
  uint     i, k, p;
  rootPtr = RF_root[treeID];
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    varCount = RF_intrPredictorSize;
    predictorPtr = RF_fobservation[treeID];
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    permuteObsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    varCount = RF_intrPredictorSize;
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
  if (!(RF_opt & OPT_VIMP_JOIN)) {
    pOuterCount = varCount;
    pInnerCount = 1;
  }
  else {
    pOuterCount = 1;
    pInnerCount = varCount;
  }
  originalVIMP = dmatrix(1, pInnerCount, 1, permuteObsSize);
  for (pOut = 1; pOut <= pOuterCount; pOut++) {
    for (pIn = 1; pIn <= pInnerCount; pIn++) {
      if (!(RF_opt & OPT_VIMP_JOIN)) {
        p = pOut;
      }
      else {
        p = pIn;
      }
      for (k=1; k<= permuteObsSize; k++) {
        originalVIMP[pIn][k] = predictorPtr[RF_intrPredictor[p]][indexVIMP[k]];
      }
      permute(treeID, permuteObsSize, permuteVIMP);
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[RF_intrPredictor[p]][indexVIMP[k]] = originalVIMP[pIn][permuteVIMP[k]];
      }
    }
    for (i=1; i <= obsSize; i++) {
      if ((membershipFlag[i] == selectionFlag)) {
        vimpMembership[pOut][i] = getProxyMember(rootPtr, predictorPtr, i);
      }
      else {
        vimpMembership[pOut][i] = NULL;
      }
    }
    for (pIn = 1; pIn <= pInnerCount; pIn++) {
      if (!(RF_opt & OPT_VIMP_JOIN)) {
        p = pOut;
      }
      else {
        p = pIn;
      }
      for (k=1; k <= permuteObsSize; k++) {
        predictorPtr[RF_intrPredictor[p]][indexVIMP[k]] = originalVIMP[pIn][k];
      }
    }
  }  
  free_dmatrix(originalVIMP, 1, pInnerCount, 1, permuteObsSize);
  free_uivector(indexVIMP, 1, permuteObsSize);
  free_uivector(permuteVIMP, 1, permuteObsSize);
} 
void getVIMPmembership (uint mode, uint treeID, Node ***vimpMembership) {
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
      getRandomMembership(mode, treeID, vimpMembership);
    }
    else {
      getPermuteMembership(mode, treeID, vimpMembership);
    }
  }
  else {
  }
}
void updateVimpEnsemble (uint       mode,
                         uint       treeID,
                         Node    ***vimpMembership) {
  uint   obsSize;
  uint   varCount;
  uint  *membershipFlag;
  char   selectionFlag;
  char   ensembleFlag;
  uint   varLoopCount;
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
    varCount = RF_intrPredictorSize;
    membershipFlag = RF_testMembershipFlag;
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    varCount = RF_intrPredictorSize;
    membershipFlag = RF_bootMembershipFlag[treeID];
    selectionFlag = FALSE;
    break;
  }
  if (!(RF_opt & OPT_VIMP_JOIN)) {
    varLoopCount = varCount;
  }
  else {
    varLoopCount = 1;
  }
  updateGenericVimpEnsemble(treeID,
                            0,
                            varLoopCount,
                            obsSize,
                            selectionFlag,
                            membershipFlag,
                            vimpMembership,
                            ensembleFlag,
                            RF_vimpOutcome, 
                            RF_sVimpEnsemble,   
                            RF_cVimpEnsemble);  
}
void updateTreeEnsemble (uint       mode,
                         uint       treeID,
                         double   **treeOutcome,
                         double ****chfTreeEnsemble,
                         double  ***mcTreeEnsemble) {
  uint   obsSize;
  uint   varCount;
  char   selectionFlag;
  char   ensembleFlag;
  uint      *membershipFlag;
  Node    ***nodeMembership;
  selectionFlag = ACTIVE;  
  ensembleFlag = FALSE;
  varCount = 1;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    membershipFlag = RF_testMembershipFlag;
    nodeMembership = RF_fnodeMembership;
    selectionFlag = ACTIVE;
    break;
  default:
    obsSize = RF_observationSize;
    membershipFlag = RF_bootMembershipFlag[treeID];
    nodeMembership = RF_nodeMembership;
    selectionFlag = FALSE;
    break;
  }
  updateGenericVimpEnsemble(treeID,
                            treeID - 1,
                            varCount,
                            obsSize,
                            selectionFlag,
                            membershipFlag,
                            nodeMembership,
                            ensembleFlag,
                            treeOutcome,
                            chfTreeEnsemble,
                            mcTreeEnsemble);
}
void updateGenericVimpEnsemble (uint       treeID,
                                uint       nodeOffset,
                                uint       varLoopCount,
                                uint       obsSize,
                                uint       selectionFlag,
                                uint      *membershipFlag,
                                Node    ***noiseMembership,
                                char       ensembleFlag,
                                double   **genOutcome,
                                double ****chfGenEnsemble,
                                double  ***mcGenEnsemble) {
  Node  *terminalNode;
  uint   i, j, k, p;
  for (p=1; p <= varLoopCount; p++) {
    for (i=1; i <= obsSize; i++) {
      if ((membershipFlag[i] == selectionFlag)) {
        terminalNode = noiseMembership[nodeOffset + p][i];
        if (!ISNA(terminalNode -> predictedOutcome)) {
          if (!ensembleFlag) {
            genOutcome[p][i] = terminalNode -> predictedOutcome;
            RF_vimpEnsembleDen[p][i] = 1;
          }
          else {
            genOutcome[p][i] += terminalNode -> predictedOutcome;
            RF_vimpEnsembleDen[p][i] ++;
          }
          if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
            if (RF_eventTypeSize == 1) {
              for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                if (!ensembleFlag) {
                  chfGenEnsemble[p][1][k][i] = (terminalNode -> nelsonAalen)[k]; 
                }
                else {
                  chfGenEnsemble[p][1][k][i] += (terminalNode -> nelsonAalen)[k];   
                }
              }
            }
            else {
              for (j=1; j <= RF_eventTypeSize; j++) {
                for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                  if (!ensembleFlag) {
                    chfGenEnsemble[p][j][k][i] = (terminalNode -> CIF)[j][k]; 
                  }
                  else {
                    chfGenEnsemble[p][j][k][i] += (terminalNode -> CIF)[j][k];                     
                  }
                }
              }
            }
          }
          else {
            if (RF_rFactorCount > 0) {
              for (j=1; j<= RF_rFactorSize[1]; j++) {
                if (!ensembleFlag) {
                  mcGenEnsemble[p][j][i] = (double) (terminalNode -> multiClassProb)[1][j] / (double) (terminalNode -> membrCount);
                }
                else {
                  mcGenEnsemble[p][j][i] += (double) (terminalNode -> multiClassProb)[1][j] / (double) (terminalNode -> membrCount);
                }
              }
            }
          }
        }
        else {
          if (!(RF_opt & OPT_OUTC_TYPE)) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  NA encountered for VIMP outcome in terminal node:  %10d", terminalNode -> leafCount);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
      else {
        if (!ensembleFlag) {
          RF_vimpEnsembleDen[p][i] = 0;
        }
      }
    }  
  }  
}
void summarizeVimpPerformance(uint       mode,
                              uint       treeID) {
  uint      obsSize;
  uint      varCount;
  double    **responsePtr;
  double    **importancePtr;
  uint      *denomPtr;
  uint      *vimpDenomPtr;
  double   ***crVimpMortality;
  double    **subVimpOutcomePtr;
  char        responseImputeFlag;  
  double    **treeOutcome;
  double  ****chfTreeEnsemble;
  double   ***mcTreeEnsemble;
  double    **crTreeMortality;
  double    **subTreeOutcomePtr;
  double    maxValue, maxClass;
  uint      normalizationFlag;
  uint i, j, p;
  if (!(RF_opt & OPT_VIMP)) {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  VIMP analysis requested while OPT bit not set.  \n");
    return;
  }
  crVimpMortality    = NULL;  
  crTreeMortality    = NULL;  
  subVimpOutcomePtr  = NULL;
  treeOutcome        = NULL;  
  chfTreeEnsemble    = NULL;  
  mcTreeEnsemble     = NULL;  
  responseImputeFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    if (RF_opt & OPT_VIMP_JOIN) {
      varCount = 1;
    }
    else {
      varCount = RF_intrPredictorSize;
    }
    if (RF_opt & OPT_VIMP_LEOB) {
      responsePtr = RF_fresponse[treeID];
      denomPtr = RF_testMembershipFlag;
    }
    else {
      responsePtr = RF_fresponseIn;
      denomPtr = RF_fullEnsembleDen;
      if (RF_fmRecordSize > 0) {
        if(RF_fmResponseFlag == TRUE) {
          responseImputeFlag = TRUE;
        }
      }
    }
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_VIMP_JOIN) {
      varCount = 1;
    }
    else {
      varCount = RF_intrPredictorSize;
    }
    if (RF_opt & OPT_VIMP_LEOB) {
      responsePtr = RF_response[treeID];
      denomPtr = RF_oobMembershipFlag[treeID];
    }
    else {
      responsePtr = RF_responseIn;
      denomPtr = RF_oobEnsembleDen;
      if (RF_mRecordSize > 0) {
        if(RF_mResponseFlag == TRUE) {
          responseImputeFlag = TRUE;
        }
      }
    }
    break;
  }
  responsePtr = stackAndImputeGenericResponse(responseImputeFlag, mode, RF_rSize, obsSize, 0, RF_forestSize, responsePtr);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      crVimpMortality = dmatrix3(1, varCount, 1, RF_eventTypeSize, 1, obsSize);
      if (RF_opt & OPT_VIMP_LEOB) {
        crTreeMortality = dmatrix(1, RF_eventTypeSize, 1, obsSize);
      }
    }
  }  
  if (RF_opt & OPT_VIMP_LEOB) {
    stackTreeEnsemble(mode, &treeOutcome, &chfTreeEnsemble, &mcTreeEnsemble);
    updateTreeEnsemble(mode, treeID, treeOutcome, chfTreeEnsemble, mcTreeEnsemble);
    importancePtr = RF_vimpLeo[treeID];
  }
  else {
    importancePtr = RF_importancePtr;
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    for (p = 1; p <= varCount; p++) {
      vimpDenomPtr = RF_vimpEnsembleDen[p];
      if (RF_eventTypeSize == 1) {
        getEnsembleMortality(mode, treeID, obsSize, RF_sVimpEnsemble[p], vimpDenomPtr, RF_vimpOutcome[p]);
      }
      else {
        getEnsembleMortalityCR(mode, treeID, obsSize, RF_sVimpEnsemble[p], vimpDenomPtr, crVimpMortality[p]);
      }
    }  
    if (RF_opt & OPT_VIMP_LEOB) {
      if (RF_eventTypeSize == 1) {
        getEnsembleMortality(mode, treeID, obsSize, chfTreeEnsemble[1], denomPtr, treeOutcome[1]);
      }
      else {
        getEnsembleMortalityCR(mode, treeID, obsSize, chfTreeEnsemble[1], denomPtr, crTreeMortality);
      }
    }
  }  
  else {
    if (RF_rFactorCount > 0) {
      for (p=1; p <= varCount; p++) {
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
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (RF_eventTypeSize > 1) {
        subTreeOutcomePtr = crTreeMortality;
      }
      else {
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
                   treeOutcome[1],
                   subTreeOutcomePtr,
                   denomPtr,
                   RF_perfLeo[treeID]);
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
  for (p=1; p <= varCount; p++) {
    vimpDenomPtr = RF_vimpEnsembleDen[p];
    if (normalizationFlag == TRUE) {
      for (i = 1; i <= obsSize; i++) {
        if (vimpDenomPtr[i] != 0) {
          RF_vimpOutcome[p][i] = RF_vimpOutcome[p][i] / vimpDenomPtr[i];
        }
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (RF_eventTypeSize > 1) {
        subVimpOutcomePtr = crVimpMortality[p];
      }
      else {
        subVimpOutcomePtr = NULL;
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
    getPerformance(treeID,
                   mode, 
                   obsSize, 
                   responsePtr, 
                   RF_vimpOutcome[p],
                   subVimpOutcomePtr,
                   vimpDenomPtr,
                   importancePtr[p]);
  }  
  unstackImputeResponse(responseImputeFlag, RF_rSize, obsSize, responsePtr);
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      free_dmatrix3(crVimpMortality, 1, varCount, 1, RF_eventTypeSize, 1, obsSize);
      if (RF_opt & OPT_VIMP_LEOB) {
        free_dmatrix(crTreeMortality, 1, RF_eventTypeSize, 1, obsSize);
      }
    }
  }
  if (RF_opt & OPT_VIMP_LEOB) {
    unstackTreeEnsemble(mode, treeOutcome, chfTreeEnsemble, mcTreeEnsemble);
  }
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
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    perfDimOne = RF_eventTypeSize;
  }
  else {
    if (RF_rFactorCount > 0) {
      perfDimOne = RF_rFactorSize[1] + 1;
    }
    else {
      perfDimOne = 1;
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
Node ***stackVimpMembership(uint mode) {
  Node ***membership;
  uint obsSize;
  uint varCount;
  membership = NULL;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        varCount = 1;
      }
      else {
        varCount = RF_intrPredictorSize;
      }
      break;
    default:
      obsSize = RF_observationSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        varCount = 1;
      }
      else {
        varCount = RF_intrPredictorSize;
      }
      break;
    }
    membership = (Node ***) vmatrix(1, varCount, 1, obsSize);
  }
  return membership;
}
void unstackVimpMembership(uint mode, Node ***membership) {
  uint obsSize;
  uint varCount;
  if (RF_opt & OPT_VIMP) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        varCount = 1;
      }
      else {
        varCount = RF_intrPredictorSize;
      }
      break;
    default:
      obsSize = RF_observationSize;
      if (RF_opt & OPT_VIMP_JOIN) {
        varCount = 1;
      }
      else {
        varCount = RF_intrPredictorSize;
      }
      break;
    }
    free_vmatrix((void **) membership, 1, varCount, 1, obsSize);
  }
}
void stackTreeEnsemble(uint         mode,
                       double    ***treeOutcome,
                       double  *****chfTreeEnsemble,
                       double   ****mcTreeEnsemble) {
  uint obsSize;
  uint ensbDimTwo;
  uint i,j,k;
  *chfTreeEnsemble  = NULL;
  *mcTreeEnsemble   = NULL;
  if (RF_timeIndex > 0) {
    ensbDimTwo = RF_sortedTimeInterestSize;
  }
  else {
    if (RF_rFactorCount > 0) {
      ensbDimTwo = RF_rFactorSize[1];
    }
    else {
      ensbDimTwo = 1;
    }
  }
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
    *chfTreeEnsemble  = dmatrix4(1, 1, 1, RF_eventTypeSize, 1, ensbDimTwo, 1, obsSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      for (i = 1; i <= obsSize; i++) {
        for (k = 1; k <= ensbDimTwo; k++) {
          (*chfTreeEnsemble)[1][j][k][i] = 0.0;
        }
      }
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      *mcTreeEnsemble = dmatrix3(1, 1, 1, ensbDimTwo, 1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        for (k = 1; k <= ensbDimTwo; k++) {
          (*mcTreeEnsemble)[1][k][i] = 0.0;
        }
      }
    }
  }
}
void unstackTreeEnsemble(uint        mode,
                         double    **treeOutcome,
                         double  ****chfTreeEnsemble,
                         double   ***mcTreeEnsemble) {
  uint obsSize;
  uint ensbDimTwo;
  if (RF_timeIndex > 0) {
    ensbDimTwo = RF_sortedTimeInterestSize;
  }
  else {
    if (RF_rFactorCount > 0) {
      ensbDimTwo = RF_rFactorSize[1];
    }
    else {
      ensbDimTwo = 1;
    }
  }
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
    free_dmatrix4(chfTreeEnsemble, 1, 1, 1, RF_eventTypeSize, 1, ensbDimTwo, 1, obsSize);
  }
  else {
    if (RF_rFactorCount > 0) {
      free_dmatrix3(mcTreeEnsemble, 1, 1, 1, ensbDimTwo, 1, obsSize);
    }
  }
}
