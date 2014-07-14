////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.4
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
#include          "impute.h"
#include "stackPreDefined.h"
void stackIncomingResponseArrays(uint mode) {
  uint i, j;
  RF_timeIndex = RF_statusIndex = 0;
  if (RF_rSize > 0) {
    RF_rType = (char **) new_vvector(1, RF_rSize, NRUTIL_CPTR);
    RF_yIndex = uivector(1, RF_rSize);
    j = 0;
    for (i = 1; i <= RF_rSize; i++) {
      RF_rType[i] = (char*) CHAR(STRING_ELT(AS_CHARACTER(RF_sexp_rType), i-1));
      if ((strcmp(RF_rType[i], "C") != 0) &&
          (strcmp(RF_rType[i], "I") != 0) &&
          (strcmp(RF_rType[i], "R") != 0) &&
          (strcmp(RF_rType[i], "T") != 0) &&
          (strcmp(RF_rType[i], "S") != 0)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Invalid type:  [%10d] = %2s", i, RF_rType[i]);
        Rprintf("\nRF-SRC:  Outcomes must be 'C', 'I', 'R', 'T', or 'S'.");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      RF_yIndex[i] = 0;
      if (strcmp(RF_rType[i], "T") == 0) {
        RF_timeIndex = i;
      }
      else if (strcmp(RF_rType[i], "S") == 0) {
        RF_statusIndex = i;
      }
      else {
        RF_yIndex[++j] = i;
      }
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        if (RF_rSize != RF_frSize) {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  TRAIN and TEST outcome/response matrices must be of the same dimension.  ");
          Rprintf("\nRF-SRC:  TRAIN vs TEST:  %10d vs %10d  ", RF_rSize, RF_frSize);
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
      else {
        if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB) | (RF_opt & OPT_VIMP)) {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  TEST outcome/response matrix must be present when PERF or VIMP is requested.  ");
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_ySize = 0;
      RF_rTarget = 0;
      RF_ptnCount = 0;
    }
    else {
      RF_ySize = RF_rSize - ((RF_timeIndex == 0) ? 0:1) - ((RF_statusIndex == 0) ? 0:1);
      if (RF_ySize != RF_rSize) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Responses must be [S], [C], [R], [R+], [C+], [M+].  ");
        Rprintf("\nRF-SRC:  The application will now exit.\n");
        error("\nRF-SRC:  The application will now exit.\n");
      }
      if (RF_ySize == 1) {
        RF_rTarget = 1;
      }
      else {
      }
    }
    RF_responseIn = (double **) new_vvector(1, RF_rSize, NRUTIL_DPTR);
    for (i=1; i <= RF_rSize; i++) {
      RF_responseIn[i] = (RF_rData + ((i-1) * RF_observationSize) - 1);
    }
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        RF_fresponseIn = (double **) new_vvector(1, RF_frSize, NRUTIL_DPTR);
        for (i=1; i <= RF_rSize; i++) {
          RF_fresponseIn[i] = (RF_frData + ((i-1) * RF_fobservationSize) - 1);
        }
      }
      else {
        RF_fresponseIn = NULL;
      }
    }
  }
  else {
    RF_rType      = NULL;
    RF_responseIn = NULL;
  }
}
void unstackIncomingResponseArrays(uint mode) {
  if (RF_rSize > 0) {
    free_new_vvector(RF_rType, 1, RF_rSize, NRUTIL_CPTR);
    free_uivector(RF_yIndex, 1, RF_rSize);
    free_new_vvector(RF_responseIn, 1, RF_rSize, NRUTIL_DPTR);
    if (mode == RF_PRED) {
      if (RF_frSize > 0) {
        free_new_vvector(RF_fresponseIn, 1, RF_frSize, NRUTIL_DPTR);
      }
    }
  }
}
void stackIncomingCovariateArrays(uint mode) {
  uint i;
  RF_xType = (char **) new_vvector(1, RF_xSize, NRUTIL_CPTR);
  for (i = 1; i <= RF_xSize; i++) {
    RF_xType[i] = (char*) CHAR(STRING_ELT(AS_CHARACTER(RF_sexp_xType), i-1));
    if ((strcmp(RF_xType[i], "C") != 0) &&
        (strcmp(RF_xType[i], "I") != 0) &&
        (strcmp(RF_xType[i], "R") != 0)) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Invalid type:  [%10d] = %2s", i, RF_xType[i]);
      Rprintf("\nRF-SRC:  Predictors must be 'C', 'I', or 'R'.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  RF_observationIn = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
  for (i=1; i <= RF_xSize; i++) {
    RF_observationIn[i] = (RF_xData + ((i-1) * RF_observationSize) - 1);
  }
  if (mode == RF_PRED) {
    RF_fobservationIn = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
    for (i=1; i <= RF_xSize; i++) {
      RF_fobservationIn[i] = (RF_fxData + ((i-1) * RF_fobservationSize) - 1);
    }
  }
}
void unstackIncomingCovariateArrays(uint mode) {
  free_new_vvector(RF_xType, 1, RF_xSize, NRUTIL_CPTR);
  free_new_vvector(RF_observationIn, 1, RF_xSize, NRUTIL_DPTR);
  if (mode == RF_PRED) {
    free_new_vvector(RF_fobservationIn, 1, RF_xSize, NRUTIL_DPTR);
  }
}
void stackIncomingArrays(uint mode) {
  stackIncomingResponseArrays(mode);
  stackIncomingCovariateArrays(mode);
  if (mode == RF_GROW) {
    if ((RF_timeIndex == 0) && (RF_statusIndex == 0)) {
      if ((RF_splitRule != RAND_SPLIT)  &&
          (RF_splitRule != REGR_WT_NRM) &&
          (RF_splitRule != REGR_WT_OFF) &&
          (RF_splitRule != REGR_WT_HVY) &&
          (RF_splitRule != CLAS_WT_NRM) &&
          (RF_splitRule != CLAS_WT_OFF) &&
          (RF_splitRule != CLAS_WT_HVY) &&
          (RF_splitRule != MVRG_SPLIT)  &&
          (RF_splitRule != MVCL_SPLIT)  &&
          (RF_splitRule != USPV_SPLIT)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  !SURV data and split rule specified are incompatible.");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else if ((RF_timeIndex != 0) && (RF_statusIndex != 0)) {
      if ((RF_splitRule != SURV_LGRNK)  &&
          (RF_splitRule != SURV_CR_LOG) &&
          (RF_splitRule != SURV_LRSCR)  &&
          (RF_splitRule != SURV_CR_LAU) &&
          (RF_splitRule != RAND_SPLIT)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  SURV data and split rule specified are incompatible.");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Data set contains mixed outcomes with no comatible split rule.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void unstackIncomingArrays(uint mode) {
  unstackIncomingResponseArrays(mode);
  unstackIncomingCovariateArrays(mode);
}
void stackPreDefinedCommonArrays() {
  uint i;
  RF_tNodeMembership = (Node ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_bootMembershipIndex = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  RF_bootMembershipFlag = (char **) new_vvector(1, RF_forestSize, NRUTIL_CPTR);
  RF_bootMembershipCount = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
  RF_oobMembershipFlag = (char **) new_vvector(1, RF_forestSize, NRUTIL_CPTR);
  RF_tNodeList = (Node ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_identityMembershipIndex = uivector(1, RF_observationSize);
  for (i = 1; i <= RF_observationSize; i++) {
    RF_identityMembershipIndex[i] = i;
  }
  RF_oobSize = uivector(1, RF_forestSize);
  RF_maxDepth = uivector(1, RF_forestSize);
  RF_serialTreeIndex = uivector(1, RF_forestSize);
  if (RF_timeIndex > 0) {
    RF_masterTime  = dvector(1, RF_observationSize);
    RF_masterTimeIndexIn  = uivector(1, RF_observationSize);
  }
  RF_root = (Node **) new_vvector(1, RF_forestSize, NRUTIL_NPTR);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_root[i] = NULL;
  }
  if (RF_ptnCount > 0) {
    RF_pNodeMembership = (Node ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pNodeList = (Node ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
    RF_pLeafCount = uivector(1, RF_forestSize);
  }
  RF_orderedLeafCount = uivector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_orderedLeafCount[i] = 0;
  }
}
void unstackPreDefinedCommonArrays() {
  free_new_vvector(RF_tNodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_new_vvector(RF_bootMembershipIndex, 1, RF_forestSize, NRUTIL_UPTR);
  free_new_vvector(RF_bootMembershipFlag, 1, RF_forestSize, NRUTIL_CPTR);
  free_new_vvector(RF_bootMembershipCount, 1, RF_forestSize, NRUTIL_UPTR);
  free_new_vvector(RF_oobMembershipFlag, 1, RF_forestSize, NRUTIL_CPTR);
  free_new_vvector(RF_tNodeList, 1, RF_forestSize, NRUTIL_NPTR2);
  free_uivector(RF_identityMembershipIndex, 1, RF_observationSize);
  free_uivector(RF_oobSize, 1, RF_forestSize);
  free_uivector(RF_maxDepth, 1, RF_forestSize);
  free_uivector(RF_serialTreeIndex, 1, RF_forestSize);
  if (RF_timeIndex > 0) {
    free_dvector(RF_masterTime, 1, RF_observationSize);
    free_uivector(RF_masterTimeIndexIn, 1, RF_observationSize);
  }
  free_new_vvector(RF_root, 1, RF_forestSize, NRUTIL_NPTR);
  if (RF_ptnCount > 0) {
    free_new_vvector(RF_pNodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
    free_new_vvector(RF_pNodeList, 1, RF_forestSize, NRUTIL_NPTR2);
    free_uivector(RF_pLeafCount, 1, RF_forestSize);
  }
  free_uivector(RF_orderedLeafCount, 1, RF_forestSize);
}
void stackPreDefinedGrowthArrays() {
  char uniformFlag, integerFlag;
  double meanWeight;
  uint i;
  if (RF_opt & OPT_TREE) {
    RF_nodeCount = uivector(1, RF_forestSize);
    RF_mwcpCount = uivector(1, RF_forestSize);
    for (i = 1; i <= RF_forestSize; i++) {
      RF_nodeCount[i] = RF_mwcpCount[i] = 0;
    }
  }
  if (RF_opt & OPT_VIMP) {
    RF_intrPredictor = uivector(1, RF_intrPredictorSize);
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_intrPredictor[i] = i;
    }
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = TRUE;
    }
  }
  RF_xWeightSorted = NULL;
  RF_xWeightDensity = NULL;
  RF_xWeightDensitySize = 0;
  meanWeight = getMeanValue(RF_xWeight, RF_xSize);
  uniformFlag = TRUE;
  i = 0;
  while (uniformFlag && (i < RF_xSize)) {
    ++i;
    if (fabs(RF_xWeight[i] - meanWeight) > EPSILON) {
      uniformFlag = FALSE;
    }
  }
  if (uniformFlag) {
    RF_xWeightType = RF_WGHT_UNIFORM;
  }
  else {
    integerFlag = TRUE;
    i = 0;
    while (integerFlag && (i < RF_xSize)) {
      i++;
      if (fabs(round(RF_xWeight[i]) - RF_xWeight[i]) > EPSILON) {
        integerFlag = FALSE;
      }
    }
    if(integerFlag) {
      RF_xWeightType = RF_WGHT_INTEGER;
    }
    else {
      RF_xWeightType = RF_WGHT_GENERIC;
    }
  }
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    RF_xWeightSorted = uivector(1, RF_xSize);
    indexx(RF_xSize, RF_xWeight, RF_xWeightSorted);
    RF_xWeightDensitySize = 0;
    for (i = 1; i <= RF_xSize; i++) {
      RF_xWeightDensitySize += (uint) RF_xWeight[i];
    }
    break;
  case RF_WGHT_GENERIC:
    RF_xWeightSorted = uivector(1, RF_xSize);
    indexx(RF_xSize, RF_xWeight, RF_xWeightSorted);
    break;
  }
  RF_gobservationIndv = uivector(1, RF_observationSize);
  for (i = 1; i <= RF_observationSize; i++) {
    RF_gobservationIndv[i] = i;
  }
}
void unstackPreDefinedGrowthArrays() {
  if (RF_opt & OPT_TREE) {
    free_uivector(RF_nodeCount, 1, RF_forestSize);
    free_uivector(RF_mwcpCount, 1, RF_forestSize);
  }
  if (RF_opt & OPT_VIMP) {
    free_uivector(RF_intrPredictor, 1, RF_intrPredictorSize);
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  switch (RF_xWeightType) {
  case RF_WGHT_UNIFORM:
    break;
  case RF_WGHT_INTEGER:
    free_uivector(RF_xWeightSorted, 1, RF_xSize);
    break;
  case RF_WGHT_GENERIC:
    free_uivector(RF_xWeightSorted, 1, RF_xSize);
    break;
  }
  free_uivector(RF_gobservationIndv, 1, RF_observationSize);
}
void stackPreDefinedRestoreArrays() {
  uint i;
  RF_nodeCount = uivector(1, RF_forestSize);
  RF_mwcpCount = uivector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_nodeCount[i] = RF_mwcpCount[i] = 0;
  }
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
  if(RF_sobservationSize > 0) {
    RF_soobSize = uivector(1, RF_forestSize);
  }
  else {
    RF_gobservationIndv = uivector(1, RF_observationSize);
    for (i = 1; i <= RF_observationSize; i++) {
      RF_gobservationIndv[i] = i;
    }
  }
}
void unstackPreDefinedRestoreArrays() {
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  if(RF_sobservationSize > 0) {
    free_uivector(RF_soobSize, 1, RF_forestSize);
  }
  else {
    free_uivector(RF_gobservationIndv, 1, RF_observationSize);
  }
}
void stackPreDefinedPredictArrays() {
  uint i;
  RF_ftNodeMembership = (Node ***) new_vvector(1, RF_forestSize, NRUTIL_NPTR2);
  RF_testMembershipFlag = cvector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_testMembershipFlag[i] = ACTIVE;
  }
  RF_nodeCount = uivector(1, RF_forestSize);
  RF_mwcpCount = uivector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_nodeCount[i] = RF_mwcpCount[i] = 0;
  }
  if (RF_opt & OPT_VIMP) {
    checkInteraction();
    RF_importanceFlag = cvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      RF_importanceFlag[i] = FALSE;
    }
    for (i = 1; i <= RF_intrPredictorSize; i++) {
      RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
    }
  }
  RF_gobservationIndv = uivector(1, RF_fobservationSize);
  for (i = 1; i <= RF_fobservationSize; i++) {
    RF_gobservationIndv[i] = i;
  }
}
void unstackPreDefinedPredictArrays() {
  free_new_vvector(RF_ftNodeMembership, 1, RF_forestSize, NRUTIL_NPTR2);
  free_cvector(RF_testMembershipFlag, 1, RF_fobservationSize);
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  if (RF_opt & OPT_VIMP) {
    free_cvector(RF_importanceFlag, 1, RF_xSize);
  }
  free_uivector(RF_gobservationIndv, 1, RF_fobservationSize);
}
void checkInteraction() {
  uint leadingIndex, i;
  if((RF_intrPredictorSize <= 0) || (RF_intrPredictorSize > RF_xSize)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of predictors to be perturbed must be greater than zero and less than %10d:  %10d \n", RF_xSize, RF_intrPredictorSize);
    error("\nRF-SRC:  The application will now exit.\n");
  }
  uint *intrPredictorCopy = uivector(1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    intrPredictorCopy[i] = RF_intrPredictor[i];
  }
  hpsortui(intrPredictorCopy, RF_intrPredictorSize);
  leadingIndex = 1;
  for (i=2; i <= RF_intrPredictorSize; i++) {
    if (intrPredictorCopy[i] > intrPredictorCopy[leadingIndex]) {
      leadingIndex++;
    }
  }
  if (RF_intrPredictorSize != leadingIndex) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Interaction terms are not unique.");
    Rprintf("\nRF-SRC:  Only %10d of %10d are unique.", leadingIndex, RF_intrPredictorSize);
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_uivector(intrPredictorCopy, 1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    if (RF_intrPredictor[i] > RF_xSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Parameter verification failed.");
      Rprintf("\nRF-SRC:  Interaction terms are not coherent.");
      Rprintf("\nRF-SRC:  Predictor encountered is %10d, maximum allowable is %10d.", RF_intrPredictor[i], RF_xSize);
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
