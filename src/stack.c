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


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include       "nodeOps.h"
#include     "factorOps.h"
#include        "impute.h"
#include      "treeUtil.h"
#include         "stack.h"
void stackIncomingArrays(uint mode) {
  uint i, j;
  RF_timeIndex = RF_statusIndex = 0;
  RF_rType = (char**) vvector(1, RF_rSize);
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
      if ((RF_opt & OPT_PERF) || (RF_opt & OPT_VIMP)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  TEST outcome/response matrix must be present when PERF or VIMP is requested.  ");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    RF_ySize = 0;
  }
  else {
    RF_ySize = RF_rSize - ((RF_timeIndex == 0) ? 0:1) - ((RF_statusIndex == 0) ? 0:1);
  }
  RF_responseIn = (double **) vvector(1, RF_rSize);
  for (i=1; i <= RF_rSize; i++) {
    RF_responseIn[i] = (RF_rData + ((i-1) * RF_observationSize) - 1);
  }
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      RF_fresponseIn = (double **) vvector(1, RF_frSize);
      for (i=1; i <= RF_rSize; i++) {
        RF_fresponseIn[i] = (RF_frData + ((i-1) * RF_fobservationSize) - 1);
      }
    }
  }
  RF_xType = (char**) vvector(1, RF_xSize);
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
  RF_observationIn = (double **) vvector(1, RF_xSize);
  for (i=1; i <= RF_xSize; i++) {
    RF_observationIn[i] = (RF_xData + ((i-1) * RF_observationSize) - 1);
  }
  if (mode == RF_PRED) {
    RF_fobservationIn = (double **) vvector(1, RF_xSize);
    for (i=1; i <= RF_xSize; i++) {
      RF_fobservationIn[i] = (RF_fxData + ((i-1) * RF_fobservationSize) - 1);
    }
  }
  if (mode == RF_GROW) {
    if ((RF_timeIndex == 0) && (RF_statusIndex == 0)) {
      if ((RF_splitRule != REGR_SPLIT) &&
          (RF_splitRule != CLAS_SPLIT) &&
          (RF_splitRule != RAND_SPLIT)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  !SURV data and split rule specified are incompatible.");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    else if ((RF_timeIndex != 0) && (RF_statusIndex != 0)) {
      if ((RF_splitRule != SURV_LGRNK)  &&
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
  free_vvector(RF_rType, 1, RF_rSize);
  free_uivector(RF_yIndex, 1, RF_rSize);
  free_vvector(RF_responseIn, 1, RF_rSize);
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      free_vvector(RF_fresponseIn, 1, RF_frSize);
    }
  }
  free_vvector(RF_xType, 1, RF_xSize);
  free_vvector(RF_observationIn, 1, RF_xSize);
  if (mode == RF_PRED) {
    free_vvector(RF_fobservationIn, 1, RF_xSize);
  }
}
void checkInteraction() {
  uint leadingIndex, i;
  if((RF_intrPredictorSize <= 0) || (RF_intrPredictorSize > RF_xSize)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of predictors to be perturbed must be greater than zero and less than %10d:  %10d \n", RF_xSize, RF_intrPredictorSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
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
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_uivector(intrPredictorCopy, 1, RF_intrPredictorSize);
  for (i=1; i <= RF_intrPredictorSize; i++) {
    if (RF_intrPredictor[i] > RF_xSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Parameter verification failed.");
      Rprintf("\nRF-SRC:  Interaction terms are not coherent.");
      Rprintf("\nRF-SRC:  Predictor encountered is %10d, maximum allowable is %10d.", RF_intrPredictor[i], RF_xSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
}
void stackPreDefinedCommonArrays() {
  uint i;
  RF_nodeMembership = (Node ***) vvector(1, RF_forestSize);
  RF_bootMembershipIndex = (uint **) vvector(1, RF_forestSize);
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    RF_trivialBootMembershipIndex = uivector(1, RF_observationSize);
    for (i = 1; i <= RF_observationSize; i++) {
      RF_trivialBootMembershipIndex[i] = i;
    }
  }
  RF_bootMembershipFlag = (uint **) vvector(1, RF_forestSize);
  RF_oobMembershipFlag = (uint **) vvector(1, RF_forestSize);
  RF_terminalNode = (Node ***) vvector(1, RF_forestSize);
  RF_oobSize = uivector(1, RF_forestSize);
  RF_serialTreeIndex = uivector(1, RF_forestSize);
  if (RF_timeIndex > 0) {
    RF_masterTime  = dvector(1, RF_observationSize);
    RF_masterTimeIndexIn  = uivector(1, RF_observationSize);
  }
  RF_root = (Node **) vvector(1, RF_forestSize);
  for (i = 1; i <= RF_forestSize; i++) {
    RF_root[i] = NULL;
  }
}
void unstackPreDefinedCommonArrays() {
  free_vvector(RF_nodeMembership, 1, RF_forestSize);
  free_vvector(RF_bootMembershipIndex, 1, RF_forestSize);
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    free_uivector(RF_trivialBootMembershipIndex, 1, RF_observationSize);
  }
  free_vvector(RF_bootMembershipFlag, 1, RF_forestSize);
  free_vvector(RF_oobMembershipFlag, 1, RF_forestSize);
  free_vvector(RF_terminalNode, 1, RF_forestSize);
  free_uivector(RF_oobSize, 1, RF_forestSize);
  free_uivector(RF_serialTreeIndex, 1, RF_forestSize);
  if (RF_timeIndex > 0) {
    free_dvector(RF_masterTime, 1, RF_observationSize);
    free_uivector(RF_masterTimeIndexIn, 1, RF_observationSize);
  }
  free_vvector((Node **) RF_root, 1, RF_forestSize);
}
void stackPreDefinedGrowthArrays() {
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
    if (RF_opt & OPT_VIMP_TYPE) {
      RF_importanceFlag = cvector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        RF_importanceFlag[i] = TRUE;
      }
    }
  }
}
void unstackPreDefinedGrowthArrays() {
  if (RF_opt & OPT_TREE) {
    free_uivector(RF_nodeCount, 1, RF_forestSize);
    free_uivector(RF_mwcpCount, 1, RF_forestSize);
  }
  if (RF_opt & OPT_VIMP) {
    free_uivector(RF_intrPredictor, 1, RF_intrPredictorSize);
    if (RF_opt & OPT_VIMP_TYPE) {
      free_cvector(RF_importanceFlag, 1, RF_xSize);
    }
  }
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
    if (RF_opt & OPT_VIMP_TYPE) {
      RF_importanceFlag = cvector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        RF_importanceFlag[i] = FALSE;
      }
      for (i = 1; i <= RF_intrPredictorSize; i++) {
        RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
      }
    }
  }
}
void unstackPreDefinedRestoreArrays() {
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_TYPE) {
      free_cvector(RF_importanceFlag, 1, RF_xSize);
    }
  }
}
void stackPreDefinedPredictArrays() {
  uint i;
  RF_fnodeMembership = (Node ***) vvector(1, RF_forestSize);
  RF_testMembershipFlag = uivector(1, RF_fobservationSize);
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
    if (RF_opt & OPT_VIMP_TYPE) {
      RF_importanceFlag = cvector(1, RF_xSize);
      for (i = 1; i <= RF_xSize; i++) {
        RF_importanceFlag[i] = FALSE;
      }
      for (i = 1; i <= RF_intrPredictorSize; i++) {
        RF_importanceFlag[RF_intrPredictor[i]] = TRUE;
      }
    }
  }
}
void unstackPreDefinedPredictArrays() {
  free_vvector(RF_fnodeMembership, 1, RF_forestSize);
  free_uivector(RF_testMembershipFlag, 1, RF_fobservationSize);
  free_uivector(RF_nodeCount, 1, RF_forestSize);
  free_uivector(RF_mwcpCount, 1, RF_forestSize);
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_TYPE) {
      free_cvector(RF_importanceFlag, 1, RF_xSize);
    }
  }
}
void initializeArrays(char mode) {
  uint i, j;
  uint leadingIndex;
  if (RF_timeIndex > 0) {
    RF_masterTimeSize = 0;
    for (j = 1; j <= RF_observationSize; j++) {
      if (!ISNA(RF_responseIn[RF_timeIndex][j])) {
        RF_masterTimeSize ++;
        RF_masterTime[RF_masterTimeSize] = RF_responseIn[RF_timeIndex][j];
      }
    }
    hpsort(RF_masterTime, RF_masterTimeSize);
    leadingIndex = 1;
    for (i=2; i <= RF_masterTimeSize; i++) {
      if (RF_masterTime[i] > RF_masterTime[leadingIndex]) {
        leadingIndex++;
        RF_masterTime[leadingIndex] = RF_masterTime[i];
      }
    }
    RF_masterTimeSize = leadingIndex;
    for (i= RF_masterTimeSize + 1; i <= RF_observationSize; i++) {
      RF_masterTime[i] = 0;
    }
    if (!(RF_opt & OPT_IMPU_ONLY)) {
      hpsort(RF_timeInterest, RF_timeInterestSize);
      RF_sortedTimeInterestSize = 1;
      for (i=2; i <= RF_timeInterestSize; i++) {
        if (RF_timeInterest[i] > RF_timeInterest[RF_sortedTimeInterestSize]) {
          (RF_sortedTimeInterestSize) ++;
          RF_timeInterest[RF_sortedTimeInterestSize] = RF_timeInterest[i];
        }
      }
      if (RF_sortedTimeInterestSize != RF_timeInterestSize) {
        Rprintf("\nRFsrc:  *** WARNING *** ");
        Rprintf("\nRFsrc:  Time points of interest are not unique.");
        Rprintf("\nRFsrc:  The ensemble estimate output matrix is being");
        Rprintf("\nRFsrc:  resized as [N'] x [n], where N' is the");
        Rprintf("\nRFsrc:  unique time points of interest and n is");
        Rprintf("\nRFsrc:  number of observations in the data.");
      }
      for (i = RF_sortedTimeInterestSize + 1; i <= RF_timeInterestSize; i++) {
        RF_timeInterest[i] = 0;
      }
    }
  }
}
void stackFactorArrays() {
  stackFactorGeneric(RF_rSize, 
                     RF_rType, 
                     &RF_rFactorMap,
                     &RF_rFactorCount,
                     &RF_rFactorIndex,
                     &RF_rFactorSize);
  stackFactorGeneric(RF_xSize, 
                     RF_xType, 
                     &RF_xFactorMap,
                     &RF_xFactorCount,
                     &RF_xFactorIndex,
                     &RF_xFactorSize);
}
void stackFactorGeneric(uint    size, 
                        char  **type, 
                        uint  **p_factorMap,
                        uint   *factorCount,
                        uint  **p_factorIndex,
                        uint  **p_factorSize) {
  uint i, j;
  *p_factorMap = uivector(1, size);
  *factorCount = 0;
  for (i = 1; i <= size; i++) {
    (*p_factorMap)[i] = 0;
    if (strcmp(type[i], "C") == 0) {
      (*factorCount) ++;
      (*p_factorMap)[i] = *factorCount;
    }
  }
  if (*factorCount > 0) {
    *p_factorIndex = uivector(1, *factorCount);
    j = 0;
    for (i = 1; i <= size; i++) {
      if ((*p_factorMap)[i] > 0) {
        (*p_factorIndex)[++j] = i;
      }
    }
    *p_factorSize = uivector(1, *factorCount);
  }
}
void unstackFactorArrays() {
  uint j, k;
  free_uivector(RF_rFactorMap, 1, RF_rSize);
  if (RF_rFactorCount > 0) {
    free_uivector(RF_rFactorIndex, 1, RF_rFactorCount);
    free_uivector(RF_rFactorSize, 1, RF_rFactorCount);
  }
  free_uivector(RF_xFactorMap, 1, RF_xSize);
  if (RF_xFactorCount > 0) {
    free_uivector(RF_xFactorIndex, 1, RF_xFactorCount);
    free_uivector(RF_xFactorSize, 1, RF_xFactorCount);
  }
  if ((RF_rFactorCount + RF_xFactorCount) > 0) {
    for (j = 1; j <= RF_forestSize; j++) {
      if (RF_factorList[j] != NULL) {
        for (k = 1; k <= RF_maxFactorLevel; k++) {
          if (RF_factorList[j][k] != NULL) {
            free_Factor(RF_factorList[j][k]);
          }
        }
        free_vvector(RF_factorList[j], 1, RF_maxFactorLevel);
      }
    }
    free_vvector(RF_factorList, 1, RF_forestSize);
  }
}
char stackMissingArrays(char mode) {
  char result;
  char mFlag;
  char dualUseFlag;
  uint recordSize;
  uint vSize;
  uint i, j, k;
  result = TRUE;
  for (j = 1; j <= RF_rSize; j++) {
    if (j == RF_timeIndex) {
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_timeIndex][i])) {
          if (RF_responseIn[RF_timeIndex][i] < 0) {
            result = FALSE;
            Rprintf("\nRF-SRC:  TRAINING time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_timeIndex][i]);
          }
        }
      }
    }
    if (j == RF_statusIndex) {
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_statusIndex][i])) {
          if (RF_responseIn[RF_statusIndex][i] < 0) {
            result = FALSE;
            Rprintf("\nRF-SRC:  TRAINING status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_responseIn[RF_statusIndex][i]);
          }
        }
      }
    }
    if (j == RF_statusIndex) {
      mFlag = FALSE;
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[RF_statusIndex][i])) {
          if (RF_responseIn[RF_statusIndex][i] >= 0) {
            mFlag = TRUE;
            i = RF_observationSize;
          }
        }
      }
      if (mFlag == FALSE) {
        Rprintf("\nRF-SRC:  All TRAINING status elements are censored or missing. \n");
        result = FALSE;
      }
    }
    if ((RF_statusIndex == 0) && (RF_timeIndex == 0)) {
      mFlag = FALSE;
      for (i = 1; i <= RF_observationSize; i++) {
        if (!ISNA(RF_responseIn[j][i])) {
          mFlag = TRUE;
          i = RF_observationSize;
        }
      }
      if (mFlag == FALSE) {
        Rprintf("\nRF-SRC:  All TRAINING outcome/response elements are missing for:  %10d \n", j);
        result = FALSE;
      }
    }
    if (result == FALSE) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Missingness verification failed.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }  
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      if (RF_timeIndex > 0) {
        for (i = 1 ; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_timeIndex][i])) {
            if (RF_fresponseIn[RF_timeIndex][i] < 0) {
              result = FALSE;
              Rprintf("\nRF-SRC:  PRED time elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_timeIndex][i]);
            }
          }
        }
      }
    }
    if (RF_frSize > 0) {
      if (RF_statusIndex > 0) {
        for (i = 1 ; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_statusIndex][i])) {
            if (RF_fresponseIn[RF_statusIndex][i] < 0) {
              result = FALSE;
              Rprintf("\nRF-SRC:  PRED status elements must be greater than or equal to zero or NA:  [%10d] = %12.4f \n", i, RF_fresponseIn[RF_statusIndex][i]);
            }
          }
        }
      }
    }
    if (result == FALSE) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Missingness verification failed.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }  
  RF_response = (double ***) vvector(1, RF_forestSize);
  for (i = 1 ; i <= RF_forestSize; i++) {
    RF_response[i] = RF_responseIn;
  }
  if (RF_timeIndex > 0) {
    RF_time = (double **) vvector(1, RF_forestSize);
    RF_masterTimeIndex = (uint **) vvector(1, RF_forestSize);
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_time[i] = RF_responseIn[RF_timeIndex];
      RF_masterTimeIndex[i] = RF_masterTimeIndexIn;
    }
    updateTimeIndexArray(0, NULL, NULL, RF_observationSize, RF_responseIn[RF_timeIndex], TRUE, FALSE, RF_masterTimeIndexIn);
  }
  if (RF_statusIndex > 0) {
    RF_status = (double **) vvector(1, RF_forestSize);
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_status[i] = RF_responseIn[RF_statusIndex];
    }
  }
  RF_observation = (double ***) vvector(1, RF_forestSize);
  for (i = 1 ; i <= RF_forestSize; i++) {
    RF_observation[i] = RF_observationIn;
  }
  RF_mRecordMap = uivector(1, RF_observationSize);
  RF_mRecordSize = getRecordMap(RF_mRecordMap, 
                              RF_observationSize, 
                              RF_responseIn, 
                              RF_observationIn);
  if (RF_mRecordSize == 0) {
    RF_mStatusFlag = RF_mTimeFlag = RF_mResponseFlag = RF_mPredictorFlag = FALSE;
    if (mode == RF_GROW) {
      RF_imputeSize = 1;
    }
  }
  else {
    stackMissingSignatures(RF_observationSize,
                           RF_rSize,
                           RF_responseIn,
                           RF_observationIn,
                           RF_mRecordMap,
                           RF_mRecordSize,
                           & RF_mRecordIndex,
                           & RF_mvSize,
                           & RF_mvSign,
                           & RF_mvIndex,
                           & RF_mrFactorSize,
                           & RF_mrFactorIndex,
                           & RF_mxFactorSize,
                           & RF_mxFactorIndex,
                           & RF_mTimeFlag,
                           & RF_mStatusFlag,
                           & RF_mResponseFlag,
                           & RF_mPredictorFlag);
    if (RF_mResponseFlag == TRUE) {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_response[i] = NULL;
        if (RF_timeIndex > 0) {
          RF_time[i] = NULL;
          RF_masterTimeIndex[i] = NULL;
        }
        if (RF_statusIndex > 0) {
          RF_status[i] = NULL;
        }
      }
    }
    if (RF_mPredictorFlag == TRUE) {
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_observation[i] = NULL;
      }
    }
  }  
  if (mode == RF_PRED) {
    RF_fobservation = (double ***) vvector(1, RF_forestSize);
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_fobservation[i] = RF_fobservationIn;
    }
    RF_fmRecordMap = uivector(1, RF_fobservationSize);
    if (RF_frSize > 0) {
      RF_fresponse = (double ***) vvector(1, RF_forestSize);
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_fresponse[i] = RF_fresponseIn;
      }
      if (RF_timeIndex > 0) {
        RF_ftime = (double **) vvector(1, RF_forestSize);
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_ftime[i] = RF_fresponseIn[RF_timeIndex];
        }
      }
      if (RF_statusIndex > 0) {
        RF_fstatus = (double **) vvector(1, RF_forestSize);
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_fstatus[i] = RF_fresponseIn[RF_statusIndex];
        }
      }
    }
    else {
      RF_fresponseIn = NULL;
    }
    RF_fmRecordSize = getRecordMap(RF_fmRecordMap, 
                                 RF_fobservationSize, 
                                 RF_fresponseIn, 
                                 RF_fobservationIn);
    if (RF_fmRecordSize == 0) {
      RF_fmStatusFlag = RF_fmTimeFlag = RF_fmResponseFlag = RF_fmPredictorFlag = FALSE;
    }  
    else {
      stackMissingSignatures(RF_fobservationSize,
                             RF_frSize,
                             RF_fresponseIn,
                             RF_fobservationIn,
                             RF_fmRecordMap,
                             RF_fmRecordSize,
                             & RF_fmRecordIndex,
                             & RF_fmvSize,
                             & RF_fmvSign,
                             & RF_fmvIndex,
                             & RF_fmrFactorSize,
                             & RF_fmrFactorIndex,
                             & RF_fmxFactorSize,
                             & RF_fmxFactorIndex,
                             & RF_fmTimeFlag,
                             & RF_fmStatusFlag,
                             & RF_fmResponseFlag,
                             & RF_fmPredictorFlag);
      if (RF_frSize > 0) {
        if (RF_fmResponseFlag == TRUE) {
          for (i = 1 ; i <= RF_forestSize; i++) {
            RF_fresponse[i] = NULL;
            if (RF_timeIndex > 0) {
              RF_ftime[i] = NULL;
            }
            if (RF_statusIndex > 0) {
              RF_fstatus[i] = NULL;
            }
          }
        }
      }
      if (RF_fmPredictorFlag == TRUE) {
        for (i = 1 ; i <= RF_forestSize; i++) {
          RF_fobservation[i] = NULL;
        }
      }
    }  
  }  
  dualUseFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      recordSize = RF_fmRecordSize;
      vSize = RF_fmvSize;
      dualUseFlag = TRUE;
      mFlag = ACTIVE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      recordSize = RF_mRecordSize;
      vSize = RF_mvSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
      RF_opt = RF_opt & (~OPT_OMIS);
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    RF_dmRecordBootFlag = cmatrix(1, RF_forestSize, 1, recordSize);
    for (j = 1; j <= RF_forestSize; j++) {
      for (i = 1; i <= recordSize; i++) {
        RF_dmRecordBootFlag[j][i] = mFlag;
      }
    }
    RF_dmvImputation = dmatrix3(1, RF_forestSize, 1, recordSize, 1, vSize);
    RF_mTerminalInfo = (Terminal ***) vvector(1, RF_forestSize);
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= recordSize; j++) {
        for (k = 1; k <= vSize; k++) {
          RF_dmvImputation[i][j][k] = NA_REAL;
        }
      }
    }
  }
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    initializeFactorArrays(mode);
  }
  return result;
}
void unstackMissingArrays(char mode) {
  char dualUseFlag;
  uint recordSize;
  uint vSize;
  free_vvector(RF_response, 1, RF_forestSize);
  if (RF_timeIndex > 0) {
    free_vvector(RF_time, 1, RF_forestSize);
    free_vvector(RF_masterTimeIndex, 1, RF_forestSize);
  }
  if (RF_statusIndex > 0) {
    free_vvector(RF_status, 1, RF_forestSize);
  }
  free_vvector(RF_observation, 1, RF_forestSize);
  free_uivector(RF_mRecordMap, 1, RF_observationSize);
  if (RF_mRecordSize == 0) {
  }
  else {
    unstackMissingSignatures(RF_rSize,
                             RF_mRecordSize,
                             RF_mRecordIndex,
                             RF_mvSize,
                             RF_mvSign,
                             RF_mvIndex,
                             RF_mrFactorSize,
                             RF_mrFactorIndex,
                             RF_mxFactorSize,
                             RF_mxFactorIndex);
  }
  if (mode == RF_PRED) {
    free_vvector(RF_fobservation, 1, RF_forestSize);
    free_uivector(RF_fmRecordMap, 1, RF_fobservationSize);
    if (RF_frSize > 0) {
      free_vvector(RF_fresponse, 1, RF_forestSize);
      if (RF_timeIndex > 0) {
        free_vvector(RF_ftime, 1, RF_forestSize);
      }
      if (RF_statusIndex > 0) {
        free_vvector(RF_fstatus, 1, RF_forestSize);
      }
    }
    if (RF_fmRecordSize == 0) {
    }
    else {
      unstackMissingSignatures(RF_frSize,
                               RF_fmRecordSize,
                               RF_fmRecordIndex,
                               RF_fmvSize,
                               RF_fmvSign,
                               RF_fmvIndex,
                               RF_fmrFactorSize,
                               RF_fmrFactorIndex,
                               RF_fmxFactorSize,
                               RF_fmxFactorIndex);
    }
  }
  dualUseFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    if (RF_fmRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = RF_fmRecordSize;
      vSize = RF_fmvSize;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = RF_mRecordSize;
      vSize = RF_mvSize;
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    free_cmatrix(RF_dmRecordBootFlag, 1, RF_forestSize, 1, recordSize);
    free_dmatrix3(RF_dmvImputation, 1, RF_forestSize, 1, recordSize, 1, vSize);
    free_vvector(RF_mTerminalInfo, 1, RF_forestSize);
  }
}
void stackMissingSignatures(uint     obsSize, 
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize, 
                            uint   **p_recordIndex, 
                            uint    *p_vSize,
                            int   ***p_vSign, 
                            int    **p_vIndex,
                            uint    *pRF_mrFactorSize,
                            uint   **pRF_mrFactorIndex,
                            uint    *pRF_mxFactorSize,
                            uint   **pRF_mxFactorIndex,
                            char    *pRF_mTimeFlag,
                            char    *pRF_mStatusFlag,
                            char    *pRF_mResponseFlag,
                            char    *pRF_mPredictorFlag) {
  uint i, j, p;
  if (recordSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to allocate for missingness in its absence.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  *p_recordIndex = uivector(1, recordSize);
  i = 0;
  for (j = 1; j <= obsSize; j++) {
    if (recordMap[j] > 0) {
      i++;
      (*p_recordIndex)[i] = j;
    }
  }
  *p_vSign = imatrix(1, rspSize + RF_xSize, 1, recordSize);
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize + RF_xSize; i++) {
      (*p_vSign)[i][j] = 0;
    }
  }
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize; i++) {
      if (ISNA(responsePtr[i][(*p_recordIndex)[j]])) {
        (*p_vSign)[i][j] = 1;
      }
    }
    for (i = 1; i <= RF_xSize; i++) {
      if (ISNA(predictorPtr[i][(*p_recordIndex)[j]])) {
        (*p_vSign)[rspSize + i][j] = 1;
      }
    }
  }
  *pRF_mStatusFlag = *pRF_mTimeFlag = *pRF_mResponseFlag = *pRF_mPredictorFlag = FALSE;
  *p_vIndex = ivector(1, rspSize + RF_xSize);
  *p_vSize = 0;
  for (i = 1; i <= rspSize; i++) {
    (*p_vIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_vSign)[i][j] == 1) {
        (*p_vSize) ++;
        (*p_vIndex)[*p_vSize] = - i;
        *pRF_mResponseFlag = TRUE;
        if (i == RF_timeIndex) {
          *pRF_mTimeFlag = TRUE;
        }
        else if (i == RF_statusIndex) {
          *pRF_mStatusFlag = TRUE;
        }
        j = recordSize;
      }
    }
  }  
  for (i = rspSize + 1; i <= rspSize + RF_xSize; i++) {
    (*p_vIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_vSign)[i][j] == 1) {
        (*p_vSize) ++;
        (*p_vIndex)[*p_vSize] =  i - rspSize;
        *pRF_mPredictorFlag = TRUE;
        j = recordSize;
      }
    }
  }  
  if (rspSize > 0) {
    *pRF_mrFactorIndex = uivector(1, rspSize);
    for (p = 1; p <= rspSize; p++) {
      (*pRF_mrFactorIndex)[p] = 0;
    }
  }
  *pRF_mxFactorIndex = uivector(1, RF_xSize);
  for (p = 1; p <= RF_xSize; p++) {
    (*pRF_mxFactorIndex)[p] = 0;
  }
  *pRF_mrFactorSize = *pRF_mxFactorSize = 0;
  for (p = 1; p <= *p_vSize; p++) {
    if ((*p_vIndex)[p] < 0) {
      if (strcmp(RF_rType[(uint) abs((*p_vIndex)[p])], "C") == 0) {
        (*pRF_mrFactorSize) ++;
        (*pRF_mrFactorIndex)[*pRF_mrFactorSize] = (uint) abs((*p_vIndex)[p]);
      }
    }
    else {
      if (strcmp(RF_xType[(*p_vIndex)[p]], "C") == 0) {
        (*pRF_mxFactorSize) ++;
        (*pRF_mxFactorIndex)[*pRF_mxFactorSize] = (*p_vIndex)[p];
      }
    }
  }
}
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize, 
                              uint     *recordIndex, 
                              uint      vSize,
                              int     **vSign, 
                              int      *vIndex,
                              uint      mrFactorSize,
                              uint     *mrFactorIndex,
                              uint      mxFactorSize,
                              uint     *mxFactorIndex) {
  if (recordSize == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to deallocate for missingness in its absence.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  free_uivector(recordIndex, 1, recordSize);
  free_imatrix(vSign, 1, rspSize + RF_xSize, 1, recordSize);
  free_ivector(vIndex, 1, rspSize + RF_xSize);
  if (rspSize > 0) {
    free_uivector(mrFactorIndex, 1, rspSize);
  }
  free_uivector(mxFactorIndex, 1, RF_xSize);
}
void initializeFactorArrays(char mode) {
  uint i, j;
  uint factorLevel;
  if (!(RF_rFactorCount + RF_xFactorCount > 0)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to initialize factorness in its absence.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  RF_rMaxFactorLevel = 0;
  for (j = 1; j <= RF_rFactorCount; j++) {
    RF_rFactorSize[j] = RF_rLevels[RF_rFactorIndex[j]];
    for (i = 1; i <= RF_observationSize; i++) {
      if (!ISNA(RF_responseIn[RF_rFactorIndex[j]][i])) {
        if (RF_responseIn[RF_rFactorIndex[j]][i] >= 1) {
          factorLevel = (uint) RF_responseIn[RF_rFactorIndex[j]][i];
          if (RF_rLevels[RF_rFactorIndex[j]] < factorLevel) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Factor level in data greater inconsistent with number of levels indicated:  %10d %10.4f", factorLevel, RF_rLevels[RF_rFactorIndex[j]]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        else {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_responseIn[RF_rFactorIndex[j]][i]);
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if (RF_rMaxFactorLevel < RF_rFactorSize[j]) {
      RF_rMaxFactorLevel = RF_rFactorSize[j];
    }
  }
  RF_xMaxFactorLevel = 0;
  for (j = 1; j <= RF_xFactorCount; j++) {
    RF_xFactorSize[j] = RF_xLevels[RF_xFactorIndex[j]];
    for (i = 1; i <= RF_observationSize; i++) {
      if (!ISNA(RF_observationIn[RF_xFactorIndex[j]][i])) {
        if (RF_observationIn[RF_xFactorIndex[j]][i] >= 1) {
          factorLevel = (uint) RF_observationIn[RF_xFactorIndex[j]][i];
          if (RF_xLevels[RF_xFactorIndex[j]] < factorLevel) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Factor level in data greater inconsistent with number of levels indicated:  %10d %10.4f", factorLevel, RF_xLevels[RF_xFactorIndex[j]]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
        else {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_observationIn[RF_xFactorIndex[j]][i]);
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    if (RF_xMaxFactorLevel < RF_xFactorSize[j]) {
      RF_xMaxFactorLevel = RF_xFactorSize[j];
    }
  }
  RF_maxFactorLevel = (RF_xMaxFactorLevel > RF_rMaxFactorLevel) ? RF_xMaxFactorLevel : RF_rMaxFactorLevel;
  if (mode == RF_PRED) {
    if (RF_frSize > 0) {
      for (j = 1; j <= RF_rFactorCount; j++) {
        factorLevel = 0;
        for (i = 1; i <= RF_fobservationSize; i++) {
          if (!ISNA(RF_fresponseIn[RF_rFactorIndex[j]][i])) {
            if (RF_fresponseIn[RF_rFactorIndex[j]][i] >= 1) {
              factorLevel = (factorLevel > (uint) RF_fresponseIn[RF_rFactorIndex[j]][i]) ? factorLevel : ((uint) RF_fresponseIn[RF_rFactorIndex[j]][i]);
            }
            else {
              Rprintf("\nRF-SRC:  *** ERROR *** ");
              Rprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_fobservationIn[RF_rFactorIndex[j]][i]);
              Rprintf("\nRF-SRC:  Please Contact Technical Support.");
              error("\nRF-SRC:  The application will now exit.\n");
            }
          }
        }
        if (factorLevel > RF_rFactorSize[j]) {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, RF_rFactorSize[j]);
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
    }
    for (j = 1; j <= RF_xFactorCount; j++) {
      factorLevel = 0;
      for (i = 1; i <= RF_fobservationSize; i++) {
        if (!ISNA(RF_fobservationIn[RF_xFactorIndex[j]][i])) {
          if (RF_fobservationIn[RF_xFactorIndex[j]][i] >= 1) {
            factorLevel = (factorLevel > (uint) RF_fobservationIn[RF_xFactorIndex[j]][i]) ? factorLevel : ((uint) RF_fobservationIn[RF_xFactorIndex[j]][i]);
          }
          else {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Factor level less than one (1):  %10.4f", RF_fobservationIn[RF_xFactorIndex[j]][i]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
      if (factorLevel > RF_xFactorSize[j]) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  !GROW factor level greater than maximum GROW factor level:  %10d vs. %10d", factorLevel, RF_xFactorSize[j]);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }
  RF_factorList = (Factor ***) vvector(1, RF_forestSize);
  for (j = 1; j <= RF_forestSize; j++) {
    RF_factorList[j] = NULL;
  }
}
char stackCompetingArrays(char mode) {
  uint obsSize;
  double  *statusPtr;
  uint    *mRecordMap;
  int    **mvSign;
  char eventAnalysisFlag, eventSubsetFlag, consistencyFlag, overrideFlag;
  char statusFlag;
  uint *eventCounter;
  uint i, j, jgrow;
  if (RF_statusIndex == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to stack competing risk structures in the absence of SURV data.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  RF_eventType = uivector(1, RF_observationSize);
  overrideFlag = FALSE;
  if (mode == RF_GROW) {
    if ((RF_splitRule == SURV_LGRNK) || 
        (RF_splitRule == SURV_LRSCR)) {
      overrideFlag = TRUE;
    }
  }
  else {
    if (!(RF_opt & OPT_COMP_RISK)) {
      overrideFlag = TRUE;
    }
  }
  getEventTypeSize(RF_observationSize, 
                   RF_responseIn[RF_statusIndex], 
                   RF_mRecordMap, 
                   RF_mvSign,  
                   overrideFlag, 
                   & RF_eventTypeSize,
                   & RF_mStatusSize,
                   RF_eventType);
  if (mode == RF_GROW) {
    if (RF_splitRule == SURV_CR_LAU) {
      if (RF_eventTypeSize == 1) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Split rule specified is for Competing Risk scenarios only.");
        error("\nRF-SRC:  The data set does not contain multiple events.");
      }
      i = 0;
      for (j = 1; j <= RF_eventTypeSize; j++) {
        if(fabs(RF_crWeight[j]) <= EPSILON) {
          i ++;
        }
        else {
          if(RF_crWeight[j] < 0.0) {
            Rprintf("\nRF-SRC:  *** ERROR *** ");
            Rprintf("\nRF-SRC:  Parameter verification failed.");
            Rprintf("\nRF-SRC:  Competing risk weight elements must be greater than or equal to zero:  %12.4f \n", RF_crWeight[j]);
            Rprintf("\nRF-SRC:  Please Contact Technical Support.");
            error("\nRF-SRC:  The application will now exit.\n");
          }
        }
      }
      if (i == RF_eventTypeSize) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Parameter verification failed.");
        Rprintf("\nRF-SRC:  Competing risk weight elements are all zero. \n");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }
  else {
    if (RF_opt & OPT_COMP_RISK) {
      if (RF_eventTypeSize == 1) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  CR analysis has been specified in !GROW mode.");
        error("\nRF-SRC:  However, the GROW data set does not contain multiple events.");
      }
    }
  }
  if (RF_eventTypeSize > 1) {
    RF_eventTypeIndex  = uivector(1, RF_eventType[RF_eventTypeSize]);
    for (j = 1; j <= RF_eventType[RF_eventTypeSize]; j++) {
      RF_eventTypeIndex[j] = 0;
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eventTypeIndex[RF_eventType[j]] = j;
    }
  }
  switch (mode) { 
  case RF_PRED: 
    if ((RF_opt & OPT_PERF) || (RF_opt & OPT_VIMP)) { 
      eventAnalysisFlag = TRUE; 
    } 
    else { 
      eventAnalysisFlag = FALSE; 
    } 
    break; 
  default: 
    eventAnalysisFlag = FALSE; 
    break; 
  } 
  if (eventAnalysisFlag == TRUE) {
    uint *feventType = uivector(1, RF_fobservationSize); 
    uint feventTypeSize;
    consistencyFlag = TRUE; 
    overrideFlag = FALSE;
    if (RF_eventTypeSize == 1) { 
      overrideFlag = TRUE; 
    }
    getEventTypeSize(RF_fobservationSize, 
                     RF_fresponseIn[RF_statusIndex], 
                     RF_fmRecordMap,
                     RF_fmvSign, 
                     overrideFlag, 
                     & feventTypeSize, 
                     & RF_mStatusSize, 
                     feventType);
    if (RF_eventTypeSize > 1) { 
      for (j = 1; j <= feventTypeSize; j++) {
        for (jgrow = 1; jgrow <= RF_eventTypeSize; jgrow++) { 
          if (feventType[j] != RF_eventType[jgrow]) { 
            if (jgrow == RF_eventTypeSize) { 
              consistencyFlag = FALSE; 
            } 
          } 
          else { 
            jgrow = RF_eventTypeSize; 
          } 
        } 
      } 
    }
    free_uivector(feventType, 1, RF_fobservationSize); 
    if (consistencyFlag == FALSE) { 
      Rprintf("\nRF-SRC: *** ERROR *** ");
      Rprintf("\nRF-SRC: Unknown event type encountered in !GROW mode. "); 
      error("\nRF-SRC: Please Contact Technical Support.");
    }
  }  
  if (RF_eventTypeSize > 1) { 
    if (mode == RF_PRED) { 
      if ((RF_opt & OPT_PERF) || (RF_opt & OPT_VIMP)) { 
        eventSubsetFlag = TRUE; 
      }
      else { 
        eventSubsetFlag = FALSE; 
      } 
    } 
    else { 
      eventSubsetFlag = TRUE; 
    } 
  } 
  else { 
    eventSubsetFlag = FALSE; 
  }
  if (eventSubsetFlag == TRUE) {
    if ((mode == RF_GROW) || (mode == RF_REST)) { 
      obsSize = RF_observationSize; 
      statusPtr = RF_responseIn[RF_statusIndex]; 
      mvSign = RF_mvSign;
      mRecordMap = RF_mRecordMap; 
    } 
    else { 
      obsSize = RF_fobservationSize;
      statusPtr = RF_fresponseIn[RF_statusIndex]; 
      mvSign = RF_fmvSign; 
      mRecordMap = RF_fmRecordMap; 
    }
    RF_eIndividualSize = uivector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) { 
      RF_eIndividualSize[j] = 0; 
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) { 
        statusFlag = TRUE; 
      } 
      else { 
        if (mvSign[RF_statusIndex][mRecordMap[i]] == 0) { 
          statusFlag = TRUE; 
        } 
      }
      if (statusFlag == TRUE) { 
        if ((uint) statusPtr[i] > 0) { 
          RF_eIndividualSize[RF_eventTypeIndex[(uint) statusPtr[i]]] ++; 
        }
        else { 
          for (j=1; j <= RF_eventTypeSize; j++) {
            RF_eIndividualSize[j] ++; 
          } 
        } 
      } 
    } 
    RF_eIndividualIn = (uint **) vvector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) { 
      RF_eIndividualIn[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    eventCounter = uivector(1, RF_eventTypeSize); 
    for (j = 1; j <= RF_eventTypeSize; j++) { 
      eventCounter[j] = 0; 
    }
    for (i = 1; i <= obsSize; i++) {
      statusFlag = FALSE;
      if (mRecordMap[i] == 0) { 
        statusFlag = TRUE; 
      } 
      else { 
        if (mvSign[RF_statusIndex][mRecordMap[i]] == 0) { 
          statusFlag = TRUE; 
        } 
      }
      if (statusFlag == TRUE) { 
        if ((uint) statusPtr[i] > 0) { 
          j = RF_eventTypeIndex[(uint) statusPtr[i]]; 
          eventCounter[j] ++;
          RF_eIndividualIn[j][eventCounter[j]] = i; 
        } 
        else { 
          for (j=1; j <= RF_eventTypeSize; j++) { 
            eventCounter[j] ++;
            RF_eIndividualIn[j][eventCounter[j]] = i; 
          } 
        } 
      } 
    }
    free_uivector(eventCounter, 1, RF_eventTypeSize); 
  }  
  return TRUE;
}
void getEventTypeSize(uint obsSize, 
                      double *status, 
                      uint *mRecordMap,
                      int **mvSign, 
                      char overrideFlag, 
                      uint *eventTypeSize, 
                      uint *msize,
                      uint *eventType) {
  uint statusFlag;
  uint leadingIndex;
  uint i;
  if (RF_statusIndex == 0) { 
    Rprintf("\nRF-SRC: *** ERROR *** ");
    Rprintf("\nRF-SRC: Attempt to stack competing risk structures in the absence of SURV data."); 
    Rprintf("\nRF-SRC: Please Contact Technical Support."); 
    error("\nRF-SRC: The application will now exit.\n");
  }
  *eventTypeSize = *msize = 0;
  for (i = 1; i <= obsSize; i++) { 
    eventType[i] = 0; 
    statusFlag = FALSE;
    if (mRecordMap[i] == 0) { 
      statusFlag = TRUE; 
    }
    else { 
      if (mvSign[RF_statusIndex][mRecordMap[i]] == 0) { 
        statusFlag = TRUE; 
      } 
    }
    if (statusFlag == TRUE) { 
      if ((uint) status[i] > 0) { 
        if (overrideFlag == TRUE) { 
        } 
        else { 
          (*eventTypeSize) ++;
          eventType[*eventTypeSize] = (uint) status[i]; 
        } 
      } 
      else { 
      } 
    } 
    else { 
      (*msize) ++; 
    }
  }  
  if (overrideFlag == TRUE) { 
    *eventTypeSize = 1; 
  } 
  else { 
    if(*eventTypeSize > 0) {
      hpsortui(eventType, *eventTypeSize);
      leadingIndex = 1;
      for (i=2; i <= *eventTypeSize; i++) { 
        if (eventType[i] > eventType[leadingIndex]) { 
          leadingIndex++;
          eventType[leadingIndex] = eventType[i]; 
        } 
      } 
      *eventTypeSize = leadingIndex;
    }
    for (i= *eventTypeSize + 1; i <= obsSize; i++) { 
      eventType[i] = 0; 
    }
  }
}
void unstackCompetingArrays(char mode) {
  char eventSubsetFlag;
  uint j;
  if (RF_statusIndex == 0) { 
    Rprintf("\nRF-SRC: *** ERROR *** ");
    Rprintf("\nRF-SRC: Attempt to unstack competing risk structures in the absence of SURV data."); 
    Rprintf("\nRF-SRC: Please Contact Technical Support."); 
    error("\nRF-SRC: The application will now exit.\n"); 
  }
  if (RF_eventTypeSize > 1) {  
    free_uivector(RF_eventTypeIndex, 1, RF_eventType[RF_eventTypeSize]);
  }
  free_uivector(RF_eventType, 1, RF_observationSize);
  if (RF_eventTypeSize > 1) { 
    if (mode == RF_PRED) { 
      if ((RF_opt & OPT_PERF) || (RF_opt & OPT_VIMP)) { 
        eventSubsetFlag = TRUE; 
      }
      else { 
        eventSubsetFlag = FALSE; 
      } 
    } 
    else { 
      eventSubsetFlag = TRUE; 
    } 
  } 
  else { 
    eventSubsetFlag = FALSE; 
  }
  if (eventSubsetFlag == TRUE) {
    for (j = 1; j <= RF_eventTypeSize; j++) { 
      free_uivector(RF_eIndividualIn[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    free_vvector(RF_eIndividualIn, 1, RF_eventTypeSize);
    free_uivector(RF_eIndividualSize, 1, RF_eventTypeSize);
  }  
}
char stackClassificationArrays(char mode) {
  char classAnalysisFlag, consistencyFlag;
  uint j, k, jgrow;
  if (RF_rFactorCount == 0) { 
    Rprintf("\nRF-SRC: *** ERROR *** ");
    Rprintf("\nRF-SRC: Attempt to stack classification structures in the absence of CLAS data."); 
    Rprintf("\nRF-SRC: Please Contact Technical Support."); 
    error("\nRF-SRC: The application will now exit.\n");
  }
  RF_classLevel = (uint **) vvector(1, RF_rFactorCount);
  RF_classLevelSize = uivector(1, RF_rFactorCount);
  getClassLevelSize(RF_observationSize, 
                    RF_responseIn, 
                    RF_mRecordMap,
                    RF_mvSign, 
                    RF_classLevelSize, 
                    RF_classLevel);
  RF_classLevelIndex = (uint **) vvector(1, RF_rFactorCount); 
  for (k = 1; k <= RF_rFactorCount; k++) {
    RF_classLevelIndex[k] = uivector(1, RF_classLevel[k][RF_classLevelSize[k]]); 
    for (j = 1; j <= RF_classLevel[k][RF_classLevelSize[k]]; j++) {
      RF_classLevelIndex[k][j] = 0; 
    } 
    for (j = 1; j <= RF_classLevelSize[k]; j++) {
      RF_classLevelIndex[k][RF_classLevel[k][j]] = j; 
    }
  }  
  switch (mode) { 
  case RF_PRED: 
    if ((RF_opt & OPT_PERF) || (RF_opt & OPT_VIMP)) {
      classAnalysisFlag = TRUE; 
    } 
    else { 
      classAnalysisFlag = FALSE; 
    }
    break; 
  default: 
    classAnalysisFlag = FALSE; 
    break; 
  } 
  if (classAnalysisFlag == TRUE) {
    uint **fclassLevel = (uint **) vvector(1, RF_rFactorCount);
    uint *fclassLevelSize = uivector(1, RF_rFactorCount);
    getClassLevelSize(RF_fobservationSize, 
                      RF_fresponseIn, 
                      RF_fmRecordMap, 
                      RF_fmvSign,  
                      fclassLevelSize, 
                      fclassLevel);
    consistencyFlag = TRUE;
    for (j = 1; j <= RF_rFactorCount; j++) {
      for (k = 1; k <= fclassLevelSize[j]; k++) {
        for (jgrow = 1; jgrow <= RF_classLevelSize[j]; jgrow++) {
          if (fclassLevel[j][k] != RF_classLevel[j][jgrow]) {
            if (jgrow == RF_classLevelSize[j]) {
              consistencyFlag = FALSE;
            }
          }
          else {
            jgrow = RF_classLevelSize[j];
          }
        }
      }
    }
    for (j = 1; j <= RF_rFactorCount; j ++) {
      free_uivector(fclassLevel[j], 1, fclassLevelSize[j]);
    }
    free_vvector(fclassLevel, 1, RF_rFactorCount);
    free_uivector(fclassLevelSize, 1, RF_rFactorCount);
    if (consistencyFlag == FALSE) {
    }
  }  
  return TRUE;
}
void getClassLevelSize(uint      obsSize, 
                       double  **response, 
                       uint     *mRecordMap, 
                       int     **mvSign,  
                       uint     *classLevelSize,
                       uint    **classLevel) {
  uint *rawClassVector;
  uint classFlag;
  uint leadingIndex;
  uint i, j, k;
  if (RF_rFactorCount == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to stack classification response structures in the absence of CLAS data.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  rawClassVector = uivector(1, obsSize);
  for (j = 1; j <= RF_rFactorCount; j++) {
    classLevelSize[j] = 0;
  }
  for (j = 1; j <= RF_rFactorCount; j++) {
    for (i = 1; i <= obsSize; i++) {
      classFlag = FALSE;
      if (mRecordMap[i] == 0) { 
        classFlag = TRUE;
      }
      else {
        if (mvSign[RF_rFactorIndex[j]][mRecordMap[i]] == 0) {
          classFlag = TRUE;
        }
      }
      if (classFlag == TRUE) {
        classLevelSize[j] ++;
        rawClassVector[classLevelSize[j]] = (uint) response[RF_rFactorIndex[j]][i];
      }
      else {
      }
    }  
    hpsortui(rawClassVector, classLevelSize[j]);
    leadingIndex = 1;
    for (k=2; k <= classLevelSize[j]; k++) {
      if (rawClassVector[k] > rawClassVector[leadingIndex]) {
        leadingIndex++;
        rawClassVector[leadingIndex] = rawClassVector[k];
      }
    }
    classLevelSize[j] = leadingIndex;
    classLevel[j] = uivector(1, classLevelSize[j]);
    for (k=1; k <= classLevelSize[j]; k++) {
      classLevel[j][k] = rawClassVector[k];
    }
  } 
  free_uivector(rawClassVector, 1, obsSize);
}
void unstackClassificationArrays(char mode) {
  uint j;
  if (RF_rFactorCount == 0) { 
    Rprintf("\nRF-SRC: *** ERROR *** ");
    Rprintf("\nRF-SRC: Attempt to unstack classification structures in the absence of CLAS data."); 
    Rprintf("\nRF-SRC: Please Contact Technical Support."); 
    error("\nRF-SRC: The application will now exit.\n");
  }
  for (j = 1; j <= RF_rFactorCount; j++) {
    free_uivector(RF_classLevelIndex[j], 1, RF_classLevel[j][RF_classLevelSize[j]]); 
  }
  free_vvector(RF_classLevelIndex, 1, RF_rFactorCount); 
  for (j = 1; j <= RF_rFactorCount; j ++) {
    free_uivector(RF_classLevel[j], 1, RF_classLevelSize[j]);
  }
  free_vvector(RF_classLevel, 1, RF_rFactorCount);
  free_uivector(RF_classLevelSize, 1, RF_rFactorCount);
}
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               double  **pRF_oobEnsemble,
                               double  **pRF_fullEnsemble,
                               double  **p_performance,
                               uint    **pRF_leafCount,
                               uint    **pRF_proximity,
                               double  **pRF_importance,
                               int     **pRF_seed,
                               double  **pRF_oobImputation,
                               double  **p_imputation,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               double ***pRF_sOOBImputeResponsePtr,
                               double ***pRF_sOOBImputePredictorPtr,
                               uint    **pRF_varUsed,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth,
                               double  **pRF_oobEnsembleCIF,
                               double  **pRF_fullEnsembleCIF,
                               double  **pRF_oobEnsembleSRV,
                               double  **pRF_fullEnsembleSRV,
                               double  **pRF_oobEnsembleMRT,
                               double  **pRF_fullEnsembleMRT,
                               uint    **pRF_terminalNodeMembership,
                               uint    **pRF_bootstrapMembership,
                               uint     *stackCount,
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint ensembleSize;
  uint performanceSize;
  uint proximitySize;
  uint imputationSize;
  uint importanceSize;
  uint xVimpSize;
  uint varUsedSize;
  uint splitDepthSize;
  uint  obsSize;
  uint  mRecordSize;
  uint *mRecordIndex;
  double **responsePtr;
  double **predictorPtr;
  uint     rspSize;
  uint     perfDim;
  uint     ensbDimOne;
  uint     ensbDimTwo;
  uint     dpthDimOne;
  uint i,j,k,p;
  sexpIndex      = 0;  
  ensembleSize   = 0;  
  performanceSize= 0;  
  proximitySize  = 0;  
  imputationSize = 0;  
  importanceSize = 0;  
  xVimpSize      = 0;  
  varUsedSize    = 0;  
  splitDepthSize = 0;  
  dpthDimOne     = 0;  
  obsSize        = 0;  
  mRecordSize    = 0;  
  responsePtr    = NULL;  
  predictorPtr   = NULL;  
  mRecordIndex   = NULL;  
  rspSize        = 0;     
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      perfDim = ensbDimOne = RF_eventTypeSize;
    }
    else {
      perfDim = ensbDimOne = 1;
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      perfDim = RF_rFactorSize[1] + 1;
      ensbDimOne = 1;
    }
    else {
      perfDim = 1;
      ensbDimOne = 1;
    }
  }
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
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    if (RF_opt & OPT_SPLDPTH_F) {
      dpthDimOne = 1;
    }
    else {
      dpthDimOne = RF_forestSize;
    }
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    mRecordSize = RF_fmRecordSize;
    rspSize = RF_frSize;
    responsePtr  = RF_fresponseIn;
    predictorPtr = RF_fobservationIn;
    mRecordIndex = RF_fmRecordIndex;
    *stackCount = 2;
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_eventTypeSize == 1) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if (RF_opt & OPT_PERF) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2; 
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_MISS) {
      imputationSize = (RF_xSize + rspSize + 1) * mRecordSize;
      (*stackCount) += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      splitDepthSize = dpthDimOne * RF_xSize * RF_observationSize;
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_VIMP) {
      (*stackCount) += 1;
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
        importanceSize = perfDim;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
        importanceSize = perfDim * RF_intrPredictorSize;
      }
    }
    if (RF_opt & OPT_MEMB) {
      (*stackCount) += 2;
    }
    break;
  default:
    obsSize = RF_observationSize;
    mRecordSize = RF_mRecordSize;
    rspSize = RF_rSize;
    responsePtr  = RF_responseIn;
    predictorPtr = RF_observationIn;
    mRecordIndex = RF_mRecordIndex;
    *stackCount = 0;
    if (RF_opt & OPT_LEAF) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_FENS) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_OENS) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_PERF) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_eventTypeSize == 1) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if (RF_opt & OPT_OENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_eventTypeSize == 1) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2; 
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_SEED) {
      if (RF_opt & OPT_TREE) {
        (*stackCount) += 7;
      }
      else {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  SEXP TREE output request inconsistent.");
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
    if ((RF_opt & OPT_MISS) | (RF_opt & OPT_OMIS)) {
      imputationSize = (RF_xSize + rspSize + 1) * mRecordSize;
      if (RF_opt & OPT_MISS) {
        (*stackCount) += 1;
      }
      if (RF_opt & OPT_OMIS) {
        (*stackCount) += 1;
      }
    }
    if (RF_opt & OPT_VUSE) {
      if (RF_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
        varUsedSize = RF_forestSize;
      }
      else {
        varUsedSize = 1;
      }
      (*stackCount) += 1;
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      splitDepthSize = dpthDimOne * RF_xSize * RF_observationSize;
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_VIMP) {
      (*stackCount) += 1;
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
        importanceSize = perfDim;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
        importanceSize = perfDim * RF_intrPredictorSize;
      }
    }
    if (RF_opt & OPT_MEMB) {
      (*stackCount) += 2;
    }
    break;
  }
  performanceSize = perfDim * RF_forestSize;
  ensembleSize = ensbDimOne * ensbDimTwo * obsSize;
  PROTECT(sexpVector[RF_OUTP_ID] = allocVector(VECSXP, *stackCount));
  PROTECT(sexpVector[RF_STRG_ID] = allocVector(STRSXP, *stackCount));
  setAttrib(sexpVector[RF_OUTP_ID], R_NamesSymbol, sexpVector[RF_STRG_ID]);
  sexpIndex = 0;
  if (RF_opt & OPT_FENS) {
    PROTECT(sexpVector[RF_FENS_ID] = NEW_NUMERIC(ensembleSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_FENS_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_FENS_ID]));
    *pRF_fullEnsemble = NUMERIC_POINTER(sexpVector[RF_FENS_ID]);
    sexpIndex ++;
    RF_fullEnsemblePtr = (double ***) vvector(1, ensbDimOne);
    RF_fullEnsembleDen = uivector(1, obsSize);
    for (j = 1; j <= ensbDimOne; j++) {
      RF_fullEnsemblePtr[j] = (double **) vvector(1, ensbDimTwo);
      for (k = 1; k <= ensbDimTwo; k++) {
        RF_fullEnsemblePtr[j][k]  = (*pRF_fullEnsemble) + ((j-1) * ensbDimTwo * obsSize) + ((k-1) * obsSize) - 1;
      }
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= ensbDimOne; j++) {
        for (k = 1; k <= ensbDimTwo; k++) {
          RF_fullEnsemblePtr[j][k][i] = 0.0;
        }
      }
      RF_fullEnsembleDen[i] = 0;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      PROTECT(sexpVector[RF_FMRT_ID] = NEW_NUMERIC(ensbDimOne * obsSize));
      SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_FMRT_ID]);
      SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_FMRT_ID]));
      *pRF_fullEnsembleMRT = NUMERIC_POINTER(sexpVector[RF_FMRT_ID]);
      sexpIndex ++;
      RF_fullMRTPtr = (double **) vvector(1, ensbDimOne);
      for (j = 1; j <= ensbDimOne; j++) {
        RF_fullMRTPtr[j] = (*pRF_fullEnsembleMRT) + ((j-1) * obsSize) - 1;
      }
      for (j = 1; j <= ensbDimOne; j++) {
        for (i = 1; i <= obsSize; i++) {
          RF_fullMRTPtr[j][i] = 0.0;
        }
      }
      if (RF_eventTypeSize == 1) {
        PROTECT(sexpVector[RF_FSRV_ID] = NEW_NUMERIC(ensbDimTwo * obsSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_FSRV_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_FSRV_ID]));
        *pRF_fullEnsembleSRV = NUMERIC_POINTER(sexpVector[RF_FSRV_ID]);
        sexpIndex ++;
        RF_fullSRVPtr = (double **) vvector(1, ensbDimTwo);
        for (j = 1; j <= ensbDimTwo; j++) {
          RF_fullSRVPtr[j]  = (*pRF_fullEnsembleSRV) + ((j-1) * obsSize) - 1;
        }
        for (j = 1; j <= ensbDimTwo; j++) {
          for (i = 1; i <= obsSize; i++) {
            RF_fullSRVPtr[j][i]  = 0.0;
          }
        }
      }  
      else {
        PROTECT(sexpVector[RF_FCIF_ID] = NEW_NUMERIC(ensembleSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_FCIF_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_FCIF_ID]));
        *pRF_fullEnsembleCIF = NUMERIC_POINTER(sexpVector[RF_FCIF_ID]);
        sexpIndex ++;
        RF_fullCIFPtr = (double ***) vvector(1, ensbDimOne);
        for (j = 1; j <= ensbDimOne; j++) {
          RF_fullCIFPtr[j] = (double **) vvector(1, ensbDimTwo);
          for (k = 1; k <= ensbDimTwo; k++) {
            RF_fullCIFPtr[j][k]  = (*pRF_fullEnsembleCIF) + ((j-1) * ensbDimTwo * obsSize) + ((k-1) * obsSize) - 1;
          }
        }
        for (j = 1; j <= ensbDimOne; j++) {
          for (k = 1; k <= ensbDimTwo; k++) {
            for (i = 1; i <= obsSize; i++) {
              RF_fullCIFPtr[j][k][i] = 0.0;
            }
          }
        }
      }
    }
  }
  if (RF_opt & OPT_OENS) {
    PROTECT(sexpVector[RF_OENS_ID] = NEW_NUMERIC(ensembleSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OENS_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OENS_ID]));
    *pRF_oobEnsemble = NUMERIC_POINTER(sexpVector[RF_OENS_ID]);
    sexpIndex ++;
    RF_oobEnsemblePtr  = (double ***) vvector(1, ensbDimOne);
    RF_oobEnsembleDen  = uivector(1, obsSize);
    for (j = 1; j <= ensbDimOne; j++) {
      RF_oobEnsemblePtr[j] = (double **) vvector(1, ensbDimTwo);
      for (k = 1; k <= ensbDimTwo; k++) {
        RF_oobEnsemblePtr[j][k]  = (*pRF_oobEnsemble) + ((j-1) * ensbDimTwo * obsSize) + ((k-1) * obsSize) - 1;
      }
    }
    for (i = 1; i <= obsSize; i++) {
      for (j = 1; j <= ensbDimOne; j++) {
        for (k = 1; k <= ensbDimTwo; k++) {
          RF_oobEnsemblePtr[j][k][i] = 0.0;
        }
      }
      RF_oobEnsembleDen[i] = 0;
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      PROTECT(sexpVector[RF_OMRT_ID] = NEW_NUMERIC(ensbDimOne * obsSize));
      SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OMRT_ID]);
      SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OMRT_ID]));
      *pRF_oobEnsembleMRT = NUMERIC_POINTER(sexpVector[RF_OMRT_ID]);
      sexpIndex ++;
      RF_oobMRTPtr = (double **) vvector(1, ensbDimOne);
      for (j = 1; j <= ensbDimOne; j++) {
        RF_oobMRTPtr[j] = (*pRF_oobEnsembleMRT) + ((j-1) * obsSize) - 1;
      }
      for (j = 1; j <= ensbDimOne; j++) {
        for (i = 1; i <= obsSize; i++) {
          RF_oobMRTPtr[j][i] = 0.0;
        }
      }
      if (RF_eventTypeSize == 1) {
        PROTECT(sexpVector[RF_OSRV_ID] = NEW_NUMERIC(ensbDimTwo * obsSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OSRV_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OSRV_ID]));
        *pRF_oobEnsembleSRV = NUMERIC_POINTER(sexpVector[RF_OSRV_ID]);
        sexpIndex ++;
        RF_oobSRVPtr = (double **) vvector(1, ensbDimTwo);
        for (j = 1; j <= ensbDimTwo; j++) {
          RF_oobSRVPtr[j]  = (*pRF_oobEnsembleSRV) + ((j-1) * obsSize) - 1;
        }
        for (j = 1; j <= ensbDimTwo; j++) {
          for (i = 1; i <= obsSize; i++) {
            RF_oobSRVPtr[j][i]  = 0.0;
          }
        }
      } 
      else {
        PROTECT(sexpVector[RF_OCIF_ID] = NEW_NUMERIC(ensembleSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OCIF_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OCIF_ID]));
        *pRF_oobEnsembleCIF = NUMERIC_POINTER(sexpVector[RF_OCIF_ID]);
        sexpIndex ++;
        RF_oobCIFPtr = (double ***) vvector(1, ensbDimOne);
        for (j = 1; j <= ensbDimOne; j++) {
          RF_oobCIFPtr[j] = (double **) vvector(1, ensbDimTwo);
          for (k = 1; k <= ensbDimTwo; k++) {
            RF_oobCIFPtr[j][k]  = (*pRF_oobEnsembleCIF) + ((j-1) * ensbDimTwo * obsSize) + ((k-1) * obsSize) - 1;
          }
        }
        for (j = 1; j <= ensbDimOne; j++) {
          for (k = 1; k <= ensbDimTwo; k++) {
            for (i = 1; i <= obsSize; i++) {
              RF_oobCIFPtr[j][k][i] = 0.0;
            }
          }
        }
      }     
    }
  }
  if (RF_opt & OPT_PERF) {
    PROTECT(sexpVector[RF_PERF_ID] = NEW_NUMERIC(performanceSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PERF_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_PERF_ID]));
    *p_performance = NUMERIC_POINTER(sexpVector[RF_PERF_ID]);
    sexpIndex ++;
    RF_performancePtr = (double **) vvector(1, RF_forestSize);
    for (i = 1; i <= RF_forestSize; i++) {
      RF_performancePtr[i]  = (*p_performance)  + ((i-1) * perfDim) - 1;
    }
    for (j = 1; j <=  RF_forestSize; j++) {
      for (k = 1; k <= perfDim; k++) {
        RF_performancePtr[j][k] = NA_REAL;
      }
    }
  }
  if (RF_opt & OPT_PROX) {
    PROTECT(sexpVector[RF_PROX_ID] = NEW_INTEGER(proximitySize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_PROX_ID]));
    *pRF_proximity = (uint*) INTEGER_POINTER(sexpVector[RF_PROX_ID]);
    sexpIndex ++;
    (*pRF_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*pRF_proximity)[i] = 0;
    }
  }
  if (RF_opt & OPT_LEAF) {
    PROTECT(sexpVector[RF_LEAF_ID] = NEW_INTEGER(RF_forestSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_LEAF_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_LEAF_ID]));
    *pRF_leafCount    = (uint*) INTEGER_POINTER(sexpVector[RF_LEAF_ID]);
    sexpIndex ++;
    (*pRF_leafCount) --;
    for (i = 1; i <= RF_forestSize; i++) {
      (*pRF_leafCount)[i] = 0;
    }
    RF_leafCount = *pRF_leafCount;
  }
  if (RF_opt & OPT_SEED) {
    PROTECT(sexpVector[RF_SEED_ID] = NEW_INTEGER(RF_forestSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_SEED_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_SEED_ID]));
    *pRF_seed = (int*) INTEGER_POINTER(sexpVector[RF_SEED_ID]);
    sexpIndex ++;
    (*pRF_seed) --;
    for (i = 1; i <= RF_forestSize; i++) {
      (*pRF_seed)[i] = -1;
    }
  }
  if (RF_opt & OPT_MISS) {
    PROTECT(sexpVector[RF_MISS_ID] = NEW_NUMERIC(imputationSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_MISS_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_MISS_ID]));
    *p_imputation = NUMERIC_POINTER(sexpVector[RF_MISS_ID]);
    sexpIndex ++;
    if (rspSize > 0) {
      *pRF_sImputeResponsePtr = (double **) vvector(1, rspSize);
      for (i = 1; i <= rspSize; i++) {
        (*pRF_sImputeResponsePtr)[i]  = (*p_imputation)  + (i * mRecordSize) - 1;
      }
    }
    *pRF_sImputePredictorPtr = (double **) vvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      (*pRF_sImputePredictorPtr)[i]  = (*p_imputation)  + ((rspSize + i) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*p_imputation)[i-1] = (double) mRecordIndex[i];
      if (rspSize > 0) {
        for (j = 1; j <= rspSize; j++) {
          (*pRF_sImputeResponsePtr)[j][i] = responsePtr[j][mRecordIndex[i]];
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_sImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
  }
  if (RF_opt & OPT_OMIS) {
    PROTECT(sexpVector[RF_OMIS_ID] = NEW_NUMERIC(imputationSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OMIS_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OMIS_ID]));
    *pRF_oobImputation = NUMERIC_POINTER(sexpVector[RF_OMIS_ID]);
    sexpIndex ++;
    if (rspSize > 0) {
      *pRF_sOOBImputeResponsePtr = (double **) vvector(1, rspSize);
      for (i = 1; i <= rspSize; i++) {
        (*pRF_sOOBImputeResponsePtr)[i] = (*pRF_oobImputation)  + (i * mRecordSize) - 1;
      }
    }
    *pRF_sOOBImputePredictorPtr = (double **) vvector(1, RF_xSize);
    for (i = 1; i <= RF_xSize; i++) {
      (*pRF_sOOBImputePredictorPtr)[i] = (*pRF_oobImputation)  + ((rspSize + i) * mRecordSize) - 1;
    }
    for (i = 1; i <= mRecordSize; i++) {
      (*pRF_oobImputation)[i-1] = (double) mRecordIndex[i];
      if (rspSize > 0) {
        for (j = 1; j <= rspSize; j++) {
          (*pRF_sOOBImputeResponsePtr)[j][i] = responsePtr[j][mRecordIndex[i]];
        }
      }
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_sOOBImputePredictorPtr)[j][i] = predictorPtr[j][mRecordIndex[i]];
      }
    }
  }
  if (RF_opt & OPT_VIMP) {
    PROTECT(sexpVector[RF_VIMP_ID] = NEW_NUMERIC(importanceSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_VIMP_ID]));
    *pRF_importance = NUMERIC_POINTER(sexpVector[RF_VIMP_ID]);
    sexpIndex ++;
    RF_importancePtr = (double **) vvector(1, xVimpSize);
    for (k = 1; k <= xVimpSize; k++) {
      RF_importancePtr[k]  = (*pRF_importance)  + ((k-1) * perfDim) - 1;
    }
    for (k = 1; k <= xVimpSize; k++) {
      for (j = 1; j <= perfDim; j++) {
        RF_importancePtr[k][j] = NA_REAL;
      }
    }
    RF_vimpEnsembleDen  = (uint **) vvector(1, xVimpSize);
    for (j = 1; j <= xVimpSize; j++) {
      RF_vimpEnsembleDen[j] = uivector(1, obsSize);
      for (i = 1; i <= obsSize; i++) {
        RF_vimpEnsembleDen[j][i] = 0;
      }
    }
    if(RF_opt & OPT_VIMP_LEOB) { 
      RF_vimpLeo = dmatrix3(1, RF_forestSize, 1, xVimpSize, 1, perfDim);
      RF_perfLeo = dmatrix(1, RF_forestSize, 1,  perfDim);
      for (i = 1; i <= RF_forestSize; i++) {
        for (k = 1; k <= perfDim; k++) {
          RF_perfLeo[i][k] = NA_REAL;
          for (j = 1; j <= xVimpSize; j++) {
            RF_vimpLeo[i][j][k] = NA_REAL;
          }
        }
      }
    }
    RF_vimpOutcome = dmatrix(1, xVimpSize, 1, obsSize);
    for (k=1; k <= xVimpSize; k++) {
      for (i = 1; i <= obsSize; i++) {
        RF_vimpOutcome[k][i] = 0.0;
      }
    }
    RF_cVimpEnsemble  = NULL;
    RF_sVimpEnsemble = NULL;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_sVimpEnsemble = dmatrix4(1, xVimpSize, 1, ensbDimOne, 1, ensbDimTwo, 1, obsSize);
      for (p = 1; p <= xVimpSize; p++) {
        for (j = 1; j <= ensbDimOne; j++) {
          for (i = 1; i <= obsSize; i++) {
            for (k = 1; k <= ensbDimTwo; k++) {
              RF_sVimpEnsemble[p][j][k][i] = 0.0;
            }
          }
        }
      }
    }
    else {
      if (RF_rFactorCount > 0) {
        RF_cVimpEnsemble = dmatrix3(1, xVimpSize, 1, ensbDimTwo, 1, obsSize);
        for (p = 1; p <= xVimpSize; p++) {
          for (i = 1; i <= obsSize; i++) {
            for (k = 1; k <= ensbDimTwo; k++) {
              RF_cVimpEnsemble[p][k][i] = 0.0;
            }
          }
        }
      }
    }
  }  
  if (RF_opt & OPT_VUSE) {
    PROTECT(sexpVector[RF_VUSE_ID] = NEW_INTEGER(varUsedSize * RF_xSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_VUSE_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_VUSE_ID]));
    *pRF_varUsed = (uint*) INTEGER_POINTER(sexpVector[RF_VUSE_ID]);
    sexpIndex ++;
    if (RF_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      *pRF_varUsedPtr = (uint **) vvector(1, RF_forestSize);
      for (i = 1; i <= RF_forestSize; i++) {
        (*pRF_varUsedPtr)[i] = (*pRF_varUsed) + ((i-1)*(RF_xSize)) - 1;
      }
    }
    else {
      *pRF_varUsedPtr = uimatrix(1, RF_forestSize, 1, RF_xSize);
    }
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= RF_xSize; j++) {
        (*pRF_varUsedPtr)[i][j] = 0;
      }
    }
    (*pRF_varUsed) --;
  }
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    PROTECT(sexpVector[RF_DPTH_ID] = NEW_NUMERIC(splitDepthSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_DPTH_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_DPTH_ID]));
    *p_splitDepth = NUMERIC_POINTER(sexpVector[RF_DPTH_ID]);
    sexpIndex ++;
    RF_splitDepthPtr  = (double ***) vvector(1, dpthDimOne);
    for (j = 1; j <= dpthDimOne; j++) {
      RF_splitDepthPtr[j] = (double **) vvector(1, RF_xSize);
      for (k = 1; k <= RF_xSize; k++) {
        RF_splitDepthPtr[j][k]  = (*p_splitDepth) + ((j-1) * RF_xSize * RF_observationSize) + ((k-1) * RF_observationSize) - 1;
      }
    }
    for (i = 1; i <= RF_observationSize; i++) {
      for (j = 1; j <= dpthDimOne; j++) {
        for (k = 1; k <= RF_xSize; k++) {
          RF_splitDepthPtr[j][k][i] = 0;
        }
      }
    }
  }
  if (RF_opt & OPT_MEMB) {
    PROTECT(sexpVector[RF_NMBR_ID] = NEW_INTEGER(RF_forestSize * obsSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_NMBR_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_NMBR_ID]));
    *pRF_terminalNodeMembership = (uint*) INTEGER_POINTER(sexpVector[RF_NMBR_ID]);
    sexpIndex++;
    RF_terminalNodeMembershipPtr = (uint **) vvector(1, RF_forestSize);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_terminalNodeMembershipPtr)[i] = (*pRF_terminalNodeMembership) + ((i-1) * obsSize) - 1;
    }
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= obsSize; j++) {
        (RF_terminalNodeMembershipPtr)[i][j] = 0;
      }
    }
    (*pRF_terminalNodeMembership) --;
    PROTECT(sexpVector[RF_BMBR_ID] = NEW_INTEGER(RF_forestSize * obsSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_BMBR_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_BMBR_ID]));
    *pRF_bootstrapMembership = (uint*) INTEGER_POINTER(sexpVector[RF_BMBR_ID]);
    sexpIndex++;
    RF_bootstrapMembershipPtr = (uint **) vvector(1, RF_forestSize);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_bootstrapMembershipPtr)[i] = (*pRF_bootstrapMembership) + ((i-1) * obsSize) - 1;
    }
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= obsSize; j++) {
        (RF_bootstrapMembershipPtr)[i][j] = 0;
      }
    }
    (*pRF_bootstrapMembership) --;
  }
  return (sexpIndex);
}
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root) {
  uint obsSize;
  uint xVimpSize;
  uint rspSize;
  uint     perfDim;
  uint     ensbDimOne;
  uint     ensbDimTwo;
  uint     dpthDimOne;
  uint i, j;
  obsSize        = 0;  
  xVimpSize      = 0;  
  rspSize        = 0;  
  dpthDimOne     = 0;  
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      perfDim = ensbDimOne = RF_eventTypeSize;
    }
    else {
      perfDim = ensbDimOne = 1;
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      perfDim = RF_rFactorSize[1] + 1;
      ensbDimOne = 1;
    }
    else {
      perfDim = 1;
      ensbDimOne = 1;
    }
  }
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
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    if (RF_opt & OPT_SPLDPTH_F) {
      dpthDimOne = 1;
    }
    else {
      dpthDimOne = RF_forestSize;
    }
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    rspSize = RF_frSize;
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
      }
    }
    break;
  default:
    obsSize = RF_observationSize;
    rspSize = RF_rSize;
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_JOIN) {
        xVimpSize = 1;
      }
      else {
        xVimpSize = RF_intrPredictorSize;
      }
    }
    break;
  }
  if (RF_opt & OPT_FENS) {
    for (j = 1; j <= ensbDimOne; j++) {
      free_vvector(RF_fullEnsemblePtr[j], 1, ensbDimTwo);
    }
    free_vvector(RF_fullEnsemblePtr, 1, ensbDimOne);
    free_uivector(RF_fullEnsembleDen, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_vvector(RF_fullMRTPtr, 1, ensbDimOne);
      if (RF_eventTypeSize == 1) {
        free_vvector(RF_fullSRVPtr, 1, ensbDimTwo);
      }
      else {
        for (j = 1; j <= ensbDimOne; j++) {
          free_vvector(RF_fullCIFPtr[j], 1, ensbDimTwo);
        }
        free_vvector(RF_fullCIFPtr, 1, ensbDimOne);
      }
    }
  }
  if (RF_opt & OPT_OENS) {
    for (j = 1; j <= ensbDimOne; j++) {
      free_vvector(RF_oobEnsemblePtr[j], 1, ensbDimTwo);
    }
    free_vvector(RF_oobEnsemblePtr, 1, ensbDimOne);
    free_uivector(RF_oobEnsembleDen, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_vvector(RF_oobMRTPtr, 1, ensbDimOne);
      if (RF_eventTypeSize == 1) {
        free_vvector(RF_oobSRVPtr, 1, ensbDimTwo);
      }
      else {
        for (j = 1; j <= ensbDimOne; j++) {
          free_vvector(RF_oobCIFPtr[j], 1, ensbDimTwo);
        }
        free_vvector(RF_oobCIFPtr, 1, ensbDimOne);
      }
    }
  }
  if (RF_opt & OPT_PERF) {
    free_vvector(RF_performancePtr, 1, RF_forestSize);
  }
  if (RF_opt & OPT_TREE) {
    for (i = 1; i <= RF_forestSize; i++) {
      freeTree(0, root[i], TRUE);
    }
  }
  if (RF_opt & OPT_MISS) {
    if (rspSize > 0) {
      free_vvector(RF_sImputeResponsePtr, 1, rspSize);
    }
    free_vvector(RF_sImputePredictorPtr, 1, RF_xSize);
  }
  if (RF_opt & OPT_OMIS) {
    if (rspSize > 0) {
      free_vvector(RF_sOOBImputeResponsePtr, 1, rspSize);
    }
    free_vvector(RF_sOOBImputePredictorPtr, 1, RF_xSize);
  }
  if (RF_opt & OPT_VIMP) {
    for (j = 1; j <= xVimpSize; j++) {
      free_uivector(RF_vimpEnsembleDen[j], 1, obsSize);
    }
    free_vvector(RF_vimpEnsembleDen, 1, xVimpSize);
    free_vvector(RF_importancePtr, 1, xVimpSize);
    if(RF_opt & OPT_VIMP_LEOB) { 
      free_dmatrix3(RF_vimpLeo, 1, RF_forestSize, 1, xVimpSize, 1, perfDim);
      free_dmatrix(RF_perfLeo, 1, RF_forestSize, 1, perfDim);
    }
    free_dmatrix(RF_vimpOutcome, 1, xVimpSize, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_dmatrix4(RF_sVimpEnsemble, 1, xVimpSize, 1, ensbDimOne, 1, ensbDimTwo, 1, obsSize);
    }
    else {
      if (RF_rFactorCount > 0) {
        free_dmatrix3(RF_cVimpEnsemble, 1, xVimpSize, 1, ensbDimTwo, 1, obsSize);
      }
    }
  }
  if (RF_opt & OPT_VUSE) {
    if (RF_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      free_vvector(RF_varUsedPtr, 1, RF_forestSize);
    }
    else {
      free_uimatrix(RF_varUsedPtr, 1, RF_forestSize, 1, RF_xSize);
    }
  }  
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    for (j = 1; j <= dpthDimOne; j++) {
      free_vvector(RF_splitDepthPtr[j], 1, RF_xSize);
    }
    free_vvector(RF_splitDepthPtr, 1, dpthDimOne);
  }  
  if (RF_opt & OPT_MEMB) {
    free_vvector(RF_terminalNodeMembershipPtr, 1, RF_forestSize);
    free_vvector(RF_bootstrapMembershipPtr,    1, RF_forestSize);
  }
}
uint stackVariableOutputObjects(uint     totalNodeCount,
                                uint     totalMWCPCount,
                                uint   **pRF_treeID,
                                uint   **pRF_nodeID,
                                uint   **pRF_parmID,                                   
                                double **pRF_contPT,
                                uint   **pRF_mwcpSZ,
                                uint   **pRF_mwcpPT,
                                uint     sexpIndex,
                                char   **sexpString,
                                SEXP    *sexpVector) {
  if (RF_opt & OPT_TREE) {
    PROTECT(sexpVector[RF_TREE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RF_NODE_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RF_PARM_ID] = NEW_INTEGER(totalNodeCount));
    PROTECT(sexpVector[RF_CONT_PT] = NEW_NUMERIC(totalNodeCount));    
    PROTECT(sexpVector[RF_MWCP_SZ] = NEW_INTEGER(totalNodeCount));  
    PROTECT(sexpVector[RF_MWCP_PT] = NEW_INTEGER(totalMWCPCount));  
    *pRF_treeID = (uint*) INTEGER_POINTER(sexpVector[RF_TREE_ID]);
    *pRF_nodeID = (uint*) INTEGER_POINTER(sexpVector[RF_NODE_ID]);
    *pRF_parmID = (uint*) INTEGER_POINTER(sexpVector[RF_PARM_ID]);
    *pRF_contPT = NUMERIC_POINTER(sexpVector[RF_CONT_PT]);
    *pRF_mwcpSZ = (uint*) INTEGER_POINTER(sexpVector[RF_MWCP_SZ]);
    *pRF_mwcpPT = (uint*) INTEGER_POINTER(sexpVector[RF_MWCP_PT]);
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_TREE_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_TREE_ID]));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_NODE_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_NODE_ID]));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PARM_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_PARM_ID]));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_CONT_PT]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_CONT_PT]));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_MWCP_SZ]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_MWCP_SZ]));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_MWCP_PT]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_MWCP_PT]));
    (*pRF_treeID) --;
    (*pRF_nodeID) --;
    (*pRF_parmID) --;
    (*pRF_contPT) --;
    (*pRF_mwcpSZ) --;
    (*pRF_mwcpPT) --;
  }
  return (sexpIndex);
}
