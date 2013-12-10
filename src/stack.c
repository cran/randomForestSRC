////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.4
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
#include       "factorOps.h"
#include          "impute.h"
#include        "treeUtil.h"
#include           "stack.h"
void initializeTimeArrays(char mode) {
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
  if (size > 0) {
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
  else {
    *factorCount = 0;
  }
}
void unstackFactorArrays() {
  uint j, k;
  if (RF_rSize > 0) {
    free_uivector(RF_rFactorMap, 1, RF_rSize);
    if (RF_rFactorCount > 0) {
      free_uivector(RF_rFactorIndex, 1, RF_rFactorCount);
      free_uivector(RF_rFactorSize, 1, RF_rFactorCount);
    }
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
  uint i, j;
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
  if (RF_rSize > 0) {
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
      updateTimeIndexArray(0, 
                           NULL, 
                           RF_observationSize, 
                           RF_responseIn[RF_timeIndex], 
                           TRUE, 
                           FALSE, 
                           RF_masterTimeIndexIn);
    }
    if (RF_statusIndex > 0) {
      RF_status = (double **) vvector(1, RF_forestSize);
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_status[i] = RF_responseIn[RF_statusIndex];
      }
    }
  }
  else {
    for (i = 1 ; i <= RF_forestSize; i++) {
      RF_response[i] = NULL;
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
  }
  else {
    stackMissingSignatures(RF_observationSize,
                           RF_rSize,
                           RF_responseIn,
                           RF_observationIn,
                           RF_mRecordMap,
                           RF_mRecordSize,
                           & RF_mRecordIndex,
                           & RF_mpIndexSize,
                           & RF_mpSign,
                           & RF_mpIndex,
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
    RF_fresponse = (double ***) vvector(1, RF_forestSize);
    if (RF_frSize > 0) {
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
      for (i = 1 ; i <= RF_forestSize; i++) {
        RF_fresponse[i] = NULL;
      }
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
                             & RF_fmpIndexSize,
                             & RF_fmpSign,
                             & RF_fmpIndex,
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
      dualUseFlag = TRUE;
      mFlag = ACTIVE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
    }
    break;
  case RF_REST:
    if (RF_mRecordSize > 0) {
      recordSize = RF_mRecordSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      recordSize = RF_mRecordSize;
      dualUseFlag = TRUE;
      mFlag = FALSE;
    }
    else {
      RF_opt = RF_opt & (~OPT_MISS);
      RF_nImpute = 1;
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
    RF_mTermList = (Terminal ***) vvector(1, RF_forestSize);
    RF_mTermMembership = (Terminal ***) vvector(1, RF_forestSize);
  }
  if (RF_rFactorCount + RF_xFactorCount > 0) {
    initializeFactorArrays(mode);
  }
  return result;
}
void unstackMissingArrays(char mode) {
  char dualUseFlag;
  uint recordSize;
  free_vvector(RF_response, 1, RF_forestSize);
  if (RF_rSize > 0) {
    if (RF_timeIndex > 0) {
      free_vvector(RF_time, 1, RF_forestSize);
      free_vvector(RF_masterTimeIndex, 1, RF_forestSize);
    }
    if (RF_statusIndex > 0) {
      free_vvector(RF_status, 1, RF_forestSize);
    }
  }
  free_vvector(RF_observation, 1, RF_forestSize);
  free_uivector(RF_mRecordMap, 1, RF_observationSize);
  if (RF_mRecordSize == 0) {
  }
  else {
    unstackMissingSignatures(RF_rSize,
                             RF_mRecordSize,
                             RF_mRecordIndex,
                             RF_mpIndexSize,
                             RF_mpSign,
                             RF_mpIndex,
                             RF_mrFactorSize,
                             RF_mrFactorIndex,
                             RF_mxFactorSize,
                             RF_mxFactorIndex);
  }
  if (mode == RF_PRED) {
    free_vvector(RF_fobservation, 1, RF_forestSize);
    free_uivector(RF_fmRecordMap, 1, RF_fobservationSize);
    free_vvector(RF_fresponse, 1, RF_forestSize);
    if (RF_frSize > 0) {
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
                               RF_fmpIndexSize,
                               RF_fmpSign,
                               RF_fmpIndex,
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
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      dualUseFlag = TRUE;
      recordSize = RF_mRecordSize;
    }
    break;
  }  
  if (dualUseFlag == TRUE) {
    free_cmatrix(RF_dmRecordBootFlag, 1, RF_forestSize, 1, recordSize);
    free_vvector(RF_mTermList, 1, RF_forestSize);
    free_vvector(RF_mTermMembership, 1, RF_forestSize);
  }
}
void stackMissingSignatures(uint     obsSize, 
                            uint     rspSize,
                            double **responsePtr,
                            double **predictorPtr,
                            uint    *recordMap,
                            uint     recordSize, 
                            uint   **p_recordIndex, 
                            uint    *p_pIndexSize,
                            int   ***p_pSign, 
                            int    **p_pIndex,
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
  *p_pSign = imatrix(1, rspSize + RF_xSize, 1, recordSize);
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize + RF_xSize; i++) {
      (*p_pSign)[i][j] = 0;
    }
  }
  for (j = 1; j <= recordSize; j++) {
    for (i = 1; i <= rspSize; i++) {
      if (ISNA(responsePtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[i][j] = 1;
      }
    }
    for (i = 1; i <= RF_xSize; i++) {
      if (ISNA(predictorPtr[i][(*p_recordIndex)[j]])) {
        (*p_pSign)[rspSize + i][j] = 1;
      }
    }
  }
  *pRF_mStatusFlag = *pRF_mTimeFlag = *pRF_mResponseFlag = *pRF_mPredictorFlag = FALSE;
  *p_pIndex = ivector(1, rspSize + RF_xSize);
  *p_pIndexSize = 0;
  for (i = 1; i <= rspSize; i++) {
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] = - i;
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
    (*p_pIndex)[i] = 0;
    for (j = 1; j <= recordSize; j++) {
      if ((*p_pSign)[i][j] == 1) {
        (*p_pIndexSize) ++;
        (*p_pIndex)[*p_pIndexSize] =  i - rspSize;
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
  for (p = 1; p <= *p_pIndexSize; p++) {
    if ((*p_pIndex)[p] < 0) {
      if (strcmp(RF_rType[(uint) abs((*p_pIndex)[p])], "C") == 0) {
        (*pRF_mrFactorSize) ++;
        (*pRF_mrFactorIndex)[*pRF_mrFactorSize] = (uint) abs((*p_pIndex)[p]);
      }
    }
    else {
      if (strcmp(RF_xType[(*p_pIndex)[p]], "C") == 0) {
        (*pRF_mxFactorSize) ++;
        (*pRF_mxFactorIndex)[*pRF_mxFactorSize] = (*p_pIndex)[p];
      }
    }
  }
}
void unstackMissingSignatures(uint      rspSize,
                              uint      recordSize, 
                              uint     *recordIndex, 
                              uint      vIndexSize,
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
            Rprintf("\nRF-SRC:  Factor level in data inconsistent with number of levels indicated:  %10d %10.4f", factorLevel, RF_rLevels[RF_rFactorIndex[j]]);
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
            Rprintf("\nRF-SRC:  Factor level in data inconsistent with number of levels indicated:  %10d %10.4f", factorLevel, RF_xLevels[RF_xFactorIndex[j]]);
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
  int    **mpSign;
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
    if ((RF_splitRule == SURV_CR_LAU) || (RF_splitRule == SURV_CR_LOG)) {
      RF_opt = RF_opt | OPT_COMP_RISK;
    }
    else {
      RF_opt = RF_opt & (~OPT_COMP_RISK);
    }
  }
  else {
  }
  getEventTypeSize(RF_observationSize, 
                   RF_responseIn[RF_statusIndex], 
                   RF_mRecordMap, 
                   RF_mpSign,  
                   overrideFlag, 
                   & RF_eventTypeSize,
                   & RF_mStatusSize,
                   RF_eventType);
  if (mode == RF_GROW) {
    if (RF_splitRule == RAND_SPLIT) {
      if (RF_eventTypeSize == 1) {
      }
      else {
        RF_opt = RF_opt | OPT_COMP_RISK;
      }
    }
    if ((RF_splitRule == SURV_CR_LAU) || (RF_splitRule == SURV_CR_LOG)) {
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
    RF_eventTypeIndex  = uivector(1, RF_eventType[RF_eventTypeSize]);
    for (j = 1; j <= RF_eventType[RF_eventTypeSize]; j++) {
      RF_eventTypeIndex[j] = 0;
    }
    for (j = 1; j <= RF_eventTypeSize; j++) {
      RF_eventTypeIndex[RF_eventType[j]] = j;
    }
  switch (mode) { 
  case RF_PRED: 
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB) | (RF_opt & OPT_VIMP)) { 
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
    getEventTypeSize(RF_fobservationSize, 
                     RF_fresponseIn[RF_statusIndex], 
                     RF_fmRecordMap,
                     RF_fmpSign, 
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
      if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB) | (RF_opt & OPT_VIMP)) { 
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
      mpSign = RF_mpSign;
      mRecordMap = RF_mRecordMap; 
    } 
    else { 
      obsSize = RF_fobservationSize;
      statusPtr = RF_fresponseIn[RF_statusIndex]; 
      mpSign = RF_fmpSign; 
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
        if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) { 
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
        if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) { 
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
                      int **mpSign, 
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
      if (mpSign[RF_statusIndex][mRecordMap[i]] == 0) { 
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
    free_uivector(RF_eventTypeIndex, 1, RF_eventType[RF_eventTypeSize]);
  free_uivector(RF_eventType, 1, RF_observationSize);
  if (RF_eventTypeSize > 1) { 
    if (mode == RF_PRED) { 
      if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB) | (RF_opt & OPT_VIMP)) { 
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
                    RF_mpSign, 
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
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB) | (RF_opt & OPT_VIMP)) {
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
                      RF_fmpSign,  
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
                       int     **mpSign,  
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
        if (mpSign[RF_rFactorIndex[j]][mRecordMap[i]] == 0) {
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
