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
#include        "treeUtil.h"
#include     "stackOutput.h"
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               double  **pRF_oobEnsemble,
                               double  **pRF_fullEnsemble,
                               double  **p_performance,
                               uint    **pRF_tLeafCount,
                               double  **pRF_proximity,
                               double  **pRF_weight,
                               double  **pRF_importance,
                               int     **pRF_seed,
                               double  **p_imputation,
                               double ***pRF_sImputeResponsePtr,
                               double ***pRF_sImputePredictorPtr,
                               uint    **pRF_varUsed,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth,
                               double  **pRF_oobEnsembleCIF,
                               double  **pRF_fullEnsembleCIF,
                               double  **pRF_oobEnsembleSRV,
                               double  **pRF_fullEnsembleSRV,
                               double  **pRF_oobEnsembleMRT,
                               double  **pRF_fullEnsembleMRT,
                               uint    **pRF_tNodeMembershipIndex,
                               uint    **pRF_pNodeMembershipIndex,
                               uint    **pRF_bootstrapMembership,
                               uint     *stackCount,
                               SEXP     *sexpVector) {
  uint sexpIndex;
  uint ensembleSize;
  uint performanceSize;
  uint proximitySize;
  uint weightSize;
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
  char     perfFlag;
  uint i,j,k,p;
  sexpIndex      = 0;  
  ensembleSize   = 0;  
  performanceSize= 0;  
  proximitySize  = 0;  
  weightSize     = 0;  
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
  perfDim = ensbDimOne = ensbDimTwo = 0;
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
    if (RF_rSize > 0) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_opt & OPT_COMP_RISK) {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = RF_eventTypeSize;
          }
          else {
            ensbDimOne = RF_eventTypeSize;
            perfDim    = 0;
          }
        }
        else {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = 1;
          }
          else {
            ensbDimOne = 1;
            perfDim    = 0;
          }
        }
      }
      else {
        if ((RF_rTarget < 1) || (RF_rTarget > RF_rSize)) {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  Target response is out of range for [C+], [R+], [M+]:  %10d  ", RF_rTarget);
          Rprintf("\nRF-SRC:  The application will now exit.\n");
          error("\nRF-SRC:  The application will now exit.\n");
        }
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          if (RF_frSize > 0) {
            perfDim = RF_rFactorSize[RF_rFactorMap[RF_rTarget]] + 1;
          }
          else {
            perfDim = 0;
          }
          ensbDimOne = 1;
        }
        else {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = 1;
          }
          else {
            ensbDimOne = 1;
            perfDim    = 0;
          }
        }
      }
      if (RF_timeIndex > 0) {
        ensbDimTwo = RF_sortedTimeInterestSize;
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          ensbDimTwo = RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
        }
        else {
          ensbDimTwo = 1;
        }
      }
    }
    else {
    }
    *stackCount = 1;
    if (RF_opt & OPT_FENS) {
      (*stackCount) += 1;
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (!(RF_opt & OPT_COMP_RISK)) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)){
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2;
      (*stackCount) += 1;
    }
    if (RF_optHigh & OPT_WGHT) {
      weightSize = obsSize * RF_observationSize;
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_NODE_STAT) {
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
      if (RF_ptnCount > 0) {
        (*stackCount) += 1;
      }
    }
    break;
  default:
    obsSize = RF_observationSize;
    mRecordSize = RF_mRecordSize;
    rspSize = RF_rSize;
    responsePtr  = RF_responseIn;
    predictorPtr = RF_observationIn;
    mRecordIndex = RF_mRecordIndex;
    if (RF_rSize == 0) {
      perfFlag = FALSE;
    }
    else {
      if (mode == RF_GROW) {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          perfFlag = TRUE;
        }
        else {
          if (RF_rSize > 1) {
            perfFlag = FALSE;
          }
          else {
            perfFlag = TRUE;
          }
        }
      }
      else {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          perfFlag = TRUE;
        }
        else {
          if (RF_rSize > 0) {
            if ((RF_rTarget < 1) || (RF_rTarget > RF_rSize)) {
              perfFlag = FALSE;  
              Rprintf("\nRF-SRC:  *** ERROR *** ");
              Rprintf("\nRF-SRC:  Target response is out of range for [C+], [R+], [M+]:  %10d  ", RF_rTarget);
              Rprintf("\nRF-SRC:  The application will now exit.\n");
              error("\nRF-SRC:  The application will now exit.\n");
            }
            else {
              perfFlag = TRUE;
            }
          }
          else {
            perfFlag = FALSE;
          }
        }
      }
    }
    if (!perfFlag) {
      RF_opt                  = RF_opt & (~OPT_PERF);
      RF_opt                  = RF_opt & (~OPT_PERF_CALB);
      RF_opt                  = RF_opt & (~OPT_VIMP);
      RF_opt                  = RF_opt & (~OPT_OENS);
      RF_opt                  = RF_opt & (~OPT_FENS);
      RF_rTarget = 0;
      perfDim = ensbDimOne = ensbDimTwo = 0;
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_opt & OPT_COMP_RISK) {
          perfDim = ensbDimOne = RF_eventTypeSize;
        }
        else {
          perfDim = ensbDimOne = 1;
        }
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          perfDim = RF_rFactorSize[RF_rFactorMap[RF_rTarget]] + 1;
          ensbDimOne = 1;
        }
        else {
          perfDim = ensbDimOne = 1;
        }
      }
      if (RF_timeIndex > 0) {
        ensbDimTwo = RF_sortedTimeInterestSize;
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          ensbDimTwo = RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
        }
        else {
          ensbDimTwo = 1;
        }
      }
    }
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
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_FENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (!(RF_opt & OPT_COMP_RISK)) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if (RF_opt & OPT_OENS) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (!(RF_opt & OPT_COMP_RISK)) {
          (*stackCount) += 2;
        }
        else {
          (*stackCount) += 2;
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) >> 1;
      (*stackCount) += 1;
    }
    if (RF_optHigh & OPT_WGHT) {
      weightSize = obsSize * RF_observationSize;
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
    if (RF_opt & OPT_NODE_STAT) {
      (*stackCount) += 1;
    }
    if (RF_opt & OPT_MISS) {
      imputationSize = (RF_xSize + rspSize + 1) * mRecordSize;
      (*stackCount) += 1;
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
      if (RF_ptnCount > 0) {
        (*stackCount) += 1;
      }
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
    RF_fullEnsemblePtr = (double ***) new_vvector(1, ensbDimOne, NRUTIL_DPTR2);
    RF_fullEnsembleDen = uivector(1, obsSize);
    for (j = 1; j <= ensbDimOne; j++) {
      RF_fullEnsemblePtr[j] = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
      RF_fullMRTPtr = (double **) new_vvector(1, ensbDimOne, NRUTIL_DPTR);
      for (j = 1; j <= ensbDimOne; j++) {
        RF_fullMRTPtr[j] = (*pRF_fullEnsembleMRT) + ((j-1) * obsSize) - 1;
      }
      for (j = 1; j <= ensbDimOne; j++) {
        for (i = 1; i <= obsSize; i++) {
          RF_fullMRTPtr[j][i] = 0.0;
        }
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        PROTECT(sexpVector[RF_FSRV_ID] = NEW_NUMERIC(ensbDimTwo * obsSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_FSRV_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_FSRV_ID]));
        *pRF_fullEnsembleSRV = NUMERIC_POINTER(sexpVector[RF_FSRV_ID]);
        sexpIndex ++;
        RF_fullSRVPtr = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
        RF_fullCIFPtr = (double ***) new_vvector(1, ensbDimOne, NRUTIL_DPTR2);
        for (j = 1; j <= ensbDimOne; j++) {
          RF_fullCIFPtr[j] = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
    RF_oobEnsemblePtr  = (double ***) new_vvector(1, ensbDimOne, NRUTIL_DPTR2);
    RF_oobEnsembleDen  = uivector(1, obsSize);
    for (j = 1; j <= ensbDimOne; j++) {
      RF_oobEnsemblePtr[j] = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
      RF_oobMRTPtr = (double **) new_vvector(1, ensbDimOne, NRUTIL_DPTR);
      for (j = 1; j <= ensbDimOne; j++) {
        RF_oobMRTPtr[j] = (*pRF_oobEnsembleMRT) + ((j-1) * obsSize) - 1;
      }
      for (j = 1; j <= ensbDimOne; j++) {
        for (i = 1; i <= obsSize; i++) {
          RF_oobMRTPtr[j][i] = 0.0;
        }
      }
      if (!(RF_opt & OPT_COMP_RISK)) {
        PROTECT(sexpVector[RF_OSRV_ID] = NEW_NUMERIC(ensbDimTwo * obsSize));
        SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_OSRV_ID]);
        SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_OSRV_ID]));
        *pRF_oobEnsembleSRV = NUMERIC_POINTER(sexpVector[RF_OSRV_ID]);
        sexpIndex ++;
        RF_oobSRVPtr = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
        RF_oobCIFPtr = (double ***) new_vvector(1, ensbDimOne, NRUTIL_DPTR2);
        for (j = 1; j <= ensbDimOne; j++) {
          RF_oobCIFPtr[j] = (double **) new_vvector(1, ensbDimTwo, NRUTIL_DPTR);
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
  if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
    PROTECT(sexpVector[RF_PERF_ID] = NEW_NUMERIC(performanceSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PERF_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_PERF_ID]));
    *p_performance = NUMERIC_POINTER(sexpVector[RF_PERF_ID]);
    sexpIndex ++;
    RF_performancePtr = (double **) new_vvector(1, RF_forestSize, NRUTIL_DPTR);
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
    PROTECT(sexpVector[RF_PROX_ID] = NEW_NUMERIC(proximitySize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PROX_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_PROX_ID]));
    *pRF_proximity = NUMERIC_POINTER(sexpVector[RF_PROX_ID]);
    sexpIndex ++;
    RF_proximityDen = dvector(1, proximitySize);
    (*pRF_proximity) --;
    for (i = 1; i <= proximitySize; i++) {
      (*pRF_proximity)[i] = 0;
      RF_proximityDen[i]  = 0;
    }
  }
  if (RF_optHigh & OPT_WGHT) {
    PROTECT(sexpVector[RF_WGHT_ID] = NEW_NUMERIC(weightSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_WGHT_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_WGHT_ID]));
    *pRF_weight = NUMERIC_POINTER(sexpVector[RF_WGHT_ID]);
    sexpIndex ++;
    RF_weightPtr = (double **) new_vvector(1, obsSize, NRUTIL_DPTR);
    for (i = 1; i <= obsSize; i++) {
      RF_weightPtr[i]  = (*pRF_weight)  + ((i-1) * RF_observationSize) - 1;
    }
    for (j = 1; j <= obsSize; j++) {
      for (k = 1; k <= RF_observationSize; k++) {
        RF_weightPtr[j][k] = 0.0;
      }
    }
  }
  if (RF_opt & OPT_LEAF) {
    PROTECT(sexpVector[RF_LEAF_ID] = NEW_INTEGER(RF_forestSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_LEAF_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_LEAF_ID]));
    *pRF_tLeafCount    = (uint*) INTEGER_POINTER(sexpVector[RF_LEAF_ID]);
    sexpIndex ++;
    (*pRF_tLeafCount) --;
    for (i = 1; i <= RF_forestSize; i++) {
      (*pRF_tLeafCount)[i] = 0;
    }
    RF_tLeafCount = *pRF_tLeafCount;
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
      *pRF_sImputeResponsePtr = (double **) new_vvector(1, rspSize, NRUTIL_DPTR);
      for (i = 1; i <= rspSize; i++) {
        (*pRF_sImputeResponsePtr)[i]  = (*p_imputation)  + (i * mRecordSize) - 1;
      }
    }
    *pRF_sImputePredictorPtr = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
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
  if (RF_opt & OPT_VIMP) {
    PROTECT(sexpVector[RF_VIMP_ID] = NEW_NUMERIC(importanceSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_VIMP_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex, mkChar(sexpString[RF_VIMP_ID]));
    *pRF_importance = NUMERIC_POINTER(sexpVector[RF_VIMP_ID]);
    sexpIndex ++;
    RF_vimpMembership = (Node ****) new_vvector(1, xVimpSize, NRUTIL_NPTR3);
    for (k = 1; k <= xVimpSize; k++) {
      RF_vimpMembership[k] = (Node ***) new_vvector(1,  RF_forestSize, NRUTIL_NPTR2);
    }
    for (k = 1; k <= xVimpSize; k++) {
      for (i = 1; i <= RF_forestSize; i++) {
        RF_vimpMembership[k][i] = NULL;
      }
    }
    RF_importancePtr = (double **) new_vvector(1, xVimpSize, NRUTIL_DPTR);
    for (k = 1; k <= xVimpSize; k++) {
      RF_importancePtr[k]  = (*pRF_importance)  + ((k-1) * perfDim) - 1;
    }
    for (k = 1; k <= xVimpSize; k++) {
      for (j = 1; j <= perfDim; j++) {
        RF_importancePtr[k][j] = NA_REAL;
      }
    }
    RF_vimpEnsembleDen  = (uint **) new_vvector(1, xVimpSize, NRUTIL_UPTR);
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
    for (p=1; p <= xVimpSize; p++) {
      for (i = 1; i <= obsSize; i++) {
        RF_vimpOutcome[p][i] = 0.0;
      }
    }
    RF_cVimpEnsemble  = NULL;
    RF_sVimpOutcome   = NULL;
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      RF_sVimpOutcome = dmatrix3(1, xVimpSize, 1, ensbDimOne, 1, obsSize);
      for (p=1; p <= xVimpSize; p++) {
        for (j = 1; j <= ensbDimOne; j++) {
          for (i = 1; i <= obsSize; i++) {
            RF_sVimpOutcome[p][j][i] = 0.0;
          }
        }
      }
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
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
      *pRF_varUsedPtr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
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
    RF_splitDepthPtr  = (double ***) new_vvector(1, dpthDimOne, NRUTIL_DPTR2);
    for (j = 1; j <= dpthDimOne; j++) {
      RF_splitDepthPtr[j] = (double **) new_vvector(1, RF_xSize, NRUTIL_DPTR);
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
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_NMBR_ID]));
    *pRF_tNodeMembershipIndex = (uint*) INTEGER_POINTER(sexpVector[RF_NMBR_ID]);
    RF_tNodeMembershipIndexPtr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
    for (i = 1; i <= RF_forestSize; i++) {
      (RF_tNodeMembershipIndexPtr)[i] = (*pRF_tNodeMembershipIndex) + ((i-1) * obsSize) - 1;
    }
    for (i = 1; i <= RF_forestSize; i++) {
      for (j = 1; j <= obsSize; j++) {
        (RF_tNodeMembershipIndexPtr)[i][j] = 0;
      }
    }
    (*pRF_tNodeMembershipIndex) --;
    if (RF_ptnCount > 0) {
      PROTECT(sexpVector[RF_PMBR_ID] = NEW_INTEGER(RF_forestSize * obsSize));
      SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_PMBR_ID]);
      SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_PMBR_ID]));
      *pRF_pNodeMembershipIndex = (uint*) INTEGER_POINTER(sexpVector[RF_PMBR_ID]);
      RF_pNodeMembershipIndexPtr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
      for (i = 1; i <= RF_forestSize; i++) {
        (RF_pNodeMembershipIndexPtr)[i] = (*pRF_pNodeMembershipIndex) + ((i-1) * obsSize) - 1;
      }
      for (i = 1; i <= RF_forestSize; i++) {
        for (j = 1; j <= obsSize; j++) {
          (RF_pNodeMembershipIndexPtr)[i][j] = 0;
        }
      }
      (*pRF_pNodeMembershipIndex) --;
    }
    PROTECT(sexpVector[RF_BMBR_ID] = NEW_INTEGER(RF_forestSize * obsSize));
    SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_BMBR_ID]);
    SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_BMBR_ID]));
    *pRF_bootstrapMembership = (uint*) INTEGER_POINTER(sexpVector[RF_BMBR_ID]);
    RF_bootstrapMembershipPtr = (uint **) new_vvector(1, RF_forestSize, NRUTIL_UPTR);
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
  uint proximitySize;
  uint rspSize;
  uint     perfDim;
  uint     ensbDimOne;
  uint     ensbDimTwo;
  uint     dpthDimOne;
  char     perfFlag;
  uint i, j, k;
  obsSize        = 0;  
  xVimpSize      = 0;  
  proximitySize  = 0;  
  rspSize        = 0;  
  dpthDimOne     = 0;  
  perfDim = ensbDimOne = ensbDimTwo = 0;
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
    if (RF_rSize > 0) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_opt & OPT_COMP_RISK) {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = RF_eventTypeSize;
          }
          else {
            ensbDimOne = RF_eventTypeSize;
            perfDim    = 0;
          }
        }
        else {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = 1;
          }
          else {
            ensbDimOne = 1;
            perfDim    = 0;
          }
        }
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          if (RF_frSize > 0) {
            perfDim = RF_rFactorSize[RF_rFactorMap[RF_rTarget]] + 1;
          }
          else {
            perfDim = 0;
          }
          ensbDimOne = 1;
        }
        else {
          if (RF_frSize > 0) {
            perfDim = ensbDimOne = 1;
          }
          else {
            ensbDimOne = 1;
            perfDim    = 0;
          }
        }
      }
      if (RF_timeIndex > 0) {
        ensbDimTwo = RF_sortedTimeInterestSize;
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          ensbDimTwo = RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
        }
        else {
          ensbDimTwo = 1;
        }
      }
    }
    else {
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2;
    }
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
    if (RF_rSize == 0) {
      perfFlag = FALSE;
    }
    else {
      if (mode == RF_GROW) {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          perfFlag = TRUE;
        }
        else {
          if (RF_rSize > 1) {
            perfFlag = FALSE;
          }
          else {
            perfFlag = TRUE;
          }
        }
      }
      else {
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          perfFlag = TRUE;
        }
        else {
          if (RF_rSize > 0) {
            perfFlag = TRUE;
          }
          else {
            perfFlag = FALSE;
          }
        }
      }
    }
    if (!perfFlag) {
      perfDim = ensbDimOne = ensbDimTwo = 0;
    }
    else {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        if (RF_opt & OPT_COMP_RISK) {
          perfDim = ensbDimOne = RF_eventTypeSize;
        }
        else {
          perfDim = ensbDimOne = 1;
        }
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          perfDim = RF_rFactorSize[RF_rFactorMap[RF_rTarget]] + 1;
          ensbDimOne = 1;
        }
        else {
          perfDim = ensbDimOne = 1;
        }
      }
      if (RF_timeIndex > 0) {
        ensbDimTwo = RF_sortedTimeInterestSize;
      }
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          ensbDimTwo = RF_rFactorSize[RF_rFactorMap[RF_rTarget]];
        }
        else {
          ensbDimTwo = 1;
        }
      }
    }
    if (RF_opt & OPT_PROX) {
      proximitySize = ((obsSize + 1)  * obsSize) / 2;
    }
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
      free_new_vvector(RF_fullEnsemblePtr[j], 1, ensbDimTwo, NRUTIL_DPTR);
    }
    free_new_vvector(RF_fullEnsemblePtr, 1, ensbDimOne, NRUTIL_DPTR2);
    free_uivector(RF_fullEnsembleDen, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_new_vvector(RF_fullMRTPtr, 1, ensbDimOne, NRUTIL_DPTR);
      if (!(RF_opt & OPT_COMP_RISK)) {
        free_new_vvector(RF_fullSRVPtr, 1, ensbDimTwo, NRUTIL_DPTR);
      }
      else {
        for (j = 1; j <= ensbDimOne; j++) {
          free_new_vvector(RF_fullCIFPtr[j], 1, ensbDimTwo, NRUTIL_DPTR);
        }
        free_new_vvector(RF_fullCIFPtr, 1, ensbDimOne, NRUTIL_DPTR2);
      }
    }
  }
  if (RF_opt & OPT_OENS) {
    for (j = 1; j <= ensbDimOne; j++) {
      free_new_vvector(RF_oobEnsemblePtr[j], 1, ensbDimTwo, NRUTIL_DPTR);
    }
    free_new_vvector(RF_oobEnsemblePtr, 1, ensbDimOne, NRUTIL_DPTR2);
    free_uivector(RF_oobEnsembleDen, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_new_vvector(RF_oobMRTPtr, 1, ensbDimOne, NRUTIL_DPTR);
      if (!(RF_opt & OPT_COMP_RISK)) {
        free_new_vvector(RF_oobSRVPtr, 1, ensbDimTwo, NRUTIL_DPTR);
      }
      else {
        for (j = 1; j <= ensbDimOne; j++) {
          free_new_vvector(RF_oobCIFPtr[j], 1, ensbDimTwo, NRUTIL_DPTR);
        }
        free_new_vvector(RF_oobCIFPtr, 1, ensbDimOne, NRUTIL_DPTR2);
      }
    }
  }
  if (RF_opt & OPT_PERF) {
    free_new_vvector(RF_performancePtr, 1, RF_forestSize, NRUTIL_DPTR);
  }
  if (RF_opt & OPT_PROX) {
    free_dvector(RF_proximityDen, 1, proximitySize);
  }
  if (RF_optHigh & OPT_WGHT) {
    free_new_vvector(RF_weightPtr, 1, obsSize, NRUTIL_DPTR);
  }
  if (RF_opt & OPT_TREE) {
    for (i = 1; i <= RF_forestSize; i++) {
      freeTree(0, root[i], TRUE);
    }
  }
  if (RF_opt & OPT_MISS) {
    if (rspSize > 0) {
      free_new_vvector(RF_sImputeResponsePtr, 1, rspSize, NRUTIL_DPTR);
    }
    free_new_vvector(RF_sImputePredictorPtr, 1, RF_xSize, NRUTIL_DPTR);
  }
  if (RF_opt & OPT_VIMP) {
    for (k = 1; k <= xVimpSize; k++) {
      free_new_vvector(RF_vimpMembership[k], 1,  RF_forestSize, NRUTIL_NPTR2);
    }
    free_new_vvector(RF_vimpMembership, 1, xVimpSize, NRUTIL_NPTR3);
    for (j = 1; j <= xVimpSize; j++) {
      free_uivector(RF_vimpEnsembleDen[j], 1, obsSize);
    }
    free_new_vvector(RF_vimpEnsembleDen, 1, xVimpSize, NRUTIL_UPTR);
    free_new_vvector(RF_importancePtr, 1, xVimpSize, NRUTIL_DPTR);
    if(RF_opt & OPT_VIMP_LEOB) {
      free_dmatrix3(RF_vimpLeo, 1, RF_forestSize, 1, xVimpSize, 1, perfDim);
      free_dmatrix(RF_perfLeo, 1, RF_forestSize, 1, perfDim);
    }
    free_dmatrix(RF_vimpOutcome, 1, xVimpSize, 1, obsSize);
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      free_dmatrix3(RF_sVimpOutcome, 1, xVimpSize, 1, ensbDimOne, 1, obsSize);
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        free_dmatrix3(RF_cVimpEnsemble, 1, xVimpSize, 1, ensbDimTwo, 1, obsSize);
      }
    }
  }
  if (RF_opt & OPT_VUSE) {
    if (RF_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      free_new_vvector(RF_varUsedPtr, 1, RF_forestSize, NRUTIL_UPTR);
    }
    else {
      free_uimatrix(RF_varUsedPtr, 1, RF_forestSize, 1, RF_xSize);
    }
  }
  if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
    for (j = 1; j <= dpthDimOne; j++) {
      free_new_vvector(RF_splitDepthPtr[j], 1, RF_xSize, NRUTIL_DPTR);
    }
    free_new_vvector(RF_splitDepthPtr, 1, dpthDimOne, NRUTIL_DPTR2);
  }
  if (RF_opt & OPT_MEMB) {
    free_new_vvector(RF_tNodeMembershipIndexPtr, 1, RF_forestSize, NRUTIL_UPTR);
    free_new_vvector(RF_bootstrapMembershipPtr,  1, RF_forestSize, NRUTIL_UPTR);
    if (RF_ptnCount > 0) {
      free_new_vvector(RF_pNodeMembershipIndexPtr, 1, RF_forestSize, NRUTIL_UPTR);
    }
  }
}
uint stackVariableOutputObjects(char     mode,
                                uint     totalNodeCount,
                                uint     totalMWCPCount,
                                uint   **pRF_treeID,
                                uint   **pRF_nodeID,
                                uint   **pRF_parmID,
                                double **pRF_contPT,
                                uint   **pRF_mwcpSZ,
                                uint   **pRF_mwcpPT,
                                double **pRF_spltST,
                                double **pRF_spltVR,
                                uint     sexpIndex,
                                char   **sexpString,
                                SEXP    *sexpVector) {
  if (mode == RF_GROW) {
    if (RF_opt & OPT_TREE) {
      PROTECT(sexpVector[RF_MWCP_PT] = NEW_INTEGER(totalMWCPCount));
      *pRF_mwcpPT = (uint*) INTEGER_POINTER(sexpVector[RF_MWCP_PT]);
      SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_MWCP_PT]);
      SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_MWCP_PT]));
      (*pRF_mwcpPT) --;
      PROTECT(sexpVector[RF_TREE_ID] = NEW_INTEGER(totalNodeCount));
      PROTECT(sexpVector[RF_NODE_ID] = NEW_INTEGER(totalNodeCount));
      PROTECT(sexpVector[RF_PARM_ID] = NEW_INTEGER(totalNodeCount));
      PROTECT(sexpVector[RF_CONT_PT] = NEW_NUMERIC(totalNodeCount));
      PROTECT(sexpVector[RF_MWCP_SZ] = NEW_INTEGER(totalNodeCount));
      *pRF_treeID = (uint*) INTEGER_POINTER(sexpVector[RF_TREE_ID]);
      *pRF_nodeID = (uint*) INTEGER_POINTER(sexpVector[RF_NODE_ID]);
      *pRF_parmID = (uint*) INTEGER_POINTER(sexpVector[RF_PARM_ID]);
      *pRF_contPT = NUMERIC_POINTER(sexpVector[RF_CONT_PT]);
      *pRF_mwcpSZ = (uint*) INTEGER_POINTER(sexpVector[RF_MWCP_SZ]);
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
      (*pRF_treeID) --;
      (*pRF_nodeID) --;
      (*pRF_parmID) --;
      (*pRF_contPT) --;
      (*pRF_mwcpSZ) --;
    }
  }
  if (RF_opt & OPT_NODE_STAT) {
      PROTECT(sexpVector[RF_SPLT_ST] = NEW_NUMERIC(totalNodeCount));
      *pRF_spltST = NUMERIC_POINTER(sexpVector[RF_SPLT_ST]);
      SET_VECTOR_ELT(sexpVector[RF_OUTP_ID], sexpIndex, sexpVector[RF_SPLT_ST]);
      SET_STRING_ELT(sexpVector[RF_STRG_ID], sexpIndex++, mkChar(sexpString[RF_SPLT_ST]));
      (*pRF_spltST) --;
      *pRF_spltVR = NULL;
  }
  return (sexpIndex);
}
