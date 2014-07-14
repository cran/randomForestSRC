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
#include         "nodeOps.h"
#include        "treeUtil.h"
#include          "impute.h"
#include      "importance.h"
#include        "survival.h"
#include       "survivalE.h"
#include  "classification.h"
#include      "regression.h"
#include       "rfsrcUtil.h"
void getVariablesUsed(uint treeID, Node *parent, uint *varUsedVector) {
  if (RF_tLeafCount[treeID] > 0) {
    if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
      varUsedVector[parent -> splitParameter] ++;
      getVariablesUsed(treeID, parent ->  left, varUsedVector);
      getVariablesUsed(treeID, parent -> right, varUsedVector);
    }
  }
  return;
}
void updateTerminalNodeOutcomes (uint b) {
  if (RF_tLeafCount[b] > 0) {
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (uint leaf = 1; leaf <= RF_tLeafCount[b]; leaf++) {
        getAtRiskAndEventCounts(b, leaf);
        getLocalRatio(b, leaf);
        getLocalSurvival(b, leaf);
        if (!(RF_opt & OPT_COMP_RISK)) {
          getLocalNelsonAalen(b, leaf);
        }
        else {
          getLocalCSH(b, leaf);
          getLocalCIF(b, leaf);
        }
        unstackAtRiskAndEventCounts(RF_tNodeList[b][leaf]);
        getSurvival(b, leaf);
        if (!(RF_opt & OPT_COMP_RISK)) {
          getNelsonAalen(b, leaf);
        }
        else {
          getCSH(b, leaf);
          getCIF(b, leaf);
        }
        getMortality(b, leaf);
        freeTerminalNodeLocalSurvivalStructures(RF_tNodeList[b][leaf]);
      }  
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        getMultiClassProb(b);
      }
      else {
        getMeanResponse(b);
      }
    }
  }  
}
void updateEnsembleCalculations (char      multipleImputeFlag,
                                 uint      mode,
                                 uint      b) {
  uint      obsSize;
  double  **ensemblePtr;
  uint     *ensembleDenPtr;
  double   *outcome;
  double  **conditionalOutcome;
  double  **responsePtr;
  uint     *denominatorCopy;
  char      respImputeFlag;
  uint      thisSerialTreeCount;
  uint      j;
  respImputeFlag       = FALSE; 
  responsePtr          = NULL;  
  obsSize              = 0;     
  outcome              = NULL;  
  conditionalOutcome   = NULL;  
  ensemblePtr          = NULL;  
  ensembleDenPtr       = NULL;  
  denominatorCopy      = NULL;  
  thisSerialTreeCount  = 0;     
  updateTerminalNodeOutcomes(b);
  if (RF_tLeafCount[b] > 0) {
    switch (mode) {
    case RF_PRED:
      obsSize = RF_fobservationSize;
      if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
        ensembleDenPtr = RF_fullEnsembleDen;
        if (RF_opt & OPT_PERF_CALB) {
          if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
            ensemblePtr = RF_fullEnsemblePtr[1];
          }
        }
      }
      break;
    default:
      obsSize = RF_observationSize;
      if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
        ensembleDenPtr = RF_oobEnsembleDen;
        if (RF_opt & OPT_PERF_CALB) {
          if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
            ensemblePtr = RF_oobEnsemblePtr[1];
          }
        }
      }
      break;
    }
    outcome = dvector(1, obsSize);
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
      denominatorCopy = uivector(1, obsSize);
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (RF_opt & OPT_PERF) {
        if (RF_opt & OPT_COMP_RISK) {
          conditionalOutcome = dmatrix(1, RF_eventTypeSize, 1, obsSize);
        }
      }
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        if (RF_opt & OPT_PERF_CALB) {
          conditionalOutcome = dmatrix(1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]], 1, obsSize);
        }
      }
      else {
      }
    }
  }
#ifdef SUPPORT_OPENMP
#pragma omp critical (_update_ensemble)
#endif
  { 
    if (RF_tLeafCount[b] > 0) {
      RF_serialTreeIndex[++RF_serialTreeCount] = b;
      thisSerialTreeCount = RF_serialTreeCount;
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        updateEnsembleSurvival(mode, b);
        if (RF_opt & OPT_PERF) {
          switch (mode) {
          case RF_PRED:
            ensemblePtr = RF_fullMRTPtr;
            break;
          default:
            ensemblePtr = RF_oobMRTPtr;
            break;
          }
          if (!(RF_opt & OPT_COMP_RISK)) {
            getEnsembleMortality(mode, b, obsSize, ensemblePtr, ensembleDenPtr, outcome);
          }
          else {
            getEnsembleMortalityCR(mode, b, obsSize, ensemblePtr, ensembleDenPtr, conditionalOutcome);
          }
        }  
      }  
      else {
        if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
          updateEnsembleMultiClass(mode, b, outcome);
          if (RF_opt & OPT_PERF_CALB) {
            copyEnsemble(mode, obsSize, ensemblePtr, conditionalOutcome);
          }
        }
        else {
          updateEnsembleMean(mode, b, outcome);
        }
      }
      if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
        copyDenominator(mode, obsSize, ensembleDenPtr, denominatorCopy);
        respImputeFlag = stackAndImputePerfResponse(mode,
                                                    multipleImputeFlag,
                                                    b,
                                                    thisSerialTreeCount,
                                                    RF_rSize,
                                                    &responsePtr);
      }
      else {
        respImputeFlag = FALSE;
      }
    }  
    else {
      RF_serialTreeIndex[++RF_serialTreeCount] = b;
    }
  } 
  if (RF_tLeafCount[b] > 0) {
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
      getPerformance(b,
                     mode,
                     obsSize,
                     responsePtr,
                     outcome,
                     conditionalOutcome,
                     denominatorCopy,
                     RF_performancePtr[thisSerialTreeCount]);
      unstackImputeResponse(respImputeFlag, RF_rSize, obsSize, responsePtr);
    }
    for (j = 1; j <= RF_tLeafCount[b]; j++) {
      freeTerminalNodeSurvivalStructures(RF_tNodeList[b][j]);
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (RF_opt & OPT_PERF) {
        if (RF_opt & OPT_COMP_RISK) {
          free_dmatrix(conditionalOutcome, 1, RF_eventTypeSize, 1, obsSize);
        }
      }
    }
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        if (RF_opt & OPT_PERF_CALB) {
          free_dmatrix(conditionalOutcome, 1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]], 1, obsSize);
        }
      }
      else {
      }
    }
    free_dvector(outcome, 1, obsSize);
    if ((RF_opt & OPT_PERF) | (RF_opt & OPT_PERF_CALB)) {
      free_uivector(denominatorCopy, 1, obsSize);
    }
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_LEOB) {
        summarizeTreePerformance(mode, b);
      }
    }
  }  
}
void copyDenominator(uint mode, uint size, uint *denomPtr, uint *denominatorCopy) {
  uint i;
  for (i = 1; i <= size; i++) {
    denominatorCopy[i] = denomPtr[i];
  }
}
void copyEnsemble(uint mode, uint size, double **ensemblePtr, double **ensembleCopy) {
  uint j, k;
  for (k = 1; k <= size; k++) {
    for (j = 1; j <= RF_rFactorSize[1]; j++) {
      ensembleCopy[j][k] = ensemblePtr[j][k];
    }
  }
}
char stackAndImputePerfResponse(uint      mode,
                                char      multipleImputeFlag,
                                uint      treeID,
                                uint      serialID,
                                uint      rSize,
                                double ***responsePtr) {
  uint     obsSize;
  char     imputeFlag;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  imputeFlag = FALSE;
  if (mode != RF_PRED) {
    *responsePtr = RF_response[treeID];
    if (multipleImputeFlag == FALSE) {
      if (RF_mRecordSize > 0) {
        if(RF_mResponseFlag == TRUE) {
          imputeFlag = TRUE;
        }
      }
    }
  } 
  else {
    *responsePtr = RF_fresponse[treeID];
    if (RF_fmRecordSize > 0) {
      if(RF_fmResponseFlag == TRUE) {
        imputeFlag = TRUE;
      }
    }
  }  
  *responsePtr = stackAndImputeGenericResponse(imputeFlag, mode, rSize, obsSize, treeID, serialID, *responsePtr);
  return imputeFlag;
}
double **stackAndImputeGenericResponse(char flag, uint mode, uint rSize, uint obsSize, uint treeID, uint serialID, double **responsePtr) {
  uint i, p;
  double **mResponsePtr;
  if (flag == TRUE) {
    mResponsePtr   = dmatrix(1, rSize, 1, obsSize);
    for (p = 1; p <= rSize; p++) {
      for (i = 1; i <= obsSize; i++) {
        mResponsePtr[p][i] = responsePtr[p][i];
      }
    }
    imputeResponse(mode, serialID, mResponsePtr);
  }
  else {
    mResponsePtr = responsePtr;
  }
  return mResponsePtr;
}
void unstackImputeResponse(char flag, uint rSize, uint obsSize, double **mResponsePtr) {
  if (flag == TRUE) {
    free_dmatrix(mResponsePtr, 1, rSize, 1, obsSize);
  }
}
void getPerformance(uint      treeID,
                    uint      mode,
                    uint      obsSize,
                    double  **responsePtr,
                    double   *outcome,
                    double  **conditionalOutcome,
                    uint     *denomPtr,
                    double   *performancePtr) {
  uint      j;
  double    concordanceIndex;
  double   *condPerformanceVector;
  concordanceIndex = NA_REAL;  
  condPerformanceVector = stackCondPerformance();
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (!(RF_opt & OPT_COMP_RISK)) {
      concordanceIndex = getConcordanceIndex(1,
                                             obsSize,
                                             responsePtr[RF_timeIndex],
                                             responsePtr[RF_statusIndex],
                                             outcome,
                                             denomPtr);
    }
    else {
      getCRPerformance(mode,
                       obsSize,
                       responsePtr,
                       conditionalOutcome,
                       denomPtr,
                       condPerformanceVector);
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      if (RF_opt & OPT_PERF_CALB) {
        concordanceIndex = getBrierScore(obsSize,
                                         responsePtr[RF_rTarget],
                                         conditionalOutcome,
                                         denomPtr,
                                         condPerformanceVector);
      }
      else {
        concordanceIndex = getClassificationIndex( obsSize,
                                                   responsePtr[RF_rTarget],
                                                   outcome,
                                                   denomPtr);
        getConditionalClassificationIndex( obsSize,
                                           responsePtr[RF_rTarget],
                                           outcome,
                                           denomPtr,
                                           condPerformanceVector);
      }
    }
    else {
      concordanceIndex = getMeanSquareError( obsSize,
                                             responsePtr[RF_rTarget],
                                             outcome,
                                             denomPtr);
    }
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (!(RF_opt & OPT_COMP_RISK)) {
      performancePtr[1] = concordanceIndex;
    }
    else {
      for (j=1; j <= RF_eventTypeSize; j++) {
        performancePtr[j] = condPerformanceVector[j];
      }
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      if (RF_opt & OPT_PERF_CALB) {
        performancePtr[1] = concordanceIndex;
        for (j=1; j <=RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
          performancePtr[1+j] = condPerformanceVector[j];
        }
      }
      else {
        performancePtr[1] = concordanceIndex;
        for (j=1; j <=RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; j++) {
          performancePtr[1+j] = condPerformanceVector[j];
        }
      }
    }
    else {
      performancePtr[1] = concordanceIndex;
    }
  }
  unstackCondPerformance(condPerformanceVector);
}
double *stackCondPerformance() {
  double *cpv;
  cpv = NULL;
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_COMP_RISK) {
      cpv = dvector(1, RF_eventTypeSize);
    }
  }
  else {
    if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      cpv = dvector(1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]]);
    }
  }
  return cpv;
}
void unstackCondPerformance(double *cpv) {
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_COMP_RISK) {
      free_dvector(cpv, 1, RF_eventTypeSize);
    }
  }
  else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
      free_dvector(cpv, 1, RF_rFactorSize[RF_rFactorMap[RF_rTarget]]);
    }
  }
}
void finalizeEnsembleEstimates(uint mode) {
  char oobFlag, fullFlag;
  uint      obsSize;
  double ***ensemblePtr;
  uint     *ensembleDenPtr;
  double  **ensSRVPtr;
  double ***ensCIFPtr;
  double  **ensMRTPtr;
  uint i, j, k;
  oobFlag = fullFlag = FALSE;
  ensemblePtr     = NULL;  
  ensembleDenPtr  = NULL;  
  ensSRVPtr       = NULL;  
  ensCIFPtr       = NULL;  
  ensMRTPtr       = NULL;  
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_OENS) {
      oobFlag = TRUE;
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensemblePtr = RF_oobEnsemblePtr;
      ensembleDenPtr = RF_oobEnsembleDen;
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        ensMRTPtr = RF_oobMRTPtr;
        if (!(RF_opt & OPT_COMP_RISK)) {
          ensSRVPtr = RF_oobSRVPtr;
        }
        else {
          ensCIFPtr = RF_oobCIFPtr;
        }
      }
    }
    else {
      if (fullFlag == TRUE) {
        ensemblePtr = RF_fullEnsemblePtr;
        ensembleDenPtr = RF_fullEnsembleDen;
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          ensMRTPtr = RF_fullMRTPtr;
          if (!(RF_opt & OPT_COMP_RISK)) {
            ensSRVPtr = RF_fullSRVPtr;
          }
          else {
            ensCIFPtr = RF_fullCIFPtr;
          }
        }
      }
    }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      for (i = 1; i <= obsSize; i++) {
        if (ensembleDenPtr[i] != 0) {
          if (!(RF_opt & OPT_COMP_RISK)) {
            ensMRTPtr[1][i] = ensMRTPtr[1][i] / ensembleDenPtr[i];
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              ensemblePtr[1][k][i] = ensemblePtr[1][k][i] / ensembleDenPtr[i];
              ensSRVPtr[k][i]  = ensSRVPtr[k][i] / ensembleDenPtr[i];
            }
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              ensMRTPtr[j][i] = ensMRTPtr[j][i] / ensembleDenPtr[i];
              for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                ensemblePtr[j][k][i] = ensemblePtr[j][k][i] / ensembleDenPtr[i];
                ensCIFPtr[j][k][i] = ensCIFPtr[j][k][i] / ensembleDenPtr[i];
              }
            }
          }
        }
        else {
          if (!(RF_opt & OPT_COMP_RISK)) {
            ensMRTPtr[1][i] = NA_REAL;
            for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
              ensemblePtr[1][k][i] = NA_REAL;
              ensSRVPtr[k][i]      = NA_REAL;
            }
          }
          else {
            for(j = 1; j <= RF_eventTypeSize; j ++) {
              ensMRTPtr[j][i] = NA_REAL;
              for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                ensemblePtr[j][k][i] = NA_REAL;
                ensCIFPtr[j][k][i] = NA_REAL;
              }
            }
          }
        }
      }
    }  
    else {
      if (strcmp(RF_rType[RF_rTarget], "C") == 0) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDenPtr[i] != 0) {
            for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; k++) {
              ensemblePtr[1][k][i] = ensemblePtr[1][k][i] / ensembleDenPtr[i];
            }
          }
          else {
            for (k=1; k <= RF_rFactorSize[RF_rFactorMap[RF_rTarget]]; k++) {
              ensemblePtr[1][k][i] = NA_REAL;
            }
          }
        }
      }
      else {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDenPtr[i] != 0) {
            ensemblePtr[1][1][i] = ensemblePtr[1][1][i] / ensembleDenPtr[i];
          }
          else {
            ensemblePtr[1][1][i] = NA_REAL;
          }
        }
      }
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
