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


#ifdef SUPPORT_OPENMP
#include             <omp.h>
#endif
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
void getVariablesUsed(Node *parent, uint *varUsedVector) {
  if (((parent -> left) != NULL) && ((parent -> right) != NULL)) {
    varUsedVector[parent -> splitParameter] ++;
    getVariablesUsed(parent ->  left, varUsedVector);
    getVariablesUsed(parent -> right, varUsedVector);
  }
  return;
}
void updateEnsembleCalculations (char      multipleImputeFlag,
                                 uint      mode,
                                 Node     *rootPtr,
                                 uint      b,
                                 uint      serialID) {
  uint j;
  Node *terminalNode;
  uint      obsSize;
  double ***ensemblePtr;
  uint     *ensembleDenPtr;
  double   *outcome;
  double  **conditionalOutcome;
  double  **responsePtr;
  Node   ***vimpMembership;
  double    **treeOutcome;
  double  ****crTreeOutcome;
  double   ***mcTreeOutcome;
  uint     *denominatorCopy;
  char      respImputeFlag;
  responsePtr          = NULL;  
  obsSize              = 0;     
  vimpMembership       = NULL;  
  outcome              = NULL;  
  conditionalOutcome   = NULL;  
  treeOutcome          = NULL;  
  crTreeOutcome        = NULL;  
  mcTreeOutcome        = NULL;  
  denominatorCopy      = NULL;  
  if (RF_leafCount[b] == 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to compute performance on a rejected tree:  %10d", b);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    break;
  default:
    obsSize = RF_observationSize;
    break;
  }
  outcome = dvector(1, obsSize);
  if (RF_opt & OPT_PERF) {
    denominatorCopy = uivector(1, obsSize);
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    getAtRiskAndEventCounts(b);
    getLocalRatio(b);
    getLocalSurvival(b);
    if (RF_eventTypeSize == 1) {
      getLocalNelsonAalen(b);
    }
    else {
      getLocalCSH(b);
      getLocalCIF(b);
    }
    getSurvival(b);
    if (RF_eventTypeSize == 1) {
      getNelsonAalen(b);
    }
    else {
      getCSH(b);
      getCIF(b);
    }
    if (RF_opt & OPT_PERF) {
      if (RF_eventTypeSize > 1) {
        conditionalOutcome = dmatrix(1, RF_eventTypeSize, 1, obsSize);
      }
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      if (RF_opt & OPT_PERF_CALB) {
        conditionalOutcome = dmatrix(1, RF_rFactorSize[1], 1, obsSize);
      }
      getMultiClassProb(mode, b);
    }
    else {
      getMeanResponse(mode, b);
    }
  }
  if (RF_opt & OPT_VIMP) {
    vimpMembership = stackVimpMembership(mode);
    if (RF_opt & OPT_VIMP_LEOB) {
      stackTreeEnsemble(mode, &treeOutcome, &crTreeOutcome, &mcTreeOutcome);
    }
    getVIMPmembership(mode, b, vimpMembership);
  }
  respImputeFlag = FALSE; 
  if (RF_opt & OPT_PERF) {
    respImputeFlag = stackAndImputePerfResponse(mode, 
                                                multipleImputeFlag, 
                                                b, 
                                                serialID, 
                                                RF_rSize, 
                                                &responsePtr);
  }
#ifdef SUPPORT_OPENMP
#pragma omp critical (_update_ensemble_)
#endif
  { 
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      updateEnsembleSurvival(mode, b);
      if (RF_opt & OPT_PERF) {
        switch (mode) {
        case RF_PRED:
          obsSize = RF_fobservationSize;
          ensembleDenPtr = RF_fullEnsembleDen;
          if (RF_eventTypeSize == 1) {
            ensemblePtr = RF_fullEnsemblePtr;
          }
          else {
            ensemblePtr = RF_fullCIFPtr;
          }
          break;
        default:
          obsSize = RF_observationSize;
          ensembleDenPtr = RF_oobEnsembleDen;
          if (RF_eventTypeSize == 1) {
            ensemblePtr = RF_oobEnsemblePtr;
          }
          else {
            ensemblePtr = RF_oobCIFPtr;
          }
          break;
        }
        if (RF_eventTypeSize == 1) {
          getEnsembleMortality(mode, b, obsSize, ensemblePtr, ensembleDenPtr, outcome);
        }
        else {
          getEnsembleMortalityCR(mode, b, obsSize, ensemblePtr, ensembleDenPtr, conditionalOutcome);
        }
      }  
    }  
    else {
      if (RF_rFactorCount > 0) {
        updateEnsembleMultiClass(mode, b, outcome);
        if (RF_opt & OPT_PERF_CALB) {
          copyEnsemble(mode, conditionalOutcome);
        }
      }
      else {
        updateEnsembleMean(mode, b, outcome);
      }
    }
    if (RF_opt & OPT_VIMP) {
      updateVimpEnsemble(mode, b, vimpMembership);
      if (RF_opt & OPT_VIMP_LEOB) {
        summarizeVimpPerformance(mode, b);
      }
    }
    if (RF_opt & OPT_PERF) {
      copyDenominator(mode, denominatorCopy);
    }
  } 
  for (j = 1; j <= RF_leafCount[b]; j++) {
    terminalNode = RF_terminalNode[b][j];
    freeTerminalNodeStructures(terminalNode);
  }
  if (RF_opt & OPT_VIMP) {
    unstackVimpMembership(mode, vimpMembership);
    if (RF_opt & OPT_VIMP_LEOB) {
      unstackTreeEnsemble(mode, treeOutcome, crTreeOutcome, mcTreeOutcome);
    }
  }
  if (RF_opt & OPT_PERF) {
    getPerformance(b,
                   mode,
                   obsSize,
                   responsePtr,
                   outcome,
                   conditionalOutcome, 
                   denominatorCopy,
                   RF_performancePtr[serialID]);
  }  
  if (RF_opt & OPT_PERF) {
    unstackImputeResponse(respImputeFlag, RF_rSize, obsSize, responsePtr);
  }
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_opt & OPT_PERF) {
      if (RF_eventTypeSize > 1) {
        free_dmatrix(conditionalOutcome, 1, RF_eventTypeSize, 1, obsSize);
      }
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      if (RF_opt & OPT_PERF_CALB) {
        free_dmatrix(conditionalOutcome, 1, RF_rFactorSize[1], 1, obsSize);
      }
    }
    else {
    }
  }
  free_dvector(outcome, 1, obsSize);
  if (RF_opt & OPT_PERF) {
    free_uivector(denominatorCopy, 1, obsSize);
  }
}
void copyDenominator(uint mode, uint *denominatorCopy) {
  uint *denomPtr;
  uint obsSize;
  uint i;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    denomPtr = RF_fullEnsembleDen;
    break;
  default:
    obsSize = RF_observationSize;
    if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
      denomPtr = RF_fullEnsembleDen;
    }
    else {
      denomPtr = RF_oobEnsembleDen;
    }
    break;
  }
  for (i = 1; i <= obsSize; i++) {
    denominatorCopy[i] = denomPtr[i];
  }
}
void copyEnsemble(uint mode, double **ensembleCopy) {
  double **ensemblePtr;
  uint obsSize;
  uint j, k;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    ensemblePtr = RF_fullEnsemblePtr[1];
    break;
  default:
    obsSize = RF_observationSize;
    ensemblePtr = RF_oobEnsemblePtr[1];
    break;
  }
  for (k = 1; k <= obsSize; k++) {
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
  if ((mode == RF_GROW) || (mode == RF_REST)) {
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
    if (RF_eventTypeSize == 1) {
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
    if (RF_rFactorCount > 0) {
      if (RF_opt & OPT_PERF_CALB) {
        concordanceIndex = getBrierScore(obsSize, 
                                         responsePtr[1],
                                         conditionalOutcome,
                                         denomPtr,
                                         condPerformanceVector);
      }
      else {
        concordanceIndex = getClassificationIndex( obsSize, 
                                                   responsePtr[1],
                                                   outcome, 
                                                   denomPtr);
        getConditionalClassificationIndex( obsSize, 
                                           responsePtr[1],
                                           outcome, 
                                           denomPtr,
                                           condPerformanceVector);
      }
    }
    else {
      concordanceIndex = getMeanSquareError( obsSize, 
                                             responsePtr[1],
                                             outcome, 
                                             denomPtr);
    }
  }
    if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
      if (RF_eventTypeSize == 1) {
        performancePtr[1] = concordanceIndex;
      }
      else {
        for (j=1; j <=RF_eventTypeSize; j++) {
          performancePtr[j] = condPerformanceVector[j];
        }
      }
    }
    else {
      if (RF_rFactorCount > 0) {
        if (RF_opt & OPT_PERF_CALB) {
          performancePtr[1] = concordanceIndex;
          for (j=1; j <=RF_rFactorSize[1]; j++) {
            performancePtr[1+j] = condPerformanceVector[j];
          }
        }
        else {
          performancePtr[1] = concordanceIndex;
          for (j=1; j <=RF_rFactorSize[1]; j++) {
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
    if (RF_eventTypeSize > 1) {
      cpv = dvector(1, RF_eventTypeSize);
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      cpv = dvector(1, RF_rFactorSize[1]);
    }
  }
  return cpv;
}
void unstackCondPerformance(double *cpv) {
  if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
    if (RF_eventTypeSize > 1) {
      free_dvector(cpv, 1, RF_eventTypeSize);
    }
  }
  else {
    if (RF_rFactorCount > 0) {
      free_dvector(cpv, 1, RF_rFactorSize[1]);
    }
  }
}
void finalizeEnsembleEstimates(uint mode, uint rejectedTreeCount) {
  char finalizeFlag, oobFlag, fullFlag;
  uint      obsSize;
  double ***ensemblePtr;
  uint     *ensembleDenPtr;
  double  **ensSRVPtr;
  double ***ensCIFPtr;
  double  **ensMRTPtr;
  uint i, j, k;
  finalizeFlag = oobFlag = fullFlag = FALSE;
  ensemblePtr     = NULL;  
  ensembleDenPtr  = NULL;  
  ensSRVPtr       = NULL;  
  ensCIFPtr       = NULL;  
  ensMRTPtr       = NULL;  
  if (rejectedTreeCount < RF_forestSize) {
    finalizeFlag = TRUE;
  }
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
   if (finalizeFlag) {
     while ((oobFlag == TRUE) || (fullFlag == TRUE)) { 
      if (oobFlag == TRUE) {
        ensemblePtr = RF_oobEnsemblePtr;
        ensembleDenPtr = RF_oobEnsembleDen;
        if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
          ensMRTPtr = RF_oobMRTPtr;
          if (RF_eventTypeSize == 1) {
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
            if (RF_eventTypeSize == 1) {
              ensSRVPtr = RF_fullSRVPtr;
            }
            else {
              ensCIFPtr = RF_fullCIFPtr;
            }
          }
        }
        else {
          Rprintf("\nRF-SRC:  *** ERROR *** ");
          Rprintf("\nRF-SRC:  Unknown case in switch encountered. ");
          Rprintf("\nRF-SRC:  Please Contact Technical Support.");
          error("\nRF-SRC:  The application will now exit.\n");
        }
      }
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        for (i = 1; i <= obsSize; i++) {
          if (ensembleDenPtr[i] != 0) {
            if (RF_eventTypeSize == 1) {
              ensMRTPtr[1][i] = 0.0;
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                ensemblePtr[1][k][i] = ensemblePtr[1][k][i] / ensembleDenPtr[i];
                ensSRVPtr[k][i]  = ensSRVPtr[k][i] / ensembleDenPtr[i];   
                ensMRTPtr[1][i] += ensemblePtr[1][k][i];
              }
            }
            else {
              for(j = 1; j <= RF_eventTypeSize; j ++) {
                for (k=1; k <= RF_sortedTimeInterestSize; k++) {
                  ensemblePtr[j][k][i] = ensemblePtr[j][k][i] / ensembleDenPtr[i];
                  ensCIFPtr[j][k][i] = ensCIFPtr[j][k][i] / ensembleDenPtr[i];
                }
                ensMRTPtr[j][i] = 0.0;
                for (k=1; k <= RF_sortedTimeInterestSize - 1; k++) {
                  ensMRTPtr[j][i] += ensCIFPtr[j][k][i] * (RF_timeInterest[k+1] - RF_timeInterest[k]);;
                }
              }
            }
          }
          else {
            if (RF_eventTypeSize == 1) {
              ensMRTPtr[1][i] = NA_REAL;
              for (k = 1; k <= RF_sortedTimeInterestSize; k++) {
                ensemblePtr[1][k][i] = NA_REAL;
                ensSRVPtr[k][i]      = NA_REAL;   
                ensMRTPtr[1][i]      = NA_REAL;   
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
        if (RF_rFactorCount > 0) {
          for (i = 1; i <= obsSize; i++) {
            if (ensembleDenPtr[i] != 0) {
              for (k=1; k <= RF_rFactorSize[1]; k++) {
                ensemblePtr[1][k][i] = ensemblePtr[1][k][i] / ensembleDenPtr[i];
              }
            }
            else {
              for (k=1; k <= RF_rFactorSize[1]; k++) {
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
}
