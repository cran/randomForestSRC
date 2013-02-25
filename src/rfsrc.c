////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.1.0
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
////    email:  ubk@kogalur.com
////    URL:    http://www.kogalur.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************


#include        "global.h"
#include         "trace.h"
#include        "nrutil.h"
#include        "random.h"
#include        "factor.h"
#include         "stack.h"
#include       "nodeOps.h"
#include          "tree.h"
#include      "treeUtil.h"
#include     "bootstrap.h"
#include        "impute.h"
#include     "rfsrcUtil.h"
#include      "survival.h"
#include    "importance.h"
#include         "rfsrc.h"
uint stackCount;
char  *sexpString[RF_SEXP_CNT] = {
  "",              
  "",              
  "fullEnsemble",  
  "oobEnsemble",   
  "performance",   
  "proximity",     
  "leafCount",     
  "treeID",        
  "nodeID",        
  "parmID",        
  "contPT",        
  "mwcpSZ",        
  "mwcpPT",        
  "seed",          
  "importance",    
  "imputation",    
  "oobImputation", 
  "varUsed",       
  "splitDepth",    
  "fullCIF",       
  "oobCIF",        
  "fullSurvival",  
  "oobSurvival",   
  "fullMortality", 
  "oobMortality",  
  "nodeMembership",
  "bootMembership" 
};
SEXP sexpVector[RF_SEXP_CNT];
uint     *RF_treeID_;
uint     *RF_nodeID_;
uint     *RF_parmID_;
uint     *RF_mwcpSZ_;
double   *RF_contPT_;
uint     *RF_mwcpPT_;
uint      RF_totalNodeCount;
double   *RF_importance_;
double   *RF_imputation_;
double   *RF_oobImputation_;
uint     *RF_varUsed_;
double   *RF_splitDepth_;
int      *RF_seed_;
uint     *RF_leafCount_;
double   *RF_performance_;
uint     *RF_terminalNodeMembership_;
uint     *RF_bootstrapMembership_;
double   *RF_fullEnsemble_;
double   *RF_oobEnsemble_;
double   *RF_fullEnsembleCIF_;
double   *RF_oobEnsembleCIF_;
double   *RF_fullEnsembleSRV_;
double   *RF_oobEnsembleSRV_;
double   *RF_fullEnsembleMRT_;
double   *RF_oobEnsembleMRT_;
uint     *RF_proximity_;
uint      RF_opt;
uint      RF_splitRule;
uint      RF_splitRandomRule;
uint      RF_imputeSize;
uint      RF_forestSize;
uint      RF_minimumNodeSize;
int       RF_maximumNodeDepth;
double    RF_crWeight;
uint      RF_randomCovariateCount;
double   *RF_randomCovariateWeight;
int       RF_numThreads;
uint      RF_observationSize;
uint      RF_rSize;
double   *RF_rData;
uint      RF_xSize;
double   *RF_xData;
double  **RF_responseIn;
double  **RF_observationIn;
SEXP      RF_sexp_xType;
char    **RF_xType;
int      *RF_xLevels;
SEXP      RF_sexp_rType;
char    **RF_rType;
int      *RF_rLevels;
uint      RF_fobservationSize;
uint      RF_frSize;
double   *RF_frData;
double   *RF_fxData;
double  **RF_fresponseIn;
double  **RF_fobservationIn;
uint      RF_timeIndex;
uint      RF_statusIndex;
uint     *RF_yIndex;
uint      RF_ySize;
uint     *RF_testMembershipFlag;  
uint      RF_intrPredictorSize;
uint     *RF_intrPredictor;
uint     *RF_intrIndividual;
char     *RF_importanceFlag;   
uint      RF_eventTypeSize;
uint      RF_mStatusSize; 
uint     *RF_eventType;
uint     *RF_eventTypeIndex;
uint     *RF_eIndividualSize;
uint    **RF_eIndividualIn;
uint      *RF_classLevelSize;
uint     **RF_classLevel;
uint     **RF_classLevelIndex;
uint    ***RF_cIndividualIn;
double   *RF_timeInterest;
uint      RF_timeInterestSize;
uint      RF_sortedTimeInterestSize;
double   *RF_masterTime;
uint      RF_masterTimeSize;
uint     *RF_masterTimeIndexIn;
uint      RF_rFactorCount;
uint     *RF_rFactorMap;
uint     *RF_rFactorIndex;
uint     *RF_rFactorSize;
uint      RF_mrFactorSize;
uint      RF_fmrFactorSize;
uint     *RF_mrFactorIndex;
uint     *RF_fmrFactorIndex;
uint      RF_xFactorCount;
uint     *RF_xFactorMap;
uint     *RF_xFactorIndex;
uint     *RF_xFactorSize;
uint      RF_mxFactorSize;
uint      RF_fmxFactorSize;
uint     *RF_mxFactorIndex;
uint     *RF_fmxFactorIndex;
uint      RF_rMaxFactorLevel;
uint      RF_xMaxFactorLevel;
uint      RF_maxFactorLevel;
char      RF_mStatusFlag; 
char      RF_mTimeFlag; 
char      RF_mResponseFlag; 
char      RF_mPredictorFlag; 
char      RF_fmStatusFlag; 
char      RF_fmTimeFlag;
char      RF_fmResponseFlag; 
char      RF_fmPredictorFlag;
uint     *RF_mRecordMap;
uint     *RF_fmRecordMap;
uint      RF_mRecordSize;
uint      RF_fmRecordSize;
uint     *RF_mRecordIndex;
uint     *RF_fmRecordIndex;
uint      RF_mvSize;
uint      RF_fmvSize;
int     **RF_mvSign;
int     **RF_fmvSign;
int      *RF_mvIndex;
int      *RF_fmvIndex;
double   **RF_importancePtr;
double **RF_sImputeResponsePtr;
double **RF_sImputePredictorPtr;
double **RF_sOOBImputeResponsePtr;
double **RF_sOOBImputePredictorPtr;
uint  **RF_terminalNodeMembershipPtr;
uint  **RF_bootstrapMembershipPtr;
double ***RF_oobEnsemblePtr;
double ***RF_fullEnsemblePtr;
double ***RF_oobCIFPtr;
double ***RF_fullCIFPtr;
double  **RF_oobSRVPtr;
double  **RF_fullSRVPtr;
double  **RF_oobMRTPtr;
double  **RF_fullMRTPtr;
uint     *RF_oobEnsembleDen;
uint     *RF_fullEnsembleDen;
uint     **RF_vimpEnsembleDen;
double   **RF_vimpOutcome;
double  ***RF_sVimpOutcome;
double  ***RF_cVimpEnsemble;
double  ***RF_vimpLeo;
double   **RF_perfLeo;
double ***RF_splitDepthPtr;
uint    *RF_serialTreeIndex;  
uint     RF_serialTreeCount;  
char    **RF_dmRecordBootFlag;
double ***RF_dmvImputation;
Terminal ***RF_mTerminalInfo;
double **RF_performancePtr;
uint   **RF_varUsedPtr;
uint    *RF_oobSize;
uint    *RF_foobSize;
uint    *RF_leafCount;
uint    *RF_nodeCount;
uint    *RF_mwcpCount;
uint    *RF_maxDepth;
double  **RF_status;
double  **RF_time;
double ***RF_response;
double  **RF_ftime;
double  **RF_fstatus;
double ***RF_fresponse;
double ***RF_observation;
double ***RF_fobservation;
Node    **RF_root;
Node   ***RF_nodeMembership;
uint    **RF_bootMembershipIndex;
uint     *RF_trivialBootMembershipIndex;
uint    **RF_bootMembershipFlag;
uint    **RF_oobMembershipFlag;
Node   ***RF_terminalNode;
Node   ***RF_fnodeMembership;
uint    **RF_masterTimeIndex;
Factor ***RF_factorList;
float (*ran1) (uint);
void  (*randomSetChain) (uint, int);
int   (*randomGetChain) (uint);
float (*ran2) (uint);
void  (*randomSetUChain) (uint, int);
int   (*randomGetUChain) (uint);
float (*ran3) (uint);
void  (*randomSetUChainCov) (uint, int);
int   (*randomGetUChainCov) (uint);
SEXP rfsrc(char mode, int seedValue, uint traceFlag) {
  uint sexpIndex;
  uint **mwcpPtrPtr;
  uint  *mwcpPtr;
  uint   totalMWCPCount, rejectedTreeCount;  
  uint i, j, k, r;
  int vimpCount, intrIndex, b, p;
  if (RF_imputeSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number imputations must be greater than zero:  %10d \n", RF_forestSize);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_forestSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of bootstrap iterations must be greater than zero:  %10d \n", RF_forestSize);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_observationSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of individuals must be greater than one:  %10d \n", RF_observationSize);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
  if (RF_xSize < 1) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Parameter verification failed.");
    Rprintf("\nRF-SRC:  Number of parameters must be greater than zero:  %10d \n", RF_xSize);
    Rprintf("\nRF-SRC:  The application will now exit.\n");
    return R_NilValue;
  }
#ifdef SUPPORT_OPENMP
  if (RF_numThreads < 0) {
    RF_numThreads = (omp_get_max_threads() > 2) ? (omp_get_max_threads() - 1) : (omp_get_max_threads());
  }
  else {
    RF_numThreads = (RF_numThreads < omp_get_max_threads()) ? (RF_numThreads) : (omp_get_max_threads());
  }
#endif
  stackIncomingArrays(mode);
  stackPreDefinedCommonArrays();
  switch (mode) {
  case RF_PRED:
    stackPreDefinedPredictArrays();
    break;
  case RF_REST:
    stackPreDefinedRestoreArrays();
    break;
  default:
    stackPreDefinedGrowthArrays();
    break;
  }
  initializeArrays(mode);
  stackFactorArrays(mode);
  stackMissingArrays(mode);
  if (RF_statusIndex > 0) {
    stackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    stackClassificationArrays(mode);
  }
  sexpIndex = stackDefinedOutputObjects(mode,
                                        sexpString,
                                        & RF_root,
                                        & RF_oobEnsemble_,
                                        & RF_fullEnsemble_,
                                        & RF_performance_,
                                        & RF_leafCount_,
                                        & RF_proximity_,
                                        & RF_importance_,
                                        & RF_seed_,
                                        & RF_oobImputation_,
                                        & RF_imputation_,
                                        & RF_sImputeResponsePtr,
                                        & RF_sImputePredictorPtr,
                                        & RF_sOOBImputeResponsePtr,
                                        & RF_sOOBImputePredictorPtr,
                                        & RF_varUsed_,
                                        & RF_varUsedPtr,
                                        & RF_splitDepth_,
                                        & RF_oobEnsembleCIF_,
                                        & RF_fullEnsembleCIF_,
                                        & RF_oobEnsembleSRV_,
                                        & RF_fullEnsembleSRV_,
                                        & RF_oobEnsembleMRT_,
                                        & RF_fullEnsembleMRT_,
                                        & RF_terminalNodeMembership_,
                                        & RF_bootstrapMembership_,
                                        & stackCount,
                                        sexpVector
                                        );
#ifdef SUPPORT_OPENMP
  ran1 = &randomChainParallel;
  ran2 = &randomUChainParallel;
  ran3 = &randomUChainParallelCov;
  randomSetChain = &randomSetChainParallel;
  randomSetUChain = &randomSetUChainParallel;
  randomSetUChainCov = &randomSetUChainParallelCov;
  randomGetChain = &randomGetChainParallel;
  randomGetUChain = &randomGetUChainParallel;
  randomGetUChainCov = &randomGetUChainParallelCov;
  randomStack(RF_forestSize, RF_xSize);
  if (mode == RF_GROW) {
    randomSetChain(1, seedValue);
    for (b = 1; b <= RF_forestSize; b++) {
      randomSetChain(b, -abs(randomGetChain(1) * 113));
      ran1(1);
    }
    randomSetChain(1, -abs(randomGetChain(1) * 113));
  }
  randomSetUChain(1, seedValue);
  for (b = 1; b <= RF_forestSize; b++) {
    randomSetUChain(b, -abs(randomGetUChain(1) * 251));
    ran2(1);
    ran2(1);
  }
  randomSetUChain(1, -abs(randomGetUChain(1) * 251));
  randomSetUChainCov(1, seedValue);
  for (p = 1; p <= RF_xSize; p++) {
    randomSetUChainCov(p, -abs(randomGetUChainCov(1) * 397));
    ran3(1);
    ran3(1);
  }
  randomSetUChainCov(1, -abs(randomGetUChainCov(1) * 397));
#else
  ran1 = &randomChainSerial;
  ran2 = &randomUChainSerial;
  ran3 = &randomUChainSerialCov;
  randomSetChain = &randomSetChainSerial;
  randomSetUChain = &randomSetUChainSerial;
  randomSetUChainCov = &randomSetUChainSerialCov;
  randomGetChain = &randomGetChainSerial;
  randomGetUChain = &randomGetUChainSerial;
  randomGetUChainCov = &randomGetUChainSerialCov;
  randomStack(1, 1);
  if (mode == RF_GROW) {
    randomSetChain(1, seedValue);
    ran1(1);
    randomSetChain(1, -abs(randomGetChain(1)));
  }
  randomSetUChain(1, seedValue);
  ran2(1);
  ran2(1);
  randomSetUChain(1, -abs(randomGetUChain(1)) * 251);
  randomSetUChainCov(1, seedValue);
  ran3(1);
  ran3(1);
  randomSetUChainCov(1, -abs(randomGetUChainCov(1)) * 397);
#endif
  if (mode != RF_GROW) {
    for (b = 1; b <= RF_forestSize; b++) {
      if(RF_seed_[b] >= 0) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Parameter verification failed.");
        Rprintf("\nRF-SRC:  Forest random seed element must be less than zero:  %10d \n", RF_seed_[b]);
        Rprintf("\nRF-SRC:  The application will now exit.\n");
        return R_NilValue;
      }
    }
  }
  if (mode != RF_GROW) {
    for (i = 1; i <= RF_totalNodeCount; i++) {
      if ((RF_treeID_[i] > 0) && (RF_treeID_[i] <= RF_forestSize)) {
        RF_nodeCount[RF_treeID_[i]] ++;
        RF_mwcpCount[RF_treeID_[i]] += RF_mwcpSZ_[i];
      }
      else {
        Rprintf("\nDiagnostic Trace of Tree Record:  \n");
        Rprintf("\n    treeID     nodeID     parmID       spltPT     mwcpSZ \n");
        Rprintf("%10d %10d %10d %12.4f %10d \n", RF_treeID_[i], RF_nodeID_[i], RF_parmID_[i], RF_contPT_[i], RF_mwcpSZ_[i]);
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Invalid forest input record at line:  %10d", i);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        Rprintf("\nRF-SRC:  The application will now exit.\n");
        return R_NilValue;
      }
    }
    for (b = 1; b <= RF_forestSize; b++) {
      RF_leafCount[b] = (RF_nodeCount[b] + 1) >> 1;
    }
  }
  for (r = 1; r <= RF_imputeSize; r++) {
    if (r == RF_imputeSize) {
#ifdef SUPPORT_OPENMP
      if (mode == RF_GROW) {
        if (RF_opt & OPT_TREE) { 
          for (b = 1; b <= RF_forestSize; b++) {
            randomSetChain(b , -abs(randomGetChain(b)));
            RF_seed_[b] = randomGetChain(b);
          }
        }
      }
      else {
        for (b = 1; b <= RF_forestSize; b++) {
          randomSetChain(b , RF_seed_[b]);
        }
      }
#else
      if (mode == RF_GROW) {
        if (RF_opt & OPT_TREE) { 
          randomSetChain(1, -abs(randomGetChain(1)));
          RF_seed_[1] = randomGetChain(1);
        }
      }
      else {
        randomSetChain(1 , RF_seed_[1]);
      }
#endif
    }  
    for(k = 1; k <= RF_forestSize; k++) {
      RF_serialTreeIndex[k] = 0;
    }
    RF_serialTreeCount = 0;
    if (RF_numThreads > 0) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (b = 1; b <= RF_forestSize; b++) {
        acquireTree(mode, r, b);
      }
    }
    else {
      for (b = 1; b <= RF_forestSize; b++) {
        acquireTree(mode, r, b);
      }
    }
    if (r == RF_imputeSize) {
      if (RF_opt & OPT_PROX) {
        for (b = 1; b <= RF_forestSize; b++) {
          updateProximity(mode, b);
        }
      }
      if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
        for (b = 1; b <= RF_forestSize; b++) {
          updateSplitDepth(b, RF_root[b], RF_maxDepth[b]);
        }
      }
      if (RF_opt & OPT_VUSE) {
        for (b = 1; b <= RF_forestSize; b++) {
          getVariablesUsed(b, RF_root[b], RF_varUsedPtr[b]);
        }
      }
    }
    if (r == RF_imputeSize) {
      if (!(RF_opt & OPT_IMPU_ONLY)) {
        char multipleImputeFlag;
        multipleImputeFlag = FALSE;
        if (mode == RF_GROW) {
          if (r > 1) {
            multipleImputeFlag = TRUE;
          } 
        }
        if (RF_numThreads > 0) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
          for (b = 1; b <= RF_forestSize; b++) {
            updateEnsembleCalculations(multipleImputeFlag, mode, b);
          }
        }
        else {
          for (b = 1; b <= RF_forestSize; b++) {
            updateEnsembleCalculations(multipleImputeFlag, mode, b);
          }
        }
      }
    }  
    else {
      for (b = 1; b <= RF_forestSize; b++) {
        freeTree(b, RF_root[b], TRUE);
        free_vvector(RF_nodeMembership[b], 1, RF_observationSize);
        free_uivector(RF_bootMembershipIndex[b], 1, RF_observationSize);
        free_uivector(RF_bootMembershipFlag[b], 1, RF_observationSize);
        free_uivector(RF_oobMembershipFlag[b], 1, RF_observationSize);
        free_vvector(RF_terminalNode[b], 1, RF_leafCount[b] + 1);
        if (mode == RF_PRED) {
          free_vvector((Node **) RF_fnodeMembership[b],  1, RF_fobservationSize);
        }
        unstackShadow(mode, b);
      }
    }
    if (mode != RF_PRED) {
      if (r == 1) {
        if (!(RF_opt & (OPT_BOOT_NODE | OPT_BOOT_NONE))) {
          if (RF_opt & OPT_OMIS) {
            imputeSummary(RF_GROW, FALSE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(FALSE);
              }
            }
          }
        }
        else {
          if (RF_opt & OPT_MISS) {
            imputeSummary(RF_GROW, ACTIVE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(ACTIVE);
              }
            }
          }
        }
      }  
      else {
        if (r < RF_imputeSize) {
          if (RF_opt & OPT_MISS) {
            imputeSummary(RF_GROW, ACTIVE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(ACTIVE);
              }
            }
          }
        }
      }
      if (r == RF_imputeSize) {
        if (RF_imputeSize == 1) {
          if (RF_opt & OPT_MISS) {
            imputeSummary(RF_GROW, TRUE);
            if (RF_timeIndex > 0) {
              if (RF_mTimeFlag == TRUE) {
                imputeMultipleTime(TRUE);
              }
            }
          }
        }  
        else {
        }  
      }  
    }  
    else {
      if (r == RF_imputeSize) {
        if (RF_opt & OPT_MISS) {
          imputeSummary(RF_PRED, ACTIVE);
          if (RF_timeIndex > 0) {
            if (RF_mTimeFlag == TRUE) {
            }
          }
        }
      }
    }
    if ((RF_opt & OPT_OMIS) | (RF_opt & OPT_MISS)) {
      for (b = 1; b <= RF_forestSize; b++) {
        for (j = 1; j <= RF_leafCount[b]; j++) {
          freeTerminal(RF_mTerminalInfo[b][j]);
        }
        free_vvector(RF_mTerminalInfo[b], 1, RF_leafCount[b] + 1);
      }
    }
  }  
  if (RF_opt & OPT_VIMP) {
    if (RF_opt & OPT_VIMP_JOIN) {
      vimpCount = 1;
    }
    else {
      vimpCount = RF_intrPredictorSize;
    }
    if (RF_numThreads > 0) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
      for (intrIndex = 1; intrIndex <= vimpCount; intrIndex++) {
        for (uint bb = 1; bb <= RF_forestSize; bb++) {
          if (RF_leafCount[bb] > 0) {
            updateVimpCalculations(mode, bb, intrIndex);
          }
        }
      }
    }
    else {
      for (intrIndex = 1; intrIndex <= vimpCount; intrIndex++) {
        for (uint bb = 1; bb <= RF_forestSize; bb++) {
          if (RF_leafCount[bb] > 0) {
            updateVimpCalculations(mode, bb, intrIndex);
          }
        }
      }
    }
  }
  for (b = 1; b <= RF_forestSize; b++) {
    for (j = 1; j <= RF_leafCount[b]; j++) {
      if ((RF_timeIndex > 0) && (RF_statusIndex > 0)) {
        unstackMortality(RF_terminalNode[b][j]);
      }
      else {
        if (RF_rFactorCount > 0) {
          unstackMultiClassProb(RF_terminalNode[b][j]);
        }
      }
    }
    if (!(RF_opt & OPT_TREE)) {
      freeTree(b, RF_root[b], TRUE);
    }
    free_vvector(RF_nodeMembership[b], 1, RF_observationSize);
    free_uivector(RF_bootMembershipIndex[b], 1, RF_observationSize);
    free_uivector(RF_bootMembershipFlag[b], 1, RF_observationSize);
    free_uivector(RF_oobMembershipFlag[b], 1, RF_observationSize);
    free_vvector(RF_terminalNode[b], 1, RF_leafCount[b] + 1);
    if (mode == RF_PRED) {
      free_vvector((Node **) RF_fnodeMembership[b],  1, RF_fobservationSize);
    }
    unstackShadow(mode, b);
  }
  rejectedTreeCount = k = 0;
  for (b = 1; b <= RF_forestSize; b++) {
    if (RF_leafCount[b] == 0) {
      rejectedTreeCount ++;
    }
    if (RF_leafCount[b] == 1) {
      k ++;
    }
  }
  if (rejectedTreeCount < RF_forestSize) {
    if (RF_opt & OPT_VIMP) {
      if (RF_opt & OPT_VIMP_JOIN) {
        vimpCount = 1;
      }
      else {
        vimpCount = RF_intrPredictorSize;
      }
      if (RF_opt & OPT_VIMP_LEOB) {
      }
      else {
        if (RF_numThreads > 0) {
#ifdef SUPPORT_OPENMP
#pragma omp parallel for num_threads(RF_numThreads)
#endif
          for (p = 1; p <= vimpCount; p++) {
            summarizeVimpPerformance(mode, 0, p);
          }
        }
        else {
          for (p = 1; p <= vimpCount; p++) {
            summarizeVimpPerformance(mode, 0, p);
          }
        }
      }
      finalizeVimpPerformance(mode, rejectedTreeCount);
    }
    finalizeEnsembleEstimates(mode, rejectedTreeCount);
    if (RF_opt & OPT_VUSE) {
      if (RF_opt & (~OPT_VUSE) & OPT_VUSE_TYPE) {
      }
      else {
        for (j = 1; j <= RF_xSize; j++) {
          RF_varUsed_[j] = 0;
          for (i = 1; i <= RF_forestSize; i++) {
            RF_varUsed_[j] += RF_varUsedPtr[i][j];
          }
        }
      }
    }
    if (RF_opt & (OPT_SPLDPTH_F | OPT_SPLDPTH_T)) {
      if (RF_opt & OPT_SPLDPTH_F) {
        for (j = 1; j <= RF_xSize; j++) {
          for (i = 1; i <= RF_observationSize; i++) {
            RF_splitDepthPtr[1][j][i] = RF_splitDepthPtr[1][j][i] / (RF_forestSize - rejectedTreeCount);
          }
        }
      }
      else {
      }
    }
  }  
  else {
    Rprintf("\nRF-SRC:  *** WARNING *** ");
    Rprintf("\nRF-SRC:  Insufficient trees for analysis.  \n");
  }
  RF_totalNodeCount = totalMWCPCount = 0;
  if (mode == RF_GROW) {
    if (RF_opt & OPT_TREE) {
      for (b = 1; b <= RF_forestSize; b++) {
        if (RF_leafCount[b] > 0) {
          totalMWCPCount += RF_mwcpCount[b];
          RF_totalNodeCount += (2 * RF_leafCount[b]) - 1;
        }
        else {
          RF_totalNodeCount ++;
        }
      }
    }
  }
  sexpIndex = 
    stackVariableOutputObjects(RF_totalNodeCount,    
                               totalMWCPCount,       
                               & RF_treeID_,         
                               & RF_nodeID_,         
                               & RF_parmID_,         
                               & RF_contPT_,         
                               & RF_mwcpSZ_,         
                               & RF_mwcpPT_,         
                               sexpIndex, 
                               sexpString,
                               sexpVector);
  if (mode == RF_GROW) {
    if (RF_opt & OPT_TREE) {
      mwcpPtr = RF_mwcpPT_;
      mwcpPtrPtr = & mwcpPtr;
      RF_totalNodeCount = 1;
      for (b = 1; b <= RF_forestSize; b++) {
        saveTree(b, 
                 RF_root[b], 
                 & RF_totalNodeCount, 
                 RF_treeID_, 
                 RF_nodeID_, 
                 RF_parmID_, 
                 RF_contPT_,
                 RF_mwcpSZ_,
                 mwcpPtrPtr);
      }
      RF_totalNodeCount --;
    }  
  }  
  unstackDefinedOutputObjects(mode, RF_root);
  if (RF_statusIndex > 0) {
    unstackCompetingArrays(mode);
  }
  if (RF_rFactorCount > 0) {
    unstackClassificationArrays(mode);
  }
  unstackMissingArrays(mode);
  unstackFactorArrays();
  switch (mode) {
  case RF_PRED:
    unstackPreDefinedPredictArrays();
    break;
  case RF_REST:
    unstackPreDefinedRestoreArrays();
    break;
  default:
    unstackPreDefinedGrowthArrays();
    break;
  }
  unstackPreDefinedCommonArrays();
  unstackIncomingArrays(mode);
#ifdef SUPPORT_OPENMP
  randomUnstack(RF_forestSize, RF_xSize);
#else
  randomUnstack(1, 1);
#endif
  UNPROTECT(stackCount + 2);
  return sexpVector[RF_OUTP_ID];
}
