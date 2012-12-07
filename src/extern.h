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


#ifndef EXTERNAL_H
#define EXTERNAL_H
#include          "node.h"
#include      "terminal.h"
#include        "factor.h"
extern uint     *RF_treeID_;
extern uint     *RF_nodeID_;
extern uint     *RF_parmID_;
extern uint     *RF_mwcpSZ_;
extern double   *RF_contPT_;
extern uint     *RF_mwcpPT_;
extern uint      RF_totalNodeCount;
extern int      *RF_seed_;
extern uint     *RF_leafCount_;
extern uint     *RF_proximity_;
extern uint      RF_opt;
extern uint      RF_splitRule;
extern uint      RF_splitRandomRule;
extern uint      RF_imputeSize;
extern uint      RF_forestSize;
extern uint      RF_minimumNodeSize;
extern double   *RF_crWeight;
extern uint      RF_randomCovariateCount;
extern double   *RF_randomCovariateWeight;
extern int       RF_numThreads;
extern uint      RF_observationSize;
extern uint      RF_rSize;
extern double   *RF_rData;
extern uint      RF_xSize;
extern double   *RF_xData;
extern double  **RF_responseIn;
extern double  **RF_observationIn;
extern SEXP      RF_sexp_xType;
extern char    **RF_xType;
extern int      *RF_xLevels;
extern SEXP      RF_sexp_rType;
extern char    **RF_rType;
extern int      *RF_rLevels;
extern uint      RF_fobservationSize;
extern uint      RF_frSize;
extern double   *RF_frData;
extern double   *RF_fxData;
extern double  **RF_fresponseIn;
extern double  **RF_fobservationIn;
extern uint      RF_timeIndex;
extern uint      RF_statusIndex;
extern uint     *RF_yIndex;
extern uint      RF_ySize;
extern uint     *RF_testMembershipFlag;  
extern uint      RF_intrPredictorSize;
extern uint     *RF_intrPredictor;
extern uint     *RF_intrIndividual;
extern char     *RF_importanceFlag;   
extern uint      RF_eventTypeSize;
extern uint      RF_mStatusSize; 
extern uint     *RF_eventType;
extern uint     *RF_eventTypeIndex;
extern uint     *RF_eIndividualSize;
extern uint    **RF_eIndividualIn;
extern uint      *RF_classLevelSize;
extern uint     **RF_classLevel;
extern uint     **RF_classLevelIndex;
extern uint    ***RF_cIndividualIn;
extern double   *RF_timeInterest;
extern uint      RF_timeInterestSize;
extern uint      RF_sortedTimeInterestSize;
extern double   *RF_masterTime;
extern uint      RF_masterTimeSize;
extern uint     *RF_masterTimeIndexIn;
extern uint      RF_rFactorCount;
extern uint     *RF_rFactorMap;
extern uint     *RF_rFactorIndex;
extern uint     *RF_rFactorSize;
extern uint      RF_mrFactorSize;
extern uint      RF_fmrFactorSize;
extern uint     *RF_mrFactorIndex;
extern uint     *RF_fmrFactorIndex;
extern uint      RF_xFactorCount;
extern uint     *RF_xFactorMap;
extern uint     *RF_xFactorIndex;
extern uint     *RF_xFactorSize;
extern uint      RF_mxFactorSize;
extern uint      RF_fmxFactorSize;
extern uint     *RF_mxFactorIndex;
extern uint     *RF_fmxFactorIndex;
extern uint      RF_rMaxFactorLevel;
extern uint      RF_xMaxFactorLevel;
extern uint      RF_maxFactorLevel;
extern char      RF_mStatusFlag; 
extern char      RF_mTimeFlag; 
extern char      RF_mResponseFlag; 
extern char      RF_mPredictorFlag; 
extern char      RF_fmStatusFlag; 
extern char      RF_fmTimeFlag;
extern char      RF_fmResponseFlag; 
extern char      RF_fmPredictorFlag;
extern uint     *RF_mRecordMap;
extern uint     *RF_fmRecordMap;
extern uint      RF_mRecordSize;
extern uint      RF_fmRecordSize;
extern uint     *RF_mRecordIndex;
extern uint     *RF_fmRecordIndex;
extern uint      RF_mvSize;
extern uint      RF_fmvSize;
extern int     **RF_mvSign;
extern int     **RF_fmvSign;
extern int      *RF_mvIndex;
extern int      *RF_fmvIndex;
extern double   **RF_importancePtr;
extern double **RF_sImputeResponsePtr;
extern double **RF_sImputePredictorPtr;
extern double **RF_sOOBImputeResponsePtr;
extern double **RF_sOOBImputePredictorPtr;
extern uint  **RF_terminalNodeMembershipPtr;
extern uint  **RF_bootstrapMembershipPtr;
extern double ***RF_oobEnsemblePtr;
extern double ***RF_fullEnsemblePtr;
extern double ***RF_oobCIFPtr;
extern double ***RF_fullCIFPtr;
extern double  **RF_oobSRVPtr;
extern double  **RF_fullSRVPtr;
extern double  **RF_oobMRTPtr;
extern double  **RF_fullMRTPtr;
extern uint     *RF_oobEnsembleDen;
extern uint     *RF_fullEnsembleDen;
extern uint     **RF_vimpEnsembleDen;
extern double ****RF_sVimpEnsemble;
extern double   **RF_vimpOutcome;
extern double  ***RF_cVimpEnsemble;
extern double  ***RF_vimpLeo;
extern double   **RF_perfLeo;
extern double ***RF_splitDepthPtr;
extern uint    *RF_serialTreeIndex;  
extern uint     RF_serialTreeCount;  
extern char    **RF_dmRecordBootFlag;
extern double ***RF_dmvImputation;
extern Terminal ***RF_mTerminalInfo;
extern double **RF_performancePtr;
extern uint   **RF_varUsedPtr;
extern uint    *RF_oobSize;
extern uint    *RF_foobSize;
extern uint    *RF_leafCount;
extern uint    *RF_nodeCount;
extern uint    *RF_mwcpCount;
extern double  **RF_status;
extern double  **RF_time;
extern double ***RF_response;
extern double  **RF_ftime;
extern double  **RF_fstatus;
extern double ***RF_fresponse;
extern double ***RF_observation;
extern double ***RF_fobservation;
extern Node    **RF_root;
extern Node   ***RF_nodeMembership;
extern uint    **RF_bootMembershipIndex;
extern uint     *RF_trivialBootMembershipIndex;
extern uint    **RF_bootMembershipFlag;
extern uint    **RF_oobMembershipFlag;
extern Node   ***RF_terminalNode;
extern Node   ***RF_fnodeMembership;
extern uint    **RF_masterTimeIndex;
extern Factor ***RF_factorList;
extern float (*ran1) (uint);
extern void  (*randomSetChain) (uint, int);
extern int   (*randomGetChain) (uint);
extern float (*ran2) (uint);
extern void  (*randomSetUChain) (uint, int);
extern int   (*randomGetUChain) (uint);
#endif
