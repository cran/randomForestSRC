////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.6.0
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


#ifndef RSFSTACKOUTPUT_H
#define RSFSTACKOUTPUT_H
#include "node.h"
uint stackDefinedOutputObjects(char      mode,
                               char    **sexpString,
                               Node   ***pRF_root,
                               double  **pRF_oobEnsemble,
                               double  **pRF_fullEnsemble,
                               double  **p_performance,
                               uint    **p_leafCount,
                               double  **pRF_proximity,
                               double  **pRF_weight,
                               double  **p_varImportance,
                               int     **pRF_seed,
                               double  **p_imputation,
                               double ***p_sImputeResponsePtr,
                               double ***p_sImputePredictorPtr,
                               uint    **pRF_varUsed,
                               uint   ***pRF_varUsedPtr,
                               double  **p_splitDepth,
                               double  **pRF_oobPOEEnsemble,
                               double  **pRF_fullPOEEnsemble,
                               double  **pRF_oobEnsembleSRV,
                               double  **pRF_fullEnsembleSRV,
                               double  **pRF_oobEnsembleMRT,
                               double  **pRF_fullEnsembleMRT,
                               uint    **pRF_tTermMembershipIndex,
                               uint    **pRF_pNodeMembershipIndex,
                               uint    **pRF_bootstrapMembership,
                               uint     *stackCount,
                               SEXP     *sexpVector);
void unstackDefinedOutputObjects(char      mode,
                                 Node    **root);
uint stackVariableOutputObjects(char     mode,
                                uint     totalNodeCount,
                                uint     totalMWCPCount,
                                uint   **p_treeID,
                                uint   **pRF_nodeID,
                                uint   **pRF_parmID,
                                double **pRF_contPT,
                                uint   **pRF_mwcpSZ,
                                uint   **pRF_mwcpPT,
                                double **pRF_spltST,
                                double **pRF_spltVR,
                                uint     sexpLength,
                                char   **sexpString,
                                SEXP    *sexpVector);
uint stackVariableTerminalNodeOutputObjects(char     mode,
                                            double **pRF_TN_SURV,
                                            double **pRF_TN_MORT,
                                            double **pRF_TN_NLSN,
                                            double **pRF_TN_CSHZ,
                                            double **pRF_TN_CIFN,
                                            double **pRF_TN_REGR,
                                            uint   **pRF_TN_CLAS,
                                            uint   **pRF_TN_MCNT,
                                            uint   **pRF_TN_MEMB,
                                            uint     sexpIndex,
                                            char   **sexpString,
                                            SEXP    *sexpVector);
void stackAuxVariableTerminalNodeOutputObjects(char    mode,
                                               double *pRF_TN_SURV,
                                               double *pRF_TN_MORT,
                                               double *pRF_TN_NLSN,
                                               double *pRF_TN_CSHZ,
                                               double *pRF_TN_CIFN,
                                               double *pRF_TN_REGR,
                                               uint   *pRF_TN_CLAS,
                                               uint   *pRF_TN_MCNT,
                                               uint   *pRF_TN_MEMB);
void unstackAuxVariableTerminalNodeOutputObjects(char mode);
void saveTerminalNodeOutputObjects (uint mode);
#endif
