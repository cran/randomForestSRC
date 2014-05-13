////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.0
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


#include        "global.h"
#include        "extern.h"
#include         "trace.h"
#include        "nrutil.h"
#include        "impute.h"
#include      "survival.h"
#include      "survivalE.h"
void updateEnsembleSurvival(uint mode, uint treeID) {
  uint obsSize;
  unsigned char oobFlag, fullFlag, selectionFlag;
  double ***ensemblePtr;
  double  **ensSRVPtr;
  double ***ensCIFPtr;
  double  **ensMRTPtr;
  Node   ***nodeMembershipPtr;
  uint     *ensembleDenPtr;
  uint i, j, k;
  Node *parent;
  ensemblePtr       = NULL;  
  ensembleDenPtr    = NULL;  
  ensSRVPtr         = NULL;  
  ensCIFPtr         = NULL;  
  ensMRTPtr         = NULL;  
  oobFlag = fullFlag = FALSE;
  switch (mode) {
  case RF_PRED:
    obsSize = RF_fobservationSize;
    oobFlag = FALSE;
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    nodeMembershipPtr = RF_ftNodeMembership;
    break;
  default:
    obsSize = RF_observationSize;
    if (RF_opt & OPT_OENS) {
      if (RF_oobSize[treeID] > 0) {
        oobFlag = TRUE;
      }
    }
    if (RF_opt & OPT_FENS) {
      fullFlag = TRUE;
    }
    nodeMembershipPtr = RF_tNodeMembership;
    break;
  }
  while ((oobFlag == TRUE) || (fullFlag == TRUE)) {
    if (oobFlag == TRUE) {
      ensemblePtr    = RF_oobEnsemblePtr;
      ensembleDenPtr = RF_oobEnsembleDen;
      ensMRTPtr = RF_oobMRTPtr;
      if (!(RF_opt & OPT_COMP_RISK)) {
        ensSRVPtr = RF_oobSRVPtr;
      }
      else {
        ensCIFPtr = RF_oobCIFPtr;
      }
    }
    else {
      if (fullFlag == TRUE) {
        ensemblePtr    = RF_fullEnsemblePtr;
        ensembleDenPtr = RF_fullEnsembleDen;
        ensMRTPtr = RF_fullMRTPtr;
      if (!(RF_opt & OPT_COMP_RISK)) {
          ensSRVPtr = RF_fullSRVPtr;
        }
        else {
          ensCIFPtr = RF_fullCIFPtr;
        }
      }
    }
    for (i=1; i <= obsSize; i++) {
      selectionFlag = TRUE;
      if (oobFlag == TRUE) {
        if (RF_bootMembershipFlag[treeID][i] == FALSE) {
          selectionFlag = TRUE;
        }
        else {
          selectionFlag = FALSE;
        }
      }
      if (selectionFlag) {
        parent = nodeMembershipPtr[treeID][i];
        if (RF_opt & OPT_OUTC_TYPE) {
          if (!ISNA(parent -> predictedOutcome)) {
          }
          else {
            selectionFlag = FALSE;
          }
        }
      }
      if (selectionFlag) {
        ensembleDenPtr[i] ++;
        if (!(RF_opt & OPT_COMP_RISK)) {
          ensMRTPtr[1][i] += parent -> mortality[1];
          for (k=1; k <= RF_sortedTimeInterestSize; k++) {
            ensemblePtr[1][k][i] += parent -> nelsonAalen[k];
            ensSRVPtr[k][i] += parent -> survival[k];
          }
        }
        else {
          for (j = 1; j <= RF_eventTypeSize; j++) {
            ensMRTPtr[j][i] += parent -> mortality[j];
            for (k=1; k <= RF_sortedTimeInterestSize; k++) {
              ensemblePtr[j][k][i] += parent -> CSH[j][k];
              ensCIFPtr[j][k][i] += parent -> CIF[j][k];
            }
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
void getEnsembleMortalityCR(uint      mode,
                            uint      treeID,
                            uint      obsSize,
                            double  **ensemblePtr,
                            uint     *ensembleDenPtr,
                            double  **crMortality) {
  uint i, j;
  for (j = 1; j <= RF_eventTypeSize; j ++) {
    for (i = 1; i <= obsSize; i++) {
      if (ensembleDenPtr[i] != 0) {
        crMortality[j][i] = ensemblePtr[j][i] / ensembleDenPtr[i];
      }
      else {
        crMortality[j][i] = NA_REAL;
      }
    }
  }
}
void getEnsembleMortality(uint      mode,
                          uint      treeID,
                          uint      obsSize,
                          double  **ensemblePtr,
                          uint     *ensembleDenPtr,
                          double   *mortality) {
  uint i;
  for (i = 1; i <= obsSize; i++) {
    if (ensembleDenPtr[i] != 0) {
      mortality[i] = ensemblePtr[1][i] / ensembleDenPtr[i];
    }
    else {
      mortality[i] = NA_REAL;
    }
  }
}
void getConditionalConcordanceArrays(uint     j,
                                     double  *timePtr,
                                     double  *statusPtr,
                                     double  *mortalityPtr,
                                     uint    *genericEnsembleDenPtr,
                                     uint    *meIndividualSize,
                                     uint   **eIndividual,
                                     double  *subsettedTime,
                                     double  *subsettedStatus,
                                     double  *subsettedMortality,
                                     uint    *subsettedEnsembleDen) {
  uint i;
  if (!(RF_opt & OPT_COMP_RISK)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt to update event type subsets in a non-CR analysis.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  for (i = 1; i <= meIndividualSize[j]; i++) {
    subsettedTime[i]        = timePtr[eIndividual[j][i]];
    subsettedStatus[i]      = statusPtr[eIndividual[j][i]];
    subsettedMortality[i]   = mortalityPtr[eIndividual[j][i]];
    subsettedEnsembleDen[i] = genericEnsembleDenPtr[eIndividual[j][i]];
  }
}
double getConcordanceIndex(int     polarity,
                           uint    size,
                           double *timePtr,
                           double *statusPtr,
                           double *predictedOutcome,
                           uint   *denCount) {
  uint i,j;
  ulong concordancePairSize;
  ulong concordanceWorseCount;
  double result;
  concordancePairSize = concordanceWorseCount = 0;
  for (i=1; i < size; i++) {
    for (j=i+1; j <= size; j++) {
      if (denCount[i] != 0  && denCount[j] != 0) {
        if ( ((timePtr[i] - timePtr[j] > EPSILON) && (statusPtr[j] > 0)) ||
             ((fabs(timePtr[i] - timePtr[j]) <= EPSILON) && (statusPtr[j] > 0) && (statusPtr[i] == 0)) ) {
          concordancePairSize += 2;
          if (predictedOutcome[j] - predictedOutcome[i] > EPSILON) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[j] - predictedOutcome[i]) < EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( ((timePtr[j] - timePtr[i]) > EPSILON  && (statusPtr[i] > 0)) ||
                  ((fabs(timePtr[j] - timePtr[i]) <= EPSILON)  && (statusPtr[i] > 0) && (statusPtr[j] == 0)) ) {
          concordancePairSize += 2;
          if ( predictedOutcome[i] - predictedOutcome[j] > EPSILON ) {
            concordanceWorseCount += 2;
          }
          else if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 1;
          }
        }
        else if ( (fabs(timePtr[i]- timePtr[j]) <= EPSILON) && (statusPtr[i] > 0) && (statusPtr[j] > 0) ) {
          concordancePairSize += 2;
          if (fabs(predictedOutcome[i] - predictedOutcome[j]) < EPSILON) {
            concordanceWorseCount += 2;
          }
          else {
            concordanceWorseCount += 1;
          }
        }
      }  
    }  
  }  
  if (concordancePairSize == 0) {
    result = NA_REAL;
  }
  else {
    result = 1.0 - ((double) concordanceWorseCount / (double) concordancePairSize);
  }
  return result;
}
void getCRPerformance (uint     mode,
                       uint     obsSize,
                       double **responsePtr,
                       double **conditionalMortality,
                       uint    *ensembleDenPtr,
                       double  *performanceVector) {
  uint   mRecordSize;
  int  **mpSign;
  uint  *mRecordIndex;
  uint  *meIndividualSize;
  uint **eIndividual;
  double concordanceIndex;
  uint j;
  if (!(RF_opt & OPT_COMP_RISK)) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Attempt at conditional performance updates in a non-CR analysis.");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  if (RF_mStatusSize > 0) {
    if (mode != RF_PRED) {
      mRecordSize = RF_mRecordSize;
      mpSign = RF_mpSign;
      mRecordIndex = RF_mRecordIndex;
    }
    else {
      mRecordSize = RF_fmRecordSize;
      mpSign = RF_fmpSign;
      mRecordIndex = RF_fmRecordIndex;
    }
    meIndividualSize  = uivector(1, RF_eventTypeSize);
    eIndividual = (uint **) vvector(1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      eIndividual[j] = uivector(1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    updateEventTypeSubsets(responsePtr[RF_statusIndex], mRecordSize, mpSign, mRecordIndex, meIndividualSize, eIndividual);
  }
  else {
    meIndividualSize  = RF_eIndividualSize;
    eIndividual = RF_eIndividualIn;
  }
  double *subsettedTime      = dvector(1, obsSize);
  double *subsettedStatus    = dvector(1, obsSize);
  double *subsettedMortality = dvector(1, obsSize);
  uint *subsettedEnsembleDen = uivector(1, obsSize);
  for (j = 1; j <= RF_eventTypeSize; j++) {
    getConditionalConcordanceArrays(j,
                                    responsePtr[RF_timeIndex],
                                    responsePtr[RF_statusIndex],
                                    conditionalMortality[j],
                                    ensembleDenPtr,
                                    meIndividualSize,
                                    eIndividual,
                                    subsettedTime,
                                    subsettedStatus,
                                    subsettedMortality,
                                    subsettedEnsembleDen);
    concordanceIndex = getConcordanceIndex(1,
                                           meIndividualSize[j],
                                           subsettedTime,
                                           subsettedStatus,
                                           subsettedMortality,
                                           subsettedEnsembleDen);
    if (ISNA(concordanceIndex)) {
      performanceVector[j] = NA_REAL;
    }
    else {
      performanceVector[j] = concordanceIndex;
    }
  }
  if (RF_mStatusSize > 0) {
    free_uivector(meIndividualSize, 1, RF_eventTypeSize);
    for (j = 1; j <= RF_eventTypeSize; j++) {
      free_uivector(eIndividual[j], 1, RF_eIndividualSize[j] + RF_mStatusSize + 1);
    }
    free_vvector(eIndividual, 1, RF_eventTypeSize);
  }
  free_dvector(subsettedTime, 1, obsSize);
  free_dvector(subsettedStatus, 1, obsSize);
  free_dvector(subsettedMortality, 1, obsSize);
  free_uivector(subsettedEnsembleDen, 1, obsSize);
}
