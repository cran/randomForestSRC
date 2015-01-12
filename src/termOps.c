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


#include        <stdlib.h>
#include        "nrutil.h"
#include       "termOps.h"
#include   <R_ext/Print.h>
#include   <Rdefines.h>
extern unsigned int getTraceFlag();
extern unsigned int getNodeDefTraceFlag();
extern unsigned int getForkDefTraceFlag();
extern unsigned int getTurnOffTraceFlag();
#include <R_ext/Arith.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
Terminal *makeTerminal() {
  Terminal *parent = (Terminal*) gblock((size_t) sizeof(Terminal));
  parent -> lmiIndex      = NULL;
  parent -> lmiSize       = 0;
  parent -> lmiValue  = NULL;
  parent -> nodeID     = 0;
  parent -> predictedOutcome     = NA_REAL;
  parent -> eTypeSize            = 0;
  parent -> mTimeSize            = 0;
  parent -> eTimeSize            = 0;
  parent -> sTimeSize            = 0;
  parent -> atRiskCount          = NULL;
  parent -> eventCount           = NULL;
  parent -> eventTimeIndex       = NULL;
  parent -> localRatio           = NULL;
  parent -> localCSH             = NULL;
  parent -> localCIF             = NULL;
  parent -> localSurvival        = NULL;
  parent -> localNelsonAalen     = NULL;
  parent -> CSH                  = NULL;
  parent -> CIF                  = NULL;
  parent -> survival             = NULL;
  parent -> nelsonAalen          = NULL;
  parent -> rfCount              = 0;
  parent -> rfSize               = NULL;
  parent -> multiClassProb       = NULL;
  parent -> rnfCount             = 0;
  parent -> meanResponse         = NULL;
  parent -> weight               = 0.0;
  parent -> membrCount           = 0;
  return parent;
}
void freeTerminal(Terminal        *parent) {
  unstackTermLMIIndex(parent);
  freeTerminalNodeSurvivalStructuresNonVimp(parent);
  freeTerminalNodeSurvivalStructuresFinal(parent);
  freeTerminalNodeNonSurvivalStructures(parent);
  free_gblock(parent, sizeof(Terminal));
}
void freeTerminalNodeLocalSurvivalStructures(Terminal *tTerm) {
  unstackLocalRatio(tTerm);
  unstackLocalSurvival(tTerm);
  unstackLocalNelsonAalen(tTerm);
  if (tTerm -> eTypeSize > 1) {
    unstackLocalCSH(tTerm);
    unstackLocalCIF(tTerm);
  }
  unstackEventTimeIndex(tTerm);
}
void freeTerminalNodeSurvivalStructuresNonVimp(Terminal *tTerm) {
  unstackSurvival(tTerm);
  unstackNelsonAalen(tTerm);
  unstackCSH(tTerm);
  unstackCIF(tTerm);
}
void freeTerminalNodeSurvivalStructuresFinal(Terminal *tTerm) {
  unstackMortality(tTerm);
}
void freeTerminalNodeNonSurvivalStructures(Terminal *tTerm) {
  unstackMultiClassProb(tTerm);
  unstackMeanResponse(tTerm);
}
void stackAtRiskAndEventCounts(Terminal *tTerm, unsigned int eTypeSize, unsigned int mTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> mTimeSize > 0) {
    if (tTerm -> mTimeSize != mTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  mTimeSize has been previously defined:  %10d vs %10d", tTerm -> mTimeSize, mTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> mTimeSize = mTimeSize;
  }
  tTerm -> atRiskCount     = uivector(1, mTimeSize);
  tTerm -> eventCount      = uimatrix(1, eTypeSize, 1, mTimeSize);
}
void stackEventTimeIndex(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> eventTimeIndex  = uivector(1, eTimeSize + 1);
}
void unstackAtRiskAndEventCounts(Terminal *tTerm) {
  if (tTerm -> atRiskCount != NULL) {
    free_uivector(tTerm -> atRiskCount, 1, tTerm -> mTimeSize);
    tTerm -> atRiskCount = NULL;
  }
  if (tTerm -> eventCount != NULL) {
    free_uimatrix(tTerm -> eventCount, 1, tTerm -> eTypeSize, 1, tTerm -> mTimeSize);
    tTerm -> eventCount = NULL;
  }
}
void unstackEventTimeIndex(Terminal *tTerm) {
  if (tTerm -> eventTimeIndex != NULL) {
    free_uivector(tTerm -> eventTimeIndex, 1, tTerm -> eTimeSize + 1);
    tTerm -> eventTimeIndex = NULL;
  }
}
void stackLocalRatio(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localRatio = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalRatio(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localRatio != NULL) {
      free_dmatrix(tTerm -> localRatio, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localRatio = NULL;
    }
  }
}
void stackLocalSurvival(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localSurvival = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalSurvival(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localSurvival != NULL) {
      free_dvector(tTerm -> localSurvival, 1, tTerm -> eTimeSize);
      tTerm -> localSurvival = NULL;
    }
  }
}
void stackLocalNelsonAalen(Terminal *tTerm, unsigned int eTimeSize) {
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localNelsonAalen = dvector(1, tTerm -> eTimeSize);
}
void unstackLocalNelsonAalen(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localNelsonAalen != NULL) {
      free_dvector(tTerm -> localNelsonAalen, 1, tTerm -> eTimeSize);
      tTerm -> localNelsonAalen = NULL;
    }
  }
}
void stackLocalCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCSH = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCSH(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCSH != NULL) {
      free_dmatrix(tTerm -> localCSH, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCSH = NULL;
    }
  }
}
void stackLocalCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int eTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> eTimeSize > 0) {
    if (tTerm -> eTimeSize != eTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTimeSize has been previously defined:  %10d vs %10d", tTerm -> eTimeSize, eTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTimeSize = eTimeSize;
  }
  tTerm -> localCIF = dmatrix(1, eTypeSize, 1, tTerm -> eTimeSize);
}
void unstackLocalCIF(Terminal *tTerm) {
  if(tTerm -> eTimeSize > 0) {
    if (tTerm -> localCIF != NULL) {
      free_dmatrix(tTerm -> localCIF, 1, tTerm -> eTypeSize, 1, tTerm -> eTimeSize);
      tTerm -> localCIF = NULL;
    }
  }
}
void stackNelsonAalen(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> nelsonAalen = dvector(1, tTerm -> sTimeSize);
}
void unstackNelsonAalen(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> nelsonAalen != NULL) {
      free_dvector(tTerm -> nelsonAalen, 1, tTerm -> sTimeSize);
      tTerm -> nelsonAalen = NULL;
    }
  }
}
void stackSurvival(Terminal *tTerm, unsigned int sTimeSize) {
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> survival = dvector(1, tTerm -> sTimeSize);
}
void unstackSurvival(Terminal *tTerm) {
  if(tTerm -> sTimeSize > 0) {
    if (tTerm -> survival != NULL) {
      free_dvector(tTerm -> survival, 1, tTerm -> sTimeSize);
      tTerm -> survival = NULL;
    }
  }
}
void stackCSH(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CSH = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCSH(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CSH != NULL) {
        free_dmatrix(tTerm -> CSH, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CSH = NULL;
      }
    }
  }
}
void stackCIF(Terminal *tTerm, unsigned int eTypeSize, unsigned int sTimeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  if (tTerm -> sTimeSize > 0) {
    if (tTerm -> sTimeSize != sTimeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  sTimeSize has been previously defined:  %10d vs %10d", tTerm -> sTimeSize, sTimeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> sTimeSize = sTimeSize;
  }
  tTerm -> CIF = dmatrix(1, eTypeSize, 1, tTerm -> sTimeSize);
}
void unstackCIF(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if(tTerm -> sTimeSize > 0) {
      if (tTerm -> CIF != NULL) {
        free_dmatrix(tTerm -> CIF, 1, tTerm -> eTypeSize, 1, tTerm -> sTimeSize);
        tTerm -> CIF = NULL;
      }
    }
  }
}
void stackMortality(Terminal *tTerm, unsigned int eTypeSize) {
  if (tTerm -> eTypeSize > 0) {
    if (tTerm -> eTypeSize != eTypeSize) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  eTypeSize has been previously defined:  %10d vs %10d", tTerm -> eTypeSize, eTypeSize);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> eTypeSize = eTypeSize;
  }
  tTerm -> mortality = dvector(1, eTypeSize);
}
void unstackMortality(Terminal *tTerm) {
  if(tTerm -> eTypeSize > 0) {
    if (tTerm -> mortality != NULL) {
      free_dvector(tTerm -> mortality, 1, tTerm -> eTypeSize);
      tTerm -> mortality = NULL;
    }
  }
}
void stackMultiClassProb(Terminal *tTerm, unsigned int rfCount, unsigned int *rfSize) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfCount != rfCount) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  rfCount has been previously defined:  %10d vs %10d", tTerm -> rfCount, rfCount);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> rfCount = rfCount;
  }
  tTerm -> rfSize = uivector(1, tTerm -> rfCount);
  tTerm -> multiClassProb = (unsigned int **) new_vvector(1, tTerm -> rfCount, NRUTIL_UPTR);
  for (j = 1; j <= tTerm -> rfCount; j++) {
    (tTerm -> rfSize)[j] = rfSize[j];
    (tTerm -> multiClassProb)[j] = uivector(1, (tTerm -> rfSize)[j]);
  }
}
void unstackMultiClassProb(Terminal *tTerm) {
  unsigned int j;
  if (tTerm -> rfCount > 0) {
    if (tTerm -> rfSize != NULL) {
      if (tTerm -> multiClassProb != NULL) {
        for (j = 1; j <= tTerm -> rfCount; j++) {
          if (tTerm -> multiClassProb[j] != NULL) {
            free_uivector(tTerm -> multiClassProb[j], 1, tTerm -> rfSize[j]);
            tTerm -> multiClassProb[j] = NULL;
          }
        }
        free_new_vvector(tTerm -> multiClassProb, 1, tTerm -> rfCount, NRUTIL_UPTR);
        tTerm -> multiClassProb = NULL;
      }
      free_uivector(tTerm -> rfSize, 1, tTerm -> rfCount);
      tTerm -> rfSize = NULL;
    }
  }
}
void stackMeanResponse(Terminal *tTerm, unsigned int rnfCount) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> rnfCount != rnfCount) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  rnfCount has been previously defined:  %10d vs %10d", tTerm -> rnfCount, rnfCount);
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  else {
    tTerm -> rnfCount = rnfCount;
  }
  tTerm -> meanResponse = dvector(1, tTerm -> rnfCount);
}
void unstackMeanResponse(Terminal *tTerm) {
  if (tTerm -> rnfCount > 0) {
    if (tTerm -> meanResponse != NULL) {
      free_dvector(tTerm -> meanResponse, 1, tTerm -> rnfCount);
      tTerm -> meanResponse = NULL;
    }
  }
}
void stackTermLMIIndex(Terminal *tTerm, unsigned int size) {
  if (tTerm -> lmiSize > 0) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  lmiIndex has been previously defined:  %10d vs %10d", tTerm -> lmiSize, size);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  else {
    tTerm -> lmiSize = size;
  }
  tTerm -> lmiIndex = uivector(1, tTerm -> lmiSize);
  tTerm -> lmiValue = dvector(1, tTerm -> lmiSize);
}
void unstackTermLMIIndex(Terminal *tTerm) {
  if(tTerm -> lmiSize > 0) {
    if (tTerm -> lmiIndex != NULL) {
      free_uivector(tTerm -> lmiIndex, 1, tTerm -> lmiSize);
      tTerm -> lmiIndex = NULL;
    }
    if (tTerm -> lmiValue != NULL) {
      free_dvector(tTerm -> lmiValue, 1, tTerm -> lmiSize);
      tTerm -> lmiValue = NULL;
    }
  }
}
void getTerminalInfo(Terminal *termPtr) {
  Rprintf("\nTerminalInfo:  %20x", termPtr);
  Rprintf("\n  LeafCnt: %10d", termPtr -> nodeID);
  Rprintf("\n");
  Rprintf("\n lmiIndex            = %20x", termPtr -> lmiIndex);
  Rprintf("\n lmiSize             = %20d", termPtr -> lmiSize);
  Rprintf("\n lmiValue            = %20x", termPtr -> lmiValue);
}
