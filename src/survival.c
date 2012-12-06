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
#include        "impute.h"
#include      "survival.h"
void getAtRiskAndEventCounts(uint treeID) {
  Node *parent;
  uint leaf, i, j, k;
  uint *membershipIndex;
  char eventFlag;
  if ((RF_opt & OPT_BOOT_NODE) | (RF_opt & OPT_BOOT_NONE)) {
    membershipIndex = RF_trivialBootMembershipIndex;
  }
  else {
    membershipIndex = RF_bootMembershipIndex[treeID];
  }
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    stackAtRisk(parent, RF_eventTypeSize, RF_masterTimeSize);
    for (j = 1; j <= RF_masterTimeSize; j++) {
      (parent -> atRiskCount)[j]    = 0;
      (parent -> eventTimeIndex)[j] = 0;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        (parent -> eventCount)[k][j] = 0;
      }
    }
    parent -> eTimeSize = 0;
    parent -> membrCount = 0;
    for (i=1; i <= RF_observationSize; i++) {
      if (RF_nodeMembership[treeID][membershipIndex[i]] == parent) {
        for (j = 1; j <= RF_masterTimeIndex[treeID][membershipIndex[i]]; j++) {
          (parent -> atRiskCount)[j] ++;
        }
        if (RF_status[treeID][membershipIndex[i]] > 0) {
          if (RF_eventTypeSize > 1) {
            k = RF_eventTypeIndex[(uint) RF_status[treeID][membershipIndex[i]]];
          }
          else {
            k = 1;
          }
          (parent -> eventCount)[k][RF_masterTimeIndex[treeID][membershipIndex[i]]] ++;
        }
        (parent -> membrCount) ++;
      }
    }
    for (j = 1; j <= RF_masterTimeSize; j++) {
      eventFlag = FALSE;
      for (k = 1; k <= RF_eventTypeSize; k++) {
        if ((parent -> eventCount)[k][j] > 0) {
          eventFlag = TRUE;
          k = RF_eventTypeSize;
        }
      }
      if (eventFlag == TRUE) {
        (parent -> eventTimeIndex)[++(parent -> eTimeSize)] = j;
      }
    }
    if (parent -> membrCount > 0) {
      parent -> predictedOutcome = 0.0;
    }  
    else {
      parent -> predictedOutcome = NA_REAL;
      if (!(RF_opt & OPT_OUTC_TYPE)) {
        Rprintf("\nRF-SRC:  *** ERROR *** ");
        Rprintf("\nRF-SRC:  Zero node count encountered for Nelson-Aalen estimate in leaf:  %10d", leaf);
        Rprintf("\nRF-SRC:  Please Contact Technical Support.");
        error("\nRF-SRC:  The application will now exit.\n");
      }
    }
  }  
}
void getLocalRatio(uint treeID) {
  Node *parent;
  uint leaf, j, q; 
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    if (parent -> membrCount > 0) {
      if(parent -> eTimeSize > 0) {
        stackLocalRatio(parent, RF_eventTypeSize, parent -> eTimeSize);
        for (j=1; j <= RF_eventTypeSize; j++) {      
          for (q = 1; q <= parent -> eTimeSize; q++) {
            (parent -> localRatio)[j][q] = 0.0;
          }
        }
        for (j=1; j <= RF_eventTypeSize; j++) {      
          for (q = 1; q <= parent -> eTimeSize; q++) {
            if ((parent -> eventCount)[j][(parent -> eventTimeIndex)[q]] > 0) {
              if ((parent -> atRiskCount)[(parent -> eventTimeIndex)[q]] >= 1) {
                (parent -> localRatio)[j][q] = ((double) (parent -> eventCount)[j][(parent -> eventTimeIndex)[q]] / (parent -> atRiskCount)[(parent -> eventTimeIndex)[q]]);
              }
              else {
                Rprintf("\nRF-SRC:  *** ERROR *** ");
                Rprintf("\nRF-SRC:  Zero At Risk Count encountered in local ratio calculation for (tree, node) = (%10d, %10d)", treeID, leaf);
                Rprintf("\nRF-SRC:  Please Contact Technical Support.");
                error("\nRF-SRC:  The application will now exit.\n");
              }
            }
          }
        }
      }
    }
  }
}
void getLocalSurvival(uint treeID) {
  Node *parent;
  uint leaf, j, q;
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    if(parent -> eTimeSize > 0) {
      stackLocalSurvival(parent, parent -> eTimeSize);
      for (q = 1; q <= parent -> eTimeSize; q++) {
        (parent -> localSurvival)[q] = 0.0;
      }
      for (q = 1; q <= parent -> eTimeSize; q++) {
        for (j = 1; j <= RF_eventTypeSize; j++) {
          (parent -> localSurvival)[q] += (parent -> localRatio)[j][q];
        }
        (parent -> localSurvival)[q] = 1.0 - (parent -> localSurvival)[q];
      }  
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localSurvival)[q] *= (parent -> localSurvival)[q-1];
      }
    }
  }  
}
void getLocalNelsonAalen(uint treeID) {
  Node *parent;
  uint leaf, q;
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    if (parent -> eTimeSize > 0) {
      stackLocalNelsonAalen(parent, parent -> eTimeSize);
      for (q = 1; q <= parent -> eTimeSize; q++) {
        (parent -> localNelsonAalen)[q] = 0.0;
      }
      for (q = 1; q <= parent -> eTimeSize; q++) {
        (parent -> localNelsonAalen)[q] = (parent -> localRatio)[1][q];
      }
      for (q = 2; q <= parent -> eTimeSize; q++) {
        (parent -> localNelsonAalen)[q] += (parent -> localNelsonAalen)[q-1];
      }
    }
  }  
}
void getLocalCSH(uint treeID) {
  Node *parent;
  uint leaf, j, q;
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    if (parent -> eTimeSize > 0) {
      stackLocalCSH(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j=1; j <= RF_eventTypeSize; j++) {      
        for (q = 1; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] = 0.0;
        }
      }
      for (j=1; j <= RF_eventTypeSize; j++) {      
        for (q = 1; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] = (parent -> localRatio)[j][q];
        }
        for (q = 2; q <= parent -> eTimeSize; q++) {
          (parent -> localCSH)[j][q] += (parent -> localCSH)[j][q-1];
        }
      }
    }
  }  
}
void getLocalCIF(uint treeID) {
  Node *parent;
  uint leaf, j, q;
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    if(parent -> eTimeSize > 0) {
      stackLocalCIF(parent, RF_eventTypeSize, parent -> eTimeSize);
      for (j=1; j <= RF_eventTypeSize; j++) {      
        for (q = 1; q <= parent -> eTimeSize; q++) {
          (parent -> localCIF)[j][q] = 0.0;
        }
      }
      for (j=1; j <= RF_eventTypeSize; j++) {
        (parent -> localCIF)[j][1] = (parent -> localRatio)[j][1];
        for (q = 2; q <= parent -> eTimeSize; q++) {
          (parent -> localCIF)[j][q] = (parent -> localSurvival)[q-1] * (parent -> localRatio)[j][q];
        }
        for (q = 2; q <= parent -> eTimeSize; q++) {
          (parent -> localCIF)[j][q] += (parent -> localCIF)[j][q-1];
        }
      }
    }
  }  
}
void getNelsonAalen(uint treeID) {
  Node *parent;
  uint leaf, i, k, q;
  uint priorTimePointIndex, currentTimePointIndex;
  for (leaf = 1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    stackNelsonAalen(parent, RF_sortedTimeInterestSize);
    for (q = 1; q <= RF_sortedTimeInterestSize; q++) {
      (parent -> nelsonAalen)[q] = 0.0;
    }
    priorTimePointIndex = 0;
    currentTimePointIndex = 1;
    for (i = 1; i <= (parent -> eTimeSize); i++) {
      for (k = priorTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
        if (RF_timeInterest[k] <= RF_masterTime[(parent -> eventTimeIndex)[i]] ) {
          currentTimePointIndex = k;
        }
        else {
          k = RF_sortedTimeInterestSize;
        }
      }
      (parent -> nelsonAalen)[currentTimePointIndex] = (parent -> localNelsonAalen)[i];
      if (i > 1) {
        for(k = priorTimePointIndex + 1; k < currentTimePointIndex; k++) {
          (parent -> nelsonAalen)[k] = (parent -> nelsonAalen)[priorTimePointIndex];
        }
      }
      if (i == (parent -> eTimeSize)) {
        for(k = currentTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
          (parent -> nelsonAalen)[k] = (parent -> nelsonAalen)[currentTimePointIndex];
        }
      }
      priorTimePointIndex = currentTimePointIndex;
    }
  }  
}
void getSurvival(uint treeID) {
  Node *parent;
  uint priorTimePointIndex, currentTimePointIndex;
  uint leaf, i, k;
  for (leaf=1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    stackSurvival(parent, RF_sortedTimeInterestSize);
    for (k=1; k <= RF_sortedTimeInterestSize; k++) {
      (parent -> survival)[k] = 0.0;
    }
    if (parent -> eTimeSize > 0) {
      priorTimePointIndex = 0;
      currentTimePointIndex = 1;
      for (i = 1; i <= (parent -> eTimeSize); i++) {
        for (k = priorTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
          if (RF_timeInterest[k] <= RF_masterTime[(parent -> eventTimeIndex)[i]] ) {
            currentTimePointIndex = k;
          }
          else {
            k = RF_sortedTimeInterestSize;
          }
        }
        (parent -> survival)[currentTimePointIndex] = (parent -> localSurvival)[i];
        if (i == 1) {
          for(k = 1; k < currentTimePointIndex; k++) {
            (parent -> survival)[k] = 1.0;
          }
        }
        if (i > 1) {
          for(k = priorTimePointIndex + 1; k < currentTimePointIndex; k++) {
            (parent -> survival)[k] = (parent -> survival)[priorTimePointIndex];
          }
        }
        if (i == (parent -> eTimeSize)) {
          for(k = currentTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
            (parent -> survival)[k] = (parent -> survival)[currentTimePointIndex];
          }
        }
        priorTimePointIndex = currentTimePointIndex;
      }
    }
    else {
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> survival)[k] = 1.0;
      }
    }
  }  
}
void getCSH(uint treeID) {
  Node *parent;
  uint priorTimePointIndex, currentTimePointIndex;
  uint leaf, i, j, k;
  for (leaf=1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    stackCSH(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {      
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CSH)[j][k] = 0.0;
      }
    }
    priorTimePointIndex = 0;
    currentTimePointIndex = 1;
    for (i = 1; i <= (parent -> eTimeSize); i++) {
      for (k = priorTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
        if (RF_timeInterest[k] <= RF_masterTime[(parent -> eventTimeIndex)[i]] ) {
          currentTimePointIndex = k;
        }
        else {
          k = RF_sortedTimeInterestSize;
        }
      }
      for (j=1; j <= RF_eventTypeSize; j++) {
        (parent -> CSH)[j][currentTimePointIndex] = (parent -> localCSH)[j][i];
        if (i > 1) {
          for(k = priorTimePointIndex + 1; k < currentTimePointIndex; k++) {
            (parent -> CSH)[j][k] = (parent -> CSH)[j][priorTimePointIndex];
          }
        }
        if (i == (parent -> eTimeSize)) {
          for(k = currentTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
            (parent -> CSH)[j][k] = (parent -> CSH)[j][currentTimePointIndex];
          }
        }
      }
      priorTimePointIndex = currentTimePointIndex;
    }
  }  
}
void getCIF(uint treeID) {
  Node *parent;
  uint priorTimePointIndex, currentTimePointIndex;
  uint leaf, i, j, k;
  for (leaf=1; leaf <= RF_leafCount[treeID]; leaf++) {
    parent = RF_terminalNode[treeID][leaf];
    stackCIF(parent, RF_eventTypeSize, RF_sortedTimeInterestSize);
    for (j=1; j <= RF_eventTypeSize; j++) {      
      for (k=1; k <= RF_sortedTimeInterestSize; k++) {
        (parent -> CIF)[j][k] = 0.0;
      }
    }
    priorTimePointIndex = 0;
    currentTimePointIndex = 1;
    for (i = 1; i <= (parent -> eTimeSize); i++) {
      for (k = priorTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
        if (RF_timeInterest[k] <= RF_masterTime[(parent -> eventTimeIndex)[i]] ) {
          currentTimePointIndex = k;
        }
        else {
          k = RF_sortedTimeInterestSize;
        }
      }
      for (j=1; j <= RF_eventTypeSize; j++) {
        (parent -> CIF)[j][currentTimePointIndex] = (parent -> localCIF)[j][i];
        if (i > 1) {
          for(k = priorTimePointIndex + 1; k < currentTimePointIndex; k++) {
            (parent -> CIF)[j][k] = (parent -> CIF)[j][priorTimePointIndex];
          }
        }
        if (i == (parent -> eTimeSize)) {
          for(k = currentTimePointIndex + 1; k <= RF_sortedTimeInterestSize; k++) {
            (parent -> CIF)[j][k] = (parent -> CIF)[j][currentTimePointIndex];
          }
        }
      }
      priorTimePointIndex = currentTimePointIndex;
    }
  }  
}
