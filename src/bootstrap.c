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
#include        "random.h"
#include       "nodeOps.h"
#include        "impute.h"
#include     "bootstrap.h"
char bootstrap (uint     mode,
                uint     treeID,
                Node    *nodePtr,
                uint    *subIndex,
                uint     subsetSize,
                uint    *index) {
  char result;
  uint i,k;
  result = TRUE;
  if (RF_opt & OPT_BOOT_NONE) {
    for (i=1; i <= subsetSize; i++) {
      index[i] = subIndex[i];
    }
  }
  else {
    for (i=1; i <= subsetSize; i++) {
      k = (uint) ceil(ran1(treeID)*(subsetSize * 1.0));
      index[i] = subIndex[k];
    }
  }
  result = getNodeSign(mode, treeID, nodePtr, index, subsetSize);
  if (result == FALSE) {
  }
  if (result == TRUE) {
    if (mode == RF_PRED) {
    }
  } 
  else {
  }
  return result;
}
char getNodeSign (uint mode, uint treeID, Node *nodePtr, uint *bmIndex, uint repMembrSize) {
  int   *mvNSptr;
  int   *fmvNSptr;
  char result;
  uint i,p,q,m;
  result = TRUE;
  switch (mode) {
  case RF_PRED:
    if (RF_mRecordSize > 0) {
      stackMVSign(nodePtr, RF_mvSize);
      mvNSptr = nodePtr -> mvSign;
    }
    else {
      mvNSptr = NULL;
    }
    if (RF_fmRecordSize > 0) {
      stackFMVSign(nodePtr, RF_fmvSize);
      fmvNSptr = nodePtr -> fmvSign;
    }
    else {
      fmvNSptr = NULL;
    }
    break;
  default:
    if (RF_mRecordSize > 0) {
      stackMVSign(nodePtr, RF_mvSize);
      mvNSptr = nodePtr -> mvSign;
    }
    else {
      mvNSptr = NULL;
    }
    fmvNSptr = NULL;
    break;
  }  
  if (mvNSptr != NULL) {
    int **mvBootstrapSign = imatrix(1, RF_mvSize, 1, repMembrSize);
    for (p = 1; p <= RF_mvSize; p++) {
      for (i = 1; i <= repMembrSize; i++) {
        mvBootstrapSign[p][i] = 0;
      }
    }
    for (p = 1; p <= RF_mvSize; p++) {
      mvNSptr[p] = 0;
    }
    for (i=1; i <= repMembrSize; i++) {
      m = bmIndex[i];
      if (RF_mRecordMap[m] != 0) {
        for (p = 1; p <= RF_mvSize; p++) {
          if (RF_mvIndex[p] < 0) {
            mvBootstrapSign[p][i] = RF_mvSign[(uint) abs(RF_mvIndex[p])][RF_mRecordMap[m]];
          }
          else {
            mvBootstrapSign[p][i] = RF_mvSign[RF_rSize + (uint) RF_mvIndex[p]][RF_mRecordMap[m]];
          }
        }
      }
      else {
        for (p=1; p <= RF_mvSize; p++) {
          mvBootstrapSign[p][i] = 0;
        }
      }
      for (p = 1; p <= RF_mvSize; p++) {
        mvNSptr[p] = mvNSptr[p] + mvBootstrapSign[p][i];
      }
    }
    m = 0;
    for (p = 1; p <= RF_mvSize; p++) {
      if (mvNSptr[p] > 0) {
        if (mvNSptr[p] == repMembrSize) {
          mvNSptr[p] = -1;
        }
        else {
          mvNSptr[p] = 1;
        }
      }
      if(RF_mvIndex[p] < 0) {
        if (mvNSptr[p] == -1) result = FALSE;
      }
      else {
        if (mvNSptr[p] == -1) m ++;
      }
    }  
    if (m == RF_mvSize) {
      result = FALSE;
    }
    free_imatrix(mvBootstrapSign, 1, RF_mvSize, 1, repMembrSize);
  }
  if (fmvNSptr != NULL) {
    for (p = 1; p <= RF_fmvSize; p++) {
      fmvNSptr[p] = 1;
    }
    if (RF_mRecordSize > 0) {
      p = q = 1;
      while ((p <= RF_mvSize) && (q <= RF_fmvSize)) {
        if (RF_mvIndex[p] == RF_fmvIndex[q]) {
          if (mvNSptr[p] == -1) {
            fmvNSptr[q] = -1;
          }
          p++;
          q++;
        }
        else if (RF_fmvIndex[q] < 0) {
          if (RF_mvIndex[p] > 0) {
            q++;
          }
          else {
            if (abs(RF_fmvIndex[q]) < abs(RF_mvIndex[p])) {
              q++;
            }
            else {
              p++;
            }
          }
        }
        else {
          if (RF_fmvIndex[q] < RF_mvIndex[p]) {
            q++;
          }
          else {
            p++;
          }
        }
      }  
    }  
  }  
  return result;
}