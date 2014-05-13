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


#include  <time.h>
#include "trace.h"
#include   <R_ext/Print.h>
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
unsigned int    RF_traceFlagDiagLevel;
unsigned int    RF_traceFlagIterValue;
size_t          RF_memor_maxMemoryAllocation;
size_t          RF_memor_minMemoryAllocation;
void setTraceFlag(unsigned int traceFlag, unsigned int tree) {
  RF_traceFlagDiagLevel = traceFlag & TRACE_MASK;
  RF_traceFlagIterValue = tree;
}
unsigned int getTraceFlag(unsigned int tree) {
  unsigned int result;
  result = 0;
  if (RF_traceFlagIterValue == tree) {
    result = RF_traceFlagDiagLevel;
  }
  else {
    if (RF_traceFlagIterValue == 0) {
      result = RF_traceFlagDiagLevel;
    }
  }
  return result;
}
unsigned int updateTimeStamp(unsigned int before) {
  unsigned int stamp;
  double cpuTimeUsed;
  stamp = clock();
  cpuTimeUsed = ((double) (stamp - before)) / CLOCKS_PER_SEC;
  if (getTraceFlag(0) & SUMM_USR_TRACE) {
    Rprintf("\nRF-SRC:  CPU process time:  %20.3f \n", cpuTimeUsed);
  }
  return stamp;
}
unsigned int getNodeDefTraceFlag() {
  return(NODE_DEF_TRACE);
}
unsigned int getForkDefTraceFlag() {
  return(FORK_DEF_TRACE);
}
unsigned int getTurnOffTraceFlag() {
  return(TURN_OFF_TRACE);
}
unsigned int getTurnOnTraceFlag() {
  return(TURN_ON_TRACE);
}
unsigned int getNumrDefTraceFlag() {
  return(NUMR_DEF_TRACE);
}
unsigned int getTimeDefTraceFlag() {
  return(TIME_DEF_TRACE);
}
void setMaxMemoryAllocation(size_t value) {
  RF_memor_maxMemoryAllocation = value;
}
void setMinMemoryAllocation(size_t value) {
  RF_memor_minMemoryAllocation = value;
}
size_t getMaxMemoryAllocation() {
  return (RF_memor_maxMemoryAllocation);
}
size_t getMinMemoryAllocation() {
  return (RF_memor_minMemoryAllocation);
}
void increaseMemoryAllocation(size_t amount) {
  RF_memor_minMemoryAllocation += amount;
  if (RF_memor_minMemoryAllocation > RF_memor_maxMemoryAllocation) {
    RF_memor_maxMemoryAllocation = RF_memor_minMemoryAllocation;
  }
}
void decreaseMemoryAllocation(size_t amount) {
    RF_memor_minMemoryAllocation -= amount;  
}
