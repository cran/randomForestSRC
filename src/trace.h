////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.4
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


#ifndef TRACE_H
#define TRACE_H
#define SUMM_USR_TRACE  0x00000001
#define SUMM_LOW_TRACE  0x00000002
#define SUMM_MED_TRACE  0x00000004
#define SUMM_HGH_TRACE  0x00000008
#define SPLT_DEF_TRACE  0x00000010
#define SPLT_LOW_TRACE  0x00000020
#define SPLT_MED_TRACE  0x00000040
#define SPLT_HGH_TRACE  0x00000080
#define FORK_DEF_TRACE  0x00000100
#define MISS_LOW_TRACE  0x00000200
#define MISS_MED_TRACE  0x00000400
#define MISS_HGH_TRACE  0x00000800
#define OUTP_DEF_TRACE  0x00001000
#define NUMR_DEF_TRACE  0x00002000
#define FACT_LOW_TRACE  0x00004000
#define FACT_HGH_TRACE  0x00008000
#define ENSB_LOW_TRACE  0x00010000
#define ENSB_HGH_TRACE  0x00020000
#define BOOT_MED_TRACE  0x00040000
#define VIMP_LOW_TRACE  0x00080000
#define NODE_DEF_TRACE  0x00100000
#define TIME_DEF_TRACE  0x00200000
#define RAND_DEF_TRACE  0x00400000
#define TURN_OFF_TRACE  0x00000000
#define TURN_ON_TRACE   0x00000001
#define TRACE_MASK      0xFFFFFFFF
void setTraceFlag(unsigned int traceFlag, unsigned int tree);
unsigned int getTraceFlag(unsigned int tree);
unsigned int updateTimeStamp(unsigned int before);
unsigned int getNodeDefTraceFlag();
unsigned int getForkDefTraceFlag();
unsigned int getTurnOffTraceFlag();
unsigned int getTurnOnTraceFlag();
unsigned int getNumrDefTraceFlag();
unsigned int getSummUsrTraceFlag();
void setMaxMemoryAllocation(size_t value);
void setMinMemoryAllocation(size_t value);
size_t getMaxMemoryAllocation();
size_t getMinMemoryAllocation();
void increaseMemoryAllocation(size_t amount);
void decreaseMemoryAllocation(size_t amount);
void changeMemoryAllocation(size_t amount, int direction);
#endif
