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


#include "global.h"
#include  "trace.h"
#include "nrutil.h"
#include "random.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
int  *ran1_iy;
int **ran1_iv;
int  *ran2_iy;
int **ran2_iv;
int      *seed1Value;
int      *seed2Value;
void randomStack(uint valSize) {
  uint b;
  ran1_iy = ivector(0, valSize-1);
  ran1_iv = imatrix(0, valSize-1, 0, NTAB-1);
  ran2_iy = ivector(0, valSize-1);
  ran2_iv = imatrix(0, valSize-1, 0, NTAB-1);
  for (b = 0; b < valSize; b++) {
    ran1_iy[b] = 0;
    ran2_iy[b] = 0;
  }
  seed1Value = ivector(0, valSize-1);
  seed2Value = ivector(0, valSize-1);
}
void randomUnstack(uint valSize) {
  free_ivector(ran1_iy, 0, valSize-1);
  free_imatrix(ran1_iv, 0, valSize-1, 0, NTAB-1);
  free_ivector(ran2_iy, 0, valSize-1);
  free_imatrix(ran2_iv, 0, valSize-1, 0, NTAB-1);
  free_ivector(seed1Value, 0, valSize-1);
  free_ivector(seed2Value, 0, valSize-1);
}
void randomSetChainParallel(uint b, int value) {
  seed1Value[b-1] = value;
}
void randomSetUChainParallel(uint b, int value) {
  seed2Value[b-1] = value;
}
void randomSetChainSerial(uint b, int value) {
  seed1Value[0] = value;
}
void randomSetUChainSerial(uint b, int value) {
  seed2Value[0] = value;
}
int randomGetChainParallel(uint b) {
  return seed1Value[b-1];
}
int randomGetUChainParallel(uint b) {
  return seed2Value[b-1];
}
int randomGetChainSerial(uint b) {
  return seed1Value[0];
}
int randomGetUChainSerial(uint b) {
  return seed2Value[0];
}
float randomChainParallel(uint b) {
  return  ran1_generic(& ran1_iy[b-1], ran1_iv[b-1], & seed1Value[b-1]);
}
float randomUChainParallel(uint b) {
  return  ran1_generic(& ran2_iy[b-1], ran2_iv[b-1], & seed2Value[b-1]);
}
float randomChainSerial(uint b) {
  return  ran1_generic(& ran1_iy[0], ran1_iv[0], & seed1Value[0]);
}
float randomUChainSerial(uint b) {
  return  ran1_generic(& ran2_iy[0], ran2_iv[0], & seed2Value[0]);
}
float ran1_generic(int *iy, int *iv, int *idum) {
  int j, k;
  float temp;
  if (*idum <= 0 || !(*iy)) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    (*iy) = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0) *idum += IM;
  j = (*iy) / NDIV;
  (*iy) = iv[j];
  iv[j] = *idum;
  if ((temp = AM * (*iy)) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}
float ran1_original(int *idum) {
  int j;
  int k;
  static int iy = 0;
  static int iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) {
      *idum = 1;
    }
    else {
      *idum = -(*idum);
    }
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ;
      *idum = IA * (*idum - k * IQ) - IR * k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }
  k = (*idum) / IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0) *idum += IM;
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp = AM * iy) > RNMX) {
    return RNMX;
  }
  else {
    return temp;
  }
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
