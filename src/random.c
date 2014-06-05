////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.2
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


#include "global.h"
#include  "trace.h"
#include "nrutil.h"
#include "random.h"
#define IA      16807
#define IM      2147483647
#define AM      (1.0/IM)
#define IQ      127773
#define IR      2836
#define NTAB    32
#define NDIV    (1+(IM-1)/NTAB)
#define EPS     1.2e-7
#define RNMX    (1.0-EPS)
#define LCG_IM  714025
#define LCG_IA  1366
#define LCG_IC  150889
int  *ran1A_iy;
int **ran1A_iv;
int  *ran1B_iy;
int **ran1B_iv;
int  *ran1C_iy;
int **ran1C_iv;
int      *seed1AValue;
int      *seed1BValue;
int      *seed1CValue;
void randomStack(uint bSize, uint pSize) {
  uint b;
  ran1A_iy = ivector(1, bSize);
  ran1A_iv = imatrix(1, bSize, 1, NTAB);
  ran1B_iy = ivector(1, bSize);
  ran1B_iv = imatrix(1, bSize, 1, NTAB);
  ran1C_iy = ivector(1, bSize);
  ran1C_iv = imatrix(1, bSize, 1, NTAB);
  for (b = 1; b <= bSize; b++) {
    ran1A_iy[b] = 0;
    ran1B_iy[b] = 0;
    ran1C_iy[b] = 0;
  }
  seed1AValue = ivector(1, bSize);
  seed1BValue = ivector(1, bSize);
  seed1CValue = ivector(1, bSize);
}
void randomUnstack(uint bSize, uint pSize) {
  free_ivector(ran1A_iy, 1, bSize);
  free_imatrix(ran1A_iv, 1, bSize, 1, NTAB);
  free_ivector(ran1B_iy, 1, bSize);
  free_imatrix(ran1B_iv, 1, bSize, 1, NTAB);
  free_ivector(ran1C_iy, 1, bSize);
  free_imatrix(ran1C_iv, 1, bSize, 1, NTAB);
  free_ivector(seed1AValue, 1, bSize);
  free_ivector(seed1BValue, 1, bSize);
  free_ivector(seed1CValue, 1, bSize);
}
void randomSetChainParallel(uint b, int value) {
  seed1AValue[b] = value;
}
void randomSetUChainParallel(uint b, int value) {
  seed1BValue[b] = value;
}
void randomSetUChainParallelCov(uint b, int value) {
  seed1CValue[b] = value;
}
void randomSetChainSerial(uint b, int value) {
  seed1AValue[1] = value;
}
void randomSetUChainSerial(uint b, int value) {
  seed1BValue[1] = value;
}
void randomSetUChainSerialCov(uint b, int value) {
  seed1CValue[1] = value;
}
int randomGetChainParallel(uint b) {
  return seed1AValue[b];
}
int randomGetUChainParallel(uint b) {
  return seed1BValue[b];
}
int randomGetUChainParallelCov(uint b) {
  return seed1CValue[b];
}
int randomGetChainSerial(uint b) {
  return seed1AValue[1];
}
int randomGetUChainSerial(uint b) {
  return seed1BValue[1];
}
int randomGetUChainSerialCov(uint b) {
  return seed1CValue[1];
}
float randomChainParallel(uint b) {
  return  ran1_generic(& ran1A_iy[b], ran1A_iv[b], & seed1AValue[b]);
}
float randomUChainParallel(uint b) {
  return  ran1_generic(& ran1B_iy[b], ran1B_iv[b], & seed1BValue[b]);
}
float randomUChainParallelCov(uint b) {
  return  ran1_generic(& ran1C_iy[b], ran1C_iv[b], & seed1CValue[b]);
}
float randomChainSerial(uint b) {
  return  ran1_generic(& ran1A_iy[1], ran1A_iv[1], & seed1AValue[1]);
}
float randomUChainSerial(uint b) {
  return  ran1_generic(& ran1B_iy[1], ran1B_iv[1], & seed1BValue[1]);
}
float randomUChainSerialCov(uint b) {
  return  ran1_generic(& ran1C_iy[1], ran1C_iv[1], & seed1CValue[1]);
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
    (*iy) = iv[1];
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
void lcgenerator(unsigned int *seed, unsigned char reset) {
  if (reset) {
    if (*seed >= LCG_IM) (*seed) %= LCG_IM;
  }
  else {
    *seed = (LCG_IA * (*seed) + LCG_IC) % LCG_IM;
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
