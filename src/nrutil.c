////**********************************************************************
////**********************************************************************
////
////  RANDOM FORESTS FOR SURVIVAL, REGRESSION, AND CLASSIFICATION (RF-SRC)
////  Version 1.5.1
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


#ifdef SUPPORT_OPENMP
#include           <omp.h>
#endif
#include      <stdlib.h>
#include      "nrutil.h"
#include <R_ext/Print.h>
#include <Rdefines.h>
#include <time.h>
extern unsigned int getTraceFlag();
extern unsigned int getNumrDefTraceFlag();
extern unsigned int getTimeDefTraceFlag();
extern void increaseMemoryAllocation(size_t amount);
extern void decreaseMemoryAllocation(size_t amount);
extern size_t getMaxMemoryAllocation();
extern size_t getMinMemoryAllocation();
#ifndef TRUE
#define TRUE   0x01
#endif
#ifndef FALSE
#define FALSE  0x00
#endif
unsigned int upower (unsigned int x, unsigned int n) {
  unsigned int p;
  if ((x >= 2) & (n > (sizeof(unsigned int) * 8) - 1)) {
    nrerror("Overflow in upower(), exponent too large.");
  }
  for (p = 1; n > 0; --n) {
    p = p * x;
  }
  return p;
}
unsigned int upower2 (unsigned int n) {
  unsigned int p;
  if (n > (sizeof(unsigned int) * 8) - 1) {
    nrerror("Overflow in upower2(), exponent too large.");
  }
  p = ((unsigned int) 1) << n;
  return p;
}
unsigned int ulog2 (unsigned int n) {
  unsigned int p;
  p = 0;
  while (n > 1) {
    n = n >> 1;
    p++;
  }
  return p;
}
void hpsort(double *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  double rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
void hpsortui(unsigned int *ra, unsigned int n) {
  unsigned int i, ir, j, l;
  unsigned int rra;
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    }
    else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      }
      else {
        j = ir+1;
      }
    }
    ra[i] = rra;
  }
}
#ifdef SWAP
#undef SWAP
#endif
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50
void indexx(unsigned int n, double *arr, unsigned int *indx) {
  unsigned int i, j, k, l;
  unsigned int indxt, itemp, ir;
  unsigned int *istack, jstack;
  double a;
  if (n < 1) nrerror("\n n of zero (0) length in indexx().");
  l  = 1;
  ir = n;
  jstack = 0;
  istack = uivector(1, NSTACK);
  for (j=1; j<=n; j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j = l+1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i=j-1; i>=l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1] = indx[i];
        }
        indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir])
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1], indx[ir])
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP(indx[l], indx[l+1])
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      a = arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in indexx().");
      if (ir-i+1 >= j-l) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j-1;
      }
      else {
        istack[jstack] = j-1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free_uivector(istack, 1, NSTACK);
}
#undef SWAP
#undef M
#undef NSTACK
#define FREE_ARG char*
#define NR_END 1
void nrerror(char error_text[]) {
  Rprintf("\n");
  Rprintf("\n  *** ERROR *** ");
  Rprintf("\n  Numerical Recipes Run-Time Error:");
  Rprintf("\n  %s", error_text);
  Rprintf("\n  Please Contact Technical Support.");
  error("\n  The application will now exit.\n");
}
void *gblock(size_t size) {
  void *v = (void *) malloc(size);
  if (!v) nrerror("\n  Allocation Failure in gblock().");
  return v;
}
void free_gblock(void *v, size_t size) {
  free((FREE_ARG) (v));
}
void *gvector(unsigned long nl, unsigned long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in gvector().");
  void *v = gblock((size_t) ((nh-nl+1+NR_END) * size));
  v = v-nl+NR_END;
  return v;
}
void free_gvector(void *v, unsigned long nl, unsigned long nh, size_t size) {
  if (nh < nl) nrerror("\n  Illegal indices in free_gvector().");
  free_gblock(v+nl-NR_END, (nh-nl+1+NR_END) * size);
}
void *vvector(unsigned long nl, unsigned long nh) {
  return ((void *) gvector(nl, nh, sizeof(void*)));
}
void free_vvector(void *v, unsigned long nl, unsigned long nh) {
  free_gvector(v, nl, nh, sizeof(void*));
}
void **vmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  void **v = (void **) vvector(nrl, nrh);
  for(unsigned long i = nrl; i <= nrh; i++) {
    v[i] = vvector(ncl, nch);
  }
  return v;
}
void free_vmatrix(void **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = nrl; i <= nrh; i++) {
    free_vvector(v[i], ncl, nch);
  }
  free_vvector(v, nrl, nrh);
}
char *cvector(unsigned long nl, unsigned long nh) {
  return ((char *) gvector(nl, nh, sizeof(char)));
}
void free_cvector(char *v, unsigned long nl, unsigned long nh) {
  free_gvector(v, nl, nh, sizeof(char));
}
char **cmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  char **v = (char **) vvector(nrl, nrh);
  for(unsigned long i = nrl; i <= nrh; i++) {
    v[i] = cvector(ncl, nch);
  }
  return v;
}
void free_cmatrix(char **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = nrl; i <= nrh; i++) {
    free_cvector(v[i], ncl, nch);
  }
  free_vvector(v, nrl, nrh);
}
int *ivector(unsigned long nl, unsigned long nh) {
  return ((int *) gvector(nl, nh, sizeof(int)));
}
void free_ivector(int *v, unsigned long nl, unsigned long nh) {
  free_gvector(v, nl, nh, sizeof(int));
}
int **imatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  int **v = (int **) vvector(nrl, nrh);
  for(unsigned long i = nrl; i <= nrh; i++) {
    v[i] = ivector(ncl, nch);
  }
  return v;
}
void free_imatrix(int **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = nrl; i <= nrh; i++) {
    free_ivector(v[i], ncl, nch);
  }
  free_vvector(v, nrl, nrh);
}
unsigned int *uivector(unsigned long nl, unsigned long nh) {
  return ((unsigned int *) gvector(nl, nh, sizeof(unsigned int)));
}
void free_uivector(unsigned int *v, unsigned long nl, unsigned long nh) {
  free_gvector(v, nl, nh, sizeof(unsigned int));
}
unsigned int **uimatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  unsigned int **v = (unsigned int **) vvector(nrl, nrh);
  for(unsigned long i = nrl; i <= nrh; i++) {
    v[i] = uivector(ncl, nch);
  }
  return v;
}
void free_uimatrix(unsigned int **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = nrl; i <= nrh; i++) {
    free_uivector(v[i], ncl, nch);
  }
  free_vvector(v, nrl, nrh);
}
double *dvector(unsigned long nl, unsigned long nh) {
  return ((double *) gvector(nl, nh, sizeof(double)));
}
void free_dvector(double *v, unsigned long nl, unsigned long nh) {
  free_gvector(v, nl, nh, sizeof(double));
}
double **dmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  double **v = (double **) vvector(nrl, nrh);
  for(unsigned long i = nrl; i <= nrh; i++) {
    v[i] = dvector(ncl, nch);
  }
  return v;
}
void free_dmatrix(double **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = nrl; i <= nrh; i++) {
    free_dvector(v[i], ncl, nch);
  }
  free_vvector(v, nrl, nrh);
}
double ***dmatrix3(unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  double ***v = (double ***) vvector(n3l, n3h);
  for(unsigned long i = n3l; i <= n3h; i++) {
    v[i] = dmatrix(nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix3(double ***v, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = n3l; i <= n3h; i++) {
    free_dmatrix(v[i], nrl, nrh, ncl, nch);
  }
  free_vvector(v, n3l, n3h);
}
double ****dmatrix4(unsigned long n4l, unsigned long n4h, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  double ****v = (double ****) vvector(n4l, n4h);
  for(unsigned long i = n4l; i <= n4h; i++) {
    v[i] = dmatrix3(n3l, n3h, nrl, nrh, ncl, nch);
  }
  return v;
}
void free_dmatrix4(double ****v, unsigned long n4l, unsigned long n4h, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch) {
  for(unsigned long i = n4l; i <= n4h; i++) {
    free_dmatrix3(v[i], n3l, n3h, nrl, nrh, ncl, nch);
  }
  free_vvector(v, n4l, n4h);
}
#undef FREE_ARG
#undef NR_END
void nrCopyMatrix(unsigned int **new, unsigned int **old, unsigned int nrow, unsigned int ncol) {
  unsigned int i,j;
  for (i = 1; i <= nrow; i++) {
    for (j = 1; j <= ncol; j++) {
      new[i][j] = old[i][j];
    }
  }
}
void nrCopyVector(char *new, char *old, unsigned int ncol) {
  unsigned int j;
  for (j = 1; j <= ncol; j++) {
    new[j] = old[j];
  }
}
void testEndianness() {
  unsigned int     test = 0x12345678;
  unsigned int *testPtr = & test;
  Rprintf("\n %2x %2x %2x %2x \n",
          *((char *) testPtr),
          *((char *) testPtr + 1),
          *((char *) testPtr + 2),
          *((char *) testPtr + 3));
}
