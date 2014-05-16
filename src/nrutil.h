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


#ifndef NRUTIL_H
#define NRUTIL_H
#include     "node.h"
#include "terminal.h"
unsigned int upower (unsigned int x, unsigned int n);
unsigned int upower2 (unsigned int n);
unsigned int ulog2 (unsigned int n);
void hpsort(
  double *ra,
  unsigned int n
);
void hpsortui(
  unsigned int *ra,
  unsigned int n
);
void hpsortxx(unsigned int n, double *arr,  unsigned *indx);
void indexx(
  unsigned int n,
  double *arr, 
  unsigned int *indx
);
void nrerror(char error_text[]);
void *gblock(size_t size);
void free_gblock(void *v, size_t size);
void *gvector(unsigned long nl, unsigned long nh, size_t size);
void free_gvector(void *v, unsigned long nl, unsigned long nh, size_t size);
void *vvector(unsigned long nl, unsigned long nh);
void free_vvector(void *v, unsigned long nl, unsigned long nh);
void **vmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_vmatrix(void **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
char *cvector(unsigned long nl, unsigned long nh);
void free_cvector(char *v, unsigned long nl, unsigned long nh);
char **cmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_cmatrix(char **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
int *ivector(unsigned long nl, unsigned long nh);
void free_ivector(int *v, unsigned long nl, unsigned long nh);
int **imatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_imatrix(int **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
unsigned int *uivector(unsigned long nl, unsigned long nh);
void free_uivector(unsigned int *v, unsigned long nl, unsigned long nh);
unsigned int **uimatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_uimatrix(unsigned int **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
double *dvector(unsigned long nl, unsigned long nh);
void free_dvector(double *v, unsigned long nl, unsigned long nh);
double **dmatrix(unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_dmatrix(double **v, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
double ***dmatrix3(unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_dmatrix3(double ***v, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
double ****dmatrix4(unsigned long n4l, unsigned long n4h, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void free_dmatrix4(double ****v, unsigned long n4l, unsigned long n4h, unsigned long n3l, unsigned long n3h, unsigned long nrl, unsigned long nrh, unsigned long ncl, unsigned long nch);
void nrCopyMatrix(
  unsigned int **new,
  unsigned int **old,
  unsigned int nrow,
  unsigned int ncol
);
void nrCopyVector(
  char *new, 
  char *old, 
  unsigned int ncol
);
#endif
