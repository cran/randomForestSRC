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


#include          <stdlib.h>
#include           <stdio.h>
#include          "global.h"
#include           "trace.h"
#include          "nrutil.h"
#include      "dataParser.h"
#define SMART_BUFFER_INCR 4096
#define TOKEN_BUFFER_INCR  7
SEXP rfsrcReadMatrix(SEXP traceFlag, 
                     SEXP fName, 
                     SEXP rowType,  
                     SEXP rowCnt, 
                     SEXP tokenDelim, 
                     SEXP colHeader, 
                     SEXP rowHeader) {
  FILE  *fopen();
  FILE  *fPtr;
  char  *_fName;
  SEXP   _sexp_rowType ;
  char **_rowType;
  uint   _rowCnt;
  char  *_tokenDelim;
  char   _colHeadF;
  char   _rowHeadF;
  SmartBuffer *sb;
  uint p, i;
  uint rowCntActual;
  uint colCntActual;
  uint sbSizeActual;
  char flag;
  setTraceFlag(INTEGER(traceFlag)[0], 0);
  _fName = (char*) CHAR(STRING_ELT(AS_CHARACTER(fName), 0));
  _sexp_rowType = rowType;
  _rowCnt = INTEGER(rowCnt)[0];
  _tokenDelim = (char*) CHAR(STRING_ELT(AS_CHARACTER(tokenDelim), 0));
  _colHeadF = (INTEGER(colHeader)[0] != 0) ? TRUE : FALSE;
  _rowHeadF = (INTEGER(rowHeader)[0] != 0) ? TRUE : FALSE;
  _rowType = (char**) vvector(1, _rowCnt);
  for (p = 1; p <= _rowCnt; p++) {
    _rowType[p] = (char*) CHAR(STRING_ELT(AS_CHARACTER(_sexp_rowType), p-1));
    if ((strcmp(_rowType[p], "X") != 0) && 
        (strcmp(_rowType[p], "C") != 0) && 
        (strcmp(_rowType[p], "c") != 0) && 
        (strcmp(_rowType[p], "I") != 0) && 
        (strcmp(_rowType[p], "R") != 0)) {
      Rprintf("\nRF-SRC:  *** ERROR *** ");
      Rprintf("\nRF-SRC:  Invalid predictor type:  [%10d] = %2s", p, _rowType[p]);
      Rprintf("\nRF-SRC:  Type must be 'C', 'c', 'I', or 'R'.");
      Rprintf("\nRF-SRC:  Please Contact Technical Support.");
      error("\nRF-SRC:  The application will now exit.\n");
    }
  }
  fPtr = fopen(_fName, "r");
  sb = parseLineSB(fPtr, *_tokenDelim, 0);
  colCntActual = _rowHeadF ? (sb -> tokenCnt - 1) : (sb -> tokenCnt);
  rowCntActual = 0;
  if (_colHeadF) {
    freeSB(sb);
    sb = parseLineSB(fPtr, *_tokenDelim, 0);
    sbSizeActual = sb -> size;
    freeSB(sb);
  }
  else {
    sbSizeActual = sb -> size;
    freeSB(sb);
  }
  ++rowCntActual;
  flag = TRUE;
  while (flag) {
    sb = parseLineSB(fPtr, *_tokenDelim, sbSizeActual);
    if (sb -> tokenCnt == 0) {
      flag = FALSE;
      freeSB(sb);
    }
    else {
      ++rowCntActual;
      freeSB(sb);
    }
  }
  if (rowCntActual != _rowCnt) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Inconsistent Predictor Count.");
    Rprintf("\nRF-SRC:  (encountered, expected) =  (%10d, %10d)", rowCntActual, _rowCnt);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  fclose(fPtr);
  if (_colHeadF) {
    fPtr = fopen(_fName, "r");
    sb = parseLineSB(fPtr, *_tokenDelim, sbSizeActual);
    freeSB(sb);
  }
  char **dataMatrix = cmatrix(1, colCntActual, 1, rowCntActual);
  for (p=1; p <= rowCntActual; p++) {
    sb = parseLineSB(fPtr, *_tokenDelim, sbSizeActual);
    if (_rowHeadF) {
    }
    for (i = 1; i <= colCntActual; i++) {
      dataMatrix[i][p] = (char) strtol(sb -> token, NULL, 10);
    }
    freeSB(sb);
  }
  free_cmatrix(dataMatrix, 1, colCntActual, 1, rowCntActual);
  return R_NilValue;
}
void getDataMatrix(FILE *fPtr, 
                   uint predictorCount, 
                   uint recordCount,
                   char rowLabelFlag,
                   char colLabelFlag, 
                   char delimiter, 
                   char **dataMatrix) {
}
void printSB(SmartBuffer *sb) {
  uint j, actualSize;
  char *orgPtr;
  orgPtr = sb -> ptr;
  sb -> ptr = sb -> buff;
  for (j=1; j <= sb -> tokenCnt; j++) {
    actualSize = getNextTokenSB(sb);
    Rprintf("\n Token: %10d %10d %s ", j, actualSize, sb -> token);
  }
  Rprintf("\n");
  sb -> ptr = orgPtr;
}
uint getNextTokenSB(SmartBuffer *sb) {
  char *tokenPtr, *ptr;
  char overFlow;
  uint size;
  overFlow = FALSE;
  if (*(sb -> ptr) == '\0') {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Token requested, none available..");
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  while (' ' == *((sb -> ptr)++)) {
  }
  ptr = sb -> ptr;
  tokenPtr = sb -> token;
  size = sb -> tokenSize;
  while ((*(sb -> ptr) != '\0') && ((sb -> delim) != *(sb -> ptr))) {
    if (size > 0) {
      *(tokenPtr++) = *(sb -> ptr);
      --size;
    }
    else {
      overFlow = TRUE;
    }
    (sb -> ptr)++;
  }
  if (overFlow) {
    Rprintf("\nRF-SRC:  *** ERROR *** ");
    Rprintf("\nRF-SRC:  Token Buffer Size Too Small.");
    Rprintf("\nRF-SRC:  (encountered, expected) =  (%10d, %10d)", (sb -> ptr) - ptr, sb -> tokenSize);
    Rprintf("\nRF-SRC:  Please Contact Technical Support.");
    error("\nRF-SRC:  The application will now exit.\n");
  }
  *tokenPtr = '\0';
  return ((sb -> tokenSize) - size);
} 
SmartBuffer *parseLineSB(FILE *fPtr, char delim, size_t sbSize) {
  size_t size;
  SmartBuffer *sb = readLineSB(fPtr, sbSize);
  (sb -> delim) = delim;
  size = 0;
  while((*(sb -> ptr)) != '\0') {
    if ((*(sb -> ptr)) == delim) {
      (sb -> tokenCnt) ++;
      if (size > (sb -> tokenSize)) {
        (sb -> tokenSize) = size;
      }
      size = 0;
    }
    else {
      ++ size;
    }
    ++(sb -> ptr);
  }
  if ((sb -> tokenCnt) > 0) {
    (sb -> tokenCnt) ++;
    (sb -> token) = cvector(0, sb -> tokenSize);
  }
  sb -> ptr = sb -> buff;
  return sb;
}
SmartBuffer *readLineSB(FILE *fPtr, size_t sbSize) {
  SmartBuffer *sb;
  char c, eolF;
  size_t i;
  if (sbSize == 0) {
    sbSize = SMART_BUFFER_INCR;
  }
  sb = reallocSB(NULL, sbSize);
  eolF  = FALSE;
  i = 0;
  while (!eolF) {
    c = getc(fPtr);
    if ((c == EOF) || (c == '\n') || (c == '\r')) {
      eolF = TRUE;
    }
    else {
      *(sb -> ptr) = c;
      ++(sb -> ptr);
      ++i;
      if (i == sb -> size) {
        sb = reallocSB(sb, (sb -> size));
      }
    }
  }
  *(sb -> ptr) = '\0';
  sb -> ptr = sb -> buff;
  return sb;
}
SmartBuffer *reallocSB(SmartBuffer *sb, size_t incr) {
  SmartBuffer *sbNew;
  size_t i;
  if (sb == NULL) {
    sbNew = (SmartBuffer*) gblock((size_t) sizeof(SmartBuffer));
    (sbNew -> buff) = cvector(0, incr - 1);
    (sbNew -> size) = incr;
    sbNew -> ptr = sbNew -> buff;
  }
  else {
    sbNew = (SmartBuffer*) gblock((size_t) sizeof(SmartBuffer));
    (sbNew -> buff) = cvector(0, (sb -> size) + incr - 1);
    (sbNew -> size) = (sb -> size) + incr;
    sb -> ptr = sb -> buff;
    sbNew -> ptr = sbNew -> buff;
    for (i = 1; i <= sb -> size; i++) {
      *(sbNew -> ptr) = *(sb -> ptr);
      ++(sbNew -> ptr);
      ++(sb -> ptr);
    }
    freeSB(sb);
  }
  sbNew -> delim = ',';
  sbNew -> tokenCnt = 0;
  sbNew -> tokenSize = 0;
  sbNew -> token = NULL;
  return sbNew;
}
void freeSB(SmartBuffer *sb) {
  free_cvector(sb -> buff, 0, (sb -> size) - 1);
  if ((sb -> token) != NULL) {
    free_cvector(sb -> token, 0, sb -> size);
  }
  free_gblock(sb, sizeof(SmartBuffer));
}
