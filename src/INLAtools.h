
/* INLAtools.h
 *
 * Copyright (C) 2025-2026 Elias T Krainski
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Elias T Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#if defined(INLA_WITH_EXTERNAL_PACKAGES)
#       include <ltdl.h>
#       include <omp.h>
#else
#       include <dlfcn.h>
#       include <R.h>
#       include <Rdefines.h>
#       include <Rinternals.h>
#       include <R_ext/Rdynload.h>			       // needed to allow user interrupts
#       include <R_ext/Utils.h>				       // needed to allow user interrupts
#endif
#include "cgeneric.h"

#if !defined(Calloc)
#       define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#endif
#define SQR(x) ((x)*(x))
#define pow2(x) ((x)*(x))
#define pow3(x) (pow2(x)*(x))
#define pow4(x) (pow2(x)*pow2(x))

#define Memcopy(dest, src, n, type) memcpy((void *) (dest), (void *) (src), (size_t) (n) * sizeof(type))
#define Free(x) if (x) { free(x); x = NULL; }

#if !defined(iszero)
#       ifdef __SUPPORT_SNAN__
#              define iszero(x) (fpclassify(x) == FP_ZERO)
#       else
#              define iszero(x) (((__typeof(x))(x)) == 0)
#       endif
#endif

#if __GNUC__ > 7
typedef size_t FORTRAN_CHARLEN_T;
#else
typedef int FORTRAN_CHARLEN_T;
#endif

#define F_ONE ((FORTRAN_CHARLEN_T)1)

// print elements of a matrix
#define printMat(_M, _nr, _nc, _msg)                                    \
if(1) {                                                                 \
  int _i, _j;                                                           \
  printf("%s (%d x %d)\n", _msg, _nr, _nc);                             \
  for(_i=0; _i<_nr; _i++) {                                             \
    for(_j=0; _j<_nc; _j++) {                                           \
      printf("%5.4f ", (_M)[(_nc) * _i + _j]);                          \
    }                                                                   \
    printf("\n");                                                       \
  }                                                                     \
  printf("\n");                                                         \
}                                                                       \

void dgemm_(const char *transa, const char *transb,
            int *m, int *n, int *k, double *alpha, double *a,
            int *lda, double *b, int *ldb, double *beta,
            double *c, int *ldc, FORTRAN_CHARLEN_T);
void dgesv_(int *n, int *nrhs, double *a, int *lda,
            int *ipiv, double *b, int *ldb,
            int *info, FORTRAN_CHARLEN_T);

void uQ2Uexpand(int *nidx, int *idx, double *x, double *xx);
void MuQ2kroneckerU(int *n1, int *n2, int *M2, int *idx2full,
                    double *x1, double *x2, double *xx);
void addMuQ2kroneckerU(int *n1, int *n2, int *M2, int *idx2full,
                       double *x1, double *x2, double *xx);

#if defined(INLA_WITH_EXTERNAL_PACKAGES)
inla_cgeneric_func_tp *inla_cgeneric_mapper(char *name);
#else
SEXP inla_cgeneric_element_get(SEXP Rcmd, SEXP Stheta, SEXP Sntheta, SEXP ints, SEXP doubles, SEXP chars, SEXP mats, SEXP smats);
#endif
inla_cgeneric_func_tp inla_cgeneric_generic0;
inla_cgeneric_func_tp inla_cgeneric_kronecker;
inla_cgeneric_func_tp inla_cgeneric_wmodel;
