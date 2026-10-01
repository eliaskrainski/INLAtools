
/* INLAtools.c

 * Copyright (C) 2026-2028 Elias T Krainski
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

#include "INLAtools.h"

void uQ2Uexpand(int *nidx, int *idx, double *x, double *xx) {
  for(int i=0; i<(*nidx); i++) {
    xx[i] = x[idx[i]];
  }
}

void MuQ2kroneckerU(int *n1, int *n2, int *M2, int *idx2expanded,
                       double *x1, double *x2, double *xx) {
  // given x1 as (all) the elements of a dense (square) n1-dimensional matrix M and
  // uQ2 as the upper side (with diagonal) elements of a sparse (square) matrix Q2
  // return the upper side of
  //   Q = M (x) Q2
  // all ordered by rows (so x1[0]=M[0,0], x1[1]=M[0,1]...)!
  int k2 = 0;
  double daux;

  // complete Q2
  int M2expanded = (*M2)*2 - (*n2);
  double x2expanded[M2expanded];
  uQ2Uexpand(&M2expanded, &idx2expanded[0], &x2[0], &x2expanded[0]);

  // loop upper triangle of x1
  for (int i1=0; i1 < (*n1); i1++) {
    // column j1 = i1
    daux = x1[(*n1)*i1+i1];
    double *to = xx + k2;
    double *from  = x2;
#ifdef _OPENMP
#pragma omp simd
#endif
    for(int k = 0; k < (*M2); k++) {
      to[k] = daux * from[k];
    }
    k2 += (*M2);
    // column j1>i1
    if(i1 < ((*n1)-1)) {
      for(int j1 = i1+1; j1 < (*n1); j1++) {
        daux = x1[(*n1)*i1+j1];
        double *to = xx + k2;
        double *from  = x2expanded;
#ifdef _OPENMP
#pragma omp simd
#endif
        for(int k = 0; k < M2expanded; k++) {
          to[k] = daux * from[k];
        }
        k2 += M2expanded;
      }
    }
  }
}

void addMuQ2kroneckerU(int *n1, int *n2, int *M2, int *idx2expanded,
                       double *x1, double *x2, double *xx) {
  // given x1 as (expanded) the elements of a dense (square) n1-dimensional matrix M and
  // uQ2 as the upper side (with diagonal) elements of a sparse (square) matrix Q2
  // return the upper side of
  //   Q = M (x) Q2
  // all ordered by rows (so x1[0]=M[0,0], x1[1]=M[0,1]...)!
  int k2 = 0;
  double daux;

  // complete Q2
  int M2expanded = (*M2)*2 - (*n2);
  double x2expanded[M2expanded];
  uQ2Uexpand(&M2expanded, &idx2expanded[0], &x2[0], &x2expanded[0]);

  // loop upper triangle of x1
  for (int i1=0; i1 < (*n1); i1++) {
    // column j1 = i1
    daux = x1[(*n1)*i1+i1];
    double *to = xx + k2;
    double *from  = x2;
#ifdef _OPENMP
#pragma omp simd
#endif
    for(int k = 0; k < (*M2); k++) {
      to[k] += daux * from[k];
    }
    k2 += (*M2);
    // column j1>i1
    if(i1 < ((*n1)-1)) {
      for(int j1 = i1+1; j1 < (*n1); j1++) {
        daux = x1[(*n1)*i1+j1];
        double *to = xx + k2;
        double *from  = x2expanded;
#ifdef _OPENMP
#pragma omp simd
#endif
        for(int k = 0; k < M2expanded; k++) {
          to[k] += daux * from[k];
        }
        k2 += M2expanded;
      }
    }
  }
}
