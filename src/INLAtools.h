
/* INLAtools.h
 *
 * Copyright (C) 2025 Elias T Krainski
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
#if defined(INLA_EXTERNAL_PACKAGE)
#include <ltdl.h>
#else 
#include <dlfcn.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>	// needed to allow user interrupts
#include <R_ext/Utils.h>	// needed to allow user interrupts
#endif
#include "cgeneric.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#if !defined(Calloc)
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#endif
#define SQR(x) ((x)*(x))
#define pow2(x) ((x)*(x))
#define pow3(x) (pow2(x)*(x))
#define pow4(x) (pow2(x)*pow2(x))

#if !defined(iszero)
#ifdef __SUPPORT_SNAN__
#define iszero(x) (fpclassify(x) == FP_ZERO)
#else
#define iszero(x) (((__typeof(x))(x)) == 0)
#endif
#endif

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif

#define F_ONE ((fortran_charlen_t)1)

#if !defined(INLA_EXTERNAL_PACKAGE)
SEXP inla_cgeneric_element_get(SEXP Rcmd, SEXP Stheta, SEXP Sntheta, SEXP ints,
			       SEXP doubles, SEXP chars, SEXP mats, SEXP smats);
#endif
inla_cgeneric_func_tp inla_cgeneric_generic0;
inla_cgeneric_func_tp inla_cgeneric_kronecker;

