
/* cgeneric_kronecker.c
 *
 * Copyright (C) 2024-2026 Elias T Krainski
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

#include <assert.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include "INLAtools.h"

typedef struct {
	inla_cgeneric_data_tp *dataM1;
	inla_cgeneric_data_tp *dataM2;
#if !defined(INLA_EXTERNAL_PACKAGES)
	void *handle1;
	void *handle2;
#endif
	inla_cgeneric_func_tp *model1_func;
	inla_cgeneric_func_tp *model2_func;
	int nth1;
} cache_tp;


#if defined(INLA_EXTERNAL_PACKAGES)
// Force the compiler to keep this symbol even with aggressive LTO enabled
__attribute__((used)) __attribute__((visibility("default")))
#       if defined(__cplusplus)
extern "C"
#       endif
#endif
double *inla_cgeneric_kronecker(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// concatenated data approach of the lists

	// Merge cgeneric data from model1 (M1) and model2 (M2)
	// to define a Kronecker product model (KM) with precision as
	// Q = Q1 (x) Q2
	// where
	// Q1 is the precision for M1
	// Q2 is the precision for M2

	// cmd: length 1 string
	// theta: {theta1, theta[d12cache->nth1]}
	// data: {data1, data2}, but data1->ints[0]->ints[0...10] contains
	// n1, ni1, nd1, nc1, nm1, nsm1, m1,
	// ni2, nd2, nc2, nm2, nsm2, m2, n, M}
	// n1 is n for mode1 (M1)
	// ni1 is the number of ints for M1
	// ...
	// m1 is M for M1
	// ni2 is the number of ints for M2
	// ...
	// n is n for KM
	// M is M for KM
	// so that
	// data->ints[<ni1] contain ints for M1
	// data->doubles[<nd1] contain doubles for M1
	// data->chars[<nc1] contain chars for M1
	// data->mats[<nm1] contain mats for M1
	// data->smats[<nsm1] contain smats for M1
	// and
	// data->ints[ni1 .. < ni1+ni2] contain ints for M2
	// data->doubles[nd1 .. < nd1+nd2] contain doubles for M2
	// data->chars[nc1 .. < nc1+nc2] contain chars for M2
	// data->mats[nm1 .. < nm1+nm2] contain mats for M2
	// data->smats[nsm1 .. < nsm1+nsm2] contain smats for M2
	// additionally:
	// data->ints[ni1+ni2] contain nu1 index
	// data->ints[ni1+ni2+1] contain nu2 index
	// data->smatrices[nsm1+nsm2] contains the graph
	// where ->x is the order

	double *ret1 = NULL;				       // to store output from M1.
	double *ret2 = NULL;				       // to store output from M2.
	double *ret = NULL;				       // to return;

	int M1, M2, n, M;
	int ni1, nd1, nc1, nm1, nsm1;
	int ni2, nd2, nc2, nm2, nsm2;

	// int n1 = data->ints[0]->ints[0];
	ni1 = data->ints[0]->ints[1];
	nd1 = data->ints[0]->ints[2];
	nc1 = data->ints[0]->ints[3];
	nm1 = data->ints[0]->ints[4];
	nsm1 = data->ints[0]->ints[5];
	M1 = data->ints[0]->ints[6];

	// int n2 = data->ints[ni1]->ints[0];
	ni2 = data->ints[0]->ints[7];
	nd2 = data->ints[0]->ints[8];
	nc2 = data->ints[0]->ints[9];
	nm2 = data->ints[0]->ints[10];
	nsm2 = data->ints[0]->ints[11];
	M2 = data->ints[0]->ints[12];
	n = data->ints[0]->ints[13];
	M = data->ints[0]->ints[14];

	assert(ni1 > 1);
	assert(ni2 > 1);
	assert(nc1 > 1);
	assert(nc2 > 1);
	assert(data->n_chars > 5);

	if (!(data->cache)) {
#pragma omp critical (Name_5bd4b7198feb5550e84446518f90d47072338c18)
		if (!(data->cache)) {
			assert(!strcasecmp(data->ints[ni1 + ni2]->name, "idx1u"));
			assert(!strcasecmp(data->ints[ni1 + ni2 + 1]->name, "idx2u"));
			assert(!strcasecmp(data->smats[nsm1 + nsm2]->name, "Kgraph"));

			cache_tp *d12cache = Calloc(1, cache_tp);
			d12cache->dataM1 = Calloc(1, inla_cgeneric_data_tp);
			d12cache->dataM2 = Calloc(1, inla_cgeneric_data_tp);

			d12cache->dataM1->n_ints = ni1;
			d12cache->dataM1->ints = &data->ints[0];
			d12cache->dataM1->n_doubles = nd1;
			if (nd1 > 0) {
				d12cache->dataM1->doubles = &data->doubles[0];
			}
			d12cache->dataM1->n_chars = nc1;
			if (nc1 > 0) {
				d12cache->dataM1->chars = &data->chars[2];	// first two is for KM!
			}
			d12cache->dataM1->n_mats = nm1;
			if (nm1 > 0) {
				d12cache->dataM1->mats = &data->mats[0];
			}
			d12cache->dataM1->n_smats = nsm1;
			if (nsm1 > 0) {
				d12cache->dataM1->smats = &data->smats[0];
			}

			d12cache->dataM2->n_ints = ni2;
			d12cache->dataM2->ints = &data->ints[ni1];
			d12cache->dataM2->n_doubles = nd2;
			if (nd2 > 0) {
				d12cache->dataM2->doubles = &data->doubles[nd1];
			}
			d12cache->dataM2->n_chars = nc2;
			if (nc2 > 0) {
				d12cache->dataM2->chars = &data->chars[2 + nc1];	// first two is for KM!
			}
			d12cache->dataM2->n_mats = nm2;
			if (nm2 > 0) {
				d12cache->dataM2->mats = &data->mats[nm1];
			}
			d12cache->dataM2->n_smats = nsm2;
			if (nsm2 > 0) {
				d12cache->dataM2->smats = &data->smats[nsm1];
			}
#if !defined(INLA_EXTERNAL_PACKAGES)
			d12cache->handle1 = dlopen(&d12cache->dataM1->chars[1]->chars[0], RTLD_LAZY);
			if (!d12cache->handle1) {
				Rf_error("Failed to load shared library '%s': %s", &d12cache->dataM1->chars[1]->chars[0], dlerror());
				exit(1);
			}
			if (!strcmp(&d12cache->dataM1->chars[1]->chars[0], &d12cache->dataM2->chars[1]->chars[0])) {
				d12cache->handle2 = dlopen(&d12cache->dataM2->chars[1]->chars[0], RTLD_LAZY);
				if (!d12cache->handle2) {
					Rf_error("Failed to load shared library '%s': %s", &d12cache->dataM2->chars[0]->chars[0], dlerror());
					exit(1);
				}
			} else {
				d12cache->handle2 = d12cache->handle1;

			}
			*(void **)(&d12cache->model1_func) = dlsym(d12cache->handle1, &d12cache->dataM1->chars[0]->chars[0]);
//			const char *error = NULL;
			if (!d12cache->model1_func) {
				Rf_error("Fail to get %s\n%s\n", &d12cache->dataM1->chars[0]->chars[0], dlerror());
				exit(1);
			}
			*(void **)(&d12cache->model2_func) = dlsym(d12cache->handle2, &d12cache->dataM2->chars[0]->chars[0]);
			if (!d12cache->model2_func) {
			  Rf_error("Fail to get %s\n%s\n", &d12cache->dataM2->chars[0]->chars[0], dlerror());
				exit(1);
			}
#else
			d12cache->model1_func = (inla_cgeneric_func_tp *) inla_cgeneric_mapper(&d12cache->dataM1->chars[0]->chars[0]);
			d12cache->model2_func = (inla_cgeneric_func_tp *) inla_cgeneric_mapper(&d12cache->dataM2->chars[0]->chars[0]);
			assert(d12cache->model1_func && "model1_func not found");
			assert(d12cache->model2_func && "model2_func not found");

#endif
			double *ret = d12cache->model1_func(INLA_CGENERIC_INITIAL, NULL, d12cache->dataM1);
			d12cache->nth1 = (int) ret[0];
			Free(ret);

			data->cache = (void *) d12cache;
		}
	}

	assert(data->cache);
	cache_tp *d12cache = (cache_tp *) data->cache;

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{

		assert(M == data->smats[nsm1 + nsm2]->n);

		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		ret[0] = n;
		ret[1] = M;

		for (int i = 0; i < M; i++) {
			ret[2 + i] = data->smats[nsm1 + nsm2]->i[i];
		}
		for (int i = 0; i < M; i++) {
			ret[2 + M + i] = data->smats[nsm1 + nsm2]->j[i];
		}

		break;
	}

	case INLA_CGENERIC_Q:
	{
		ret = Calloc(2 + M, double);
		assert(ret);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;

		ret1 = d12cache->model1_func(INLA_CGENERIC_Q, &theta[0], d12cache->dataM1);
		ret2 = d12cache->model2_func(INLA_CGENERIC_Q, &theta[d12cache->nth1], d12cache->dataM2);

		int nu1 = data->ints[ni1 + ni2]->len;
		int nu2 = data->ints[ni1 + ni2 + 1]->len;

		double retE[M];
		double daux;
		int ox;

		int k = 0;
		for (int i = 0; i < M1; i++) {
			daux = ret1[2 + i];

			if (1) {
				double *to = retE + k;
				double *from = ret2;
#pragma omp simd
				for (int j = 0; j < M2; j++) {
					to[j] = daux * from[j];
				}
				k += M2;
			} else {
				for (int j = 0; j < M2; j++) {
					retE[k++] = daux * ret2[2 + j];
				}
			}
		}

		if ((nu1 > 0) & (nu2 > 0)) {
			for (int i = 0; i < nu1; i++) {
				daux = ret1[2 + data->ints[ni1 + ni2]->ints[i]];
				for (int j = 0; j < nu2; j++) {
					retE[k + j] = daux * ret2[2 + data->ints[ni1 + ni2 + 1]->ints[j]];
				}
				k += nu2;
			}
		}
		// ==============> does this works IF nu1==0 or nu2==0 ????
		assert(k == data->smats[nsm1 + nsm2]->n);

		for (k = 0; k < data->smats[nsm1 + nsm2]->n; k++) {
			ox = (int) data->smats[nsm1 + nsm2]->x[k];
			ret[2 + k] = retE[ox];
		}

		break;
	}

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
		break;
	}

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret1 = d12cache->model1_func(INLA_CGENERIC_INITIAL, NULL, d12cache->dataM1);
		ret2 = d12cache->model2_func(INLA_CGENERIC_INITIAL, NULL, d12cache->dataM2);

		int nth1 = (int) ret1[0];
		int nth2 = (int) ret2[0];
		ret = Calloc(1 + nth1 + nth2, double);
		assert(ret);
		ret[0] = nth1 + nth2;

		if (1) {
			Memcopy(ret + 1, ret1 + 1, nth1, double);
			Memcopy(ret + nth1 + 1, ret2 + 1, nth2, double);
		} else {
			for (int i = 0; i < nth1; i++) {
				ret[1 + i] = ret1[1 + i];
			}
			for (int i = 0; i < nth2; i++) {
				ret[1 + nth1 + i] = ret2[1 + i];
			}
		}

		break;
	}

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		break;
	}

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)
		ret1 = d12cache->model1_func(INLA_CGENERIC_LOG_PRIOR, &theta[0], d12cache->dataM1);
		ret2 = d12cache->model2_func(INLA_CGENERIC_LOG_PRIOR, &theta[d12cache->nth1], d12cache->dataM2);

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = ret1[0] + ret2[0];
		break;
	}

	case INLA_CGENERIC_QUIT:
	{

#if   !defined(INLA_EXTERNAL_PACKAGES)
		dlclose(d12cache->handle1);
		if (d12cache->handle1 != d12cache->handle2) {
			dlclose(d12cache->handle2);
		}
#endif
		// ==============> ?????
		// Free(d12cache);
		Free(data->cache);
	}
	default:
		break;
	}

	// strictly speaking, free(NULL) is undefined, so you have to do it properly: see macro on top
	Free(ret1);
	Free(ret2);

	return (ret);
}
