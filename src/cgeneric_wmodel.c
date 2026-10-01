/* cgeneric_wmodel.c
 */

#include <assert.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include "INLAtools.h"

typedef struct {
	inla_cgeneric_data_tp *dataMc;
#if defined(INLA_WITH_EXTERNAL_PACKAGES)
  lt_dlhandle handleMc;
#else
	void *handleMc;
#endif
	inla_cgeneric_func_tp *modelMc_func;
	int nthMc;
} cache_tp;


#if defined(INLA_WITH_EXTERNAL_PACKAGES)
// Force the compiler to keep this symbol even with aggressive LTO enabled
__attribute__((used)) __attribute__((visibility("default")))
#       if defined(__cplusplus)
extern "C"
#       endif
#endif

double *inla_cgeneric_wmodel_dev(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data) {
   return inla_cgeneric_wmodel(cmd, &theta[0], &data[0]);
}

double *inla_cgeneric_wmodel(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

  double *retMc = NULL;			// to store output from Mc.
  double *ret = NULL;				// to return;

  // cgeneric model to combine K weight vectors with
	// a cgeneric N dimension model Mc to define a W model with precision as
	//   Q = [ W (o) I_N ] (x) bdiag(Q_1, ..., Q_K) (x) [ W' (o) I_N ]
	// where
	// Q1, ..., Q_K are NxN precision matrices for Mc, each with different parameters
	// W is a KxK weights matrix
	//
	// Q_k: the precision matrix from the k-th Mc instance
	//  is build with q>=0 parameter(s)

	// constraints:
	//  c1: \sum_j W_{ij}^2 = 1
	//  c2: W_{ii}>0
	//  c3: order as exp(theta[oprm_1])>exp(theta[oprm_2])>...>exp(theta[oprm_K])

	// model parameters:
	//  theta[0, ..., K(K-1)-1, K(K-1), ..., K(K-1)q-1]
	//   K(K-1) to define W
	//   q*K to define Q_1, ..., Q_K

	// accounting for c1 and c2 we define
	//   x_{ij} = 1, if i == j, theta[K(i-1)+j] otherwise
	//   a_i = sum_j x_{ij}^2
	//   w_{ij} = x_{ij} / sqrt(a_i)

	// to have c3 we use
	// theta*[v]_K     = theta[v]_K
	// theta*[v]_{K-1} = log( exp(theta[v]_{K-1}) + exp(theta*[v]_K) )
	// theta*[v]_{K-2} = log( exp(theta[v]_{K-2}) + exp(theta*[v]_{K-1}) )
	// ...
	// theta*[v]_1 = log( exp(theta[v]_1) + exp(theta*[v]_2) )

	// cmd: length 1 string
	// theta: {nthW, theta[dCache->nthMc]}
	// data:
	// dataMc->ints[0]->ints[0...14] contains
	//  [nMc, niMc, ndMc, ncMc, nmMc, nsMc, Mc,
	//   K, niW, ndW, ncW, nmW, nsW, n, M]
	int nMc; // the size of the Mc model
	int McU;  // number of non-zeros in the upper side of Q_j
	int K;   // W model dimension
	int N;   // size of the combined model, equal K times nMc
	int M;   // n + number of non-zeros in the upper side of Q
	// so that
	//  data->ints[<niMc] contain ints for Mc
	//  data->doubles[<ndMc] contain doubles for Mc
	//  data->chars[<ncMc] contain chars for Mc
	//  data->mats[<nmMc] contain mats for Mc
	//  data->smats[<nsMc] contain smats for Mc
	nMc = data->ints[0]->ints[0];
	int niMc = data->ints[0]->ints[1];
	int ndMc = data->ints[0]->ints[2];
	int ncMc = data->ints[0]->ints[3];
	int nmMc = data->ints[0]->ints[4];
	int nsMc = data->ints[0]->ints[5];
	McU = data->ints[0]->ints[6];
	K = data->ints[0]->ints[7];  /*
	int niW = data->ints[0]->ints[8];
	int ndW = data->ints[0]->ints[9];
	int ncW = data->ints[0]->ints[10];
	int nmW = data->ints[0]->ints[11];
	int nsW = data->ints[0]->ints[12]; */
	N = data->ints[0]->ints[13];
	M = data->ints[0]->ints[14];
	assert(niMc > 1);
	assert(ncMc > 1);
	assert(data->n_chars > 3);
	int nthW = K*(K-1);

	/*
	printf("nMc: (%d) %d %d %d %d %d %d\nK  : (%d) %d %d %d %d %d N=%d M=%d\n%d\n",
        nMc, niMc, ndMc, ncMc, nmMc, nsMc, McU,
        K, niW, ndW, ncW, nmW, nsW, N, M, nthW);
	 */

	if (!(data->cache)) {
#ifdef _OPENMP
#pragma omp critical (Name_5bd4b7198feb5550e84446518f2c9f8c52c4b058)
#endif
		if (!(data->cache)) {
		  assert(!strcasecmp(data->ints[0]->name, "n"));
		  assert(!strcasecmp(data->ints[niMc]->name, "idParam"));
		  assert(!strcasecmp(data->ints[niMc+1]->name, "idxMc"));
		  assert((nMc+McU) == data->ints[niMc]->len);
		  assert(!strcasecmp(data->ints[niMc+2]->name, "ii"));
		  assert(!strcasecmp(data->ints[niMc+3]->name, "jj"));

			cache_tp *dCache = Calloc(1, cache_tp);
#if defined(INLA_WITH_EXTERNAL_PACKAGES)
			static int dCache->ltdl_cgwm = 0;
#endif
			dCache->dataMc = Calloc(1, inla_cgeneric_data_tp);

			dCache->dataMc->n_ints = niMc;
			dCache->dataMc->ints = &data->ints[0];
			dCache->dataMc->n_doubles = ndMc;
			if (ndMc > 0) {
				dCache->dataMc->doubles = &data->doubles[0];
			}
			dCache->dataMc->n_chars = ncMc;
			if (ncMc > 0) {
				dCache->dataMc->chars = &data->chars[2]; // first 2 are for Wmodel
			}
			dCache->dataMc->n_mats = nmMc;
			if (nmMc > 0) {
				dCache->dataMc->mats = &data->mats[0];
			}
			dCache->dataMc->n_smats = nsMc;
			if (nsMc > 0) {
				dCache->dataMc->smats = &data->smats[0];
			}

#if defined(INLA_WITH_EXTERNAL_PACKAGES)
			dCache->modelMc_func = (inla_cgeneric_func_tp *) inla_cgeneric_mapper(&dCache->dataMc->chars[0]->chars[0]);
			if(!dCache->modelMc_func) { // not in main INLA program, use from the shlib
			  if (dCache->ltdl_cgwm == 0) {
			    if(lt_dlinit() != 0) {
			      fprintf(stderr,"\n\n\t*** ERROR *** Failed to start libltdl:\n %s\n\n", lt_dlerror());
			      abort();
			    }
			  }
			  dCache->ltdl_cgwm = 1;
			  dCache->handleMc = lt_dlopen(&dCache->dataMc->chars[1]->chars[0]);
			  if (!dCache->handleMc) {
			    fprintf(stderr,"\n\n\t*** ERROR *** Failed to load shared library '%s':\n\n%s\n\n",
               &dCache->dataMc->chars[1]->chars[0], lt_dlerror());
			    abort();
			  }
			  *(void **)(&dCache->modelMc_func) =
			    lt_dlsym(dCache->handleMc, &dCache->dataMc->chars[0]->chars[0]);
			  if(!dCache->modelMc_func){
			    lt_dlclose(dCache->handleMc);
			  }
			}
			assert(dCache->modelMc_func && "modelMc_func not found");
#else
			if(dCache->dataMc->ints[1]->ints[0]) {
			  Rprintf("Mc shlib: %s\n", &dCache->dataMc->chars[1]->chars[0]);
			}
			dCache->handleMc = dlopen(&dCache->dataMc->chars[1]->chars[0], RTLD_LAZY);
			if (!dCache->handleMc) {
			  Rprintf("Mc shlib: %s\n", &dCache->dataMc->chars[1]->chars[0]);
			  Rf_error("Failed to load shared library '%s':\n %s\n",
              &dCache->dataMc->chars[1]->chars[0], dlerror());
				exit(1);
			} else {
			  if(dCache->dataMc->ints[1]->ints[0]) {
			    Rprintf("The shlib is loaded\n");
			  }
			}
			if(dCache->dataMc->ints[1]->ints[0]) {
			  Rprintf("Getting symbol  %s\n", &dCache->dataMc->chars[0]->chars[0]);
			}
			*(void **)(&dCache->modelMc_func) = dlsym(dCache->handleMc, &dCache->dataMc->chars[0]->chars[0]);
			//			const char *error = NULL;
			if (!dCache->modelMc_func) {
			  Rprintf("M1 symbol: %s\n", &dCache->dataMc->chars[0]->chars[0]);
			  Rprintf("M1 shlib: %s\n", &dCache->dataMc->chars[1]->chars[0]);
			  Rf_error("Fail to get %s\n%s\n", &dCache->dataMc->chars[0]->chars[0], dlerror());
			  exit(1);
			}

#endif
			// get the number of parameters of Mc
			double *ret = dCache->modelMc_func(INLA_CGENERIC_INITIAL, NULL, dCache->dataMc);
			dCache->nthMc = (int) ret[0];
			Free(ret);
			data->cache = (void *) dCache;
		}
	}

	assert(data->cache);
	cache_tp *dCache = (cache_tp *) data->cache;

	int i, j, k, kk, K2 = K*K;
	int info, ipiv[K];
	int nthMc = dCache->nthMc;
	int nparamsMc = nthMc * K;
	double thetaMc[nparamsMc+1]; // just to not have it as length 0
	double daux, xaux[K];
	double W[K2], Ww[K2], iW[K2]; // W, copy and W^{-1}
	if(theta) {
/*
	  printf("theta for W\n");
	  for(i=0; i<nthW; i++) {
	    printf("%d %2.6f \n", i, theta[i]);
	  }
	  printf("theta for Mc params\n");
	  for(i=0; i<nparamsMc; i++) {
	    printf("%d %2.6f \n", nthW+i, theta[nthW+i]);
	  } */

	  k = 0; kk = 0;
	  for(i=0; i<K; i++) {
	    xaux[i] = 1.0;
	    if(i>0) {
	      for(j=0; j<i; j++) {
	        xaux[j] = theta[k++];
	      }
	    }
	    if(i<K) {
	      for(j=i+1; j<K; j++) {
	        xaux[j] = theta[k++];
	      }
	    }
	    daux = 0.0;
	    for(j=0;j<K;j++) {
	      daux += SQR(xaux[j]);
	    }
	    daux = sqrt(daux);
	    for(j=0; j<K; j++) {
	      W[kk] = xaux[j]/daux;
	      Ww[kk] = W[kk]; // copy of W to use for W^{-1}
	      if(i==j) {
	        iW[kk++] = 1.0;
	      } else {
	        iW[kk++] = 0.0;
	      }
	    }
	  }

	  // compute W^{-1}
	  dgesv_(&K, &K, &Ww[0], &K, &ipiv[0], &iW[0], &K, &info, F_ONE);

//	  printMat(W, K, K, "W:\n");
	//  printMat(iW, K, K, "inverse of W:\n");

	  //  c3: order as exp(theta[oprm_1])>exp(theta[oprm_2])>...>exp(theta[oprm_K])
	  if(nthMc>0) {
	    daux = 0.0;
	    k = nthW + nparamsMc -1;
	    for(i=0; i<K; i++) {
	      daux += exp(theta[k]);
	      thetaMc[nparamsMc-i*nthMc-1] = log(daux);
	      k -= nthMc;
	    }
	    if(nthMc>1) {
	      k = 1;
	      for(i=0; i<K; i++){
	        for(j=1; j<nthMc; j++) {
	          thetaMc[k] = theta[nthW+k];
	          k++;
	        }
	        k++;
	      }
	    }
	  }
	  //printMat(thetaMc, K, nthMc, "thetaMc\n");

	}

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{

		assert(M == data->smats[nsMc]->n);

		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		ret[0] = N;
		ret[1] = M;

		for (i = 0; i < M; i++) {
			ret[2 + i] = data->ints[niMc+2]->ints[i];
		}
		for (i = 0; i < M; i++) {
			ret[2 + M + i] = data->ints[niMc+3]->ints[i];
		}

		break;
	}

	case INLA_CGENERIC_Q:
	{
		ret = Calloc(2 + M, double);
		assert(ret);
//		memset(ret+2, 0.0, M*sizeof(double));
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		double ret0[M];

		int i1, i2, ithetak=0;
		double iwk[K], iWk2[K*K];

		for(i=0; i<K; i++) {

		  // column iW[,k]
		  k = i;
		  for(i1=0; i1<K; i1++) {
		    iwk[i1] = iW[k];
		    k += K;
		  }
//		  printMat(iwk, 1, K, "iW[,k]:\n");

		  // compute iW[,k] iW[,k]'
		  kk = 0;
		  for(i1=0; i1<K; i1++) {
		      for(i2=0; i2<K; i2++) {
		        iWk2[kk++] = iwk[i1] * iwk[i2];
		    }
		  }
//		  printMat(iWk2, K, K, "iW[,k] iW[,k]':\n");

		  retMc = dCache->modelMc_func(INLA_CGENERIC_Q, &thetaMc[ithetak], dCache->dataMc);
		  ithetak += nthMc;

		  if(i==0) {
		    MuQ2kroneckerU(&K, &nMc, &McU, &dCache->dataMc->ints[niMc+1]->ints[0],
                     &iWk2[0], &retMc[2], &ret0[0]);
		  } else {
		    addMuQ2kroneckerU(&K, &nMc, &McU, &dCache->dataMc->ints[niMc+1]->ints[0],
                        &iWk2[0], &retMc[2], &ret0[0]);
		  }
		  Free(retMc);
		}

      for(i=0; i<M; i++) {
        ret[2+i] = ret0[dCache->dataMc->ints[niMc+4]->ints[i]];
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

		retMc = dCache->modelMc_func(INLA_CGENERIC_INITIAL, NULL, dCache->dataMc);

		ret = Calloc(1 + nthW + nparamsMc, double);
		assert(ret);
		ret[0] = nthW + nparamsMc;

		// initials for W
		for(i=0; i<nthW; i++) {
		  ret[1+i] = 1.0; // w[-i]/w[i] assuming all equal
		}

		// same initials for each instance
		if(dCache->nthMc>0) {
		  k = nthW+1;
		  for(i=0; i<K; i++) {
		    for(j=0; j<dCache->nthMc; j++)
		    ret[k++] = retMc[1+j];
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
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0.0;
		k = nthW;
		kk = 0;
		double lam, r;
		double dm = (double)(K-1); // dimension of the prior
		double w0ref = 1.0/SQR(K);
		for(i=0; i<K; i++) {
		  retMc = dCache->modelMc_func(INLA_CGENERIC_LOG_PRIOR, &thetaMc[k], dCache->dataMc);
		  k += nthMc;
		  ret[0] += retMc[0];
		  daux = 0;
		  for(j=0; j<K; j++) {
		    daux += SQR(W[kk] - w0ref);
		    kk++;
		  }
		  r = sqrt(daux);
		  lam = dCache->dataMc->doubles[ndMc]->doubles[k];
		  ret[0] += log(lam) - lam * r;
		  if(K>1) {
		    ret[0] += lgamma(dm * 0.5);
		    ret[0] -= (dm*0.5) * log(M_PI);
		    ret[0] -= (dm -1.0) * log(r);
		  }
		}
		break;
	}

	case INLA_CGENERIC_QUIT:
	{

#if defined(INLA_WITH_EXTERNAL_PACKAGES)
	  if(dCache->ltdl_cgwm>0) {
	    if(dCache->handleMc) {
	      lt_dlclose(dCache->handleMc);
	    }
	  }
#else
	  if(dCache->handleMc) {
	    dlclose(dCache->handleMc);
	  }
#endif
		// ==============> ?????
		// Free(dCache);
		Free(data->cache);
	}
	default:
		break;
	}

	// strictly speaking, free(NULL) is undefined, so you have to do it properly: see macro on top
	Free(retMc);

	return (ret);
}
