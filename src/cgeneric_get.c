
/* cgeneric_get.c
 *
 * Copyright (C) 2024 Elias T Krainski
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
 *        Elias Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include "INLAtools.h"

SEXP inla_cgeneric_element_get(SEXP Rcmd, SEXP Stheta, SEXP Sntheta, SEXP ints,
			       SEXP doubles, SEXP chars, SEXP mats, SEXP smats)
{

	int ni = 0, nd = 0, nc = 0, nm = 0, nsm = 0, nout = 0;
	int i, j, debug;
	int *ntheta = INTEGER(Sntheta);
	int nth = length(Stheta);
	if (ntheta[0] > 0) {
		nth /= ntheta[0];
	}
	// initial check: mandatory arguments
	if (isNewList(ints)) {
		ni = length(ints);
		if (ni < 2) {
			error("at least length 2 'ints' must be provided!");
		}
	} else {
		error("'ints' must be a list!");
	}

	if (isNewList(chars)) {
		nc = length(chars);
		if (nc < 2) {
			error("at least length 2 'chars' must be provided!");
		}
	} else {
		error("'chars' must be a list!");
	}

	// check the optional arguments
	if (isNull(doubles)) {
		nd = 0;
	} else {
		if (isNewList(doubles)) {
			nd = length(doubles);
		} else {
			error("'doubles' must be a list!");
		}
	}
	if (isNull(mats)) {
		nm = 0;
	} else {
		if (isNewList(mats)) {
			nm = length(mats);
		} else {
			error("'mats' must be a list");
		}
	}
	if (isNull(smats)) {
		nsm = 0;
	} else {
		if (isNewList(smats)) {
			nsm = length(smats);
		} else {
			error("'smats' must be a list");
		}
		for (i = 0; i < nsm; i++) {
			if (isNewList(VECTOR_ELT(smats, i))) {
				if (length(VECTOR_ELT(smats, i)) < 6) {
					error
					    ("'length(smats[[%d]])' must not be <6",
					     i + 1);
				}
			} else {
				error("'smats[[%d]]' must be a list", i + 1);
			}
		}
	}

	// get initial info
	char *CMD = (char *)CHAR(STRING_ELT(Rcmd, 0));
	debug = asInteger(VECTOR_ELT(ints, 1));
	if (debug > 0) {
		Rprintf("Rcmd is %s, debug = %d\n", CMD, debug);
		Rprintf("ni = %d, nd = %d, nc = %d, nm = %d, ", ni, nd, nc, nm);
		Rprintf("nsm = %d, ntheta = %d, nth = %d\n",
			nsm, ntheta[0], nth);
	}
	// check theta (can be a matrix, each column one instance)
	double *theta = NULL;
	if (!isNull(Stheta)) {
		if (!isReal(Stheta)) {
			Rprintf("Expected REAL for 'theta'!");
		}
		theta = REAL(Stheta);
		if (debug) {
			Rprintf("theta: ");
			if ((ntheta[0]) < 2) {
				for (i = 0; i < nth; i++) {
					Rprintf("%f ", theta[i]);
				}
				Rprintf("\n");
			} else {
				for (j = 0; j < (ntheta[0]); j++) {
					for (i = 0; i < nth; i++) {
						Rprintf("%f ", theta[i]);
					}
					Rprintf("\n");
				}
			}
		}
	}
	// define objects
	int ilen[ni], clen[nc];
	inla_cgeneric_data_tp *cgeneric_data = Calloc(1, inla_cgeneric_data_tp);
	char *cgeneric_shlib;
	char *cgeneric_model;
	int *iaux;
	double *daux, *ret = NULL;
	char *pcaux;
	const char *caux;

	// collect ints
	SEXP inames = PROTECT(getAttrib(ints, R_NamesSymbol));
	cgeneric_data->n_ints = ni;
	cgeneric_data->ints = Calloc(ni, inla_cgeneric_vec_tp *);
	for (i = 0; i < ni; i++) {
		ilen[i] = length(VECTOR_ELT(ints, i));
		if (debug > 0) {
			caux = CHAR(STRING_ELT(inames, i));
			Rprintf("length(ints[[%d]]), %s, is %d\n",
				i + 1, caux, ilen[i]);
		}
		cgeneric_data->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
		pcaux = (char *)CHAR(STRING_ELT(inames, i));
		cgeneric_data->ints[i]->name = pcaux;
		cgeneric_data->ints[i]->len = ilen[i];
		cgeneric_data->ints[i]->ints = INTEGER(VECTOR_ELT(ints, i));
	}
	UNPROTECT(1);

	if (nd > 0) {
		// collect doubles
		int dlen[nd];
		SEXP dnames = PROTECT(getAttrib(doubles, R_NamesSymbol));
		cgeneric_data->n_doubles = nd;
		cgeneric_data->doubles = Calloc(nd, inla_cgeneric_vec_tp *);
		for (i = 0; i < nd; i++) {
			dlen[i] = length(VECTOR_ELT(doubles, i));
			if (debug > 0) {
				caux = CHAR(STRING_ELT(dnames, i));
				Rprintf("length(doubles[[%d]]), %s, is %d\n",
					i + 1, caux, dlen[i]);
			}
			cgeneric_data->doubles[i] =
			    Calloc(1, inla_cgeneric_vec_tp);
			pcaux = (char *)CHAR(STRING_ELT(dnames, i));
			cgeneric_data->doubles[i]->name = pcaux;
			cgeneric_data->doubles[i]->len = dlen[i];
			cgeneric_data->doubles[i]->doubles =
			    REAL(VECTOR_ELT(doubles, i));
		}
		UNPROTECT(1);
	}
	// collect characters
	SEXP cnames = PROTECT(getAttrib(chars, R_NamesSymbol));
	cgeneric_data->n_chars = nc;
	cgeneric_data->chars = Calloc(nc, inla_cgeneric_vec_tp *);
	for (i = 0; i < nc; i++) {
		clen[i] = length(VECTOR_ELT(chars, i));
		caux = CHAR(STRING_ELT(cnames, i));
		if (debug > 0) {
			Rprintf("length(chars[[%d]]), %s, is %d\n",
				i + 1, caux, clen[i]);
		}
		cgeneric_data->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
		pcaux = (char *)CHAR(STRING_ELT(cnames, i));
		cgeneric_data->chars[i]->name = pcaux;
		cgeneric_data->chars[i]->len = clen[i];
		cgeneric_data->chars[i]->chars = Calloc(clen[i] + 1L, char);
		if (debug > 0) {
			Rprintf("%d: length(%s) is %d:\n", i + 1,
				cgeneric_data->chars[i]->name, clen[i]);
		}
		pcaux = (char *)CHAR(STRING_ELT(VECTOR_ELT(chars, i), 0));
		cgeneric_data->chars[i]->chars = pcaux;
		if (debug > 0) {
			Rprintf("%s \n", cgeneric_data->chars[i]->chars);
		}
	}
	UNPROTECT(1);

	// check the mandatory strings
	if (strcmp(cgeneric_data->chars[0]->name, "model") != 0)
		error("'chars[[1]]' name is not equal 'model'");
	cgeneric_model = cgeneric_data->chars[0]->chars;
	if (strcmp(cgeneric_data->chars[1]->name, "shlib") != 0)
		error("'chars[[2]]' name is not equal 'shlib'");
	cgeneric_shlib = cgeneric_data->chars[1]->chars;

	if (nm > 0) {
		// collect matrices
		int mnr[nm], mnc[nm];
		SEXP mnames = PROTECT(getAttrib(mats, R_NamesSymbol));
		cgeneric_data->n_mats = nm;
		cgeneric_data->mats = Calloc(nm, inla_cgeneric_mat_tp *);
		pcaux = (char *)CHAR(STRING_ELT(mnames, 0));
		for (i = 0; i < nm; i++) {
			daux = REAL(VECTOR_ELT(mats, i));
			mnr[i] = (int)daux[0];
			mnc[i] = (int)daux[1];
			if (debug > 0) {
				caux = CHAR(STRING_ELT(mnames, i));
				Rprintf("dim(mats[[%d]]), %s, is %d %d\n",
					i + 1, caux, mnr[i], mnc[i]);
			}
			cgeneric_data->mats[i] =
			    Calloc(1, inla_cgeneric_mat_tp);
			pcaux = (char *)CHAR(STRING_ELT(mnames, i));
			cgeneric_data->mats[i]->name = pcaux;
			cgeneric_data->mats[i]->nrow = mnr[i];
			cgeneric_data->mats[i]->ncol = mnc[i];
			cgeneric_data->mats[i]->x = &daux[2];
		}
		UNPROTECT(1);
	}

	if (nsm > 0) {
		// collect smatrices
		int smnr[nsm], smnc[nsm], smn[nsm];
		SEXP smnames = PROTECT(getAttrib(smats, R_NamesSymbol));
		cgeneric_data->n_smats = nsm;
		cgeneric_data->smats = Calloc(nsm, inla_cgeneric_smat_tp *);
		pcaux = (char *)CHAR(STRING_ELT(smnames, 0));
		for (i = 0; i < nsm; i++) {
			iaux = INTEGER(VECTOR_ELT(VECTOR_ELT(smats, i), 0));
			smnr[i] = iaux[0];
			iaux = INTEGER(VECTOR_ELT(VECTOR_ELT(smats, i), 1));
			smnc[i] = iaux[0];
			iaux = INTEGER(VECTOR_ELT(VECTOR_ELT(smats, i), 2));
			smn[i] = iaux[0];
			if (debug > 0) {
				caux = CHAR(STRING_ELT(smnames, i));
				Rprintf("smats[[%d]] %s: %d %d %d\n",
					i + 1, caux, smnr[i], smnc[i], smn[i]);
			}
			cgeneric_data->smats[i] =
			    Calloc(1, inla_cgeneric_smat_tp);
			pcaux = (char *)CHAR(STRING_ELT(smnames, i));
			cgeneric_data->smats[i]->name = pcaux;
			cgeneric_data->smats[i]->nrow = smnr[i];
			cgeneric_data->smats[i]->ncol = smnc[i];
			cgeneric_data->smats[i]->n = smn[i];
			cgeneric_data->smats[i]->i =
			    INTEGER(VECTOR_ELT(VECTOR_ELT(smats, i), 3));
			cgeneric_data->smats[i]->j =
			    INTEGER(VECTOR_ELT(VECTOR_ELT(smats, i), 4));
			cgeneric_data->smats[i]->x =
			    REAL(VECTOR_ELT(VECTOR_ELT(smats, i), 5));
		}
		UNPROTECT(1);
	}
	// load shlib
	void *handle;
	handle = dlopen(cgeneric_shlib, RTLD_LAZY);
	if (!handle) {
		Rf_error("Failed to load shared library '%s': %s",
			 cgeneric_shlib, dlerror());
	}
	inla_cgeneric_func_tp *model_func = NULL;
	*(void **)(&model_func) = dlsym(handle, cgeneric_model);
	SEXP Rret = R_NilValue;

	if (strcmp(CMD, "graph") == 0) {
		ret = model_func(INLA_CGENERIC_GRAPH, theta, cgeneric_data);
		nout = (int)ret[1];
		SEXP ii = PROTECT(allocVector(INTSXP, nout));
		iaux = INTEGER(ii);
		for (i = 0; i < nout; i++) {
			iaux[i] = ret[2 + i];
		}
		SEXP jj = PROTECT(allocVector(INTSXP, nout));
		iaux = INTEGER(jj);
		for (i = 0; i < nout; i++) {
			iaux[i] = ret[2 + nout + i];
		}
		Rret = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(Rret, 0, ii);
		SET_VECTOR_ELT(Rret, 1, jj);
		if (debug > 0) {
			Rprintf("graph with n = %d and %d nz\n",
				(int)ret[0], nout);
		}
		UNPROTECT(3);
	}

	if (strcmp(CMD, "Q") == 0) {
		ret = model_func(INLA_CGENERIC_Q, theta, cgeneric_data);
		nout = (int)ret[1];
		Rret = PROTECT(allocVector(REALSXP, nout));
		daux = REAL(Rret);
		for (i = 0; i < nout; i++) {
			daux[i] = ret[2 + i];
		}
		if (debug > 0) {
			Rprintf("Q with %d nz\n", nout);
		}
		UNPROTECT(1);
	}

	if (strcmp(CMD, "mu") == 0) {
		ret = model_func(INLA_CGENERIC_MU, theta, cgeneric_data);
		Rret = PROTECT(allocVector(REALSXP, 1));
		REAL(Rret)[0] = ret[0];
		UNPROTECT(1);
	}

	if (strcmp(CMD, "initial") == 0) {
		ret = model_func(INLA_CGENERIC_INITIAL, theta, cgeneric_data);
		nout = (int)ret[0];
		if (nout > 0) {
		  Rret = PROTECT(allocVector(REALSXP, nout));
		  daux = REAL(Rret);
		  for (i = 0; i < nout; i++) {
		    daux[i] = ret[1 + i];
		  }
		  if (debug > 0) {
		    Rprintf("intial: %d elements (%f)\n", nout, ret[0]);
				for (i = 0; i < nout; i++)
					Rprintf("theta[%d] = %f\n", i, ret[1 + i]);
		  }
			UNPROTECT(1);
		} else {
		  if (debug > 0) {
		    Rprintf("no intial\n");
		  }
		}
	}

	if (strcmp(CMD, "log_prior") == 0) {
		Rret = PROTECT(allocVector(REALSXP, ntheta[0]));
		for (j = 0; j < (ntheta[0]); j++) {
			ret = model_func(INLA_CGENERIC_LOG_PRIOR,
					 &theta[j * nth], cgeneric_data);
			REAL(Rret)[j] = ret[0];
		}
		UNPROTECT(1);
	}

	dlclose(handle);
	free(cgeneric_data);
	free(ret);

	return Rret;
}
