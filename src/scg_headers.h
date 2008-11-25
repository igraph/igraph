/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *	as described in the paper: Shrinking Matrices while preserving their
 *	eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *	Preprint available at <http://people.epfl.ch/david.morton>
 *  
 *	Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *	-----------
 *    This file contains the headers of the library SCGlib.
 *	  For use with R software <http://www.r-project.org/> define
 *	  the constant R_COMPIL and refer to the R documentation to compile
 *	  a dynamic library. The scg_r_wrapper function should be useful. 
 */

#ifndef SCG_HEADERS_H
#define SCG_HEADERS_H

#define R_COMPIL

#ifdef R_COMPIL
	#include <R.h>
	#define CALLOC calloc
	#define FREE free
	/*This does not work. Standard memory allocation used instead
	#define CALLOC(n,t) (t *) R_chk_calloc( (size_t) (n), sizeof(t) )
	#define FREE(p) (R_chk_free( (void *)(p) ), (p) = NULL)*/
	#define INT int
	#define UINT unsigned int
	#define REAL double
	#define FIRST_GROUP_NB 1
#else
	#include <stdio.h>
	#include <stdlib.h>
	#include <float.h> //for LDBL_MAX in kmeans.c
	#define CALLOC calloc
	#define FREE free
	#define INT long int
	#define UINT unsigned long int
	#define REAL double
	#define FIRST_GROUP_NB 0
#endif

struct ind_val{
	UINT ind;
	REAL val;
};
#define INDVAL struct ind_val
int compare_ind_val(const void *a, const void *b);

struct groups{
	UINT ind;
	UINT n;
	UINT* gr;
};
#define GROUPS struct groups

/*-------------------------------------------------
------------DEFINED IN scg_approximate_methods.c---
---------------------------------------------------*/	
void breaks_computation(const REAL *v,const UINT n, REAL *breaks,
						const UINT nb,const UINT method);
INT intervals_plus_kmeans(const REAL *v, UINT *gr, const UINT n,
							const UINT n_interv, const UINT maxiter);						
void intervals_method(const REAL *v, UINT *gr, const UINT n, const UINT n_interv);
/*-------------------------------------------------
------------DEFINED IN scg_optimal_method.c--------
---------------------------------------------------*/	
void cost_matrix(REAL *Cv, const INDVAL *vs, const UINT n, const UINT matrix, const REAL *ps);
REAL optimal_partition(const REAL *v, UINT *gr,const UINT n,
						const UINT nt,const UINT matrix, const REAL *p);
/*-------------------------------------------------
------------DEFINED IN scg_grouping.c--------------
---------------------------------------------------*/							
void grouping(REAL **v, UINT *gr, const UINT n, const UINT *nt, const UINT nev,
			const UINT matrix, const REAL *p, const UINT algo, const UINT maxiter);
/*-------------------------------------------------
------------DEFINED IN scg_kmeans.c----------------
---------------------------------------------------*/
INT kmeans_Lloyd(const REAL *x, const UINT n, const UINT p, REAL *cen,
						const UINT k, INT *cl, const UINT maxiter);					
/*-------------------------------------------------
------------DEFINED IN scg_exact_scg.c-------------
---------------------------------------------------*/
void exact_coarse_graining(const REAL *v, UINT *gr, const UINT n);				
/*-------------------------------------------------
------------DEFINED IN scg_utils.c-----------------
---------------------------------------------------*/	
int compare_groups(const void *a,const void *b);
int compare_real(const void *a, const void *b);
int compare_int(const void *a, const void *b);

REAL *real_sym_matrix(const UINT size);
#define real_sym_mat_get(S,i,j) S[i+j*(j+1)/2]
#define real_sym_mat_set(S,i,j,val) S[i+j*(j+1)/2] = val
#define free_real_sym_matrix(S) FREE(S)

REAL **real_matrix(const UINT nrow, const UINT ncol);
void free_real_matrix(REAL **M,const UINT nrow);

UINT **uint_matrix(const UINT nrow, const UINT ncol);
void free_uint_matrix(UINT **M, const UINT nrow);

REAL *real_vector(const UINT n);
REAL min_real_vector(const REAL *v, const UINT n);
REAL max_real_vector(const REAL *v, const UINT n);
//for unity
#define free_real_vector(v) FREE(v)

UINT *uint_vector(const UINT n);
UINT min_uint_vector(const UINT *v, const UINT n);
UINT max_uint_vector(const UINT *v, const UINT n);
//for unity
#define free_uint_vector(v) FREE(v)

#ifndef R_COMPIL
	void error(const char error_text[]);
	void warning(const char error_text[]);
#endif

#ifdef R_COMPIL
	void scg_r_wrapper(double *v, int *gr, int *n, int *nt, int *nev,
						int *nmatrix, int *nalgo, double *p,int *maxiter);
#endif

#endif
