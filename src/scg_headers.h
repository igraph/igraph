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

#include <stdio.h>
#include <stdlib.h>
#include <float.h> //for LDBL_MAX in kmeans.c
#define CALLOC calloc
#define FREE free
#define INT int
#define UINT int
#define REAL double
#define FIRST_GROUP_NB 0

struct ind_val{
	UINT ind;
	REAL val;
};
#define INDVAL struct ind_val
int igraph_i_compare_ind_val(const void *a, const void *b);

struct groups{
	UINT ind;
	UINT n;
	UINT* gr;
};
#define GROUPS struct groups

/*-------------------------------------------------
------------DEFINED IN scg_approximate_methods.c---
---------------------------------------------------*/	
int igraph_i_breaks_computation(const REAL *v,const UINT n, REAL *breaks,
				 const UINT nb,const UINT method);
INT igraph_i_intervals_plus_kmeans(const REAL *v, UINT *gr, const UINT n,
				   const UINT n_interv, const UINT maxiter);						
int igraph_i_intervals_method(const REAL *v, UINT *gr, const UINT n, const UINT n_interv);
/*-------------------------------------------------
------------DEFINED IN scg_optimal_method.c--------
---------------------------------------------------*/	
void igraph_i_cost_matrix(REAL *Cv, const INDVAL *vs, const UINT n, const UINT matrix, const REAL *ps);
REAL igraph_i_optimal_partition(const REAL *v, UINT *gr,const UINT n,
						const UINT nt,const UINT matrix, const REAL *p);

/*-------------------------------------------------
------------DEFINED IN scg_kmeans.c----------------
---------------------------------------------------*/
INT igraph_i_kmeans_Lloyd(const REAL *x, const UINT n, const UINT p, REAL *cen,
						const UINT k, INT *cl, const UINT maxiter);					
/*-------------------------------------------------
------------DEFINED IN scg_exact_scg.c-------------
---------------------------------------------------*/
int igraph_i_exact_coarse_graining(const REAL *v, UINT *gr, const UINT n);				
/*-------------------------------------------------
------------DEFINED IN scg_utils.c-----------------
---------------------------------------------------*/	
int igraph_i_compare_groups(const void *a,const void *b);
int igraph_i_compare_real(const void *a, const void *b);
int igraph_i_compare_int(const void *a, const void *b);

REAL *igraph_i_real_sym_matrix(const UINT size);
#define igraph_i_real_sym_mat_get(S,i,j) S[i+j*(j+1)/2]
#define igraph_i_real_sym_mat_set(S,i,j,val) S[i+j*(j+1)/2] = val
#define igraph_i_free_real_sym_matrix(S) FREE(S)

REAL **igraph_i_real_matrix(const UINT nrow, const UINT ncol);
void igraph_i_free_real_matrix(REAL **M,const UINT nrow);

UINT **igraph_i_uint_matrix(const UINT nrow, const UINT ncol);
void igraph_i_free_uint_matrix(UINT **M, const UINT nrow);

REAL *igraph_i_real_vector(const UINT n);
REAL igraph_i_min_real_vector(const REAL *v, const UINT n);
REAL igraph_i_max_real_vector(const REAL *v, const UINT n);
//for unity
#define igraph_i_free_real_vector(v) FREE(v)

UINT *igraph_i_uint_vector(const UINT n);
UINT igraph_i_min_uint_vector(const UINT *v, const UINT n);
UINT igraph_i_max_uint_vector(const UINT *v, const UINT n);
//for unity
#define igraph_i_free_uint_vector(v) FREE(v)

#endif
