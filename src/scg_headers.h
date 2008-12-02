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

#include "igraph.h"
#include "memory.h"

#ifdef USING_R
#  define FIRST_GROUP_NB 1
#else
#  define FIRST_GROUP_NB 0
#endif

typedef struct ind_val {
	unsigned int ind;
	igraph_real_t val;
} igraph_i_scg_indval_t;

int igraph_i_scg_compare_ind_val(const void *a, const void *b);

typedef struct groups {
	unsigned int ind;
	unsigned int n;
	unsigned int* gr;
} igraph_i_scg_groups_t;

/*-------------------------------------------------
------------DEFINED IN scg_approximate_methods.c---
---------------------------------------------------*/	
int igraph_i_scg_breaks_computation(const igraph_real_t *v,const unsigned int n, igraph_real_t *breaks,
				    const unsigned int nb,const unsigned int method);
int igraph_i_scg_intervals_plus_kmeans(const igraph_real_t *v, unsigned int *gr, const unsigned int n,
			  const unsigned int n_interv, const unsigned int maxiter);						
void igraph_i_scg_intervals_method(const igraph_real_t *v, unsigned int *gr, const unsigned int n, const unsigned int n_interv);
/*-------------------------------------------------
------------DEFINED IN scg_optimal_method.c--------
---------------------------------------------------*/	
void igraph_i_scg_cost_matrix(igraph_real_t *Cv, const igraph_i_scg_indval_t *vs, const unsigned int n, const unsigned int matrix, const igraph_real_t *ps);
igraph_real_t igraph_i_scg_optimal_partition(const igraph_real_t *v, unsigned int *gr,const unsigned int n,
					     const unsigned int nt,const unsigned int matrix, const igraph_real_t *p);
/*-------------------------------------------------
------------DEFINED IN scg_kmeans.c----------------
---------------------------------------------------*/
int igraph_i_scg_kmeans_Lloyd(const igraph_real_t *x, const unsigned int n, const unsigned int p, igraph_real_t *cen,
			      const unsigned int k, int *cl, const unsigned int maxiter);					
/*-------------------------------------------------
------------DEFINED IN scg_exact_scg.c-------------
---------------------------------------------------*/
void igraph_i_scg_exact_coarse_graining(const igraph_real_t *v, unsigned int *gr, const unsigned int n);				
/*-------------------------------------------------
------------DEFINED IN scg_utils.c-----------------
---------------------------------------------------*/	
int igraph_i_scg_compare_groups(const void *a,const void *b);
int igraph_i_scg_compare_real(const void *a, const void *b);
int igraph_i_scg_compare_int(const void *a, const void *b);

igraph_real_t *igraph_real_sym_matrix(const unsigned int size);
#define igraph_real_sym_mat_get(S,i,j) S[i+j*(j+1)/2]
#define igraph_real_sym_mat_set(S,i,j,val) S[i+j*(j+1)/2] = val
#define igraph_free_real_sym_matrix(S) igraph_Free(S)

igraph_real_t **igraph_real_matrix(const unsigned int nrow, const unsigned int ncol);
void igraph_free_real_matrix(igraph_real_t **M,const unsigned int nrow);

unsigned int **igraph_uint_matrix(const unsigned int nrow, const unsigned int ncol);
void igraph_free_uint_matrix(unsigned int **M, const unsigned int nrow);

igraph_real_t *igraph_real_vector(const unsigned int n);
igraph_real_t igraph_min_real_vector(const igraph_real_t *v, const unsigned int n);
igraph_real_t igraph_max_real_vector(const igraph_real_t *v, const unsigned int n);
//for unity
#define igraph_free_real_vector(v) igraph_Free(v)

#endif
