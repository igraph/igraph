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

struct ind_val{
	unsigned int ind;
	igraph_real_t val;
};
#define INDVAL struct ind_val
int compare_ind_val(const void *a, const void *b);

struct groups{
	unsigned int ind;
	unsigned int n;
	unsigned int* gr;
};
#define GROUPS struct groups

/*-------------------------------------------------
------------DEFINED IN scg_approximate_methods.c---
---------------------------------------------------*/	
int breaks_computation(const igraph_real_t *v,const unsigned int n, igraph_real_t *breaks,
			const unsigned int nb,const unsigned int method);
int intervals_plus_kmeans(const igraph_real_t *v, unsigned int *gr, const unsigned int n,
			  const unsigned int n_interv, const unsigned int maxiter);						
void intervals_method(const igraph_real_t *v, unsigned int *gr, const unsigned int n, const unsigned int n_interv);
/*-------------------------------------------------
------------DEFINED IN scg_optimal_method.c--------
---------------------------------------------------*/	
void cost_matrix(igraph_real_t *Cv, const INDVAL *vs, const unsigned int n, const unsigned int matrix, const igraph_real_t *ps);
igraph_real_t optimal_partition(const igraph_real_t *v, unsigned int *gr,const unsigned int n,
				const unsigned int nt,const unsigned int matrix, const igraph_real_t *p);
/*-------------------------------------------------
------------DEFINED IN scg_grouping.c--------------
---------------------------------------------------*/							
int grouping(igraph_real_t **v, unsigned int *gr, const unsigned int n, const unsigned int *nt, const unsigned int nev,
	     const unsigned int matrix, const igraph_real_t *p, const unsigned int algo, const unsigned int maxiter);
/*-------------------------------------------------
------------DEFINED IN scg_kmeans.c----------------
---------------------------------------------------*/
int kmeans_Lloyd(const igraph_real_t *x, const unsigned int n, const unsigned int p, igraph_real_t *cen,
		 const unsigned int k, int *cl, const unsigned int maxiter);					
/*-------------------------------------------------
------------DEFINED IN scg_exact_scg.c-------------
---------------------------------------------------*/
void exact_coarse_graining(const igraph_real_t *v, unsigned int *gr, const unsigned int n);				
/*-------------------------------------------------
------------DEFINED IN scg_utils.c-----------------
---------------------------------------------------*/	
int compare_groups(const void *a,const void *b);
int compare_real(const void *a, const void *b);
int compare_int(const void *a, const void *b);

igraph_real_t *real_sym_matrix(const unsigned int size);
#define real_sym_mat_get(S,i,j) S[i+j*(j+1)/2]
#define real_sym_mat_set(S,i,j,val) S[i+j*(j+1)/2] = val
#define free_real_sym_matrix(S) igraph_Free(S)

igraph_real_t **real_matrix(const unsigned int nrow, const unsigned int ncol);
void free_real_matrix(igraph_real_t **M,const unsigned int nrow);

unsigned int **uint_matrix(const unsigned int nrow, const unsigned int ncol);
void free_uint_matrix(unsigned int **M, const unsigned int nrow);

igraph_real_t *real_vector(const unsigned int n);
igraph_real_t min_real_vector(const igraph_real_t *v, const unsigned int n);
igraph_real_t max_real_vector(const igraph_real_t *v, const unsigned int n);
//for unity
#define free_real_vector(v) igraph_Free(v)

unsigned int *uint_vector(const unsigned int n);
unsigned int min_uint_vector(const unsigned int *v, const unsigned int n);
unsigned int max_uint_vector(const unsigned int *v, const unsigned int n);
//for unity
#define free_uint_vector(v) igraph_Free(v)

#endif
