/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *  as described in the paper: Shrinking Matrices while preserving their
 *  eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *  Preprint available at <http://people.epfl.ch/david.morton>
 *
 *  Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
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
 *  -----------
 *    This file contains the headers of the library SCGlib.
 *    For use with R software <http://www.r-project.org/> define
 *    the constant R_COMPIL and refer to the R documentation to compile
 *    a dynamic library. The scg_r_wrapper function should be useful.
 */

#ifndef SCG_HEADERS_H
#define SCG_HEADERS_H

#include "igraph_types.h"
#include "igraph_vector.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct ind_val {
    igraph_integer_t ind;
    igraph_real_t val;
} igraph_i_scg_indval_t;

igraph_integer_t igraph_i_compare_ind_val(const void *a, const void *b);

typedef struct groups {
    igraph_integer_t ind;
    igraph_integer_t n;
    igraph_integer_t* gr;
} igraph_i_scg_groups_t;

/*-------------------------------------------------
------------DEFINED IN scg_approximate_methods.c---
---------------------------------------------------*/

igraph_integer_t igraph_i_breaks_computation(const igraph_vector_t *v,
                                igraph_vector_t *breaks, igraph_integer_t nb,
                                igraph_integer_t method);
igraph_integer_t igraph_i_intervals_plus_kmeans(const igraph_vector_t *v, igraph_integer_t *gr,
                                   igraph_integer_t n, igraph_integer_t n_interv,
                                   igraph_integer_t maxiter);
igraph_integer_t igraph_i_intervals_method(const igraph_vector_t *v, igraph_integer_t *gr,
                              igraph_integer_t n, igraph_integer_t n_interv);

/*-------------------------------------------------
------------DEFINED IN scg_optimal_method.c--------
---------------------------------------------------*/

igraph_integer_t igraph_i_cost_matrix(igraph_real_t *Cv, const igraph_i_scg_indval_t *vs,
                         igraph_integer_t n, igraph_integer_t matrix, const igraph_vector_t *ps);
igraph_integer_t igraph_i_optimal_partition(const igraph_real_t *v, igraph_integer_t *gr, igraph_integer_t n, igraph_integer_t nt,
                               igraph_integer_t matrix, const igraph_real_t *p,
                               igraph_real_t *value);

/*-------------------------------------------------
------------DEFINED IN scg_kmeans.c----------------
---------------------------------------------------*/

igraph_integer_t igraph_i_kmeans_Lloyd(const igraph_vector_t *x, igraph_integer_t n,
                          igraph_integer_t p, igraph_vector_t *centers,
                          igraph_integer_t k, igraph_integer_t *cl, igraph_integer_t maxiter);

/*-------------------------------------------------
------------DEFINED IN scg_exact_scg.c-------------
---------------------------------------------------*/

igraph_integer_t igraph_i_exact_coarse_graining(const igraph_real_t *v, igraph_integer_t *gr,
                                   igraph_integer_t n);

/*-------------------------------------------------
------------DEFINED IN scg_utils.c-----------------
---------------------------------------------------*/

igraph_integer_t igraph_i_compare_groups(const void *a, const void *b);
igraph_integer_t igraph_i_compare_real(const void *a, const void *b);
igraph_integer_t igraph_i_compare_int(const void *a, const void *b);

igraph_real_t *igraph_i_real_sym_matrix(igraph_integer_t size);
#define igraph_i_real_sym_mat_get(S,i,j) S[i+j*(j+1)/2]
#define igraph_i_real_sym_mat_set(S,i,j,val) S[i+j*(j+1)/2] = val
#define igraph_i_free_real_sym_matrix(S) igraph_Free(S)

#endif
