/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_arpack.h"
#include "igraph_lapack.h"
#include "igraph_sparsemat.h"

#ifndef IGRAPH_EIGEN_H
#define IGRAPH_EIGEN_H

#include "igraph_decls.h"

__BEGIN_DECLS

typedef enum { IGRAPH_EIGEN_AUTO = 0,
               IGRAPH_EIGEN_LAPACK,
               IGRAPH_EIGEN_ARPACK,
               IGRAPH_EIGEN_COMP_AUTO,
               IGRAPH_EIGEN_COMP_LAPACK,
               IGRAPH_EIGEN_COMP_ARPACK
             } igraph_eigen_algorithm_t;

typedef enum { IGRAPH_EIGEN_LM = 0,
               IGRAPH_EIGEN_SM, /* 1 */
               IGRAPH_EIGEN_LA, /* 2 */
               IGRAPH_EIGEN_SA, /* 3 */
               IGRAPH_EIGEN_BE, /* 4 */
               IGRAPH_EIGEN_LR, /* 5 */
               IGRAPH_EIGEN_SR, /* 6 */
               IGRAPH_EIGEN_LI, /* 7 */
               IGRAPH_EIGEN_SI, /* 8 */
               IGRAPH_EIGEN_ALL, /* 9 */
               IGRAPH_EIGEN_INTERVAL, /* 10 */
               IGRAPH_EIGEN_SELECT
             }  /* 11 */
igraph_eigen_which_position_t;

typedef struct igraph_eigen_which_t {
    igraph_eigen_which_position_t pos;
    int howmany;
    int il, iu;
    igraph_real_t vl, vu;
    int vestimate;
    igraph_lapack_dgeevx_balance_t balance;
} igraph_eigen_which_t;

DECLDIR int igraph_eigen_matrix_symmetric(const igraph_matrix_t *A,
        const igraph_sparsemat_t *sA,
        igraph_arpack_function_t *fun, int n,
        void *extra,
        igraph_eigen_algorithm_t algorithm,
        const igraph_eigen_which_t *which,
        igraph_arpack_options_t *options,
        igraph_arpack_storage_t *storage,
        igraph_vector_t *values,
        igraph_matrix_t *vectors);

DECLDIR int igraph_eigen_matrix(const igraph_matrix_t *A,
                                const igraph_sparsemat_t *sA,
                                igraph_arpack_function_t *fun, int n,
                                void *extra,
                                igraph_eigen_algorithm_t algorithm,
                                const igraph_eigen_which_t *which,
                                igraph_arpack_options_t *options,
                                igraph_arpack_storage_t *storage,
                                igraph_vector_complex_t *values,
                                igraph_matrix_complex_t *vectors);

DECLDIR int igraph_eigen_adjacency(const igraph_t *graph,
                                   igraph_eigen_algorithm_t algorithm,
                                   const igraph_eigen_which_t *which,
                                   igraph_arpack_options_t *options,
                                   igraph_arpack_storage_t *storage,
                                   igraph_vector_t *values,
                                   igraph_matrix_t *vectors,
                                   igraph_vector_complex_t *cmplxvalues,
                                   igraph_matrix_complex_t *cmplxvectors);

DECLDIR int igraph_eigen_laplacian(const igraph_t *graph,
                                   igraph_eigen_algorithm_t algorithm,
                                   const igraph_eigen_which_t *which,
                                   igraph_arpack_options_t *options,
                                   igraph_arpack_storage_t *storage,
                                   igraph_vector_t *values,
                                   igraph_matrix_t *vectors,
                                   igraph_vector_complex_t *cmplxvalues,
                                   igraph_matrix_complex_t *cmplxvectors);


__END_DECLS

#endif
