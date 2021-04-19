/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#ifndef IGRAPH_EMBEDDING_H
#define IGRAPH_EMBEDDING_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_arpack.h"
#include "igraph_eigen.h"
#include "igraph_constants.h"

__BEGIN_DECLS

IGRAPH_EXPORT int igraph_adjacency_spectral_embedding(const igraph_t *graph,
                                                      igraph_integer_t no,
                                                      const igraph_vector_t *weights,
                                                      igraph_eigen_which_position_t which,
                                                      igraph_bool_t scaled,
                                                      igraph_matrix_t *X,
                                                      igraph_matrix_t *Y,
                                                      igraph_vector_t *D,
                                                      const igraph_vector_t *cvec,
                                                      igraph_arpack_options_t *options);

typedef enum {
    IGRAPH_EMBEDDING_D_A = 0,
    IGRAPH_EMBEDDING_I_DAD,
    IGRAPH_EMBEDDING_DAD,
    IGRAPH_EMBEDDING_OAP
} igraph_laplacian_spectral_embedding_type_t;

IGRAPH_EXPORT int igraph_laplacian_spectral_embedding(const igraph_t *graph,
                                                      igraph_integer_t no,
                                                      const igraph_vector_t *weights,
                                                      igraph_eigen_which_position_t which,
                                                      igraph_laplacian_spectral_embedding_type_t type,
                                                      igraph_bool_t scaled,
                                                      igraph_matrix_t *X,
                                                      igraph_matrix_t *Y,
                                                      igraph_vector_t *D,
                                                      igraph_arpack_options_t *options);

IGRAPH_EXPORT int igraph_dim_select(const igraph_vector_t *sv, igraph_integer_t *dim);

__END_DECLS

#endif
