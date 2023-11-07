/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_MIXING_H
#define IGRAPH_MIXING_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

__BEGIN_DECLS

IGRAPH_EXPORT igraph_error_t igraph_assortativity_nominal(const igraph_t *graph,
                                               const igraph_vector_int_t *types,
                                               igraph_real_t *res,
                                               igraph_bool_t directed,
                                               igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_assortativity(const igraph_t *graph,
                                       const igraph_vector_t *values,
                                       const igraph_vector_t *values_in,
                                       igraph_real_t *res,
                                       igraph_bool_t directed,
                                       igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_assortativity_degree(const igraph_t *graph,
                                              igraph_real_t *res,
                                              igraph_bool_t directed);

IGRAPH_EXPORT igraph_error_t igraph_joint_degree_matrix(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *jdm,
        igraph_integer_t dout, igraph_integer_t din);

IGRAPH_EXPORT igraph_error_t igraph_joint_degree_distribution(
        const igraph_t *graph, const igraph_vector_t *weights, igraph_matrix_t *p,
        igraph_neimode_t from_mode, igraph_neimode_t to_mode,
        igraph_bool_t directed_neighbors,
        igraph_bool_t normalized,
        igraph_integer_t max_from_degree, igraph_integer_t max_to_degree);

IGRAPH_EXPORT igraph_error_t igraph_joint_type_distribution(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *p,
        const igraph_vector_int_t *from_types, const igraph_vector_int_t *to_types,
        igraph_bool_t directed, igraph_bool_t normalized);

__END_DECLS

#endif
