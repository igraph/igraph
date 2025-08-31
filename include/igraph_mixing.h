/*
   igraph library.
   Copyright (C) 2009-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_MIXING_H
#define IGRAPH_MIXING_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

IGRAPH_BEGIN_C_DECLS

IGRAPH_EXPORT igraph_error_t igraph_assortativity_nominal(
        const igraph_t *graph, const igraph_vector_t *weights,
        const igraph_vector_int_t *types,
        igraph_real_t *res,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_assortativity(
        const igraph_t *graph, const igraph_vector_t *weights,
        const igraph_vector_t *values, const igraph_vector_t *values_in,
        igraph_real_t *res,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_EXPORT igraph_error_t igraph_assortativity_degree(const igraph_t *graph,
                                              igraph_real_t *res,
                                              igraph_bool_t directed);

IGRAPH_EXPORT igraph_error_t igraph_joint_degree_matrix(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *jdm,
        igraph_int_t dout, igraph_int_t din);

IGRAPH_EXPORT igraph_error_t igraph_joint_degree_distribution(
        const igraph_t *graph, const igraph_vector_t *weights, igraph_matrix_t *p,
        igraph_neimode_t from_mode, igraph_neimode_t to_mode,
        igraph_bool_t directed_neighbors,
        igraph_bool_t normalized,
        igraph_int_t max_from_degree, igraph_int_t max_to_degree);

IGRAPH_EXPORT igraph_error_t igraph_joint_type_distribution(
        const igraph_t *graph, const igraph_vector_t *weights,
        igraph_matrix_t *p,
        const igraph_vector_int_t *from_types, const igraph_vector_int_t *to_types,
        igraph_bool_t directed, igraph_bool_t normalized);

IGRAPH_END_C_DECLS

#endif
