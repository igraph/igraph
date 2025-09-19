/*
   igraph library.
   Copyright (C) 2009-2024  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CYCLES_FEEDBACK_SETS_H
#define IGRAPH_CYCLES_FEEDBACK_SETS_H

#include "igraph_decls.h"
#include "igraph_datatype.h"
#include "igraph_vector.h"

IGRAPH_BEGIN_C_DECLS

igraph_error_t igraph_i_feedback_arc_set_eades(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights, igraph_vector_int_t *layering
);
igraph_error_t igraph_i_feedback_arc_set_ip_ti(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights);
igraph_error_t igraph_i_feedback_arc_set_ip_cg(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights);
igraph_error_t igraph_i_feedback_arc_set_ip_cb(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights);
igraph_error_t igraph_i_feedback_arc_set_undirected(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights, igraph_vector_int_t *layering
);

igraph_error_t igraph_i_feedback_vertex_set_ip_cg(
        const igraph_t *graph, igraph_vector_int_t *result,
        const igraph_vector_t *weights);

IGRAPH_END_C_DECLS

#endif /* IGRAPH_CYCLES_FEEDBACK_SETS_H */
