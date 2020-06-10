/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constants.h"
#include "igraph_datatype.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_glpk_support.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_structural.h"
#include "igraph_types.h"
#include "igraph_visitor.h"

/**
 * \ingroup structural
 * \function igraph_feedback_vertex_set
 * \brief Calculates a feedback vertex set of the graph.
 *
 * </para><para>
 * A feedback vertex set is a set of vertices whose removal makes the graph acyclic.
 *
 * </para><para>
 * The problem is NP complete except in undirected graphs of maximum degree three.
 *
 * \param graph  The graph object.
 * \param result An initialized vector, the result will be returned here.
 * \param weights Weight vector or NULL if no weights are specified.
 *
 * \return Error code:
 *         \c IGRAPH_EINVAL if an unknown method was specified or the weight vector
 *            is invalid.
 *
 */
int igraph_feedback_vertex_set(const igraph_t *graph,
		igraph_vector_t *result,
                const igraph_vector_t *weights) {

    if (!igraph_is_directed(graph)) {
        return igraph_i_feedback_vertex_set_undirected(graph, result, weights);
    }

    IGRAPH_ERROR("Feedback vertex set for directed graphs is not implemented",
		 IGRAPH_EINVAL);

}


int igraph_i_feedback_vertex_set_undirected(const igraph_t *graph,
	igraph_vector_t *result,
        const igraph_vector_t *weights) {

    long int i;
    long int no_of_nodes = igraph_vcount(graph);

    /* TODO: implement */
    /* For now, return all vertices */
    igraph_vector_clear(vertices);
    IGRAPH_CHECK(igraph_vector_reserve(vertices, no_of_nodes));

    igraph_vs_t vertices_all = igraph_vss_all();
    for (i = 0; i < no_of_nodes; i++) {
	    result[i] = vertices_all[i];
    }

    return IGRAPH_SUCCESS;
}
