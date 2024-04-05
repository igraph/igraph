/*
   IGraph library.
   Copyright (C) 2022 The igraph development team

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

#include "igraph_operators.h"

#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_datatype.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_vector.h"

#include "graph/attributes.h"
#include "graph/internal.h"

/**
 * \function igraph_reverse_edges
 * \brief Reverses some edges of a directed graph.
 *
 * This function reverses some edges of a directed graph. The modification is done in place.
 * All attributes, as well as the ordering of edges and vertices are preserved.
 *
 * </para><para>
 * Note that is rarely necessary to reverse \em all edges, as almost all functions that
 * handle directed graphs take a \c mode argument that can be set to \c IGRAPH_IN to
 * effectively treat edges as reversed.
 *
 * \param graph The graph whose edges will be reversed.
 * \param es    The edges to be reversed.
 *              Pass <code>igraph_ess_all(IGRAPH_EDGEORDER_ID)</code> to reverse all edges.
 * \return Error code.
 *
 * Time complexity: O(1) if all edges are reversed, otherwise
 * O(|E|) where |E| is the number of edges in the graph.
 */
igraph_error_t igraph_reverse_edges(igraph_t *graph, const igraph_es_t eids) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t edges;
    igraph_eit_t eit;
    igraph_t new_graph;

    /* Nothing to do on undirected graph. */
    if (! igraph_is_directed(graph)) {
        return IGRAPH_SUCCESS;
    }

    /* Use fast method when all edges are to be reversed. */
    if (igraph_es_is_all(&eids)) {
        return igraph_i_reverse(graph);
    }

    /* Convert graph to edge list. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 2*no_of_edges);
    IGRAPH_CHECK(igraph_get_edgelist(graph, &edges, /* bycol= */ false));

    /* Reverse the edges. */

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        igraph_integer_t eid = IGRAPH_EIT_GET(eit);
        igraph_integer_t tmp = VECTOR(edges)[2*eid];
        VECTOR(edges)[2*eid] = VECTOR(edges)[2*eid + 1];
        VECTOR(edges)[2*eid + 1] = tmp;
    }

    /* Re-create graph from edge list and transfer attributes. */
    IGRAPH_CHECK(igraph_create(&new_graph, &edges, no_of_nodes, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_destroy, &new_graph);

    IGRAPH_I_ATTRIBUTE_DESTROY(&new_graph);
    IGRAPH_I_ATTRIBUTE_COPY(&new_graph, graph, true, true, true); /* does IGRAPH_CHECK */

    igraph_eit_destroy(&eit);
    igraph_vector_int_destroy(&edges);
    igraph_destroy(graph);
    IGRAPH_FINALLY_CLEAN(3);

    *graph = new_graph;

    return IGRAPH_SUCCESS;
}
