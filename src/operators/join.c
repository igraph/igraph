/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "igraph_operators.h"
#include "igraph_interface.h"

#include "math/safe_intop.h"

/**
 * \function igraph_join
 * \brief Creates the join of two disjoint graphs.
 *
 * \experimental
 *
 * First the vertices of the second graph will be relabeled with new
 * vertex IDs to have two disjoint sets of vertex IDs, then the union
 * of the two graphs will be formed. Finally, the vertces from the
 * first graph will have edges added to each vertex from the second.
 * If the two graphs have |V1| and |V2| vertices and |E1| and |E2|
 * edges respectively then the new graph will have |V1|+|V2| vertices
 * and |E1|+|E2|+|V1|*|V2| edges.
 *
 * </para><para>
 * The vertex ordering of the graphs will be preserved.
 * In other words, the vertex IDs of the first graph map to
 * identical values in the new graph, while the vertex  IDs
 * of the second graph map to IDs incremented by the vertex
 * count of the first graph. The new edges will be grouped with the
 * other edges that share a from vertex.
 *
 * </para><para>
 * Both graphs need to have the same directedness, i.e. either both
 * directed or both undirected. If both graphs are directed, then for each
 * vertex v, u in graphs G1, G2 we add edges (v, u), (u, v) to maintain
 * completeness.
 *
 * </para><para>
 * The current version of this function cannot handle graph, vertex
 * and edge attributes, they will be lost.
 *
 * \param res  Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param left The first graph.
 * \param right The second graph.
 * \return Error code.
 *
 * Time complexity: O(|V1|*|V2|+|E1|+|E2|).
 *
 */
igraph_error_t igraph_join(igraph_t *res,
                           const igraph_t *left,
                           const igraph_t *right) {

    igraph_integer_t no_of_nodes_left = igraph_vcount(left);
    igraph_integer_t no_of_nodes_right = igraph_vcount(right);
    igraph_integer_t no_of_new_edges;
    igraph_vector_int_t new_edges;
    igraph_bool_t directed_left = igraph_is_directed(left);
    igraph_integer_t i;
    igraph_integer_t j;

    if (directed_left != igraph_is_directed(right)) {
        IGRAPH_ERROR("Cannot create join of directed and undirected graphs.",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_disjoint_union(res,left,right));
    IGRAPH_SAFE_MULT(no_of_nodes_left, no_of_nodes_right ,&no_of_new_edges);
    IGRAPH_SAFE_MULT(no_of_new_edges, 2 ,&no_of_new_edges);
    if (directed_left) {
        IGRAPH_SAFE_MULT(no_of_new_edges, 2 ,&no_of_new_edges);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&new_edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&new_edges, no_of_new_edges));

    for(i = 0; i < no_of_nodes_left; i++) {
        for(j = 0; j < no_of_nodes_right; j++) {
            igraph_vector_int_push_back(&new_edges, i);  /* reserved */
            igraph_vector_int_push_back(&new_edges, j + no_of_nodes_left);  /* reserved */
            if (directed_left) {
                igraph_vector_int_push_back(&new_edges, j + no_of_nodes_left);  /* reserved */
                igraph_vector_int_push_back(&new_edges, i);  /* reserved */
            }
        }
    }

    IGRAPH_CHECK(igraph_add_edges(res, &new_edges, NULL));

    igraph_vector_int_destroy(&new_edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
