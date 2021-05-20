/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_topology.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"

/**
 * \function igraph_simplify_and_colorize
 * \brief Simplify the graph and compute self-loop and edge multiplicities.
 *
 * </para><para>
 * This function creates a vertex and edge colored simple graph from the input
 * graph. The vertex colors are computed as the number of incident self-loops
 * to each vertex in the input graph. The edge colors are computed as the number of
 * parallel edges in the input graph that were merged to create each edge
 * in the simple graph.
 *
 * </para><para>
 * The resulting colored simple graph is suitable for use by isomorphism checking
 * algorithms such as VF2, which only support simple graphs, but can consider
 * vertex and edge colors.
 *
 * \param graph The graph object, typically having self-loops or multi-edges.
 * \param res An uninitialized graph object. The result will be stored here
 * \param vertex_color Computed vertex colors corresponding to self-loop multiplicities.
 * \param edge_color Computed edge colors corresponding to edge multiplicities
 * \return Error code.
 *
 * \sa \ref igraph_simplify(), \ref igraph_isomorphic_vf2(), \ref igraph_subisomorphic_vf2()
 *
 */
int igraph_simplify_and_colorize(
    const igraph_t *graph, igraph_t *res,
    igraph_vector_int_t *vertex_color, igraph_vector_int_t *edge_color) {
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_vector_t edges;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int pto = -1, pfrom = -1;
    long int i;

    IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_FROM));
    IGRAPH_FINALLY(igraph_es_destroy, &es);
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    IGRAPH_CHECK(igraph_vector_int_resize(vertex_color, no_of_nodes));
    igraph_vector_int_null(vertex_color);

    IGRAPH_CHECK(igraph_vector_int_resize(edge_color, no_of_edges));
    igraph_vector_int_null(edge_color);

    i = -1;
    for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        long int edge = IGRAPH_EIT_GET(eit);
        long int from = IGRAPH_FROM(graph, edge);
        long int to   = IGRAPH_TO(graph, edge);

        if (to == from) {
            VECTOR(*vertex_color)[to]++;
            continue;
        }

        if (to == pto && from == pfrom) {
            VECTOR(*edge_color)[i]++;
        } else {
            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);
            i++;
            VECTOR(*edge_color)[i] = 1;
        }

        pfrom = from; pto = to;
    }

    igraph_vector_int_resize(edge_color, i + 1);

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
