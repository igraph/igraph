/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "core/interruption.h"

#include "igraph_bipartite.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_operators.h"
#include "igraph_structural.h"
#include "igraph_topology.h"

/**
 * \function igraph_is_perfect
 * \brief Check if the graph is perfect.
 *
 * </para><para>
 *
 * </para><para>
 * A perfect graph is a graph in which the chromatic number of every induced
 * subgraph equals the order of the largest clique of that subgraph.
 * The chromatic number of a graph G is the smallest number of colors needed to
 * color the vertices of G so that no two adjacent vertices share the same color
 *
 * </para><para>
 * This implementation is based on the strong perfect graph theorem which was
 * conjectured by Claude Berge and proved by Maria Chudnovsky, Neil Robertson,
 * Paul Seymour, and Robin Thomas.
 *
 * \param graph The input graph. The current implementation doesn't work with directed graphs
 *     or graphs with multi - edges.
 * \param perfect Pointer to an integer, if not \c NULL then the result
 *     will be stored here.
 * \return Error code.
 *
 */
int igraph_is_perfect(const igraph_t *graph, igraph_bool_t *perfect) {

    igraph_bool_t is_bipartite, is_chordal, iso, is_simple;
    igraph_integer_t girth, comp_girth, num_of_vertices = igraph_vcount(graph);
    igraph_integer_t start;
    long int i;
    igraph_t comp_graph, cycle;

    if (!perfect) {
        return IGRAPH_SUCCESS;
    }

    // If the graph is directed return error
    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("perfect graph function doesn't support directed graphs", IGRAPH_EINVAL);
    }

    //If the graph isn't simple then return an error
    IGRAPH_CHECK(igraph_is_simple(graph, &is_simple));
    if (!is_simple) {
        IGRAPH_ERROR("perfect graph function doesn't support graphs with multi edges", IGRAPH_EINVAL);
    }


    IGRAPH_CHECK(igraph_is_bipartite(graph, &is_bipartite, NULL));
    IGRAPH_CHECK(igraph_is_chordal(graph, NULL, NULL, &is_chordal, NULL, NULL));
    // chordal and bipartite graph types are perfect.
    // possibly more optimizations found here: http://www.or.uni-bonn.de/~hougardy/paper/ClassesOfPerfectGraphs.pdf

    if (is_chordal || is_bipartite) {
        *perfect = 1;
        return IGRAPH_SUCCESS;
    }

    // The weak perfect graph theorem - a graph is perfect iff its complement is perfect
    IGRAPH_CHECK(igraph_complementer(&comp_graph, graph, 0));
    IGRAPH_FINALLY(igraph_destroy, &comp_graph);
    IGRAPH_CHECK(igraph_is_bipartite(&comp_graph, &is_bipartite, NULL));
    IGRAPH_CHECK(igraph_is_chordal(&comp_graph, NULL, NULL, &is_chordal, NULL, NULL));
    if (is_chordal || is_bipartite) {
        *perfect = 1;
        igraph_destroy(&comp_graph);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    // If the girth (or the smallest circle in the graph) is bigger than 3 and have odd number of vertices then
    // the graph isn't perfect
    IGRAPH_CHECK(igraph_girth(graph, &girth, NULL));
    IGRAPH_CHECK(igraph_girth(&comp_graph, &comp_girth, NULL));
    if ((girth > 3) && (girth % 2 == 1)) {
        *perfect = 0;
        igraph_destroy(&comp_graph);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }
    if ((comp_girth > 3) && (comp_girth % 2 == 1)) {
        *perfect = 0;
        igraph_destroy(&comp_graph);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    }

    // since igraph_is_bipartite also catches trees, at this point girth and comp_girth are both at least 3.
    // for trees, their value would have been 0

    // strong perfect graph theorem
    // a graph is perfect iff neither it or its complement contains an induced odd cycle of length >= 5
    start = girth > comp_girth ? girth : comp_girth;
    start = start % 2 == 0 ? start + 1 : start + 2;
    for (i = start; i <= num_of_vertices; i += 2) {

        IGRAPH_CHECK(igraph_ring(&cycle, i, 0, 0, 1));
        IGRAPH_FINALLY(igraph_destroy, &cycle);
        IGRAPH_CHECK(igraph_subisomorphic_lad(&cycle, graph, NULL, &iso, NULL, NULL, /* induced */ 1, 0));
        if ((i > girth) && (iso)) {
            *perfect = 0;
            igraph_destroy(&cycle);
            igraph_destroy(&comp_graph);
            IGRAPH_FINALLY_CLEAN(2);
            return IGRAPH_SUCCESS;
        }

        IGRAPH_CHECK(igraph_subisomorphic_lad(&cycle, &comp_graph, NULL, &iso, NULL, NULL, /* induced */ 1, 0));
        if ((i > girth) && iso) {
            *perfect = 0;
            igraph_destroy(&cycle);
            igraph_destroy(&comp_graph);
            IGRAPH_FINALLY_CLEAN(2);
            return IGRAPH_SUCCESS;
        }
        igraph_destroy(&cycle);
        IGRAPH_FINALLY_CLEAN(1);
    }
    IGRAPH_ALLOW_INTERRUPTION();
    *perfect = 1;
    igraph_destroy(&comp_graph);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}
