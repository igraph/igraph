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

#include "igraph_structural.h"

#include "igraph_bipartite.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_operators.h"
#include "igraph_topology.h"

#include "core/interruption.h"

/**
 * \function igraph_is_perfect
 * \brief Checks if the graph is perfect.
 *
 * A perfect graph is an undirected graph in which the chromatic number of every induced
 * subgraph equals the order of the largest clique of that subgraph.
 * The chromatic number of a graph G is the smallest number of colors needed to
 * color the vertices of G so that no two adjacent vertices share the same color.
 *
 * </para><para>
 * Warning: This function may create the complement of the graph internally,
 * which consumes a lot of memory. For moderately sized graphs, consider
 * decomposing them into biconnected components and running the check separately
 * on each component.
 *
 * </para><para>
 * This implementation is based on the strong perfect graph theorem which was
 * conjectured by Claude Berge and proved by Maria Chudnovsky, Neil Robertson,
 * Paul Seymour, and Robin Thomas.
 *
 * \param graph The input graph. It is expected to be undirected and simple.
 * \param perfect Pointer to an integer, the result will be stored here.
 * \return Error code.
 *
 * Time complexity: worst case exponenital, often faster in practice.
 */
igraph_error_t igraph_is_perfect(const igraph_t *graph, igraph_bool_t *perfect) {

    igraph_bool_t is_bipartite, is_chordal, iso, is_simple;
    igraph_real_t girth, comp_girth;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t start;
    igraph_integer_t cycle_len;
    igraph_t comp_graph, cycle;

    // If the graph is directed return error.
    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("The concept of perfect graphs is only defined for undirected graphs.", IGRAPH_EINVAL);
    }

    // If the graph isn't simple then return an error.
    IGRAPH_CHECK(igraph_is_simple(graph, &is_simple));
    if (!is_simple) {
        IGRAPH_ERROR("Perfect graph testing is implemented for simple graphs only. Simplify the graph.", IGRAPH_EINVAL);
    }

    // All graphs with less than 5 vertices are perfect.
    if (no_of_nodes < 5) {
        *perfect = true;
        return IGRAPH_SUCCESS;
    }

    // Graphs with less than 5 edges or a complement with less than 5 edges
    // are also perfect. The following check handles most 5-vertex graphs,
    // but its usefulness quickly diminishes, with only 0.3% of unlabelled
    // 8-vertex graphs handled.
    // In order to avoid bad results due to integer overflow with large graphs,
    // we limit this check for small graphs only.
    if ( no_of_nodes < 10000 &&
         (no_of_edges < 5 || no_of_edges > (no_of_nodes - 1) * no_of_nodes / 2 - 5)) {
        *perfect = true;
        return IGRAPH_SUCCESS;
    }

    // Chordal and bipartite graph types are perfect.
    // Possibly more optimizations found here: http://www.or.uni-bonn.de/~hougardy/paper/ClassesOfPerfectGraphs.pdf
    IGRAPH_CHECK(igraph_is_bipartite(graph, &is_bipartite, NULL));
    if (is_bipartite) {
        *perfect = true;
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(igraph_is_chordal(graph, NULL, NULL, &is_chordal, NULL, NULL));
    if (is_chordal) {
        *perfect = true;
        return IGRAPH_SUCCESS;
    }

    // The weak perfect graph theorem:
    // A graph is perfect iff its complement is perfect.
    IGRAPH_CHECK(igraph_complementer(&comp_graph, graph, 0));
    IGRAPH_FINALLY(igraph_destroy, &comp_graph);

    IGRAPH_CHECK(igraph_is_bipartite(&comp_graph, &is_bipartite, NULL));
    if (is_bipartite) {
        *perfect = true;
        goto clean1;
    }

    IGRAPH_CHECK(igraph_is_chordal(&comp_graph, NULL, NULL, &is_chordal, NULL, NULL));
    if (is_chordal) {
        *perfect = true;
        goto clean1;
    }

    // Since igraph_is_bipartite also catches trees, at this point the girth
    // of the graph and its complementer (to be stored in girth and comp_girth)
    // are both guaranteed to be finite.

    // If the girth (or the smallest circle in the graph) is bigger than 3 and have odd number of vertices then
    // the graph isn't perfect.
    IGRAPH_CHECK(igraph_girth(graph, &girth, NULL));
    if ((girth > 3) && (((igraph_integer_t)girth) % 2 == 1)) {
        *perfect = false;
        goto clean1;
    }

    IGRAPH_CHECK(igraph_girth(&comp_graph, &comp_girth, NULL));
    if ((comp_girth > 3) && (((igraph_integer_t)comp_girth) % 2 == 1)) {
        *perfect = false;
        goto clean1;
    }

    // At this point girth and comp_girth are both at least 3.

    // Strong perfect graph theorem:
    // A graph is perfect iff neither it or its complement contains an induced odd cycle of length >= 5
    // (i.e. an odd hole). TODO: Find a more efficient way to check for odd holes.
    start = (igraph_integer_t) (girth < comp_girth ? girth : comp_girth);
    start = start % 2 == 0 ? start + 1 : start + 2;
    for (cycle_len = start; cycle_len <= no_of_nodes ; cycle_len += 2) {

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_ring(&cycle, cycle_len, IGRAPH_UNDIRECTED, /* mutual */ 0, /* circular */ 1));
        IGRAPH_FINALLY(igraph_destroy, &cycle);

        if (cycle_len > girth) {
            IGRAPH_CHECK(igraph_subisomorphic_lad(&cycle, graph, NULL, &iso, NULL, NULL, /* induced */ 1, 0));
            if (iso) {
                *perfect = false;
                goto clean2;
            }
        }

        if (cycle_len > comp_girth) {
            IGRAPH_CHECK(igraph_subisomorphic_lad(&cycle, &comp_graph, NULL, &iso, NULL, NULL, /* induced */ 1, 0));
            if (iso) {
                *perfect = false;
                goto clean2;
            }
        }

        igraph_destroy(&cycle);
        IGRAPH_FINALLY_CLEAN(1);
    }

    *perfect = true;

clean1:
    /* normal exit route */
    igraph_destroy(&comp_graph);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;

clean2:
    /* exit route if we also have a cycle to destroy */
    igraph_destroy(&cycle);
    IGRAPH_FINALLY_CLEAN(1);
    goto clean1;
}
