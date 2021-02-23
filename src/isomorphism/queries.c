/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020 The igraph development team

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

#include "igraph_interface.h"
#include "igraph_structural.h"

/**
 * \section about_graph_isomorphism
 *
 * <para>igraph provides four set of functions to deal with graph
 * isomorphism problems.</para>
 *
 * <para>The \ref igraph_isomorphic() and \ref igraph_subisomorphic()
 * functions make up the first set (in addition with the \ref
 * igraph_permute_vertices() function). These functions choose the
 * algorithm which is best for the supplied input graph. (The choice is
 * not very sophisticated though, see their documentation for
 * details.)</para>
 *
 * <para>The VF2 graph (and subgraph) isomorphism algorithm is implemented in
 * igraph, these functions are the second set. See \ref
 * igraph_isomorphic_vf2() and \ref igraph_subisomorphic_vf2() for
 * starters.</para>
 *
 * <para>Functions for the Bliss algorithm constitute the third set,
 * see \ref igraph_isomorphic_bliss().</para>
 *
 * <para>Finally, the isomorphism classes of all graphs with three and
 * four vertices are precomputed and stored in igraph, so for these
 * small graphs there is a very simple fast way to decide isomorphism.
 * See \ref igraph_isomorphic_34().
 * </para>
 */

/**
 * \function igraph_isomorphic
 * \brief Decides whether two graphs are isomorphic
 *
 * </para><para>
 * In simple terms, two graphs are isomorphic if they become indistinguishable
 * from each other once their vertex labels are removed (rendering the vertices
 * within each graph indistiguishable). More precisely, two graphs are isomorphic
 * if there is a one-to-one mapping from the vertices of the first one
 * to the vertices of the second such that it transforms the edge set of the
 * first graph into the edge set of the second. This mapping is called
 * an \em isomorphism.
 *
 * </para><para>Currently, this function supports simple graphs and graphs
 * with self-loops, but does not support multigraphs.
 *
 * </para><para>This function decides which graph isomorphism algorithm to be
 * used based on the input graphs. Right now it does the following:
 * \olist
 * \oli If one graph is directed and the other undirected then an
 *    error is triggered.
 * \oli If one of the graphs has multi-edges then an error is triggered.
 * \oli If the two graphs does not have the same number of vertices
 *    and edges it returns with \c FALSE.
 * \oli Otherwise, if the graphs have three or four vertices then an O(1)
 *    algorithm is used with precomputed data.
 * \oli Otherwise Bliss is used, see \ref igraph_isomorphic_bliss().
 * \endolist
 *
 * </para><para>Please call the VF2 and Bliss functions directly if you need
 * something more sophisticated, e.g. you need the isomorphic mapping.
 *
 * \param graph1 The first graph.
 * \param graph2 The second graph.
 * \param iso Pointer to a logical variable, will be set to TRUE (1)
 *        if the two graphs are isomorphic, and FALSE (0) otherwise.
 * \return Error code.
 * \sa \ref igraph_isoclass(), \ref igraph_isoclass_subgraph(),
 * \ref igraph_isoclass_create().
 *
 * Time complexity: exponential.
 */
int igraph_isomorphic(const igraph_t *graph1, const igraph_t *graph2,
                      igraph_bool_t *iso) {

    long int nodes1 = igraph_vcount(graph1), nodes2 = igraph_vcount(graph2);
    long int edges1 = igraph_ecount(graph1), edges2 = igraph_ecount(graph2);
    igraph_bool_t dir1 = igraph_is_directed(graph1), dir2 = igraph_is_directed(graph2);
    igraph_bool_t loop1, loop2, multi1, multi2;

    IGRAPH_CHECK(igraph_has_multiple(graph1, &multi1));
    IGRAPH_CHECK(igraph_has_multiple(graph2, &multi2));

    if (multi1 || multi2) {
        IGRAPH_ERROR("Isomorphism testing is not implemented for multigraphs", IGRAPH_UNIMPLEMENTED);
    }

    if (dir1 != dir2) {
        IGRAPH_ERROR("Cannot compare directed and undirected graphs", IGRAPH_EINVAL);
    } else if (nodes1 != nodes2 || edges1 != edges2) {
        *iso = 0;
    } else if (nodes1 == 3 || nodes1 == 4) {
        IGRAPH_CHECK(igraph_has_loop(graph1, &loop1));
        IGRAPH_CHECK(igraph_has_loop(graph2, &loop2));
        if (!loop1 && !loop2) {
            IGRAPH_CHECK(igraph_isomorphic_34(graph1, graph2, iso));
        } else {
            IGRAPH_CHECK(igraph_isomorphic_bliss(graph1, graph2, NULL, NULL, iso,
                                                 0, 0, /*sh=*/ IGRAPH_BLISS_FL, 0, 0));
        }
    } else {
        IGRAPH_CHECK(igraph_isomorphic_bliss(graph1, graph2, NULL, NULL, iso,
                                             0, 0, /*sh=*/ IGRAPH_BLISS_FL, 0, 0));
    }

    return 0;
}

/**
 * \function igraph_isomorphic_34
 * Graph isomorphism for 3-4 vertices
 *
 * This function uses precomputed indices to decide isomorphism
 * problems for graphs with only 3 or 4 vertices. Multi-edges
 * and self-loops are ignored by this function.
 * \param graph1 The first input graph.
 * \param graph2 The second input graph. Must have the same
 *   directedness as \p graph1.
 * \param iso Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 */
int igraph_isomorphic_34(const igraph_t *graph1, const igraph_t *graph2,
                         igraph_bool_t *iso) {

    igraph_integer_t class1, class2;
    IGRAPH_CHECK(igraph_isoclass(graph1, &class1));
    IGRAPH_CHECK(igraph_isoclass(graph2, &class2));
    *iso = (class1 == class2);
    return 0;
}

/**
 * \function igraph_subisomorphic
 * \brief Decide subgraph isomorphism.
 *
 * Check whether \p graph2 is isomorphic to a subgraph of \p graph1.
 * Currently this function just calls \ref igraph_subisomorphic_vf2()
 * for all graphs.
 *
 * </para><para>
 * Currently this function does not support non-simple graphs.
 *
 * \param graph1 The first input graph, may be directed or
 *   undirected. This is supposed to be the bigger graph.
 * \param graph2 The second input graph, it must have the same
 *   directedness as \p graph2, or an error is triggered. This is
 *   supposed to be the smaller graph.
 * \param iso Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: exponential.
 */
int igraph_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
                         igraph_bool_t *iso) {

    return igraph_subisomorphic_vf2(graph1, graph2, NULL, NULL, NULL, NULL, iso, NULL, NULL, NULL, NULL, NULL);
}
