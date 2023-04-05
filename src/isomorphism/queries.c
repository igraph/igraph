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
 * <para>Finally, the isomorphism classes of all directed graphs with three and
 * four vertices and all undirected graphs with 3-6 vertices are precomputed
 * and stored in igraph, so for these small graphs there is a separate fast
 * path in the code that does not use more complex, generic isomorphism
 * algorithms.</para>
 */

static igraph_error_t igraph_i_isomorphic_small(
    const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *iso
);

/**
 * \function igraph_isomorphic
 * \brief Are two graphs isomorphic?
 *
 * In simple terms, two graphs are isomorphic if they become indistinguishable
 * from each other once their vertex labels are removed (rendering the vertices
 * within each graph indistiguishable). More precisely, two graphs are isomorphic
 * if there is a one-to-one mapping from the vertices of the first one
 * to the vertices of the second such that it transforms the edge set of the
 * first graph into the edge set of the second. This mapping is called
 * an \em isomorphism.
 *
 * </para><para>This function decides which graph isomorphism algorithm to be
 * used based on the input graphs. Right now it does the following:
 * \olist
 * \oli If one graph is directed and the other undirected then an
 *    error is triggered.
 * \oli If one of the graphs has multi-edges then both graphs are
 *    simplified and colorized using \ref igraph_simplify_and_colorize() and sent to VF2.
 * \oli If the two graphs does not have the same number of vertices
 *    and edges it returns with \c false.
 * \oli Otherwise, if the \ref igraph_isoclass() function supports both
 *    graphs (which is true for directed graphs with 3 and 4 vertices, and
 *    undirected graphs with 3-6 vertices), an O(1) algorithm is used with
 *    precomputed data.
 * \oli Otherwise Bliss is used, see \ref igraph_isomorphic_bliss().
 * \endolist
 *
 * </para><para>Please call the VF2 and Bliss functions directly if you need
 * something more sophisticated, e.g. you need the isomorphic mapping.
 *
 * \param graph1 The first graph.
 * \param graph2 The second graph.
 * \param iso Pointer to a logical variable, will be set to \c true
 *        if the two graphs are isomorphic, and \c false otherwise.
 * \return Error code.
 * \sa \ref igraph_isoclass(), \ref igraph_isoclass_subgraph(),
 * \ref igraph_isoclass_create().
 *
 * Time complexity: exponential.
 */
igraph_error_t igraph_isomorphic(const igraph_t *graph1, const igraph_t *graph2,
                      igraph_bool_t *iso) {

    igraph_integer_t nodes1 = igraph_vcount(graph1), nodes2 = igraph_vcount(graph2);
    igraph_integer_t edges1 = igraph_ecount(graph1), edges2 = igraph_ecount(graph2);
    igraph_bool_t dir1 = igraph_is_directed(graph1), dir2 = igraph_is_directed(graph2);
    igraph_bool_t loop1, loop2, multi1, multi2;

    if (dir1 != dir2) {
        IGRAPH_ERROR("Cannot compare directed and undirected graphs for isomorphism.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_has_multiple(graph1, &multi1));
    IGRAPH_CHECK(igraph_has_multiple(graph2, &multi2));

    if (multi1 || multi2) {
        igraph_t r1;
        igraph_t r2;
        igraph_vector_int_t vc1;
        igraph_vector_int_t vc2;
        igraph_vector_int_t ec1;
        igraph_vector_int_t ec2;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&vc1, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vc2, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&ec1, 0);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&ec2, 0);
        IGRAPH_CHECK(igraph_simplify_and_colorize(graph1, &r1, &vc1, &ec1));
        IGRAPH_FINALLY(igraph_destroy, &r1);
        IGRAPH_CHECK(igraph_simplify_and_colorize(graph2, &r2, &vc2, &ec2));
        IGRAPH_FINALLY(igraph_destroy, &r2);
        IGRAPH_CHECK(igraph_isomorphic_vf2(&r1, &r2, &vc1, &vc2, &ec1, &ec2, iso,
                NULL, NULL, NULL, NULL, NULL));
        igraph_destroy(&r2);
        igraph_destroy(&r1);
        igraph_vector_int_destroy(&ec2);
        igraph_vector_int_destroy(&ec1);
        igraph_vector_int_destroy(&vc2);
        igraph_vector_int_destroy(&vc1);
        IGRAPH_FINALLY_CLEAN(6);

        return IGRAPH_SUCCESS;
    }

    if (nodes1 != nodes2 || edges1 != edges2) {
        *iso = false;
    } else if (nodes1 >= 3 && nodes1 <= (dir1 ? 4 : 6)) {
        IGRAPH_CHECK(igraph_has_loop(graph1, &loop1));
        IGRAPH_CHECK(igraph_has_loop(graph2, &loop2));
        if (!loop1 && !loop2) {
            IGRAPH_CHECK(igraph_i_isomorphic_small(graph1, graph2, iso));
        } else {
            IGRAPH_CHECK(igraph_isomorphic_bliss(graph1, graph2, NULL, NULL, iso,
                                                 NULL, NULL, /*sh=*/ IGRAPH_BLISS_FL, NULL, NULL));
        }
    } else {
        IGRAPH_CHECK(igraph_isomorphic_bliss(graph1, graph2, NULL, NULL, iso,
                                             NULL, NULL, /*sh=*/ IGRAPH_BLISS_FL, NULL, NULL));
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_isomorphic_34
 * \brief Graph isomorphism for 3-4 vertices (deprecated).
 *
 * \deprecated-by igraph_isomorphic 0.10.0
 *
 * If you really care about performance and you \em know for sure that your
 * input graphs are simple and have either 3 or 4 vertices for directed graphs,
 * or 3-6 vertices for undirected graphs, you can compare their isomorphism
 * classes obtained from \ref igraph_isoclass() directly instead of calling
 * \ref igraph_isomorphic(); this saves the cost of checking whether the graphs
 * do not contain multiple edges or self-loops.
 *
 * \param graph1 The first input graph.
 * \param graph2 The second input graph. Must have the same
 *   directedness as \p graph1.
 * \param iso Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_isomorphic_34(
    const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *iso
) {
    return igraph_i_isomorphic_small(graph1, graph2, iso);
}

/**
 * \function igraph_i_isomorphic_small
 * \brief Graph isomorphism for small graphs.
 *
 * This function uses precomputed indices to decide isomorphism
 * problems for directed graphs with only 3 or 4 vertices, or for undirected
 * graphs with 3, 4, 5 or 6 vertices. Multi-edges and self-loops are ignored by
 * this function.
 *
 * \param graph1 The first input graph.
 * \param graph2 The second input graph. Must have the same
 *   directedness as \p graph1.
 * \param iso Pointer to a boolean, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(1).
 */
igraph_error_t igraph_i_isomorphic_small(
    const igraph_t *graph1, const igraph_t *graph2, igraph_bool_t *iso
) {
    igraph_integer_t class1, class2;
    IGRAPH_CHECK(igraph_isoclass(graph1, &class1));
    IGRAPH_CHECK(igraph_isoclass(graph2, &class2));
    *iso = (class1 == class2);
    return IGRAPH_SUCCESS;
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
igraph_error_t igraph_subisomorphic(const igraph_t *graph1, const igraph_t *graph2,
                         igraph_bool_t *iso) {

    return igraph_subisomorphic_vf2(graph1, graph2, NULL, NULL, NULL, NULL, iso, NULL, NULL, NULL, NULL, NULL);
}
