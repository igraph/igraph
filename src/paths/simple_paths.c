/*
   IGraph library.
   Copyright (C) 2014-2022  The igraph development team <igraph@igraph.org>

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


#include "igraph_paths.h"

#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_adjlist.h"

#include "core/interruption.h"

/**
 * \function igraph_get_all_simple_paths
 * \brief List all simple paths from one source.
 *
 * A path is simple if its vertices are unique, i.e. no vertex
 * is visited more than once.
 *
 * </para><para>
 * Note that potentially there are exponentially many
 * paths between two vertices of a graph, and you may
 * run out of memory when using this function when the
 * graph has many cycles. Consider using the \p cutoff
 * parameter when you do not need long paths.
 *
 * \param graph The input graph.
 * \param res Initialized integer vector. The paths are
 *        returned here in terms of their vertices, separated
 *        by <code>-1</code> markers. The paths are included in arbitrary
 *        order, as they are found.
 * \param from The start vertex.
 * \param to The target vertices.
 * \param cutoff Maximum length of path that is considered. If
 *        negative, paths of all lengths are considered.
 * \param mode The type of the paths to consider, it is ignored
 *        for undirected graphs.
 * \return Error code.
 *
 * \sa \ref igraph_get_k_shortest_paths()
 *
 * Time complexity: O(n!) in the worst case, n is the number of
 * vertices.
 */

igraph_error_t igraph_get_all_simple_paths(const igraph_t *graph,
                                igraph_vector_int_t *res,
                                igraph_integer_t from,
                                const igraph_vs_t to,
                                igraph_integer_t cutoff,
                                igraph_neimode_t mode) {

    igraph_integer_t no_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    igraph_bool_t toall = igraph_vs_is_all(&to);
    igraph_lazy_adjlist_t adjlist;
    igraph_vector_int_t stack, dist; /* used as a stack, but represented as a vector,
                                        in order to be appendable to other vectors */
    igraph_vector_bool_t markto, added;
    igraph_vector_int_t nptr;
    int iter = 0;

    if (from < 0 || from >= no_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }

    if (!toall) {
        IGRAPH_VECTOR_BOOL_INIT_FINALLY(&markto, no_nodes);
        IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
            VECTOR(markto)[ IGRAPH_VIT_GET(vit) ] = true;
        }
        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&added, no_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&stack, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&dist, 100);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(
        graph, &adjlist, mode, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE
    ));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nptr, no_nodes);

    igraph_vector_int_clear(res);

    igraph_vector_int_clear(&stack);
    igraph_vector_int_clear(&dist);
    igraph_vector_int_push_back(&stack, from);
    igraph_vector_int_push_back(&dist, 0);
    VECTOR(added)[from] = true;
    while (!igraph_vector_int_empty(&stack)) {
        igraph_integer_t act = igraph_vector_int_tail(&stack);
        igraph_integer_t curdist = igraph_vector_int_tail(&dist);
        igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, act);
        igraph_integer_t n;
        igraph_integer_t *ptr = igraph_vector_int_get_ptr(&nptr, act);
        igraph_bool_t any;
        igraph_bool_t within_dist;
        igraph_integer_t nei;

        IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");

        n = igraph_vector_int_size(neis);

        within_dist = (curdist < cutoff || cutoff < 0);
        if (within_dist) {
            /* Search for a neighbor that was not yet visited */
            any = false;
            while (!any && (*ptr) < n) {
                nei = VECTOR(*neis)[(*ptr)];
                any = !VECTOR(added)[nei];
                (*ptr) ++;
            }
        }
        if (within_dist && any) {
            /* There is such a neighbor, add it */
            IGRAPH_CHECK(igraph_vector_int_push_back(&stack, nei));
            IGRAPH_CHECK(igraph_vector_int_push_back(&dist, curdist + 1));
            VECTOR(added)[nei] = true;
            /* Add to results */
            if (toall || VECTOR(markto)[nei]) {
                IGRAPH_CHECK(igraph_vector_int_append(res, &stack));
                IGRAPH_CHECK(igraph_vector_int_push_back(res, -1));
            }
        } else {
            /* There is no such neighbor, finished with the subtree */
            igraph_integer_t up = igraph_vector_int_pop_back(&stack);
            igraph_vector_int_pop_back(&dist);
            VECTOR(added)[up] = false;
            VECTOR(nptr)[up] = 0;
        }

        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 13);
    }

    igraph_vector_int_destroy(&nptr);
    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&dist);
    igraph_vector_int_destroy(&stack);
    igraph_vector_bool_destroy(&added);
    IGRAPH_FINALLY_CLEAN(5);

    if (!toall) {
        igraph_vector_bool_destroy(&markto);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
