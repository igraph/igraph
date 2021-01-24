/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_paths.h"
#include "igraph_interface.h"
#include "core/interruption.h"
#include "igraph_vector_ptr.h"
#include "igraph_iterators.h"
#include "igraph_adjlist.h"
#include "igraph_stack.h"

/**
 * \function igraph_get_all_simple_paths
 * List all simple paths from one source
 *
 * A path is simple, if its vertices are unique, no vertex
 * is visited more than once.
 *
 * </para><para>
 * Note that potentially there are exponentially many
 * paths between two vertices of a graph, and you may
 * run out of memory when using this function, if your
 * graph is lattice-like.
 *
 * </para><para>
 * This function currently ignored multiple and loop edges.
 * \param graph The input graph.
 * \param res Initialized integer vector, all paths are
 *        returned here, separated by -1 markers. The paths
 *        are included in arbitrary order, as they are found.
 * \param from The start vertex.
 * \param to The target vertices.
 * \param cutoff Maximum length of path that is considered. If
 *        negative, paths of all lengths are considered.
 * \param mode The type of the paths to consider, it is ignored
 *        for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(n!) in the worst case, n is the number of
 * vertices.
 */
igraph_error_t igraph_get_all_simple_paths(const igraph_t *graph,
                                igraph_vector_long_t *res,
                                igraph_long_t from,
                                const igraph_vs_t to,
                                igraph_long_t cutoff,
                                igraph_neimode_t mode) {

    igraph_long_t no_nodes = igraph_vcount(graph);
    igraph_vit_t vit;
    igraph_bool_t toall = igraph_vs_is_all(&to);
    igraph_vector_char_t markto;
    igraph_lazy_adjlist_t adjlist;
    igraph_vector_long_t stack, dist;
    igraph_vector_char_t added;
    igraph_vector_long_t nptr;
    igraph_long_t iteration = 0;

    if (from < 0 || from >= no_nodes) {
        IGRAPH_ERROR("Invalid starting vertex", IGRAPH_EINVAL);
    }

    if (!toall) {
        igraph_vector_char_init(&markto, no_nodes);
        IGRAPH_FINALLY(igraph_vector_char_destroy, &markto);
        IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
            VECTOR(markto)[ IGRAPH_VIT_GET(vit) ] = 1;
        }
        igraph_vit_destroy(&vit);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_vector_char_init(&added, no_nodes));
    IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
    IGRAPH_CHECK(igraph_vector_long_init(&stack, 100));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &stack);
    IGRAPH_CHECK(igraph_vector_long_init(&dist, 100));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &dist);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode,
                                          /*simplify=*/ 1));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_CHECK(igraph_vector_long_init(&nptr, no_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &nptr);

    igraph_vector_long_clear(res);

    igraph_vector_long_clear(&stack);
    igraph_vector_long_clear(&dist);
    igraph_vector_long_push_back(&stack, from);
    igraph_vector_long_push_back(&dist, 0);
    VECTOR(added)[from] = 1;
    while (!igraph_vector_long_empty(&stack)) {
        igraph_long_t act = igraph_vector_long_tail(&stack);
        igraph_long_t curdist = igraph_vector_long_tail(&dist);
        igraph_vector_t *neis = igraph_lazy_adjlist_get(&adjlist, act);
        igraph_long_t n = igraph_vector_size(neis);
        igraph_long_t *ptr = igraph_vector_long_e_ptr(&nptr, act);
        igraph_bool_t any;
        igraph_bool_t within_dist;
        igraph_long_t nei;

        if (iteration == 0) {
            IGRAPH_ALLOW_INTERRUPTION();
        }

        within_dist = (curdist < cutoff || cutoff < 0);
        if (within_dist) {
            /* Search for a neighbor that was not yet visited */
            any = 0;
            while (!any && (*ptr) < n) {
                nei = (igraph_long_t) VECTOR(*neis)[(*ptr)];
                any = !VECTOR(added)[nei];
                (*ptr) ++;
            }
        }
        if (within_dist && any) {
            /* There is such a neighbor, add it */
            IGRAPH_CHECK(igraph_vector_long_push_back(&stack, nei));
            IGRAPH_CHECK(igraph_vector_long_push_back(&dist, curdist + 1));
            VECTOR(added)[nei] = 1;
            /* Add to results */
            if (toall || VECTOR(markto)[nei]) {
                IGRAPH_CHECK(igraph_vector_long_append(res, &stack));
                IGRAPH_CHECK(igraph_vector_long_push_back(res, -1));
            }
        } else {
            /* There is no such neighbor, finished with the subtree */
            igraph_long_t up = igraph_vector_long_pop_back(&stack);
            igraph_vector_long_pop_back(&dist);
            VECTOR(added)[up] = 0;
            VECTOR(nptr)[up] = 0;
        }

        iteration++;
        if (iteration >= 10000) {
            iteration = 0;
        }
    }

    igraph_vector_long_destroy(&nptr);
    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_vector_long_destroy(&dist);
    igraph_vector_long_destroy(&stack);
    igraph_vector_char_destroy(&added);
    IGRAPH_FINALLY_CLEAN(5);

    if (!toall) {
        igraph_vector_char_destroy(&markto);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}
