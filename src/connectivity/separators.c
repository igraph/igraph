/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_separators.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_flow.h"
#include "igraph_interface.h"
#include "igraph_operators.h"
#include "igraph_structural.h"
#include "igraph_vector.h"

#include "core/interruption.h"

static igraph_error_t igraph_i_is_separator(const igraph_t *graph,
                                 igraph_vit_t *vit,
                                 igraph_integer_t except,
                                 igraph_bool_t *res,
                                 igraph_vector_bool_t *removed,
                                 igraph_dqueue_int_t *Q,
                                 igraph_vector_int_t *neis,
                                 igraph_integer_t no_of_nodes) {

    igraph_integer_t start = 0;

    if (IGRAPH_VIT_SIZE(*vit) >= no_of_nodes - 1) {
        /* Just need to check that we really have at least n-1 vertices in it */
        igraph_vector_bool_t hit;
        igraph_integer_t nohit = 0;
        IGRAPH_CHECK(igraph_vector_bool_init(&hit, no_of_nodes));
        IGRAPH_FINALLY(igraph_vector_bool_destroy, &hit);
        for (IGRAPH_VIT_RESET(*vit);
             !IGRAPH_VIT_END(*vit);
             IGRAPH_VIT_NEXT(*vit)) {
            igraph_integer_t v = IGRAPH_VIT_GET(*vit);
            if (!VECTOR(hit)[v]) {
                nohit++;
                VECTOR(hit)[v] = 1;
            }
        }
        igraph_vector_bool_destroy(&hit);
        IGRAPH_FINALLY_CLEAN(1);
        if (nohit >= no_of_nodes - 1) {
            *res = 0;
            return IGRAPH_SUCCESS;
        }
    }

    /* Remove the given vertices from the graph, do a breadth-first
       search and check the number of components */

    if (except < 0) {
        for (IGRAPH_VIT_RESET(*vit);
             !IGRAPH_VIT_END(*vit);
             IGRAPH_VIT_NEXT(*vit)) {
            VECTOR(*removed)[ IGRAPH_VIT_GET(*vit) ] = 1;
        }
    } else {
        /* There is an exception */
        igraph_integer_t i;
        for (i = 0, IGRAPH_VIT_RESET(*vit);
             i < except;
             i++, IGRAPH_VIT_NEXT(*vit)) {
            VECTOR(*removed)[ IGRAPH_VIT_GET(*vit) ] = 1;
        }
        for (IGRAPH_VIT_NEXT(*vit);
             !IGRAPH_VIT_END(*vit);
             IGRAPH_VIT_NEXT(*vit)) {
            VECTOR(*removed)[ IGRAPH_VIT_GET(*vit) ] = 1;
        }
    }

    /* Look for the first node that is not removed */
    while (start < no_of_nodes && VECTOR(*removed)[start]) {
        start++;
    }

    if (start == no_of_nodes) {
        IGRAPH_ERROR("All vertices are included in the separator",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_dqueue_int_push(Q, start));
    VECTOR(*removed)[start] = 1;
    while (!igraph_dqueue_int_empty(Q)) {
        igraph_integer_t node = igraph_dqueue_int_pop(Q);
        igraph_integer_t j, n;
        IGRAPH_CHECK(igraph_neighbors(graph, neis, node, IGRAPH_ALL));
        n = igraph_vector_int_size(neis);
        for (j = 0; j < n; j++) {
            igraph_integer_t nei = VECTOR(*neis)[j];
            if (!VECTOR(*removed)[nei]) {
                IGRAPH_CHECK(igraph_dqueue_int_push(Q, nei));
                VECTOR(*removed)[nei] = 1;
            }
        }
    }

    /* Look for the next node that was neighter removed, not visited */
    while (start < no_of_nodes && VECTOR(*removed)[start]) {
        start++;
    }

    /* If there is another component, then we have a separator */
    *res = (start < no_of_nodes);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_separator
 * \brief Would removing this set of vertices disconnect the graph?
 *
 * \param graph The input graph. It may be directed, but edge
 *        directions are ignored.
 * \param candidate The candidate separator. It must not contain all
 *        vertices.
 * \param res Pointer to a Boolean variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number vertices and edges.
 *
 * \example examples/simple/igraph_is_separator.c
 */

igraph_error_t igraph_is_separator(const igraph_t *graph,
                        const igraph_vs_t candidate,
                        igraph_bool_t *res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_bool_t removed;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t neis;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, candidate, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vector_bool_init(&removed, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    IGRAPH_CHECK(igraph_i_is_separator(graph, &vit, -1, res, &removed,
                                       &Q, &neis, no_of_nodes));

    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&Q);
    igraph_vector_bool_destroy(&removed);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_minimal_separator
 * \brief Decides whether a set of vertices is a minimal separator.
 *
 * A set of vertices is a minimal separator, if the removal of the
 * vertices disconnects the graph, and this is not true for any subset
 * of the set.
 *
 * </para><para>This implementation first checks that the given
 * candidate is a separator, by calling \ref
 * igraph_is_separator(). If it is a separator, then it checks that
 * each subset of size n-1, where n is the size of the candidate, is
 * not a separator.
 *
 * \param graph The input graph. It may be directed, but edge
 *        directions are ignored.
 * \param candidate The candidate minimal separators.
 * \param res Pointer to a boolean variable, the result is stored
 *        here.
 * \return Error code.
 *
 * Time complexity: O(n(|V|+|E|)), |V| is the number of vertices, |E|
 * is the number of edges, n is the number vertices in the candidate
 * separator.
 *
 * \example examples/simple/igraph_is_minimal_separator.c
 */

igraph_error_t igraph_is_minimal_separator(const igraph_t *graph,
                                const igraph_vs_t candidate,
                                igraph_bool_t *res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_bool_t removed;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t neis;
    igraph_integer_t candsize;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_vit_create(graph, candidate, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    candsize = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_bool_init(&removed, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* Is it a separator at all? */
    IGRAPH_CHECK(igraph_i_is_separator(graph, &vit, -1, res, &removed,
                                       &Q, &neis, no_of_nodes));
    if (!(*res)) {
        /* Not a separator at all, nothing to do, *res is already set */
    } else if (candsize == 0) {
        /* Nothing to do, minimal, *res is already set */
    } else {
        /* General case, we need to remove each vertex from 'candidate'
         * and check whether the remainder is a separator. If this is
         * false for all vertices, then 'candidate' is a minimal
         * separator.
         */
        igraph_integer_t i;
        for (i = 0, *res = 0; i < candsize && (!*res); i++) {
            igraph_vector_bool_null(&removed);
            IGRAPH_CHECK(igraph_i_is_separator(graph, &vit, i, res, &removed,
                                               &Q, &neis, no_of_nodes));
        }
        (*res) = (*res) ? 0 : 1;    /* opposite */
    }

    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&Q);
    igraph_vector_bool_destroy(&removed);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/* --------------------------------------------------------------------*/

#define UPDATEMARK() do {                              \
        (*mark)++;                         \
        if (!(*mark)) {                    \
            igraph_vector_int_null(leaveout);                \
            (*mark)=1;                       \
        }                                                  \
    } while (0)

static igraph_error_t igraph_i_connected_components_leaveout(const igraph_adjlist_t *adjlist,
                                      igraph_vector_int_t *components,
                                      igraph_vector_int_t *leaveout,
                                      igraph_integer_t *mark,
                                      igraph_dqueue_int_t *Q) {

    /* Another trick: we use the same 'leaveout' vector to mark the
     * vertices that were already found in the BFS
     */

    igraph_integer_t i, no_of_nodes = igraph_adjlist_size(adjlist);

    igraph_dqueue_int_clear(Q);
    igraph_vector_int_clear(components);

    for (i = 0; i < no_of_nodes; i++) {

        if (VECTOR(*leaveout)[i] == *mark) {
            continue;
        }

        VECTOR(*leaveout)[i] = *mark;
        IGRAPH_CHECK(igraph_dqueue_int_push(Q, i));
        IGRAPH_CHECK(igraph_vector_int_push_back(components, i));

        while (!igraph_dqueue_int_empty(Q)) {
            igraph_integer_t act_node = igraph_dqueue_int_pop(Q);
            igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, act_node);
            igraph_integer_t j, n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t nei = VECTOR(*neis)[j];
                if (VECTOR(*leaveout)[nei] == *mark) {
                    continue;
                }
                IGRAPH_CHECK(igraph_dqueue_int_push(Q, nei));
                VECTOR(*leaveout)[nei] = *mark;
                IGRAPH_CHECK(igraph_vector_int_push_back(components, nei));
            }
        }

        IGRAPH_CHECK(igraph_vector_int_push_back(components, -1));
    }

    UPDATEMARK();

    return IGRAPH_SUCCESS;
}

static igraph_bool_t igraph_i_separators_is_not_seen_yet(
    const igraph_vector_int_list_t *comps, const igraph_vector_int_t *newc
) {

    igraph_integer_t co, nocomps = igraph_vector_int_list_size(comps);

    for (co = 0; co < nocomps; co++) {
        igraph_vector_int_t *act = igraph_vector_int_list_get_ptr(comps, co);
        if (igraph_vector_int_all_e(act, newc)) {
            return false;
        }
    }

    /* If not found, then it is new */
    return true;
}

static igraph_error_t igraph_i_separators_store(igraph_vector_int_list_t *separators,
                                     const igraph_adjlist_t *adjlist,
                                     igraph_vector_int_t *components,
                                     igraph_vector_int_t *leaveout,
                                     igraph_integer_t *mark,
                                     igraph_vector_int_t *sorter) {

    /* We need to store N(C), the neighborhood of C, but only if it is
     * not already stored among the separators.
     */

    igraph_integer_t cptr = 0, next, complen = igraph_vector_int_size(components);

    while (cptr < complen) {
        igraph_integer_t saved = cptr;
        igraph_vector_int_clear(sorter);

        /* Calculate N(C) for the next C */

        while ( (next = VECTOR(*components)[cptr++]) != -1) {
            VECTOR(*leaveout)[next] = *mark;
        }
        cptr = saved;

        while ( (next = VECTOR(*components)[cptr++]) != -1) {
            igraph_vector_int_t *neis = igraph_adjlist_get(adjlist, next);
            igraph_integer_t j, nn = igraph_vector_int_size(neis);
            for (j = 0; j < nn; j++) {
                igraph_integer_t nei = VECTOR(*neis)[j];
                if (VECTOR(*leaveout)[nei] != *mark) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(sorter, nei));
                    VECTOR(*leaveout)[nei] = *mark;
                }
            }
        }
        igraph_vector_int_sort(sorter);

        UPDATEMARK();

        /* Add it to the list of separators, if it is new */
        if (igraph_i_separators_is_not_seen_yet(separators, sorter)) {
            IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(separators, sorter));
        }
    } /* while cptr < complen */

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_all_minimal_st_separators
 * \brief List all vertex sets that are minimal (s,t) separators for some s and t.
 *
 * This function lists all vertex sets that are minimal (s,t)
 * separators for some (s,t) vertex pair.
 *
 * </para><para>
 * Note that some vertex sets returned by this function may not be minimal
 * with respect to disconnecting the graph (or increasing the number of
 * connected components). Take for example the 5-vertex graph with edges
 * <code>0-1-2-3-4-1</code>. This function returns the vertex sets
 * <code>{1}</code>, <code>{2,4}</code> and <code>{1,3}</code>.
 * Notice that <code>{1,3}</code> is not minimal with respect to disconnecting
 * the graph, as <code>{1}</code> would be sufficient for that. However, it is
 * minimal with respect to separating vertices \c 2 and \c 4.
 *
 * </para><para>
 * See more about the implemented algorithm in
 * Anne Berry, Jean-Paul Bordat and Olivier Cogis: Generating All the
 * Minimal Separators of a Graph, In: Peter Widmayer, Gabriele Neyer
 * and Stephan Eidenbenz (editors): Graph-theoretic concepts in
 * computer science, 1665, 167--172, 1999. Springer.
 * https://doi.org/10.1007/3-540-46784-X_17
 *
 * \param graph The input graph. It may be directed, but edge
 *        directions are ignored.
 * \param separators An initialized pointer vector, the separators
 *        are stored here. It is a list of pointers to <type>igraph_vector_int_t</type>
 *        objects. Each vector will contain the ids of the vertices in
 *        the separator.
 *        To free all memory allocated for \p separators, you need call
 *        \ref igraph_vector_destroy() and then \ref igraph_free() on
 *        each element, before destroying the pointer vector itself.
 * \return Error code.
 *
 * \sa \ref igraph_minimum_size_separators()
 *
 * Time complexity: O(n|V|^3), |V| is the number of vertices, n is the
 * number of separators.
 *
 * \example examples/simple/igraph_minimal_separators.c
 */

igraph_error_t igraph_all_minimal_st_separators(
    const igraph_t *graph, igraph_vector_int_list_t *separators
) {

    /*
     * Some notes about the tricks used here. For finding the components
     * of the graph after removing some vertices, we do the
     * following. First we mark the vertices with the actual mark stamp
     * (mark), then run breadth-first search on the graph, but not
     * considering the marked vertices. Then we increase the mark. If
     * there is integer overflow here, then we zero out the mark and set
     * it to one. (We might as well just always zero it out.)
     *
     * For each separator the vertices are stored in vertex ID order.
     * This facilitates the comparison of the separators when we find a
     * potential new candidate.
     *
     * To keep track of which separator we already used as a basis, we
     * keep a boolean vector (already_tried). The try_next pointer show
     * the next separator to try as a basis.
     */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t leaveout;
    igraph_vector_bool_t already_tried;
    igraph_integer_t try_next = 0;
    igraph_integer_t mark = 1;
    igraph_integer_t v;

    igraph_adjlist_t adjlist;
    igraph_vector_int_t components;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t sorter;

    igraph_vector_int_list_clear(separators);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&leaveout, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_bool_init(&already_tried, 0));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &already_tried);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&components, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&components, no_of_nodes * 2));
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    IGRAPH_CHECK(igraph_dqueue_int_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &Q);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sorter, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&sorter, no_of_nodes));

    /* ---------------------------------------------------------------
     * INITIALIZATION, we check whether the neighborhoods of the
     * vertices separate the graph. The ones that do will form the
     * initial basis.
     */

    for (v = 0; v < no_of_nodes; v++) {

        /* Mark v and its neighbors */
        igraph_vector_int_t *neis = igraph_adjlist_get(&adjlist, v);
        igraph_integer_t i, n = igraph_vector_int_size(neis);
        VECTOR(leaveout)[v] = mark;
        for (i = 0; i < n; i++) {
            igraph_integer_t nei = VECTOR(*neis)[i];
            VECTOR(leaveout)[nei] = mark;
        }

        /* Find the components */
        IGRAPH_CHECK(igraph_i_connected_components_leaveout(
            &adjlist, &components, &leaveout, &mark, &Q));

        /* Store the corresponding separators, N(C) for each component C */
        IGRAPH_CHECK(igraph_i_separators_store(separators, &adjlist, &components,
                                               &leaveout, &mark, &sorter));

    }

    /* ---------------------------------------------------------------
     * GENERATION, we need to use all already found separators as
     * basis and see if they generate more separators
     */

    while (try_next < igraph_vector_int_list_size(separators)) {
        /* copy "basis" out of the vector_list because we are going to
         * mutate the vector_list later, and this can potentially invalidate
         * the pointer */
        igraph_vector_int_t basis = *(igraph_vector_int_list_get_ptr(separators, try_next));
        igraph_integer_t b, basislen = igraph_vector_int_size(&basis);
        for (b = 0; b < basislen; b++) {

            /* Remove N(x) U basis */
            igraph_integer_t x = VECTOR(basis)[b];
            igraph_vector_int_t *neis = igraph_adjlist_get(&adjlist, x);
            igraph_integer_t i, n = igraph_vector_int_size(neis);
            for (i = 0; i < basislen; i++) {
                igraph_integer_t sn = VECTOR(basis)[i];
                VECTOR(leaveout)[sn] = mark;
            }
            for (i = 0; i < n; i++) {
                igraph_integer_t nei = VECTOR(*neis)[i];
                VECTOR(leaveout)[nei] = mark;
            }

            /* Find the components */
            IGRAPH_CHECK(igraph_i_connected_components_leaveout(
                &adjlist, &components, &leaveout, &mark, &Q));

            /* Store the corresponding separators, N(C) for each component C */
            IGRAPH_CHECK(igraph_i_separators_store(separators, &adjlist,
                                                   &components, &leaveout, &mark,
                                                   &sorter));
        }

        try_next++;
    }

    /* --------------------------------------------------------------- */

    igraph_vector_int_destroy(&sorter);
    igraph_dqueue_int_destroy(&Q);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&components);
    igraph_vector_bool_destroy(&already_tried);
    igraph_vector_int_destroy(&leaveout);
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}

#undef UPDATEMARK

static igraph_error_t igraph_i_minimum_size_separators_append(
    igraph_vector_int_list_t *old, igraph_vector_int_list_t *new
) {

    igraph_integer_t olen = igraph_vector_int_list_size(old);
    igraph_integer_t j;

    while (!igraph_vector_int_list_empty(new)) {
        igraph_vector_int_t *newvec = igraph_vector_int_list_tail_ptr(new);

        /* Check whether the separator is already in `old' */
        for (j = 0; j < olen; j++) {
            igraph_vector_int_t *oldvec = igraph_vector_int_list_get_ptr(old, j);
            if (igraph_vector_int_all_e(oldvec, newvec)) {
                break;
            }
        }

        if (j == olen) {
            /* We have found a new separator, append it to `old' */
            /* TODO: we should have a more efficient method for moving a vector
             * from one vector_list to another */
            IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(old, newvec));
            olen++;
        }

        igraph_vector_int_list_discard_back(new);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_minimum_size_separators_topkdeg(
    const igraph_t *graph, igraph_vector_int_t *res, igraph_integer_t k
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t deg, order;
    igraph_integer_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_ALL,
                               /*loops=*/ 0));

    IGRAPH_CHECK(igraph_vector_int_order1(&deg, &order, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_int_resize(res, k));
    for (i = 0; i < k; i++) {
        VECTOR(*res)[i] = VECTOR(order)[no_of_nodes - 1 - i];
    }

    igraph_vector_int_destroy(&order);
    igraph_vector_int_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_minimum_size_separators
 * \brief Find all minimum size separating vertex sets.
 *
 * This function lists all separator vertex sets of minimum size.
 * A vertex set is a separator if its removal disconnects the graph.
 *
 * </para><para>
 * The implementation is based on the following paper:
 * Arkady Kanevsky: Finding all minimum-size separating vertex sets in
 * a graph, Networks 23, 533--541, 1993.
 *
 * \param graph The input graph, which must be undirected.
 * \param separators An initialized list of integer vectors, the separators
 *        are stored here. It is a list of pointers to igraph_vector_int_t
 *        objects. Each vector will contain the IDs of the vertices in
 *        the separator. The separators are returned in an arbitrary order.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_minimum_size_separators.c
 */

igraph_error_t igraph_minimum_size_separators(
    const igraph_t *graph, igraph_vector_int_list_t *separators
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t conn;
    igraph_vector_int_t X;
    igraph_integer_t i, j, k, n;
    igraph_bool_t issepX;
    igraph_t Gbar;
    igraph_vector_t phi;
    igraph_t graph_copy;
    igraph_vector_t capacity;
    igraph_maxflow_stats_t stats;

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Minimum size separators currently only works on undirected graphs",
                     IGRAPH_EINVAL);
    }

    igraph_vector_int_list_clear(separators);

    /* ---------------------------------------------------------------- */
    /* 1 Find the vertex connectivity of 'graph' */
    IGRAPH_CHECK(igraph_vertex_connectivity(graph, &conn, /* checks= */ true));
    k = conn;

    /* Special cases for low connectivity, two exits here! */
    if (conn == 0) {
        /* Nothing to do */
        return IGRAPH_SUCCESS;
    } else if (conn == 1) {
        igraph_vector_int_t ap;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&ap, 0);
        IGRAPH_CHECK(igraph_articulation_points(graph, &ap));
        n = igraph_vector_int_size(&ap);
        IGRAPH_CHECK(igraph_vector_int_list_resize(separators, n));
        for (i = 0; i < n; i++) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(separators, i);
            IGRAPH_CHECK(igraph_vector_int_push_back(v, VECTOR(ap)[i]));
        }
        igraph_vector_int_destroy(&ap);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    } else if (conn == no_of_nodes - 1) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(separators, no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(separators, i);
            IGRAPH_CHECK(igraph_vector_int_resize(v, no_of_nodes - 1));
            for (j = 0, k = 0; j < no_of_nodes; j++) {
                if (j != i) {
                    VECTOR(*v)[k++] = j;
                }
            }
        }
        return IGRAPH_SUCCESS;
    }

    /* Work on a copy of 'graph' */
    IGRAPH_CHECK(igraph_copy(&graph_copy, graph));
    IGRAPH_FINALLY(igraph_destroy, &graph_copy);
    IGRAPH_CHECK(igraph_simplify(&graph_copy, /* multiple */ true, /* loops */ true, NULL));

    /* ---------------------------------------------------------------- */
    /* 2 Find k vertices with the largest degrees (x1;..,xk). Check
       if these k vertices form a separating k-set of G */
    IGRAPH_CHECK(igraph_vector_int_init(&X, conn));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &X);
    IGRAPH_CHECK(igraph_i_minimum_size_separators_topkdeg(&graph_copy, &X, k));
    IGRAPH_CHECK(igraph_is_separator(&graph_copy, igraph_vss_vector(&X),
                                     &issepX));
    if (issepX) {
        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(separators, &X));
    }

    /* Create Gbar, the Even-Tarjan reduction of graph */
    IGRAPH_VECTOR_INIT_FINALLY(&capacity, 0);
    IGRAPH_CHECK(igraph_even_tarjan_reduction(&graph_copy, &Gbar, &capacity));
    IGRAPH_FINALLY(igraph_destroy, &Gbar);

    IGRAPH_VECTOR_INIT_FINALLY(&phi, no_of_edges);

    /* ---------------------------------------------------------------- */
    /* 3 If v[j] != x[i] and v[j] is not adjacent to x[i] then */
    for (i = 0; i < k; i++) {

        IGRAPH_ALLOW_INTERRUPTION();

        for (j = 0; j < no_of_nodes; j++) {
            igraph_integer_t ii = VECTOR(X)[i];
            igraph_real_t phivalue;
            igraph_bool_t conn;

            if (ii == j) {
                continue;    /* the same vertex */
            }
            IGRAPH_CHECK(igraph_are_connected(&graph_copy, ii, j, &conn));
            if (conn) {
                continue;    /* they are connected */
            }

            /* --------------------------------------------------------------- */
            /* 4 Compute a maximum flow phi in Gbar from x[i] to v[j].
            If |phi|=k, then */
            IGRAPH_CHECK(igraph_maxflow(&Gbar, &phivalue, &phi, /*cut=*/ 0,
                                        /*partition=*/ 0, /*partition2=*/ 0,
                                        /* source= */ ii + no_of_nodes,
                                        /* target= */ j,
                                        &capacity, &stats));

            if (phivalue == k) {

                /* ------------------------------------------------------------- */
                /* 5-6-7. Find all k-sets separating x[i] and v[j]. */
                igraph_vector_int_list_t stcuts;
                IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&stcuts, 0);
                IGRAPH_CHECK(igraph_all_st_mincuts(&Gbar, /*value=*/ 0,
                                                   /*cuts=*/ &stcuts,
                                                   /*partition1s=*/ 0,
                                                   /*source=*/ ii + no_of_nodes,
                                                   /*target=*/ j,
                                                   /*capacity=*/ &capacity));

                IGRAPH_CHECK(igraph_i_minimum_size_separators_append(separators, &stcuts));
                igraph_vector_int_list_destroy(&stcuts);
                IGRAPH_FINALLY_CLEAN(1);

            } /* if phivalue == k */

            /* --------------------------------------------------------------- */
            /* 8 Add edge (x[i],v[j]) to G. */
            IGRAPH_CHECK(igraph_add_edge(&graph_copy, ii, j));
            IGRAPH_CHECK(igraph_add_edge(&Gbar, ii + no_of_nodes, j));
            IGRAPH_CHECK(igraph_add_edge(&Gbar, j + no_of_nodes, ii));
            IGRAPH_CHECK(igraph_vector_push_back(&capacity, no_of_nodes));
            IGRAPH_CHECK(igraph_vector_push_back(&capacity, no_of_nodes));

        } /* for j<no_of_nodes */
    } /* for i<k */

    igraph_vector_destroy(&phi);
    igraph_destroy(&Gbar);
    igraph_vector_destroy(&capacity);
    igraph_vector_int_destroy(&X);
    igraph_destroy(&graph_copy);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}
