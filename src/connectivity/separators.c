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


/**
 * \function igraph_i_is_separator
 *
 * Checks if vertex set \c S is a separator. Optionally, also checks if it
 * is a _minimal_ separator.
 *
 * \param graph The graph, edge directions are ignored.
 * \param S Candidate vertex set.
 * \param is_separator Pointer to boolean, is S a separator?
 * \param is_minimal If not \c NULL, it is also checked whether the separator
 *    is minimal. This takes additional time.  If S is not a separator, this is
 *    set to false
 */
static igraph_error_t igraph_i_is_separator(
    const igraph_t *graph,
    igraph_vs_t S,
    igraph_bool_t *is_separator,
    igraph_bool_t *is_minimal
) {

    const igraph_integer_t vcount = igraph_vcount(graph);

    /* mark[v] means:
     *
     * bit 0: Was v visited?
     * bit 1: Is v in S?
     * bit 2: Used to keep track of which vertices were reachable in S
     *        from its neighbourhood in the minimal separator test.
     */
    igraph_vector_char_t mark;

    igraph_vector_int_t neis;
    igraph_dqueue_int_t Q;
    igraph_vit_t vit;
    igraph_integer_t S_size = 0, S_visited_count = 0;

    IGRAPH_CHECK(igraph_vit_create(graph, S, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&mark, vcount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 10);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, 100);

#define VISITED(x) (VECTOR(mark)[x] & 1) /* Was x visited? */
#define SET_VISITED(x) (VECTOR(mark)[x] |= 1) /* Mark x as visited. */
#define IN_S(x) (VECTOR(mark)[x] & 2) /* Is x in S? */
#define SET_IN_S(x) (VECTOR(mark)[x] |= 2) /* Mark x as being in S. */

    /* Mark and count vertices in S, taking care not to double-count
     * when duplicate vertices were passed in. */
    for (; ! IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        const igraph_integer_t u = IGRAPH_VIT_GET(vit);
        if (! IN_S(u)) {
            SET_IN_S(u);
            S_size++;
        }
    }

    /* If S contains more than |V| - 2 vertices, it is not a separator. */
    if (S_size > vcount-2) {
        *is_separator = false;
        goto done;
    }

    /* Assume that S is not a separator until proven otherwise. */
    *is_separator = false;

    for (IGRAPH_VIT_RESET(vit); ! IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {

        const igraph_integer_t u = IGRAPH_VIT_GET(vit);
        if (VISITED(u)) {
            continue;
        }

        /* For the sake of the following discussion, and counting of degrees,
         * let us consider each undirected edge as a pair of directed ones.
         */

        /* We start at a vertex in S, and traverse the graph with the following
         * restriction:
         *
         * Once we exited S through an edge, we may not exit through
         * any other edge again (but we may re-enter S). With this traversal,
         * we can determine: once we find ourselves outside of S, are there
         * any vertices which we can only reach by passing through S again?
         *
         * Theorem: The visited vertices in S constitute a separating set if and
         * only if when the traversal is complete, some of them still have
         * untraversed out-edges.
         *
         * Note that this refers only to the subset of S visited during this
         * traversal, not the rest of S. If there are remaining vertices in S,
         * they may constitute a separating set as well. We cannot conclude that
         * S is NOT a separator until all its vertices have been considered, hence
         * the outer for-loop.
         */

        /* Sum of in-degrees of visited vertices in S. */
        igraph_integer_t degsum = 0;

        /* Number of in-edges of vertices in S that were traversed. */
        igraph_integer_t edgecount = 0;

        /* Have we already traversed an edge leaving S? */
        igraph_bool_t exited = false;

        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, u));

        while (! igraph_dqueue_int_empty(&Q)) {
            const igraph_integer_t v = igraph_dqueue_int_pop(&Q);
            if (VISITED(v)) {
                continue;
            }

            SET_VISITED(v);
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, IGRAPH_ALL));

            const igraph_integer_t dv = igraph_vector_int_size(&neis);

            if (IN_S(v)) {
                degsum += dv;
                S_visited_count++;
            }

            for (igraph_integer_t i=0; i < dv; i++) {
                const igraph_integer_t w = VECTOR(neis)[i];

                /* Decide whether to traverse the v -> w edge. */
                if (!exited || !IN_S(v) || IN_S(w)) {
                    if (IN_S(w)) {
                        edgecount++;
                    }

                    if (!VISITED(w)) {
                        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, w));
                    }

                    if (!exited && IN_S(v) && !IN_S(w)) {
                        exited = true;
                    }
                }
            }
        }

        /* If some incident edges of visited vertices in S are
         * still untraversed, then S is a separator. We are done. */
        if (degsum > edgecount) {
            *is_separator = true;
            break;
        }
    }

done:

    /* Optionally, check if S is also a _minimal_ separator. */
    if (is_minimal != NULL) {
        /* Be optimistic, and assume that if S was found to be a separator,
         * it is a minimal separator. */
        *is_minimal = *is_separator;

        /* If S was proven to be a separator before visiting all of its
         * vertices, it is not minimal. We are done. */
        if (*is_minimal && S_visited_count < S_size) {
            *is_minimal = false;
        }

        /* The subset of S with untraversed out-edges forms a separator.
         * Therefore, if S contains vertices with no untraversed out-edges,
         * S is not minimal.
         *
         * Suppose it doesn't contain any.
         *
         * Then for S to be minimal, each of its vertices should be reachable
         * from any vertex in the unvisited neighbourhood of S. In order to
         * separate a vertex in this neighbourhood, we need to cut precisely
         * those vertices in S which are reachable from it.
         *
         * The check below verifies BOTH of the above conditions by traversing
         * the graph from each unvisited neighbour of S with the constraint
         * of never entering S.
         */
        if (*is_minimal) {
            igraph_vector_int_t Sneis;
            IGRAPH_VECTOR_INT_INIT_FINALLY(&Sneis, 10);

            /* If the 2nd bit of mark[v] is the same as 'bit', it indicates
             * that v in S was reached in the current testing round.
             * We flip 'bit' between testing rounds, assuming that the previous
             * round reached all of S. */
            igraph_bool_t bit = true;

#define REACHED(x) (!(VECTOR(mark)[x] & 4) == !bit) /* Was x in S reached already? */
#define FLIP_REACHED(x) (VECTOR(mark)[x] ^= 4) /* Flip the reachability status of x in S. */

            for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
                const igraph_integer_t u = IGRAPH_VIT_GET(vit);
                IGRAPH_CHECK(igraph_neighbors(graph, &Sneis, u, IGRAPH_ALL));
                const igraph_integer_t du = igraph_vector_int_size(&Sneis);
                for (igraph_integer_t i=0; i < du; i++) {

                    igraph_integer_t v = VECTOR(Sneis)[i];
                    if (VISITED(v)) {
                        continue;
                    }

                    /* How many vertices in S were reachable from u? */
                    igraph_integer_t S_reached = 0;
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, v));
                    while (! igraph_dqueue_int_empty(&Q)) {
                        v = igraph_dqueue_int_pop(&Q);
                        if (VISITED(v)) {
                            continue;
                        }

                        SET_VISITED(v);
                        IGRAPH_CHECK(igraph_neighbors(graph, &neis, v, IGRAPH_ALL));

                        const igraph_integer_t dv = igraph_vector_int_size(&neis);

                        for (igraph_integer_t j=0; j < dv; j++) {
                            const igraph_integer_t w = VECTOR(neis)[j];

                            if (! VISITED(w)) {
                                IGRAPH_CHECK(igraph_dqueue_int_push(&Q, w));
                            } else if (IN_S(w) && !REACHED(w)) {
                                S_reached++;
                                FLIP_REACHED(w); /* set as reachable */
                            }
                        }
                    }

                    bit = !bit;

                    if (S_reached < S_size) {
                        *is_minimal = false;
                        break;
                    }
                }

                if (! *is_minimal) {
                    break;
                }
            }

            igraph_vector_int_destroy(&Sneis);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

#undef REACHED
#undef FLIP_REACHED
#undef VISITED
#undef SET_VISITED
#undef IN_S
#undef SET_IN_S


    igraph_dqueue_int_destroy(&Q);
    igraph_vector_int_destroy(&neis);
    igraph_vector_char_destroy(&mark);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_is_separator
 * \brief Would removing this set of vertices disconnect the graph?
 *
 * A vertex set \c S is a separator if there are vertices \c u and \c v
 * in the graph such that all paths between \c u and \c v pass through
 * some vertices in \c S.
 *
 * \param graph The input graph. It may be directed, but edge
 *        directions are ignored.
 * \param candidate The candidate separator.
 * \param res Pointer to a boolean variable, the result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number vertices and edges.
 *
 * \example examples/simple/igraph_is_separator.c
 */

igraph_error_t igraph_is_separator(const igraph_t *graph,
                                   const igraph_vs_t candidate,
                                   igraph_bool_t *res) {

    return igraph_i_is_separator(graph, candidate, res, NULL);
}

/**
 * \function igraph_is_minimal_separator
 * \brief Decides whether a set of vertices is a minimal separator.
 *
 * A vertex separator \c S is minimal is no proper subset of \c S
 * is also a separator.
 *
 * \param graph The input graph. It may be directed, but edge
 *        directions are ignored.
 * \param candidate The candidate minimal separators.
 * \param res Pointer to a boolean variable, the result is stored
 *        here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number vertices and edges.
 *
 * \example examples/simple/igraph_is_minimal_separator.c
 */

igraph_error_t igraph_is_minimal_separator(const igraph_t *graph,
                                           const igraph_vs_t candidate,
                                           igraph_bool_t *res) {

    igraph_bool_t is_separator;
    return igraph_i_is_separator(graph, candidate, &is_separator, res);
}

/* --------------------------------------------------------------------*/

#define UPDATEMARK() do { \
        (*mark)++; \
        if (!(*mark)) { \
            igraph_vector_int_null(leaveout); \
            (*mark)=1; \
        } \
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
        /* TODO: Is there a cleaner way to avoid empty separators,
         * or is this an inherent limitation of the algorithm?
         * See https://github.com/igraph/igraph/issues/2517 */
        if (
            igraph_vector_int_size(sorter) > 0 &&
            igraph_i_separators_is_not_seen_yet(separators, sorter)
        ) {
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
 *    directions are ignored.
 * \param separators Pointer to a list of integer vectors, the separators
 *    will be stored here.
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
     * The try_next pointer show the next separator to try as a basis.
     */

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t leaveout;
    igraph_integer_t try_next = 0;
    igraph_integer_t mark = 1;
    igraph_integer_t v;

    igraph_adjlist_t adjlist;
    igraph_vector_int_t components;
    igraph_dqueue_int_t Q;
    igraph_vector_int_t sorter;

    igraph_vector_int_list_clear(separators);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&leaveout, no_of_nodes);
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
    igraph_vector_int_destroy(&leaveout);
    IGRAPH_FINALLY_CLEAN(5);

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
            /* We have found a new separator, append it to `old'. We do it by
             * extending it with an empty vector and then swapping it with
             * the new vector to be appended */
            igraph_vector_int_t *oldvec;
            IGRAPH_CHECK(igraph_vector_int_list_push_back_new(old, &oldvec));
            igraph_vector_int_swap(oldvec, newvec);
            olen++;
        }

        igraph_vector_int_list_discard_back(new);
    }

    return IGRAPH_SUCCESS;
}

/**
 * Finds the k largest degree vertices.
 */
static igraph_error_t igraph_i_minimum_size_separators_topkdeg(
    const igraph_t *graph, igraph_vector_int_t *res, const igraph_integer_t k
) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t deg, order;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);

    /* It is assumed that this function receives only simple graphs, so we can use the
     * faster IGRAPH_LOOPS here instead of the slower IGRAPH_NO_LOOPS. */
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS));

    IGRAPH_CHECK(igraph_vector_int_order1(&deg, &order, no_of_nodes));
    IGRAPH_CHECK(igraph_vector_int_resize(res, k));
    for (igraph_integer_t i = 0; i < k; i++) {
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
 * https://doi.org/10.1002/net.3230230604
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
    igraph_integer_t k, n;
    igraph_bool_t issepX;
    igraph_t Gbar;
    igraph_vector_t phi;
    igraph_t graph_copy;
    igraph_vector_t capacity;
    igraph_maxflow_stats_t stats;

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("Minimum size separators currently only works on undirected graphs.",
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
        for (igraph_integer_t i = 0; i < n; i++) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(separators, i);
            IGRAPH_CHECK(igraph_vector_int_push_back(v, VECTOR(ap)[i]));
        }
        igraph_vector_int_destroy(&ap);
        IGRAPH_FINALLY_CLEAN(1);
        return IGRAPH_SUCCESS;
    } else if (conn == no_of_nodes - 1) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(separators, no_of_nodes));
        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(separators, i);
            IGRAPH_CHECK(igraph_vector_int_resize(v, no_of_nodes - 1));
            for (igraph_integer_t j = 0, k = 0; j < no_of_nodes; j++) {
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
    for (igraph_integer_t i = 0; i < k; i++) {

        IGRAPH_ALLOW_INTERRUPTION();

        for (igraph_integer_t j = 0; j < no_of_nodes; j++) {
            igraph_integer_t xi = VECTOR(X)[i];
            igraph_real_t phivalue;
            igraph_bool_t conn;

            if (xi == j) {
                continue;    /* the same vertex */
            }
            IGRAPH_CHECK(igraph_are_adjacent(&graph_copy, xi, j, &conn));
            if (conn) {
                continue;    /* they are connected */
            }

            /* --------------------------------------------------------------- */
            /* 4 Compute a maximum flow phi in Gbar from x[i] to v[j].
            If |phi|=k, then */
            IGRAPH_CHECK(igraph_maxflow(&Gbar, &phivalue, &phi, /*cut=*/ NULL,
                                        /*partition=*/ NULL, /*partition2=*/ NULL,
                                        /* source= */ xi + no_of_nodes,
                                        /* target= */ j,
                                        &capacity, &stats));

            if (phivalue == k) {

                /* ------------------------------------------------------------- */
                /* 5-6-7. Find all k-sets separating x[i] and v[j]. */
                igraph_vector_int_list_t stcuts;
                IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&stcuts, 0);
                IGRAPH_CHECK(igraph_all_st_mincuts(&Gbar, /*value=*/ NULL,
                                                   /*cuts=*/ &stcuts,
                                                   /*partition1s=*/ NULL,
                                                   /*source=*/ xi + no_of_nodes,
                                                   /*target=*/ j,
                                                   /*capacity=*/ &capacity));

                IGRAPH_CHECK(igraph_i_minimum_size_separators_append(separators, &stcuts));
                igraph_vector_int_list_destroy(&stcuts);
                IGRAPH_FINALLY_CLEAN(1);

            } /* if phivalue == k */

            /* --------------------------------------------------------------- */
            /* 8 Add edge (x[i],v[j]) to G. */
            IGRAPH_CHECK(igraph_add_edge(&graph_copy, xi, j));
            IGRAPH_CHECK(igraph_add_edge(&Gbar, xi + no_of_nodes, j));
            IGRAPH_CHECK(igraph_add_edge(&Gbar, j + no_of_nodes, xi));
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
