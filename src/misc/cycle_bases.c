/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_cycles.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_error.h"
#include "igraph_interface.h"

#include "core/interruption.h"
#include "misc/order_cycle.h"

/**** Fundamental cycles *****/

/* Computes fundamental cycles for the connected component containing 'start_vid',
 * and appends them to 'result'.
 *
 * 'visited' must be a vector of length igraph_vcount(graph).
 * visited[u] will be set to mark+1 or mark+2 for each visited vertex 'u'.
 * No elements of 'visited' must have these values when calling this function.
 * 'mark' can be specified in order to be able to re-use a 'visited' vector
 * multiple times without having to re-set all its elements.
 *
 * During the operation of the function, mark+1 indicates that a vertex has been
 * queued for processing, but not processed yet. mark+2 indicates that it has
 * been processed.
 */
static igraph_error_t
igraph_i_fundamental_cycles_bfs(
        const igraph_t *graph,
        igraph_vector_int_list_t *result,
        igraph_integer_t start_vid,
        igraph_integer_t bfs_cutoff,
        const igraph_inclist_t *inclist,
        igraph_vector_int_t *visited,
        igraph_integer_t mark /* mark used in 'visited' */) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vector_int_t pred_edge;
    igraph_vector_int_t u_back, v_back;

    if (start_vid < 0 || start_vid >= no_of_nodes) {
        IGRAPH_ERROR("Invalid starting vertex id.", IGRAPH_EINVAL);
    }

    if (mark > IGRAPH_INTEGER_MAX - 2) {
        IGRAPH_ERROR("Graph too large for cycle basis.", IGRAPH_EOVERFLOW);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&pred_edge, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&u_back, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&v_back, 0);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 0);

    IGRAPH_CHECK(igraph_dqueue_int_push(&q, start_vid)); /* vertex id */
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0)); /* distance from start_vid*/
    VECTOR(*visited)[start_vid] = mark + 1; /* mark as seen */
    VECTOR(pred_edge)[start_vid] = -1; /* non-valid predecessor edge id for root vertex */

    while (! igraph_dqueue_int_empty(&q)) {
        igraph_integer_t v = igraph_dqueue_int_pop(&q);
        igraph_integer_t vdist = igraph_dqueue_int_pop(&q);

        igraph_vector_int_t *incs = igraph_inclist_get(inclist, v);
        igraph_integer_t n = igraph_vector_int_size(incs);
        igraph_integer_t i, j;

        IGRAPH_ALLOW_INTERRUPTION();

        for (i=0; i < n; ++i) {
            igraph_integer_t e = VECTOR(*incs)[i];
            igraph_integer_t u = IGRAPH_OTHER(graph, e, v);

            if (e == VECTOR(pred_edge)[v]) {
                /* do not follow the edge through which we came to v */
                continue;
            }

            if (VECTOR(*visited)[u] == mark + 2) {
                /* u has already been processed */
                continue;
            } else if (VECTOR(*visited)[u] == mark + 1) {
                /* u has been seen but not yet processed */

                /* Found cycle edge u-v. Now we walk back up the BFS tree
                 * in order to find the common ancestor of u and v. We exploit
                 * that the distance of u from the start vertex is either the
                 * same as that of v, or one greater. */

                igraph_integer_t up = u, vp = v;
                igraph_integer_t u_back_len, v_back_len;
                igraph_vector_int_t cycle;

                IGRAPH_CHECK(igraph_vector_int_push_back(&v_back, e));
                for (;;) {
                    igraph_integer_t upe, vpe;

                    if (up == vp) {
                        break;
                    }

                    upe = VECTOR(pred_edge)[up];
                    IGRAPH_CHECK(igraph_vector_int_push_back(&u_back, upe));
                    up = IGRAPH_OTHER(graph, upe, up);

                    if (up == vp) {
                        break;
                    }

                    vpe = VECTOR(pred_edge)[vp];
                    IGRAPH_CHECK(igraph_vector_int_push_back(&v_back, vpe));
                    vp = IGRAPH_OTHER(graph, vpe, vp);
                }

                u_back_len = igraph_vector_int_size(&u_back);
                v_back_len = igraph_vector_int_size(&v_back);
                IGRAPH_VECTOR_INT_INIT_FINALLY(&cycle, u_back_len + v_back_len);

                for (j=0; j < v_back_len; ++j) {
                    VECTOR(cycle)[j] = VECTOR(v_back)[j];
                }
                for (j=0; j < u_back_len; ++j) {
                    VECTOR(cycle)[v_back_len + j] = VECTOR(u_back)[u_back_len - j - 1];
                }

                igraph_vector_int_clear(&v_back);
                igraph_vector_int_clear(&u_back);

                IGRAPH_CHECK(igraph_vector_int_list_push_back(result, &cycle));
                IGRAPH_FINALLY_CLEAN(1); /* pass ownership of 'cycle' to 'result' */
            } else {
                /* encountering u for the first time, queue it for processing */

                /* Only queue vertices with distance at most 'bfs_cutoff' from the root. */
                /* Negative 'bfs_cutoff' indicates no cutoff. */
                if (bfs_cutoff < 0 || vdist < bfs_cutoff) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, u));
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, vdist + 1));
                    VECTOR(*visited)[u] = mark + 1;
                    VECTOR(pred_edge)[u] = e;
                }
            }
        }

        VECTOR(*visited)[v] = mark + 2; /* mark v as processed */
    } /* ! igraph_dqueue_int_empty(&q) */

    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&v_back);
    igraph_vector_int_destroy(&u_back);
    igraph_vector_int_destroy(&pred_edge);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_fundamental_cycles
 * \brief Finds a fundamental cycle basis.
 *
 * \experimental
 *
 * This function computes a fundamental cycle basis associated with a breadth-first
 * search tree of the graph.
 *
 * </para><para>
 * Edge directions are ignored. Multi-edges and self-loops are supported.
 *
 * \param graph The graph object.
 * \param result An initialized integer vector list. The result will be stored here,
 *   each vector containing the edge IDs of a basis element.
 * \param start_vid If negative, a complete fundamental cycle basis is returned.
 *   If a vertex ID, the fundamental cycles associated with the BFS tree rooted
 *   in that vertex will be returned, only for the weakly connected component
 *   containing that vertex.
 * \param bfs_cutoff If negative, a complete cycle basis is returned. Otherwise, only
 *   cycles of length <code>2*bfs_cutoff + 1</code> or shorter are included. \p bfs_cutoff
 *   is used to limit the depth of the BFS tree when searching for cycle edges.
 * \param weights Currently unused.
 * \return Error code.
 *
 * Time complexity: O(|V| + |E|).
 */
igraph_error_t igraph_fundamental_cycles(const igraph_t *graph,
                                         igraph_vector_int_list_t *result,
                                         igraph_integer_t start_vid,
                                         igraph_integer_t bfs_cutoff,
                                         const igraph_vector_t *weights) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t estimated_rank;
    igraph_integer_t i;
    igraph_inclist_t inclist;
    igraph_vector_int_t visited; /* see comments before igraph_i_fundamental_cycles_bfs() */

    IGRAPH_UNUSED(weights);

    if (start_vid >= no_of_nodes) {
        IGRAPH_ERROR("Vertex id out of range.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&visited, no_of_nodes);

    /* Compute cycle rank assuming that the graph is connected. */
    estimated_rank = no_of_edges - no_of_nodes + 1;
    estimated_rank = estimated_rank < 0 ? 0 : estimated_rank;

    igraph_vector_int_list_clear(result);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(result, estimated_rank));

    if (start_vid < 0) {
        for (i=0; i < no_of_nodes; ++i) {
            if (! VECTOR(visited)[i]) {
                IGRAPH_CHECK(igraph_i_fundamental_cycles_bfs(graph, result, i, bfs_cutoff, &inclist,
                                                             &visited, /* mark */ 0));
            }
        }
    } else {
        IGRAPH_CHECK(igraph_i_fundamental_cycles_bfs(graph, result, start_vid, bfs_cutoff, &inclist,
                                                     &visited, /* mark */ 0));
    }

    igraph_vector_int_destroy(&visited);
    igraph_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/***** Minimum weight cycle basis *****/

/* In this implementation, the cycle vectors (basis elements) are stored as a sparse representation:
 * a sorted list of edge indices. */


/* qsort-compatible comparison for sparse cycle vectors: shorter ones come first, use lexicographic
 * order for equal length ones. Lexicographic order helps keep row insertion into the reduced matrix
 * efficient during Gaussian elimination, by ensuring that insertions usually happen near the end. */
static int cycle_cmp(const igraph_vector_int_t *v1, const igraph_vector_int_t *v2) {
    igraph_integer_t n1 = igraph_vector_int_size(v1), n2 = igraph_vector_int_size(v2);

    if (n1 < n2) {
        return -1;
    } else if (n1 > n2) {
        return 1;
    } else {
        return igraph_vector_int_lex_cmp(v1, v2);
    }
}

/* Adding cycle vectors produces the symmetric difference of the corresponding edge sets. */
static igraph_error_t cycle_add(const igraph_vector_int_t *a, const igraph_vector_int_t *b, igraph_vector_int_t *res) {
    igraph_integer_t na = igraph_vector_int_size(a), nb = igraph_vector_int_size(b);
    const igraph_integer_t *pa = VECTOR(*a), *pb = VECTOR(*b);
    const igraph_integer_t *pa_end = pa + na, *pb_end = pb + nb;

    igraph_vector_int_clear(res);
    for (;;) {
        while (pa != pa_end && pb != pb_end && *pa < *pb) {
            IGRAPH_CHECK(igraph_vector_int_push_back(res, *pa));
            pa++;
        }
        while (pa != pa_end && pb != pb_end && *pa == *pb) {
            pa++; pb++;
        }
        while (pa != pa_end && pb != pb_end && *pb < *pa) {
            IGRAPH_CHECK(igraph_vector_int_push_back(res, *pb));
            pb++;
        }
        if (pa == pa_end) {
            while (pb != pb_end) {
                IGRAPH_CHECK(igraph_vector_int_push_back(res, *pb));
                pb++;
            }
            break;
        }
        if (pb == pb_end) {
            while (pa != pa_end) {
                IGRAPH_CHECK(igraph_vector_int_push_back(res, *pa));
                pa++;
            }
            break;
        }
    }

    return IGRAPH_SUCCESS;
}


#define MATROW(m, i) (&VECTOR(m)[i])
#define MATEL(m, i, j) VECTOR(*MATROW(m, i))[j]

/* Gaussian elimination for sparse cycle vectors. 'reduced_matrix' is always maintained
 * in row-echelon form. This function decides if 'cycle' is linearly independent of this
 * matrix, and if not, it adds it to the matrix. */
static igraph_error_t gaussian_elimination(igraph_vector_int_list_t *reduced_matrix,
                                           const igraph_vector_int_t *cycle,
                                           igraph_bool_t *independent) {

    const igraph_integer_t nrow = igraph_vector_int_list_size(reduced_matrix);
    igraph_integer_t i;

    igraph_vector_int_t work, tmp;

    IGRAPH_CHECK(igraph_vector_int_init_copy(&work, cycle));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &work);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);

    for (i=0; i < nrow; ++i) {
        igraph_vector_int_t *row = MATROW(*reduced_matrix, i);

        if ( VECTOR(*row)[0] < VECTOR(work)[0] ) {
            continue;
        } else if ( VECTOR(*row)[0] == VECTOR(work)[0] ) {
            IGRAPH_CHECK(cycle_add(row, &work, &tmp));
            if (igraph_vector_int_empty(&tmp)) {
                *independent = false;
                igraph_vector_int_destroy(&work);
                igraph_vector_int_destroy(&tmp);
                IGRAPH_FINALLY_CLEAN(2);
                return IGRAPH_SUCCESS;
            }
            IGRAPH_CHECK(igraph_vector_int_swap(&work, &tmp));
        } else { /* VECTOR(*row)[0] > VECTOR(work)[0] */
            break;
        }
    }

    /* 'cycle' was linearly independent, insert new row into matrix */
    *independent = true;
    IGRAPH_CHECK(igraph_vector_int_list_insert(reduced_matrix, i, &work)); /* transfers ownership */

    igraph_vector_int_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(2); /* +1, transferring ownership of 'work' to 'reduced_matrix' */

    return IGRAPH_SUCCESS;
}

#undef MATEL
#undef MATROW


/**
 * \function igraph_minimum_cycle_basis
 * \brief Computes a minimum weight cycle basis.
 *
 * \experimental
 *
 * This function computes a minimum weight cycle basis of a graph. Currently,
 * a modified version of Horton's algorithm is used that allows for cutoffs.
 *
 * </para><para>
 * Edge directions are ignored. Multi-edges and self-loops are supported.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Horton, J. D. (1987)
 * A polynomial-time algorithm to find the shortest cycle basis of a graph,
 * SIAM Journal on Computing, 16 (2): 358–366.
 * https://doi.org/10.1137%2F0216026
 *
 * \param graph The graph object.
 * \param result An initialized integer vector list, the elements of the cycle
 *   basis will be stored here as vectors of edge IDs.
 * \param bfs_cutoff If negative, an exact minimum cycle basis is returned. Otherwise
 *   only those cycles in the result will be part of some minimum cycle basis which
 *   are of size <code>2*bfs_cutoff + 1</code> or smaller. Cycles longer than this limit
 *   may not be of the smallest possible size.
 *   \p bfs_cutoff is used to limit the depth of the BFS tree when computing candidate
 *   cycles. Specifying a bfs_cutoff can speed up the computation substantially.
 * \param complete Boolean value. Used only when \p bfs_cutoff was given.
 *   If true, a complete basis is returned. If false, only cycles not greater
 *   than <code>2*bfs_cutoff + 1</code> are returned. This may save computation
 *   time, however, the result will not span the entire cycle space.
 * \param use_cycle_order If true, each cycle is returned in natural order:
 *   the edge IDs will appear ordered along the cycle. This comes at a small
 *   performance cost. If false, no guarantees are given about the ordering
 *   of edge IDs within cycles. This parameter exists solely to control
 *   performance tradeoffs.
 * \param weights Currently unused.
 * \return Error code.
 *
 * Time complexity: TODO.
 */
igraph_error_t igraph_minimum_cycle_basis(const igraph_t *graph,
                                          igraph_vector_int_list_t *result,
                                          igraph_integer_t bfs_cutoff,
                                          igraph_bool_t complete,
                                          igraph_bool_t use_cycle_order,
                                          const igraph_vector_t *weights) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    const igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t rank;
    igraph_vector_int_list_t candidates;

    IGRAPH_UNUSED(weights);

    /* Compute candidate elements for the minimum weight basis. */
    {
        igraph_inclist_t inclist;
        igraph_vector_int_t visited; /* visited[v] % 3 is zero for unvisited vertices, see igraph_i_fundamental_cycles_bfs() */
        igraph_vector_int_t degrees;
        igraph_integer_t no_of_comps;
        igraph_integer_t mark;

        /* We use the degrees to avoid doing a BFS from vertices with d < 3, except in special cases.
         * Degrees cannot be computed from the inclist because there we use IGRAPH_LOOPS_ONCE. */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&degrees, no_of_nodes);
        IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS));

        IGRAPH_CHECK(igraph_connected_components(graph, NULL, NULL, &no_of_comps, IGRAPH_WEAK));
        rank = no_of_edges - no_of_nodes + no_of_comps;

        IGRAPH_VECTOR_INT_INIT_FINALLY(&visited, no_of_nodes);

        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_ONCE));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

         /* TODO: estimate space to reserve. 'rank' is a lower bound only. */
        IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&candidates, 0);
        IGRAPH_CHECK(igraph_vector_int_list_reserve(&candidates, rank));

        mark = 0;
        for (igraph_integer_t i=0; i < no_of_nodes; ++i) {
            igraph_integer_t degree = VECTOR(degrees)[i];
            igraph_bool_t vis = VECTOR(visited)[i] % 3 != 0; /* was vertex i visited already? */

            /* Generally, we only need to run a BFS from vertices of degree 3 or greater.
             * The exception is a connected component which is itself a cycle, and therefore
             * only contains vertices of degree 2. Thus from unvisited vertices we always run
             * a full BFS while from already visited ones only if their degree is at least 3. */

            /* TODO: mark entire component as visited, not just vertex. */
            if (degree <= 1 || (vis && degree < 3)) {
                continue;
            }

            /* TODO: BFS is only necessary from a feedback vertex set, find fast FVS approximation algorithm. */

            IGRAPH_CHECK(igraph_i_fundamental_cycles_bfs(
                    graph, &candidates, i, (vis || !complete) ? bfs_cutoff : -1, &inclist, &visited, mark));
            mark += 3;
        }

        igraph_inclist_destroy(&inclist);
        igraph_vector_int_destroy(&visited);
        igraph_vector_int_destroy(&degrees);
        IGRAPH_FINALLY_CLEAN(3);
    }

    /* Sort candidates by size (= weight) and remove duplicates. */
    {
        igraph_integer_t cand_count = igraph_vector_int_list_size(&candidates);

        for (igraph_integer_t i=0; i < cand_count; ++i) {
            igraph_vector_int_sort(igraph_vector_int_list_get_ptr(&candidates, i));
        }
        igraph_vector_int_list_sort(&candidates, &cycle_cmp);
        igraph_vector_int_list_remove_consecutive_duplicates(&candidates, igraph_vector_int_all_e);
    }

    igraph_vector_int_list_clear(result);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(result, rank));

    /* Find a complete basis, starting with smallest elements. */
    /* This is typically the slowest part of the algorithm. */
    {
        igraph_integer_t cand_len = igraph_vector_int_list_size(&candidates);
        igraph_vector_int_list_t reduced_matrix;
        igraph_bool_t independent;

        IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&reduced_matrix, 0);

        for (igraph_integer_t i=0; i < cand_len; ++i) {
            const igraph_vector_int_t *cycle = igraph_vector_int_list_get_ptr(&candidates, i);

            IGRAPH_ALLOW_INTERRUPTION();

            IGRAPH_CHECK(gaussian_elimination(&reduced_matrix, cycle, &independent));
            if (independent) {
                IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(result, cycle));
            }

            if (igraph_vector_int_list_size(&reduced_matrix) == rank) {
                /* We have a complete basis. */
                break;
            }
        }

        igraph_vector_int_list_destroy(&reduced_matrix);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_list_destroy(&candidates);
    IGRAPH_FINALLY_CLEAN(1);

    if (use_cycle_order) {
        igraph_integer_t result_size = igraph_vector_int_list_size(result);
        igraph_vector_int_t tmp;
        IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
        for (igraph_integer_t i=0; i < result_size; ++i) {
            igraph_vector_int_t *cycle = igraph_vector_int_list_get_ptr(result, i);
            IGRAPH_CHECK(igraph_vector_int_update(&tmp, cycle));
            IGRAPH_CHECK(igraph_i_order_cycle(graph, &tmp, cycle));
        }
        igraph_vector_int_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}
