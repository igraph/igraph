/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#include "igraph_centrality.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_structural.h"

#include "core/interruption.h"

/**
 * \function igraph_cd_index
 * \brief The CD index, I-index and mCD-index of citation disruption.
 *
 * \experimental
 *
 * These measures characterize how a focal vertex in a timestamped, directed
 * citation network is used by subsequent work. An edge \c u -&gt; \c v is
 * interpreted as "\c u cites \c v", i.e. \c v is the earlier, cited work.
 *
 * </para><para>
 * For a focal vertex \c f with timestamp <code>t_f</code> and a time window
 * \c W, let \c C be the set of vertices timestamped in the interval
 * <code>(t_f, t_f + W]</code> that either cite \c f directly, or cite some
 * vertex that \c f itself cites. For each <code>c \in C</code>, let
 * <code>fbit(c) = 1</code> if \c c cites \c f directly, and let
 * <code>bbit(c) = 1</code> if \c c also cites at least one of \c f's own
 * references (i.e. \c c "backs up" through \c f rather than displacing it).
 * The CD index of \c f is then
 * <code>sum_{c \in C} (fbit(c) - 2 fbit(c) bbit(c)) / |C|</code>: each
 * citing vertex contributes +1 if it displaces \c f's references
 * (disruptive), -1 if it cites both \c f and its references
 * (consolidating), and 0 if it does not cite \c f directly. The CD index
 * is \c NaN for vertices with an empty \c C (no relevant future
 * citations).
 *
 * </para><para>
 * The I-index of \c f is simply the number of vertices citing \c f directly
 * that are timestamped at most <code>t_f + W</code> (note that unlike \c C,
 * this is not bounded below by <code>t_f</code>, matching the original
 * definition). The mCD-index is the product of the CD index and the
 * I-index, a magnitude-weighted variant of the CD index.
 *
 * </para><para>
 * This is a port of the algorithm implemented by the fast-cdindex library
 * (itself based on Funk &amp; Owen-Smith, 2017), reimplemented using igraph's
 * own data structures. When igraph is compiled with OpenMP support, the
 * computation is parallelized across the focal vertices in \p vids; the
 * number of threads used is controlled through the standard OpenMP
 * mechanisms (e.g. the \c OMP_NUM_THREADS environment variable).
 *
 * \param graph The input graph. It must be directed and must not contain
 *    self-loops.
 * \param timestamps A vector containing a timestamp for each vertex, in
 *    vertex ID order. Its length must equal the number of vertices.
 * \param res An initialized vector. The CD index of each vertex in \p vids
 *    will be stored here, in the order given by \p vids.
 * \param i_index An initialized vector, or \c NULL. If not \c NULL, the
 *    I-index of each vertex in \p vids will be stored here.
 * \param mcd_index An initialized vector, or \c NULL. If not \c NULL, the
 *    mCD-index of each vertex in \p vids will be stored here.
 * \param vids The vertices to compute the indices for.
 * \param time_window The width of the time window following each vertex's
 *    own timestamp within which citing vertices are taken into account.
 *    Must not be negative.
 * \return Error code.
 *
 * \sa \ref igraph_convergence_degree() for a related, non-temporal
 *    edge-based disruption-like measure.
 *
 * Time complexity: O(|V| + |E| + sum of local neighborhood products around
 * the vertices in \p vids).
 */
igraph_error_t igraph_cd_index(
        const igraph_t *graph,
        const igraph_vector_int_t *timestamps,
        igraph_vector_t *res,
        igraph_vector_t *i_index,
        igraph_vector_t *mcd_index,
        const igraph_vs_t vids,
        igraph_int_t time_window) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_adjlist_t out_adj, in_adj;
    igraph_vit_t vit;
    igraph_vector_int_t focal_ids;
    igraph_int_t no_of_focal;
    igraph_bool_t has_loops;
    igraph_bool_t oom_occurred = false;
    igraph_int_t idx;

    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("The CD index is only defined for directed graphs.", IGRAPH_EINVAL);
    }

    if (!timestamps) {
        IGRAPH_ERROR("Vertex timestamps must be given.", IGRAPH_EINVAL);
    }
    if (igraph_vector_int_size(timestamps) != no_of_nodes) {
        IGRAPH_ERROR("Invalid timestamp vector length.", IGRAPH_EINVAL);
    }

    if (time_window < 0) {
        IGRAPH_ERROR("Time window must not be negative.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_has_loop(graph, &has_loops));
    if (has_loops) {
        IGRAPH_ERROR("The CD index does not support graphs with self-loops.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &out_adj, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &out_adj);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &in_adj, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &in_adj);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    no_of_focal = IGRAPH_VIT_SIZE(vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&focal_ids, no_of_focal);
    for (idx = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), idx++) {
        VECTOR(focal_ids)[idx] = IGRAPH_VIT_GET(vit);
    }
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_vector_resize(res, no_of_focal));
    if (i_index) {
        IGRAPH_CHECK(igraph_vector_resize(i_index, no_of_focal));
    }
    if (mcd_index) {
        IGRAPH_CHECK(igraph_vector_resize(mcd_index, no_of_focal));
    }

    IGRAPH_ALLOW_INTERRUPTION();

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        /* marker[v] == 0 means "v not yet seen for the current focal vertex".
         * A nonzero value packs a per-focal tag in the upper bits with two
         * set-membership flags in the low bits (bit 0: v is a citer of focal
         * itself, i.e. Set A below; bit 1: v is a citer of one of focal's
         * references, i.e. Set B below). Calloc's zero-fill already means
         * "untouched", so -- unlike a sentinel of -1 -- no per-thread
         * initialization pass over the array is needed at all. */
        igraph_uint_t *marker = IGRAPH_CALLOC(no_of_nodes, igraph_uint_t);
        if (marker == NULL) {
            #ifdef _OPENMP
            #pragma omp atomic write
            #endif
            oom_occurred = true;
        } else {
            #ifdef _OPENMP
            #pragma omp for
            #endif
            for (idx = 0; idx < no_of_focal; idx++) {
                igraph_bool_t skip;
                #ifdef _OPENMP
                #pragma omp atomic read
                #endif
                skip = oom_occurred;
                if (skip) {
                    continue;
                }

                igraph_int_t focal = VECTOR(focal_ids)[idx];
                igraph_int_t t_focal = VECTOR(*timestamps)[focal];
                igraph_int_t window_end = t_focal + time_window;
                const igraph_vector_int_t *focal_out = igraph_adjlist_get(&out_adj, focal);
                const igraph_vector_int_t *focal_in = igraph_adjlist_get(&in_adj, focal);
                igraph_int_t n_focal_out = igraph_vector_int_size(focal_out);
                igraph_int_t n_focal_in = igraph_vector_int_size(focal_in);
                igraph_uint_t tag = ((igraph_uint_t) idx + 1) << 2;
                igraph_int_t count_a = 0, count_b = 0, n_distinct = 0;
                igraph_real_t sum;
                igraph_int_t i_count = 0;
                igraph_int_t r, j;

                /* Sixt & Pasin's (2024) reformulation of the CD index: rather than
                 * building one candidate set and, for each member, cross-checking
                 * whether it cites both focal and a reference (the fbit/bbit
                 * approach), decompose it into two independent set-membership
                 * counts that never need to look at each other:
                 *   Set A = vertices citing focal directly;
                 *   Set B = vertices citing any of focal's own references.
                 * CD = (-|A| - 2|B|) / |A union B| + 2, which is algebraically
                 * identical to the classic per-candidate sum/|C| formulation
                 * (s'(c) = -1 for c in A, s''(c) = -2 for c in B, and summing
                 * s'+s'' over the union telescopes to -|A|-2|B| regardless of
                 * overlap, since s' is 0 outside A and s'' is 0 outside B).
                 * This removes the need to ever look up a candidate's own
                 * out-edges, which is where the earlier fbit/bbit version spent
                 * most of its time. */

                /* Set A: focal's own citers. */
                for (j = 0; j < n_focal_in; j++) {
                    igraph_int_t cand = VECTOR(*focal_in)[j];
                    igraph_int_t t_cand = VECTOR(*timestamps)[cand];

                    /* I-index: bounded above only, by construction (matches
                     * both the original cdindex and fast-cdindex, neither of
                     * which excludes in-edges timestamped at or before t_focal). */
                    if (t_cand <= window_end) {
                        i_count++;
                    }

                    /* Set A additionally requires t_cand to be strictly later
                     * than focal's own timestamp. */
                    if (t_cand > t_focal && t_cand <= window_end) {
                        igraph_uint_t val = marker[cand];
                        if ((val & ~(igraph_uint_t) 3) != tag) {
                            val = tag;
                        }
                        if (!(val & 1)) {
                            if ((val & 3) == 0) {
                                n_distinct++;
                            }
                            marker[cand] = val | 1;
                            count_a++;
                        }
                    }
                }

                /* Set B: citers of any of focal's own references. */
                for (r = 0; r < n_focal_out; r++) {
                    igraph_int_t ref = VECTOR(*focal_out)[r];
                    const igraph_vector_int_t *ref_in = igraph_adjlist_get(&in_adj, ref);
                    igraph_int_t n_ref_in = igraph_vector_int_size(ref_in);
                    for (j = 0; j < n_ref_in; j++) {
                        igraph_int_t cand = VECTOR(*ref_in)[j];
                        igraph_int_t t_cand = VECTOR(*timestamps)[cand];
                        if (t_cand > t_focal && t_cand <= window_end) {
                            igraph_uint_t val = marker[cand];
                            if ((val & ~(igraph_uint_t) 3) != tag) {
                                val = tag;
                            }
                            if (!(val & 2)) {
                                if ((val & 3) == 0) {
                                    n_distinct++;
                                }
                                marker[cand] = val | 2;
                                count_b++;
                            }
                        }
                    }
                }

                sum = -(igraph_real_t) count_a - 2.0 * count_b;
                VECTOR(*res)[idx] = n_distinct > 0 ? sum / n_distinct + 2.0 : IGRAPH_NAN;
                if (i_index) {
                    VECTOR(*i_index)[idx] = i_count;
                }
                if (mcd_index) {
                    VECTOR(*mcd_index)[idx] = VECTOR(*res)[idx] * i_count;
                }
            }

            IGRAPH_FREE(marker);
        }
    }

    igraph_vector_int_destroy(&focal_ids);
    igraph_adjlist_destroy(&in_adj);
    igraph_adjlist_destroy(&out_adj);
    IGRAPH_FINALLY_CLEAN(3);

    if (oom_occurred) {
        IGRAPH_ERROR("Cannot calculate CD index.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }

    return IGRAPH_SUCCESS;
}
