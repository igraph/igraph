/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_transitivity.h"

#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_adjlist.h"

#include "core/interruption.h"

/* Computes the size of the intersection of two sorted vectors, treated as sets.
 * It is assumed that the vectors contain no duplicates.
 *
 * We rely on (lazy_)adjlist_get() producing sorted neighbor lists and
 * (lazy_)adjlist_init() being called with IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE
 * to prevent duplicate entries.
 */
static igraph_integer_t vector_int_intersection_size_sorted(
        const igraph_vector_int_t *v1, const igraph_vector_int_t *v2) {
    igraph_integer_t n1 = igraph_vector_int_size(v1), n2 = igraph_vector_int_size(v2);
    igraph_integer_t i1 = 0, i2 = 0;
    igraph_integer_t count = 0;

    while (i1 < n1 && i2 < n2) {
        igraph_integer_t e1 = VECTOR(*v1)[i1], e2 = VECTOR(*v2)[i2];
        if (e1 < e2) {
            i1++;
        } else if (e1 == e2) {
            count++;
            i1++; i2++;
        } else { /* e2 > e1 */
            i2++;
        }
    }

    return count;
}


/* Optimized for the case when computing ECC for all edges. */
static igraph_error_t igraph_i_ecc3_1(
        const igraph_t *graph, igraph_vector_t *res, const igraph_es_t eids,
        igraph_bool_t offset, igraph_bool_t normalize) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t degree;
    igraph_adjlist_t al;
    igraph_eit_t eit;
    const igraph_real_t c = offset ? 1.0 : 0.0;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS));

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (igraph_integer_t i=0;
         ! IGRAPH_EIT_END(eit);
         IGRAPH_EIT_NEXT(eit), i++) {

        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t v1 = IGRAPH_FROM(graph, edge), v2 = IGRAPH_TO(graph, edge);

        igraph_real_t z; /* number of triangles the edge participates in */
        igraph_real_t s; /* max number of triangles the edge could be part of */

        IGRAPH_ALLOW_INTERRUPTION();

        if (v1 == v2) {
            /* A self-loop isn't, and cannot be part of any triangles. */
            z = 0.0;
            s = 0.0;
        } else {
            const igraph_vector_int_t *a1 = igraph_adjlist_get(&al, v1), *a2 = igraph_adjlist_get(&al, v2);
            igraph_integer_t d1 = VECTOR(degree)[v1], d2 = VECTOR(degree)[v2];

            z = vector_int_intersection_size_sorted(a1, a2);
            s = (d1 < d2 ? d1 : d2) - 1.0;
        }

        VECTOR(*res)[i] = z + c;
        if (normalize) VECTOR(*res)[i] /= s;
    }

    igraph_eit_destroy(&eit);
    igraph_vector_int_destroy(&degree);
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Optimized for computing ECC for a small subset of edges. */
static igraph_error_t igraph_i_ecc3_2(
        const igraph_t *graph, igraph_vector_t *res,
        const igraph_es_t eids, igraph_bool_t offset, igraph_bool_t normalize) {

    igraph_lazy_adjlist_t al;
    igraph_eit_t eit;
    const igraph_real_t c = offset ? 1.0 : 0.0;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (igraph_integer_t i=0;
         ! IGRAPH_EIT_END(eit);
         IGRAPH_EIT_NEXT(eit), i++) {

        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t v1 = IGRAPH_FROM(graph, edge), v2 = IGRAPH_TO(graph, edge);

        igraph_real_t z; /* number of triangles the edge participates in */
        igraph_real_t s; /* max number of triangles the edge could be part of */

        IGRAPH_ALLOW_INTERRUPTION();

        if (v1 == v2) {
            /* A self-loop isn't, and cannot be part of any triangles. */
            z = 0.0;
            s = 0.0;
        } else {
            igraph_vector_int_t *a1 = igraph_lazy_adjlist_get(&al, v1);
            igraph_vector_int_t *a2 = igraph_lazy_adjlist_get(&al, v2);

            igraph_integer_t d1, d2;
            IGRAPH_CHECK(igraph_degree_1(graph, &d1, v1, IGRAPH_ALL, IGRAPH_LOOPS));
            IGRAPH_CHECK(igraph_degree_1(graph, &d2, v2, IGRAPH_ALL, IGRAPH_LOOPS));

            z = vector_int_intersection_size_sorted(a1, a2);
            s = (d1 < d2 ? d1 : d2) - 1.0;
        }

        VECTOR(*res)[i] = z + c;
        if (normalize) VECTOR(*res)[i] /= s;
    }

    igraph_eit_destroy(&eit);
    igraph_lazy_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/* Optimized for the case when computing ECC for all edges. */
static igraph_error_t igraph_i_ecc4_1(
        const igraph_t *graph, igraph_vector_t *res,
        const igraph_es_t eids, igraph_bool_t offset, igraph_bool_t normalize) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t degree;
    igraph_adjlist_t al;
    igraph_eit_t eit;
    igraph_real_t c = offset ? 1.0 : 0.0;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&degree, no_of_nodes);
    IGRAPH_CHECK(igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS));

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (igraph_integer_t i=0;
         ! IGRAPH_EIT_END(eit);
         IGRAPH_EIT_NEXT(eit), i++) {

        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t v1 = IGRAPH_FROM(graph, edge), v2 = IGRAPH_TO(graph, edge);

        igraph_real_t z; /* number of 4-cycles the edge participates in */
        igraph_real_t s; /* max number of 4-cycles the edge could be part of */

        IGRAPH_ALLOW_INTERRUPTION();

        if (v1 == v2) {
            z = 0.0;
            s = 0.0;
        } else {
            /* ensure that v1 is the vertex with the smaller degree */
            if (VECTOR(degree)[v1] > VECTOR(degree)[v2]) {
                igraph_integer_t tmp = v1;
                v1 = v2;
                v2 = tmp;
            }

            z = 0.0;
            const igraph_vector_int_t *a1 = igraph_adjlist_get(&al, v1);
            const igraph_integer_t n = igraph_vector_int_size(a1);
            for (igraph_integer_t j=0; j < n; j++) {
                igraph_integer_t v3 = VECTOR(*a1)[j];

                /* It is not possible that v3 == v1 because self-loops have been removed from the adjlist. */

                if (v3 == v2) continue;

                const igraph_vector_int_t *a2 = igraph_adjlist_get(&al, v2), *a3 = igraph_adjlist_get(&al, v3);

                z += vector_int_intersection_size_sorted(a2, a3) - 1.0;
            }

            igraph_integer_t d1 = VECTOR(degree)[v1], d2 = VECTOR(degree)[v2];
            s = (d1 - 1.0) * (d2 - 1.0);
        }

        VECTOR(*res)[i] = z + c;
        if (normalize) VECTOR(*res)[i] /= s;
    }

    igraph_eit_destroy(&eit);
    igraph_vector_int_destroy(&degree);
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}


/* Optimized for computing ECC for a small subset of edges. */
static igraph_error_t igraph_i_ecc4_2(
        const igraph_t *graph, igraph_vector_t *res,
        const igraph_es_t eids, igraph_bool_t offset, igraph_bool_t normalize) {

    igraph_lazy_adjlist_t al;
    igraph_eit_t eit;
    igraph_real_t c = offset ? 1.0 : 0.0;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

    IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_EIT_SIZE(eit)));

    for (igraph_integer_t i=0;
         ! IGRAPH_EIT_END(eit);
         IGRAPH_EIT_NEXT(eit), i++) {

        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t v1 = IGRAPH_FROM(graph, edge), v2 = IGRAPH_TO(graph, edge);

        igraph_real_t z; /* number of 4-cycles the edge participates in */
        igraph_real_t s; /* max number of 4-cycles the edge could be part of */

        IGRAPH_ALLOW_INTERRUPTION();

        igraph_integer_t d1, d2;
        IGRAPH_CHECK(igraph_degree_1(graph, &d1, v1, IGRAPH_ALL, IGRAPH_LOOPS));
        IGRAPH_CHECK(igraph_degree_1(graph, &d2, v2, IGRAPH_ALL, IGRAPH_LOOPS));

        if (v1 == v2) {
            z = 0.0;
            s = 0.0;
        } else {
            /* ensure that v1 is the vertex with the smaller degree */
            if (d1 > d2) {
                igraph_integer_t tmp = v1;
                v1 = v2;
                v2 = tmp;

                tmp = d1;
                d1 = d2;
                d2 = tmp;
            }

            z = 0.0;

            igraph_vector_int_t *a1 = igraph_lazy_adjlist_get(&al, v1);

            const igraph_integer_t n = igraph_vector_int_size(a1);
            for (igraph_integer_t j=0; j < n; j++) {
                igraph_integer_t v3 = VECTOR(*a1)[j];

                /* It is not possible that v3 == v1 because self-loops have been removed from the adjlist. */

                if (v3 == v2) continue;

                igraph_vector_int_t *a2 = igraph_lazy_adjlist_get(&al, v2);
                igraph_vector_int_t *a3 = igraph_lazy_adjlist_get(&al, v3);

                z += vector_int_intersection_size_sorted(a2, a3) - 1.0;
            }

            s = (d1 - 1.0) * (d2 - 1.0);
        }

        VECTOR(*res)[i] = z + c;
        if (normalize) VECTOR(*res)[i] /= s;
    }

    igraph_eit_destroy(&eit);
    igraph_lazy_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_ecc
 * \brief Edge clustering coefficient of some edges.
 *
 * \experimental
 *
 * The edge clustering coefficient <code>C^(k)_ij</code> of an edge (i, j)
 * is defined based on the number of k-cycles the edge participates in,
 * <code>z^(k)_ij</code>, and the largest number of such cycles it could
 * participate in given the degrees of its endpoints, <code>s^(k)_ij</code>.
 * The original definition given in the reference below is:
 *
 * </para><para>
 * <code>C^(k)_ij = (z^(k)_ij + 1) / s^(k)_ij</code>
 *
 * </para><para>
 * For <code>k=3</code>, <code>s^(k)_ij = min(d_i - 1, d_j - 1)</code>,
 * where \c d_i and \c d_j are the edge endpoint degrees.
 * For <code>k=4</code>, <code>s^(k)_ij = (d_i - 1) (d_j - 1)</code>.
 *
 * </para><para>
 * The \p normalize and \p offset parameters allow for skipping normalization
 * by <code>s^(k)</code> and offsetting the cycle count <code>z^(k)</code>
 * by one in the numerator of <code>C^(k)</code>. Set both to \c true to
 * compute the original definition of this metric.
 *
 * </para><para>
 * This function ignores edge multiplicities when listing k-cycles
 * (i.e. <code>z^(k)</code>), but not when computing the maximum number of
 * cycles an edge can participate in (<code>s^(k)</code>).
 *
 * </para><para>
 * Reference:
 *
 * </para><para>
 * F. Radicchi, C. Castellano, F. Cecconi, V. Loreto, and D. Parisi,
 * PNAS 101, 2658 (2004).
 * https://doi.org/10.1073/pnas.0400054101
 *
 * \param graph The input graph.
 * \param res Initialized vector, the result will be stored here.
 * \param eids The edges for which the edge clustering coefficient will be computed.
 * \param k Size of cycles to use in calculation. Must be at least 3. Currently
 *   only values of 3 and 4 are supported.
 * \param offset Boolean, whether to add one to cycle counts. When \c false,
 *   <code>z^(k)</code> is used instead of <code>z^(k) + 1</code>. In this case
 *   the maximum value of the normalized metric is 1. For <code>k=3</code> this
 *   is achieved for all edges in a complete graph.
 * \param normalize Boolean, whether to normalize cycle counts by the maximum
 *   possible count <code>s^(k)</code> given the degrees.
 * \return Error code.
 *
 * Time complexity: When \p k is 3, O(|V| d log d + |E| d).
 * When \p k is 4, O(|V| d log d + |E| d^2). d denotes the degree of vertices.
 */
igraph_error_t igraph_ecc(const igraph_t *graph, igraph_vector_t *res,
                          const igraph_es_t eids, igraph_integer_t k,
                          igraph_bool_t offset, igraph_bool_t normalize) {

    if (k < 3) {
        IGRAPH_ERRORF("Cycle size for edge clustering coefficient must be at least 3, got %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, k);
    }

    switch (k) {
    case 3:
        if (igraph_es_is_all(&eids)) {
            return igraph_i_ecc3_1(graph, res, eids, offset, normalize);
        } else {
            return igraph_i_ecc3_2(graph, res, eids, offset, normalize);
        }
    case 4:
        if (igraph_es_is_all(&eids)) {
            return igraph_i_ecc4_1(graph, res, eids, offset, normalize);
        } else {
            return igraph_i_ecc4_2(graph, res, eids, offset, normalize);
        }
    default:
        IGRAPH_ERROR("Edge clustering coefficient calculation is only implemented for cycle sizes 3 and 4.",
                     IGRAPH_UNIMPLEMENTED);
    }
}
