/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2020  The igraph development team <igraph@igraph.org>

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

#include "igraph_graphicality.h"

#define IGRAPH_I_MULTI_EDGES_SW 0x02 /* 010, more than one edge allowed between distinct vertices */
#define IGRAPH_I_MULTI_LOOPS_SW 0x04 /* 100, more than one self-loop allowed on the same vertex   */

static igraph_error_t igraph_i_is_graphical_undirected_multi_loops(const igraph_vector_int_t *degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_undirected_loopless_multi(const igraph_vector_int_t *degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_undirected_loopy_simple(const igraph_vector_int_t *degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_undirected_simple(const igraph_vector_int_t *degrees, igraph_bool_t *res);

static igraph_error_t igraph_i_is_graphical_directed_loopy_multi(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_directed_loopless_multi(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_directed_loopy_simple(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res);
static igraph_error_t igraph_i_is_graphical_directed_simple(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res);

static igraph_error_t igraph_i_is_bigraphical_multi(const igraph_vector_int_t *degrees1, const igraph_vector_int_t *degrees2, igraph_bool_t *res);
static igraph_error_t igraph_i_is_bigraphical_simple(const igraph_vector_int_t *degrees1, const igraph_vector_int_t *degrees2, igraph_bool_t *res);


/**
 * \function igraph_is_graphical
 * \brief Is there a graph with the given degree sequence?
 *
 * Determines whether a sequence of integers can be the degree sequence of some graph.
 * The classical concept of graphicality assumes simple graphs. This function can perform
 * the check also when either self-loops, multi-edge, or both are allowed in the graph.
 *
 * </para><para>
 * For simple undirected graphs, the Erdős-Gallai conditions are checked using the linear-time
 * algorithm of Cloteaux. If both self-loops and multi-edges are allowed,
 * it is sufficient to chek that that sum of degrees is even. If only multi-edges are allowed, but
 * not self-loops, there is an additional condition that the sum of degrees be no smaller than twice
 * the maximum degree. If at most one self-loop is allowed per vertex, but no multi-edges, a modified
 * version of the Erdős-Gallai conditions are used (see Cairns &amp; Mendan).
 *
 * </para><para>
 * For simple directed graphs, the Fulkerson-Chen-Anstee theorem is used with the relaxation by Berger.
 * If both self-loops and multi-edges are allowed, then it is sufficient to check that the sum of
 * in- and out-degrees is the same. If only multi-edges are allowed, but not self loops, there is an
 * additional condition that the sum of out-degrees (or equivalently, in-degrees) is no smaller than
 * the maximum total degree. If single self-loops are allowed, but not multi-edges, the problem is equivalent
 * to realizability as a simple bipartite graph, thus the Gale-Ryser theorem can be used; see
 * \ref igraph_is_bigraphical() for more information.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * P. Erdős and T. Gallai, Gráfok előírt fokú pontokkal, Matematikai Lapok 11, pp. 264–274 (1960).
 * https://users.renyi.hu/~p_erdos/1961-05.pdf
 *
 * </para><para>
 * Z Király, Recognizing graphic degree sequences and generating all realizations.
 * TR-2011-11, Egerváry Research Group, H-1117, Budapest, Hungary. ISSN 1587-4451 (2012).
 * http://bolyai.cs.elte.hu/egres/tr/egres-11-11.pdf
 *
 * </para><para>
 * B. Cloteaux, Is This for Real? Fast Graphicality Testing, Comput. Sci. Eng. 17, 91 (2015).
 * https://dx.doi.org/10.1109/MCSE.2015.125
 *
 * </para><para>
 * A. Berger, A note on the characterization of digraphic sequences, Discrete Math. 314, 38 (2014).
 * https://dx.doi.org/10.1016/j.disc.2013.09.010
 *
 * </para><para>
 * G. Cairns and S. Mendan, Degree Sequence for Graphs with Loops (2013).
 * https://arxiv.org/abs/1303.2145v1
 *
 * \param out_degrees A vector of integers specifying the degree sequence for
 *     undirected graphs or the out-degree sequence for directed graphs.
 * \param in_degrees A vector of integers specifying the in-degree sequence for
 *     directed graphs. For undirected graphs, it must be \c NULL.
 * \param allowed_edge_types The types of edges to allow in the graph:
 *     \clist
 *     \cli IGRAPH_SIMPLE_SW
 *       simple graphs (i.e. no self-loops or multi-edges allowed).
 *     \cli IGRAPH_LOOPS_SW
 *       single self-loops are allowed, but not multi-edges.
 *     \cli IGRAPH_MULTI_SW
 *       multi-edges are allowed, but not self-loops.
 *     \cli IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW
 *       both self-loops and multi-edges are allowed.
 *     \endclist
 * \param res Pointer to a Boolean. The result will be stored here.
 *
 * \return Error code.
 *
 * \sa \ref igraph_is_bigraphical() to check if a bi-degree-sequence can be realized as a bipartite graph;
 * \ref igraph_realize_degree_sequence() to construct a graph with a given degree sequence.
 *
 * Time complexity: O(n log n) for directed graphs with at most one self-loop per vertex,
 * and O(n) for all other cases, where n is the length of the degree sequence(s).
 */
igraph_error_t igraph_is_graphical(const igraph_vector_int_t *out_degrees,
                        const igraph_vector_int_t *in_degrees,
                        const igraph_edge_type_sw_t allowed_edge_types,
                        igraph_bool_t *res)
{
    /* Undirected case: */
    if (in_degrees == NULL)
    {
        if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW )) {
            /* Typically this case is used when multiple edges are allowed both as self-loops and
             * between distinct vertices. However, the conditions are the same even if multi-edges
             * are not allowed between distinct vertices (only as self-loops). Therefore, we
             * do not test IGRAPH_I_MULTI_EDGES_SW in the if (...). */
            return igraph_i_is_graphical_undirected_multi_loops(out_degrees, res);
        }
        else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_undirected_loopless_multi(out_degrees, res);
        }
        else if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_undirected_loopy_simple(out_degrees, res);
        }
        else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_undirected_simple(out_degrees, res);
        } else {
            /* Remaining case:
             *  - At most one self-loop per vertex but multi-edges between distinct vertices allowed.
             * These cases cannot currently be requested through the documented API,
             * so no explanatory error message for now. */
            return IGRAPH_UNIMPLEMENTED;
        }
    }
    /* Directed case: */
    else
    {
        if (igraph_vector_int_size(in_degrees) != igraph_vector_int_size(out_degrees)) {
            IGRAPH_ERROR("The length of out- and in-degree sequences must be the same.", IGRAPH_EINVAL);
        }

        if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) && (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW ) ) {
            return igraph_i_is_graphical_directed_loopy_multi(out_degrees, in_degrees, res);
        }
        else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_directed_loopless_multi(out_degrees, in_degrees, res);
        }
        else if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_directed_loopy_simple(out_degrees, in_degrees, res);
        }
        else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) ) {
            return igraph_i_is_graphical_directed_simple(out_degrees, in_degrees, res);
        } else {
            /* Remaining cases:
             *  - At most one self-loop per vertex but multi-edges between distinct vertices allowed.
             *  - At most one edge between distinct vertices but multi-self-loops allowed.
             * These cases cannot currently be requested through the documented API,
             * so no explanatory error message for now. */
            return IGRAPH_UNIMPLEMENTED;
        }
    }

    /* can't reach here */
}

/**
 * \function igraph_is_bigraphical
 * \brief Is there a bipartite graph with the given bi-degree-sequence?
 *
 * Determines whether two sequences of integers can be the degree sequences of
 * a bipartite graph. Such a pair of degree sequence is called \em bigraphical.
 *
 * </para><para>
 * When multi-edges are allowed, it is sufficient to check that the sum of degrees is the
 * same in the two partitions. For simple graphs, the Gale-Ryser theorem is used
 * with Berger's relaxation.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * H. J. Ryser, Combinatorial Properties of Matrices of Zeros and Ones, Can. J. Math. 9, 371 (1957).
 * https://dx.doi.org/10.4153/cjm-1957-044-3
 *
 * </para><para>
 * D. Gale, A theorem on flows in networks, Pacific J. Math. 7, 1073 (1957).
 * https://dx.doi.org/10.2140/pjm.1957.7.1073
 *
 * </para><para>
 * A. Berger, A note on the characterization of digraphic sequences, Discrete Math. 314, 38 (2014).
 * https://dx.doi.org/10.1016/j.disc.2013.09.010
 *
 * \param degrees1 A vector of integers specifying the degrees in the first partition
 * \param degrees2 A vector of integers specifying the degrees in the second partition
 * \param allowed_edge_types The types of edges to allow in the graph:
 *     \clist
 *     \cli IGRAPH_SIMPLE_SW
 *       simple graphs (i.e. no multi-edges allowed).
 *     \cli IGRAPH_MULTI_SW
 *       multi-edges are allowed.
 *     \endclist
 * \param res Pointer to a Boolean. The result will be stored here.
 *
 * \return Error code.
 *
 * \sa \ref igraph_is_graphical()
 *
 * Time complexity: O(n log n) for simple graphs, O(n) for multigraphs,
 * where n is the length of the larger degree sequence.
 */
igraph_error_t igraph_is_bigraphical(const igraph_vector_int_t *degrees1,
                          const igraph_vector_int_t *degrees2,
                          const igraph_edge_type_sw_t allowed_edge_types,
                          igraph_bool_t *res)
{
    /* Note: Bipartite graphs can't have self-loops so we ignore the IGRAPH_LOOPS_SW bit. */
    if (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) {
        return igraph_i_is_bigraphical_multi(degrees1, degrees2, res);
    } else {
        return igraph_i_is_bigraphical_simple(degrees1, degrees2, res);
    }
}


/***** Undirected case *****/

/* Undirected graph with multi-self-loops:
 *  - Degrees must be non-negative.
 *  - The sum of degrees must be even.
 *
 * These conditions are valid regardless of whether multi-edges are allowed between distinct vertices.
 */
static igraph_error_t igraph_i_is_graphical_undirected_multi_loops(const igraph_vector_int_t *degrees, igraph_bool_t *res) {
    igraph_integer_t sum_parity = 0; /* 0 if the degree sum is even, 1 if it is odd */
    igraph_integer_t n = igraph_vector_int_size(degrees);
    igraph_integer_t i;

    for (i = 0; i < n; ++i) {
        igraph_integer_t d = VECTOR(*degrees)[i];

        if (d < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }
        sum_parity = (sum_parity + d) & 1;
    }

    *res = (sum_parity == 0);

    return IGRAPH_SUCCESS;
}


/* Undirected loopless multigraph:
 *  - Degrees must be non-negative.
 *  - The sum of degrees must be even.
 *  - The sum of degrees must be no smaller than 2*d_max.
 */
static igraph_error_t igraph_i_is_graphical_undirected_loopless_multi(const igraph_vector_int_t *degrees, igraph_bool_t *res) {
    igraph_integer_t i;
    igraph_integer_t n = igraph_vector_int_size(degrees);
    igraph_integer_t dsum, dmax;

    /* Zero-length sequences are considered graphical. */
    if (n == 0) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    dsum = 0; dmax = 0;
    for (i = 0; i < n; ++i) {
        igraph_integer_t d = VECTOR(*degrees)[i];

        if (d < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }
        dsum += d;
        if (d > dmax) {
            dmax = d;
        }
    }

    *res = (dsum % 2 == 0) && (dsum >= 2*dmax);

    return IGRAPH_SUCCESS;
}


/* Undirected graph with no multi-edges and at most one self-loop per vertex:
 *  - Degrees must be non-negative.
 *  - The sum of degrees must be even.
 *  - Use the modification of the Erdős-Gallai theorem due to Cairns and Mendan.
 */
static igraph_error_t igraph_i_is_graphical_undirected_loopy_simple(const igraph_vector_int_t *degrees, igraph_bool_t *res) {
    igraph_vector_int_t work;
    igraph_integer_t w, b, s, c, n, k;

    n = igraph_vector_int_size(degrees);

    /* Zero-length sequences are considered graphical. */
    if (n == 0) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    /* The conditions from the loopy multigraph case are necessary here as well. */
    IGRAPH_CHECK(igraph_i_is_graphical_undirected_multi_loops(degrees, res));
    if (! *res) {
        return IGRAPH_SUCCESS;
    }

    /*
     * We follow this paper:
     *
     * G. Cairns & S. Mendan: Degree Sequences for Graphs with Loops, 2013
     * https://arxiv.org/abs/1303.2145v1
     *
     * They give the following modification of the Erdős-Gallai theorem:
     *
     * A non-increasing degree sequence d_1 >= ... >= d_n has a realization as
     * a simple graph with loops (i.e. at most one self-loop allowed on each vertex)
     * iff
     *
     * \sum_{i=1}^k d_i <= k(k+1) + \sum_{i=k+1}^{n} min(d_i, k)
     *
     * for each k=1..n
     *
     * The difference from Erdős-Gallai is that here we have the term
     * k(k+1) instead of k(k-1).
     *
     * The implementation is analogous to igraph_i_is_graphical_undirected_simple(),
     * which in turn is based on Király 2012. See comments in that function for details.
     * w and k are zero-based here, unlike in the statement of the theorem above.
     */

    IGRAPH_CHECK(igraph_vector_int_init_copy(&work, degrees));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &work);

    igraph_vector_int_reverse_sort(&work);

    *res = true;
    w = n - 1; b = 0; s = 0; c = 0;
    for (k = 0; k < n; k++) {
        b += VECTOR(work)[k];
        c += w;
        while (w > k && VECTOR(work)[w] <= k + 1) {
            s += VECTOR(work)[w];
            c -= (k + 1);
            w--;
        }
        if (b > c + s + 2*(k + 1)) {
            *res = false;
            break;
        }
        if (w == k) {
            break;
        }
    }

    igraph_vector_int_destroy(&work);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Undirected simple graph:
 *  - Degrees must be non-negative.
 *  - The sum of degrees must be even.
 *  - Use the Erdős-Gallai theorem.
 */
static igraph_error_t igraph_i_is_graphical_undirected_simple(const igraph_vector_int_t *degrees, igraph_bool_t *res) {
    igraph_vector_int_t num_degs; /* num_degs[d] is the # of vertices with degree d */
    const igraph_integer_t p = igraph_vector_int_size(degrees);
    igraph_integer_t dmin, dmax, dsum;
    igraph_integer_t n; /* number of non-zero degrees */
    igraph_integer_t k, sum_deg, sum_ni, sum_ini;
    igraph_integer_t i, dk;
    igraph_integer_t zverovich_bound;

    if (p == 0) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    /* The following implementation of the Erdős-Gallai test
     * is mostly a direct translation of the Python code given in
     *
     * Brian Cloteaux, Is This for Real? Fast Graphicality Testing,
     * Computing Prescriptions, pp. 91-95, vol. 17 (2015)
     * https://dx.doi.org/10.1109/MCSE.2015.125
     *
     * It uses counting sort to achieve linear runtime.
     */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&num_degs, p);

    dmin = p; dmax = 0; dsum = 0; n = 0;
    for (i = 0; i < p; ++i) {
        igraph_integer_t d = VECTOR(*degrees)[i];

        if (d < 0 || d >= p) {
            *res = false;
            goto finish;
        }

        if (d > 0) {
            dmax = d > dmax ? d : dmax;
            dmin = d < dmin ? d : dmin;
            dsum += d;
            n++;
            VECTOR(num_degs)[d] += 1;
        }
    }

    if (dsum % 2 != 0) {
        *res = false;
        goto finish;
    }

    if (n == 0) {
        *res = true;
        goto finish; /* all degrees are zero => graphical */
    }

    /* According to:
     *
     * G. Cairns, S. Mendan, and Y. Nikolayevsky, A sharp refinement of a result of Zverovich-Zverovich,
     * Discrete Math. 338, 1085 (2015).
     * https://dx.doi.org/10.1016/j.disc.2015.02.001
     *
     * a sufficient but not necessary condition of graphicality for a sequence of
     * n strictly positive integers is that
     *
     * dmin * n >= floor( (dmax + dmin + 1)^2 / 4 ) - 1
     * if dmin is odd or (dmax + dmin) mod 4 == 1
     *
     * or
     *
     * dmin * n >= floor( (dmax + dmin + 1)^2 / 4 )
     * otherwise.
     */

    zverovich_bound = ((dmax + dmin + 1) * (dmax + dmin + 1)) / 4;
    if (dmin % 2 == 1 || (dmax + dmin) % 4 == 1) {
        zverovich_bound -= 1;
    }

    if (dmin*n >= zverovich_bound) {
        *res = true;
        goto finish;
    }

    k = 0; sum_deg = 0; sum_ni = 0; sum_ini = 0;
    for (dk = dmax; dk >= dmin; --dk) {
        igraph_integer_t run_size, v;

        if (dk < k+1) {
            *res = true;
            goto finish;
        }

        run_size = VECTOR(num_degs)[dk];
        if (run_size > 0) {
            if (dk < k + run_size) {
                run_size = dk - k;
            }
            sum_deg += run_size * dk;
            for (v=0; v < run_size; ++v) {
                sum_ni += VECTOR(num_degs)[k+v];
                sum_ini += (k+v) * VECTOR(num_degs)[k+v];
            }
            k += run_size;
            if (sum_deg > k*(n-1) - k*sum_ni + sum_ini) {
                *res = false;
                goto finish;
            }
        }
    }

    *res = true;

finish:
    igraph_vector_int_destroy(&num_degs);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/***** Directed case *****/

/* Directed loopy multigraph:
 *  - Degrees must be non-negative.
 *  - The sum of in- and out-degrees must be the same.
 */
static igraph_error_t igraph_i_is_graphical_directed_loopy_multi(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res) {
    igraph_integer_t sumdiff; /* difference between sum of in- and out-degrees */
    igraph_integer_t n = igraph_vector_int_size(out_degrees);
    igraph_integer_t i;

    IGRAPH_ASSERT(igraph_vector_int_size(in_degrees) == n);

    sumdiff = 0;
    for (i = 0; i < n; ++i) {
        igraph_integer_t dout = VECTOR(*out_degrees)[i];
        igraph_integer_t din  = VECTOR(*in_degrees)[i];

        if (dout < 0 || din < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }

        sumdiff += din - dout;
    }

    *res = sumdiff == 0;

    return IGRAPH_SUCCESS;
}


/* Directed loopless multigraph:
 *  - Degrees must be non-negative.
 *  - The sum of in- and out-degrees must be the same.
 *  - The sum of out-degrees must be no smaller than d_max,
 *    where d_max is the largest total degree.
 */
static igraph_error_t igraph_i_is_graphical_directed_loopless_multi(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res) {
    igraph_integer_t i, sumin, sumout, dmax;
    igraph_integer_t n = igraph_vector_int_size(out_degrees);

    IGRAPH_ASSERT(igraph_vector_int_size(in_degrees) == n);

    sumin = 0; sumout = 0;
    dmax = 0;
    for (i = 0; i < n; ++i) {
        igraph_integer_t dout = VECTOR(*out_degrees)[i];
        igraph_integer_t din  = VECTOR(*in_degrees)[i];
        igraph_integer_t d = dout + din;

        if (dout < 0 || din < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }

        sumin += din; sumout += dout;

        if (d > dmax) {
            dmax = d;
        }
    }

    *res = (sumin == sumout) && (sumout >= dmax);

    return IGRAPH_SUCCESS;
}


/* Directed graph with no multi-edges and at most one self-loop per vertex:
 *  - Degrees must be non-negative.
 *  - Equivalent to bipartite simple graph.
 */
static igraph_error_t igraph_i_is_graphical_directed_loopy_simple(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res) {
    return igraph_i_is_bigraphical_simple(out_degrees, in_degrees, res);
}


/* Directed simple graph:
 *  - Degrees must be non-negative.
 *  - The sum of in- and out-degrees must be the same.
 *  - Use the Fulkerson-Chen-Anstee theorem
 */
static igraph_error_t igraph_i_is_graphical_directed_simple(const igraph_vector_int_t *out_degrees, const igraph_vector_int_t *in_degrees, igraph_bool_t *res) {
    igraph_vector_int_t in_degree_cumcounts, in_degree_counts;
    igraph_vector_int_t sorted_in_degrees, sorted_out_degrees;
    igraph_vector_int_t left_pq, right_pq;
    igraph_integer_t lhs, rhs, left_pq_size, right_pq_size, left_i, right_i, left_sum, right_sum;

    /* The conditions from the loopy multigraph case are necessary here as well. */
    IGRAPH_CHECK(igraph_i_is_graphical_directed_loopy_multi(out_degrees, in_degrees, res));
    if (! *res) {
        return IGRAPH_SUCCESS;
    }

    const igraph_integer_t vcount = igraph_vector_int_size(out_degrees);
    if (vcount == 0) {
        *res = true;
        return IGRAPH_SUCCESS;
    }


    IGRAPH_VECTOR_INT_INIT_FINALLY(&in_degree_cumcounts, vcount+1);

    /* Compute in_degree_cumcounts[d+1] to be the no. of in-degrees == d */
    for (igraph_integer_t v = 0; v < vcount; v++) {
        igraph_integer_t indeg = VECTOR(*in_degrees)[v];
        igraph_integer_t outdeg = VECTOR(*out_degrees)[v];
        if (indeg >= vcount || outdeg >= vcount) {
            *res = false;
            igraph_vector_int_destroy(&in_degree_cumcounts);
            IGRAPH_FINALLY_CLEAN(1);
            return IGRAPH_SUCCESS;
        }
        VECTOR(in_degree_cumcounts)[indeg + 1]++;
    }

    /* Compute in_degree_cumcounts[d] to be the no. of in-degrees < d */
    for (igraph_integer_t indeg = 0; indeg < vcount; indeg++) {
        VECTOR(in_degree_cumcounts)[indeg+1] += VECTOR(in_degree_cumcounts)[indeg];
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&sorted_out_degrees, vcount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&sorted_in_degrees, vcount);

    /* In the following loop, in_degree_counts[d] keeps track of the number of vertices
     * with in-degree d that were already placed. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&in_degree_counts, vcount);

    for (igraph_integer_t v = 0; v < vcount; v++) {
        igraph_integer_t outdeg = VECTOR(*out_degrees)[v];
        igraph_integer_t indeg  = VECTOR(*in_degrees)[v];
        igraph_integer_t idx = VECTOR(in_degree_cumcounts)[indeg] + VECTOR(in_degree_counts)[indeg];
        VECTOR(sorted_out_degrees)[vcount - idx - 1] = outdeg;
        VECTOR(sorted_in_degrees)[vcount - idx - 1] = indeg;
        VECTOR(in_degree_counts)[indeg]++;
    }

    igraph_vector_int_destroy(&in_degree_counts);
    igraph_vector_int_destroy(&in_degree_cumcounts);
    IGRAPH_FINALLY_CLEAN(2);

    /* Be optimistic, then check whether the Fulkerson–Chen–Anstee condition
     * holds for every k. In particular, for every k in [0; n), it must be true
     * that:
     *
     * \sum_{i=0}^k indegree[i] <=
     *     \sum_{i=0}^k min(outdegree[i], k) +
     *     \sum_{i=k+1}^{n-1} min(outdegree[i], k + 1)
     */

#define INDEGREE(x) (VECTOR(sorted_in_degrees)[x])
#define OUTDEGREE(x) (VECTOR(sorted_out_degrees)[x])

    IGRAPH_VECTOR_INT_INIT_FINALLY(&left_pq, vcount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&right_pq, vcount);

    left_pq_size = 0;
    right_pq_size = vcount;
    left_i = 0;
    right_i = 0;
    left_sum = 0;
    right_sum = 0;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        VECTOR(right_pq)[OUTDEGREE(i)]++;
    }

    *res = true;
    lhs = 0;
    rhs = 0;
    for (igraph_integer_t i = 0; i < vcount; i++) {
        lhs += INDEGREE(i);

        /* It is enough to check for indexes where the in-degree is about to
         * decrease in the next step; see "Stronger condition" in the Wikipedia
         * entry for the Fulkerson-Chen-Anstee condition. However, this does not
         * provide any noticeable benefits for the current implementation. */

        if (OUTDEGREE(i) < i) {
            left_sum += OUTDEGREE(i);
        }
        else {
            VECTOR(left_pq)[OUTDEGREE(i)]++;
            left_pq_size++;
        }
        while (left_i < i) {
            while (VECTOR(left_pq)[left_i] > 0) {
                VECTOR(left_pq)[left_i]--;
                left_pq_size--;
                left_sum += left_i;
            }
            left_i++;
        }

        while (right_i < i + 1) {
            while (VECTOR(right_pq)[right_i] > 0) {
                VECTOR(right_pq)[right_i]--;
                right_pq_size--;
                right_sum += right_i;
            }
            right_i++;
        }
        if (OUTDEGREE(i) < i + 1) {
            right_sum -= OUTDEGREE(i);
        }
        else {
            VECTOR(right_pq)[OUTDEGREE(i)]--;
            right_pq_size--;
        }

        rhs = left_sum + i * left_pq_size + right_sum + (i + 1) * right_pq_size;
        if (lhs > rhs) {
            *res = false;
            break;
        }
    }

#undef INDEGREE
#undef OUTDEGREE

    igraph_vector_int_destroy(&sorted_in_degrees);
    igraph_vector_int_destroy(&sorted_out_degrees);
    igraph_vector_int_destroy(&left_pq);
    igraph_vector_int_destroy(&right_pq);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}



/***** Bipartite case *****/

/* Bipartite graph with multi-edges:
 *  - Degrees must be non-negative.
 *  - Sum of degrees must be the same in the two partitions.
 */
static igraph_error_t igraph_i_is_bigraphical_multi(const igraph_vector_int_t *degrees1, const igraph_vector_int_t *degrees2, igraph_bool_t *res) {
    igraph_integer_t i;
    igraph_integer_t sum1, sum2;
    igraph_integer_t n1 = igraph_vector_int_size(degrees1), n2 = igraph_vector_int_size(degrees2);

    sum1 = 0;
    for (i = 0; i < n1; ++i) {
        igraph_integer_t d = VECTOR(*degrees1)[i];

        if (d < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }

        sum1 += d;
    }

    sum2 = 0;
    for (i = 0; i < n2; ++i) {
        igraph_integer_t d = VECTOR(*degrees2)[i];

        if (d < 0) {
            *res = false;
            return IGRAPH_SUCCESS;
        }

        sum2 += d;
    }

    *res = (sum1 == sum2);

    return IGRAPH_SUCCESS;
}


/* Bipartite simple graph:
 *  - Degrees must be non-negative.
 *  - Sum of degrees must be the same in the two partitions.
 *  - Use the Gale-Ryser theorem.
 */
static igraph_error_t igraph_i_is_bigraphical_simple(const igraph_vector_int_t *degrees1, const igraph_vector_int_t *degrees2, igraph_bool_t *res) {
    igraph_vector_int_t sorted_deg1, sorted_deg2;
    igraph_integer_t n1 = igraph_vector_int_size(degrees1), n2 = igraph_vector_int_size(degrees2);
    igraph_integer_t i, k;
    igraph_integer_t lhs_sum, partial_rhs_sum;

    if (n1 == 0 && n2 == 0) {
        *res = true;
        return IGRAPH_SUCCESS;
    }

    /* The conditions from the multigraph case are necessary here as well. */
    IGRAPH_CHECK(igraph_i_is_bigraphical_multi(degrees1, degrees2, res));
    if (! *res) {
        return IGRAPH_SUCCESS;
    }

    /* Ensure that degrees1 is the shorter vector as a minor optimization: */
    if (n2 < n1) {
        const igraph_vector_int_t *tmp;
        igraph_integer_t n;

        tmp = degrees1;
        degrees1 = degrees2;
        degrees2 = tmp;

        n = n1;
        n1 = n2;
        n2 = n;
    }

    /* Copy and sort both vectors: */

    IGRAPH_CHECK(igraph_vector_int_init_copy(&sorted_deg1, degrees1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &sorted_deg1);
    igraph_vector_int_reverse_sort(&sorted_deg1); /* decreasing sort */

    IGRAPH_CHECK(igraph_vector_int_init_copy(&sorted_deg2, degrees2));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &sorted_deg2);
    igraph_vector_int_sort(&sorted_deg2); /* increasing sort */

    /*
     * We follow the description of the Gale-Ryser theorem in:
     *
     * A. Berger, A note on the characterization of digraphic sequences, Discrete Math. 314, 38 (2014).
     * https://doi.org/10.1016/j.disc.2013.09.010
     *
     * Gale-Ryser condition with 0-based indexing:
     *
     * a_i and b_i denote the degree sequences of the two partitions.
     *
     * Assuming that a_0 >= a_1 >= ... >= a_{n_1 - 1},
     *
     * \sum_{i=0}^k a_i <= \sum_{j=0}^{n_2} min(b_i, k+1)
     *
     * for all 0 <= k < n_1
     */

    /* While this formulation does not require sorting degree2,
     * doing so allows for a linear-time incremental computation
     * of the inequality's right-hand-side.
     */

    *res = true; /* be optimistic */
    lhs_sum = 0;
    partial_rhs_sum = 0; /* the sum of those elements in sorted_deg2 which are <= (k+1) */
    i = 0; /* points past the first element of sorted_deg2 which > (k+1) */
    for (k = 0; k < n1; ++k) {
        lhs_sum += VECTOR(sorted_deg1)[k];

        /* Based on Theorem 3 in [Berger 2014], it is sufficient to do the check
         * for k such that a_k > a_{k+1} and for k=(n_1-1).
         */
        if (k < n1-1 && VECTOR(sorted_deg1)[k] == VECTOR(sorted_deg1)[k+1])
            continue;

        while (i < n2 && VECTOR(sorted_deg2)[i] <= k+1) {
            partial_rhs_sum += VECTOR(sorted_deg2)[i];
            i++;
        }

        /* rhs_sum for a given k is partial_rhs_sum + (n2 - i) * (k+1) */
        if (lhs_sum > partial_rhs_sum + (n2 - i) * (k+1) ) {
            *res = false;
            break;
        }
    }

    igraph_vector_int_destroy(&sorted_deg2);
    igraph_vector_int_destroy(&sorted_deg1);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
