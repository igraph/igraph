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

#include "igraph_qsort.h"

#define IGRAPH_I_MULTI_EDGES_SW 0x02 /* 010, more than one edge allowed between distinct vertices */
#define IGRAPH_I_MULTI_LOOPS_SW 0x04 /* 100, more than one self-loop allowed on the same vertex   */

static int igraph_i_is_graphical_undirected_multi_loops(const igraph_vector_t *degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_undirected_loopless_multi(const igraph_vector_t *degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_undirected_loopy_simple(const igraph_vector_t *degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_undirected_simple(const igraph_vector_t *degrees, igraph_bool_t *res);

static int igraph_i_is_graphical_directed_loopy_multi(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_directed_loopless_multi(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_directed_loopy_simple(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res);
static int igraph_i_is_graphical_directed_simple(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res);

static int igraph_i_is_bigraphical_multi(const igraph_vector_t *degrees1, const igraph_vector_t *degrees2, igraph_bool_t *res);
static int igraph_i_is_bigraphical_simple(const igraph_vector_t *degrees1, const igraph_vector_t *degrees2, igraph_bool_t *res);


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
 * Time complexity: O(n^2) for simple directed graphs, O(n log n) for graphs with self-loops,
 * and O(n) for all other cases, where n is the length of the degree sequence(s).
 */
int igraph_is_graphical(const igraph_vector_t *out_degrees,
                        const igraph_vector_t *in_degrees,
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
            /* Remainig case:
             *  - At most one self-loop per vertex but multi-edges between distinct vertices allowed.
             * These cases cannot currently be requested through the documented API,
             * so no explanatory error message for now. */
            return IGRAPH_UNIMPLEMENTED;
        }
    }
    /* Directed case: */
    else
    {
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
            /* Remainig cases:
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
int igraph_is_bigraphical(const igraph_vector_t *degrees1,
                          const igraph_vector_t *degrees2,
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
static int igraph_i_is_graphical_undirected_multi_loops(const igraph_vector_t *degrees, igraph_bool_t *res) {
    long int sum_parity = 0; /* 0 if the degree sum is even, 1 if it is odd */
    long int n = igraph_vector_size(degrees);
    long int i;

    for (i=0; i < n; ++i) {
        long int d = VECTOR(*degrees)[i];

        if (d < 0) {
            *res = 0;
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
static int igraph_i_is_graphical_undirected_loopless_multi(const igraph_vector_t *degrees, igraph_bool_t *res) {
    long int i;
    long int n = igraph_vector_size(degrees);
    long int dsum, dmax;

    /* Zero-length sequences are considered graphical. */
    if (n == 0) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }

    dsum = 0; dmax = 0;
    for (i=0; i < n; ++i) {
        long int d = VECTOR(*degrees)[i];

        if (d < 0) {
            *res = 0;
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
static int igraph_i_is_graphical_undirected_loopy_simple(const igraph_vector_t *degrees, igraph_bool_t *res) {
    igraph_vector_t work;
    long int w, b, s, c, n, k;

    n = igraph_vector_size(degrees);

    /* Zero-length sequences are considered graphical. */
    if (n == 0) {
        *res = 1;
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

    IGRAPH_CHECK(igraph_vector_copy(&work, degrees));
    IGRAPH_FINALLY(igraph_vector_destroy, &work);

    igraph_vector_reverse_sort(&work);

    *res = 1;
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
            *res = 0;
            break;
        }
        if (w == k) {
            break;
        }
    }

    igraph_vector_destroy(&work);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/* Undirected simple graph:
 *  - Degrees must be non-negative.
 *  - The sum of degrees must be even.
 *  - Use the Erdős-Gallai theorem.
 */
static int igraph_i_is_graphical_undirected_simple(const igraph_vector_t *degrees, igraph_bool_t *res) {
    igraph_vector_int_t num_degs; /* num_degs[d] is the # of vertices with degree d */
    const long int p = igraph_vector_size(degrees);
    long int dmin, dmax, dsum;
    long int n; /* number of non-zero degrees */
    long int k, sum_deg, sum_ni, sum_ini;
    long int i, dk;
    long int zverovich_bound;

    if (p == 0) {
        *res = 1;
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
    for (i=0; i < p; ++i) {
        long int d = VECTOR(*degrees)[i];

        if (d < 0 || d >= p) {
            *res = 0;
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
        *res = 0;
        goto finish;
    }

    if (n == 0) {
        *res = 1;
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
        *res = 1;
        goto finish;
    }

    k = 0; sum_deg = 0; sum_ni = 0; sum_ini = 0;
    for (dk = dmax; dk >= dmin; --dk) {
        long int run_size, v;

        if (dk < k+1) {
            *res = 1;
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
                *res = 0;
                goto finish;
            }
        }
    }

    *res = 1;

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
static int igraph_i_is_graphical_directed_loopy_multi(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    long int sumdiff; /* difference between sum of in- and out-degrees */
    long int n = igraph_vector_size(out_degrees);
    long int i;

    if (igraph_vector_size(in_degrees) != n) {
        IGRAPH_ERROR("The length of out- and in-degree sequences must be the same.", IGRAPH_EINVAL);
    }

    sumdiff = 0;
    for (i=0; i < n; ++i) {
        long int dout = VECTOR(*out_degrees)[i];
        long int din  = VECTOR(*in_degrees)[i];

        if (dout < 0 || din < 0) {
            *res = 0;
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
static int igraph_i_is_graphical_directed_loopless_multi(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    long int i, sumin, sumout, dmax;
    long int n = igraph_vector_size(out_degrees);

    if (igraph_vector_size(in_degrees) != n) {
        IGRAPH_ERROR("The length of out- and in-degree sequences must be the same.", IGRAPH_EINVAL);
    }

    sumin = 0; sumout = 0;
    dmax = 0;
    for (i=0; i < n; ++i) {
        long int dout = VECTOR(*out_degrees)[i];
        long int din  = VECTOR(*in_degrees)[i];
        long int d = dout + din;

        if (dout < 0 || din < 0) {
            *res = 0;
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
static int igraph_i_is_graphical_directed_loopy_simple(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    long int n = igraph_vector_size(out_degrees);

    if (igraph_vector_size(in_degrees) != n) {
        IGRAPH_ERROR("The length of out- and in-degree sequences must be the same.", IGRAPH_EINVAL);
    }

    return igraph_i_is_bigraphical_simple(out_degrees, in_degrees, res);
}


/* Directed simple graph:
 *  - Degrees must be non-negative.
 *  - The sum of in- and out-degrees must be the same.
 *  - Use the Fulkerson-Chen-Anstee theorem
 */

typedef struct {
    const igraph_vector_t* first;
    const igraph_vector_t* second;
} igraph_i_qsort_dual_vector_cmp_data_t;

static int igraph_i_qsort_dual_vector_cmp_desc(void* data, const void *p1, const void *p2) {
    igraph_i_qsort_dual_vector_cmp_data_t* sort_data =
        (igraph_i_qsort_dual_vector_cmp_data_t*)data;
    long int index1 = *((long int*)p1);
    long int index2 = *((long int*)p2);
    if (VECTOR(*sort_data->first)[index1] < VECTOR(*sort_data->first)[index2]) {
        return 1;
    }
    if (VECTOR(*sort_data->first)[index1] > VECTOR(*sort_data->first)[index2]) {
        return -1;
    }
    if (VECTOR(*sort_data->second)[index1] < VECTOR(*sort_data->second)[index2]) {
        return 1;
    }
    if (VECTOR(*sort_data->second)[index1] > VECTOR(*sort_data->second)[index2]) {
        return -1;
    }
    return 0;
}

static int igraph_i_is_graphical_directed_simple(const igraph_vector_t *out_degrees, const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    igraph_vector_long_t index_array;
    long int i, j, vcount, lhs, rhs;
    igraph_i_qsort_dual_vector_cmp_data_t sort_data;

    /* The conditions from the loopy multigraph case are necessary here as well. */
    IGRAPH_CHECK(igraph_i_is_graphical_directed_loopy_multi(out_degrees, in_degrees, res));
    if (! *res) {
        return IGRAPH_SUCCESS;
    }

    vcount = igraph_vector_size(out_degrees);
    if (vcount == 0) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }

    /* Create an index vector that sorts the vertices by decreasing in-degree */    
    IGRAPH_CHECK(igraph_vector_long_init_seq(&index_array, 0, vcount - 1));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &index_array);

    /* Set up the auxiliary struct for sorting */
    sort_data.first  = in_degrees;
    sort_data.second = out_degrees;

    /* Sort the index vector */
    igraph_qsort_r(VECTOR(index_array), vcount, sizeof(long int), &sort_data,
                   igraph_i_qsort_dual_vector_cmp_desc);

    /* Be optimistic, then check whether the Fulkerson–Chen–Anstee condition
     * holds for every k. In particular, for every k in [0; n), it must be true
     * that:
     *
     * \sum_{i=0}^k indegree[i] <=
     *     \sum_{i=0}^k min(outdegree[i], k) +
     *     \sum_{i=k+1}^{n-1} min(outdegree[i], k + 1)
     */

#define INDEGREE(x) (VECTOR(*in_degrees)[VECTOR(index_array)[x]])
#define OUTDEGREE(x) (VECTOR(*out_degrees)[VECTOR(index_array)[x]])

    *res = 1;
    lhs = 0;
    for (i = 0; i < vcount; i++) {
        lhs += INDEGREE(i);

        /* It is enough to check for indexes where the in-degree is about to
         * decrease in the next step; see "Stronger condition" in the Wikipedia
         * entry for the Fulkerson-Chen-Anstee condition */
        if (i != vcount - 1 && INDEGREE(i) == INDEGREE(i + 1)) {
            continue;
        }

        rhs = 0;
        for (j = 0; j <= i; j++) {
            rhs += OUTDEGREE(j) < i ? OUTDEGREE(j) : i;
        }
        for (j = i + 1; j < vcount; j++) {
            rhs += OUTDEGREE(j) < (i + 1) ? OUTDEGREE(j) : (i + 1);
        }

        if (lhs > rhs) {
            *res = 0;
            break;
        }
    }

#undef INDEGREE
#undef OUTDEGREE

    igraph_vector_long_destroy(&index_array);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}



/***** Bipartite case *****/

/* Bipartite graph with multi-eges:
 *  - Degrees must be non-negative.
 *  - Sum of degrees must be the same in the two partitions.
 */
static int igraph_i_is_bigraphical_multi(const igraph_vector_t *degrees1, const igraph_vector_t *degrees2, igraph_bool_t *res) {
    long int i;
    long int sum1, sum2;
    long int n1 = igraph_vector_size(degrees1), n2 = igraph_vector_size(degrees2);

    sum1 = 0;
    for (i=0; i < n1; ++i) {
        long int d = VECTOR(*degrees1)[i];

        if (d < 0) {
            *res = 0;
            return IGRAPH_SUCCESS;
        }

        sum1 += d;
    }

    sum2 = 0;
    for (i=0; i < n2; ++i) {
        long int d = VECTOR(*degrees2)[i];

        if (d < 0) {
            *res = 0;
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
static int igraph_i_is_bigraphical_simple(const igraph_vector_t *degrees1, const igraph_vector_t *degrees2, igraph_bool_t *res) {
    igraph_vector_t sorted_deg1, sorted_deg2;
    long int n1 = igraph_vector_size(degrees1), n2 = igraph_vector_size(degrees2);
    long int i, k;
    long lhs_sum, partial_rhs_sum;

    if (n1 == 0 && n2 == 0) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }

    /* The conditions from the multigraph case are necessary here as well. */
    IGRAPH_CHECK(igraph_i_is_bigraphical_multi(degrees1, degrees2, res));
    if (! *res) {
        return IGRAPH_SUCCESS;
    }

    /* Ensure that degrees1 is the shorter vector as a minor optimization: */
    if (n2 < n1) {
        const igraph_vector_t *tmp;
        long int n;

        tmp = degrees1;
        degrees1 = degrees2;
        degrees2 = tmp;

        n = n1;
        n1 = n2;
        n2 = n;
    }

    /* Copy and sort both vectors: */

    IGRAPH_CHECK(igraph_vector_copy(&sorted_deg1, degrees1));
    IGRAPH_FINALLY(igraph_vector_destroy, &sorted_deg1);
    igraph_vector_reverse_sort(&sorted_deg1); /* decreasing sort */

    IGRAPH_CHECK(igraph_vector_copy(&sorted_deg2, degrees2));
    IGRAPH_FINALLY(igraph_vector_destroy, &sorted_deg2);
    igraph_vector_sort(&sorted_deg2); /* increasing sort */

    /*
     * We follow the description of the Gale-Ryser theorem in:
     *
     * A. Berger, A note on the characterization of digraphic sequences, Discrete Math. 314, 38 (2014).
     * http://dx.doi.org/10.1016/j.disc.2013.09.010
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

    *res = 1; /* be optimistic */
    lhs_sum = 0;
    partial_rhs_sum = 0; /* the sum of those elements in sorted_deg2 which are <= (k+1) */
    i = 0; /* points past the first element of sorted_deg2 which > (k+1) */
    for (k=0; k < n1; ++k) {
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
            *res = 0;
            break;
        }
    }

    igraph_vector_destroy(&sorted_deg2);
    igraph_vector_destroy(&sorted_deg1);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}


/***** Legacy functions *****/

#define SUCCEED {   \
        if (res) {        \
            *res = 1;       \
        }                 \
        return IGRAPH_SUCCESS; \
    }

#define FAIL {   \
        if (res) {     \
            *res = 0;    \
        }              \
        return IGRAPH_SUCCESS; \
    }

/**
 * \function igraph_is_degree_sequence
 * \brief Determines whether a degree sequence is valid.
 *
 * \deprecated-by igraph_is_graphical 0.9
 *
 * </para><para>
 * A sequence of n integers is a valid degree sequence if there exists some
 * graph where the degree of the i-th vertex is equal to the i-th element of the
 * sequence. Note that the graph may contain multiple or loop edges; if you are
 * interested in whether the degrees of some \em simple graph may realize the
 * given sequence, use \ref igraph_is_graphical_degree_sequence.
 *
 * </para><para>
 * In particular, the function checks whether all the degrees are non-negative.
 * For undirected graphs, it also checks whether the sum of degrees is even.
 * For directed graphs, the function checks whether the lengths of the two
 * degree vectors are equal and whether their sums are also equal. These are
 * known sufficient and necessary conditions for a degree sequence to be
 * valid.
 *
 * \param out_degrees  an integer vector specifying the degree sequence for
 *     undirected graphs or the out-degree sequence for directed graphs.
 * \param in_degrees   an integer vector specifying the in-degrees of the
 *     vertices for directed graphs. For undirected graphs, this must be null.
 * \param res  pointer to a boolean variable, the result will be stored here
 * \return Error code.
 *
 * Time complexity: O(n), where n is the length of the degree sequence.
 */
int igraph_is_degree_sequence(const igraph_vector_t *out_degrees,
                              const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    IGRAPH_WARNING("igraph_is_degree_sequence is deprecated, use igraph_is_graphical.");

    /* degrees must be non-negative */
    if (igraph_vector_any_smaller(out_degrees, 0)) {
        FAIL;
    }
    if (in_degrees && igraph_vector_any_smaller(in_degrees, 0)) {
        FAIL;
    }

    if (in_degrees == 0) {
        /* sum of degrees must be even */
        if (((long int)igraph_vector_sum(out_degrees) % 2) != 0) {
            FAIL;
        }
    } else {
        /* length of the two degree vectors must be equal */
        if (igraph_vector_size(out_degrees) != igraph_vector_size(in_degrees)) {
            FAIL;
        }
        /* sum of in-degrees must be equal to sum of out-degrees */
        if (igraph_vector_sum(out_degrees) != igraph_vector_sum(in_degrees)) {
            FAIL;
        }
    }

    SUCCEED;
}

#undef SUCCEED
#undef FAIL


/**
 * \function igraph_is_graphical_degree_sequence
 * \brief Determines whether a sequence of integers can be the degree sequence of some simple graph.
 *
 * \deprecated-by igraph_is_graphical 0.9
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Hakimi SL: On the realizability of a set of integers as degrees of the
 * vertices of a simple graph. J SIAM Appl Math 10:496-506, 1962.
 *
 * </para><para>
 * PL Erdős, I Miklós and Z Toroczkai: A simple Havel-Hakimi type algorithm
 * to realize graphical degree sequences of directed graphs.
 * The Electronic Journal of Combinatorics 17(1):R66, 2010.
 * https://dx.doi.org/10.1017/S0963548317000499
 *
 * </para><para>
 * Z Kiraly: Recognizing graphic degree sequences and generating all
 * realizations. TR-2011-11, Egervary Research Group, H-1117, Budapest,
 * Hungary. ISSN 1587-4451, 2012.
 * https://www.cs.elte.hu/egres/tr/egres-11-11.pdf
 *
 * \param out_degrees  an integer vector specifying the degree sequence for
 *     undirected graphs or the out-degree sequence for directed graphs.
 * \param in_degrees   an integer vector specifying the in-degrees of the
 *     vertices for directed graphs. For undirected graphs, this must be null.
 * \param res  pointer to a boolean variable, the result will be stored here
 * \return Error code.
 *
 * Time complexity: O(n log n) for undirected graphs, O(n^2) for directed
 *                  graphs, where n is the length of the degree sequence.
 */
int igraph_is_graphical_degree_sequence(const igraph_vector_t *out_degrees,
                                        const igraph_vector_t *in_degrees, igraph_bool_t *res) {
    IGRAPH_WARNING("igraph_is_graphical_degree_sequence is deprecated, use igraph_is_graphical.");
    return igraph_is_graphical(out_degrees, in_degrees, IGRAPH_SIMPLE_SW, res);
}
