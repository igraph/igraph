/*
   IGraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_community.h"

#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_iterators.h"
#include "igraph_nongraph.h"
#include "igraph_paths.h"
#include "igraph_structural.h"
#include "igraph_transitivity.h"

#include "core/indheap.h"

/**
 * Unweighted local relative density for some vertices.
 *
 * This function ignores self-loops and edge multiplicities.
 * For isolated vertices, zero is returned.
 *
 * \param graph The input graph.
 * \param res Pointer to a vector, the result will be stored here.
 * \param vs Vertex selector, the vertices for which to perform the calculation.
 * \return Error code.
 *
 * Time complexity: TODO.
 */
static igraph_error_t igraph_i_local_relative_density(const igraph_t *graph, igraph_vector_t *res, igraph_vs_t vs) {
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t vs_size;
    igraph_vector_int_t nei_mask; /* which nodes are in the local neighbourhood? */
    igraph_vector_int_t nei_done; /* which local nodes have already been processed? -- avoids duplicate processing in multigraphs */
    igraph_lazy_adjlist_t al;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, IGRAPH_ALL, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&nei_mask, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nei_done, no_of_nodes);

    IGRAPH_CHECK(igraph_vit_create(graph, vs, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    vs_size = IGRAPH_VIT_SIZE(vit);

    IGRAPH_CHECK(igraph_vector_resize(res, vs_size));

    for (igraph_integer_t i=0; ! IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t w = IGRAPH_VIT_GET(vit);
        igraph_integer_t int_count = 0, ext_count = 0;

        igraph_vector_int_t *w_neis = igraph_lazy_adjlist_get(&al, w);
        IGRAPH_CHECK_OOM(w_neis, "Cannot calculate local relative density.");

        igraph_integer_t dw = igraph_vector_int_size(w_neis);

        /* mark neighbours of w, as well as w itself */
        for (igraph_integer_t j=0; j < dw; ++j) {
            VECTOR(nei_mask)[ VECTOR(*w_neis)[j] ] = i + 1;
        }
        VECTOR(nei_mask)[w] = i + 1;

        /* all incident edges of w are internal */
        int_count += dw;
        VECTOR(nei_done)[w] = i + 1;

        for (igraph_integer_t j=0; j < dw; ++j) {
            igraph_integer_t v = VECTOR(*w_neis)[j];

            if (VECTOR(nei_done)[v] == i + 1) {
                continue;
            } else {
                VECTOR(nei_done)[v] = i + 1;
            }

            igraph_vector_int_t *v_neis = igraph_lazy_adjlist_get(&al, v);
            IGRAPH_CHECK_OOM(v_neis, "Cannot calculate local relative density.");

            igraph_integer_t dv = igraph_vector_int_size(v_neis);

            for (igraph_integer_t k=0; k < dv; ++k) {
                igraph_integer_t u = VECTOR(*v_neis)[k];

                if (VECTOR(nei_mask)[u] == i + 1) {
                    int_count += 1;
                } else {
                    ext_count += 1;
                }
            }
        }

        IGRAPH_ASSERT(int_count % 2 == 0);
        int_count /= 2;

        VECTOR(*res)[i] = int_count == 0 ? 0.0 : (igraph_real_t) int_count / (igraph_real_t) (int_count + ext_count);
    }

    igraph_vit_destroy(&vit);
    igraph_vector_int_destroy(&nei_done);
    igraph_vector_int_destroy(&nei_mask);
    igraph_lazy_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/* Weighted local density: we simply multiply the unweighted local relative density with the undirected strength. */
static igraph_error_t weighted_local_density(const igraph_t *graph, igraph_vector_t *res, const igraph_vector_t *weights) {
    igraph_vector_t str;

    IGRAPH_CHECK(igraph_i_local_relative_density(graph, res, igraph_vss_all()));

    IGRAPH_VECTOR_INIT_FINALLY(&str, igraph_vcount(graph));

    IGRAPH_CHECK(igraph_strength(graph, &str, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS, weights));

    igraph_vector_mul(res, &str);

    igraph_vector_destroy(&str);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * Chooses the generated points for the Voronoi partitioning.
 *
 * Each generator has the highest local density within a radius \p r around it.
 *
 * Additionally, if rmax != NULL, the longest distance reached will be stored here.
 * This may be smaller than \p r. This feature is used to determine the largest r
 * value worth considering, through calling this function with r = INFINITY.
 */
static igraph_error_t choose_generators(
        const igraph_t *graph,
        igraph_vector_int_t *generators,
        igraph_real_t *rmax,
        const igraph_vector_t *local_rel_dens,
        const igraph_vector_t *lengths,
        igraph_neimode_t mode,
        igraph_real_t r) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t ord;
    igraph_vector_bool_t excluded;
    igraph_integer_t excluded_count;
    igraph_inclist_t il;
    igraph_2wheap_t q;
    igraph_real_t radius_max;

    /* ord[i] is the index of the ith largest element of local_rel_dens */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ord, 0);
    IGRAPH_CHECK(igraph_vector_qsort_ind(local_rel_dens, &ord, IGRAPH_DESCENDING));

    /* If excluded[v] is true, then v is closer to some already chosen generator than r */
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&excluded, no_of_nodes);
    excluded_count = 0;

    /* The input graph is expected to be simple, but we still set IGRAPH_LOOPS,
     * as inclist_init() performs better this way. */
    IGRAPH_CHECK(igraph_inclist_init(graph, &il, mode, IGRAPH_LOOPS));
    IGRAPH_FINALLY(igraph_inclist_destroy, &il);

    IGRAPH_CHECK(igraph_2wheap_init(&q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &q);

    radius_max = -IGRAPH_INFINITY;
    igraph_vector_int_clear(generators);
    for (igraph_integer_t i=0; i < no_of_nodes; i++) {
        igraph_integer_t g = VECTOR(ord)[i];

        if (VECTOR(excluded)[g]) continue;

        IGRAPH_CHECK(igraph_vector_int_push_back(generators, g));

        igraph_2wheap_clear(&q);
        IGRAPH_CHECK(igraph_2wheap_push_with_index(&q, g, -0.0));
        while (!igraph_2wheap_empty(&q)) {
            igraph_integer_t vid = igraph_2wheap_max_index(&q);
            igraph_real_t mindist = -igraph_2wheap_deactivate_max(&q);

            /* Exceeded cutoff distance, do not search further along this path. */
            if (mindist > r) continue;

            /* Note: We cannot stop the search after hitting an excluded vertex
             * because it is possible that another non-excluded one is reachable only
             * through this one. */
            if (! VECTOR(excluded)[vid]) {
                VECTOR(excluded)[vid] = true;
                excluded_count++;
            }

            if (mindist > radius_max) {
                radius_max = mindist;
            }

            igraph_vector_int_t *inc_edges = igraph_inclist_get(&il, vid);
            igraph_integer_t inc_count = igraph_vector_int_size(inc_edges);
            for (igraph_integer_t j=0; j < inc_count; j++) {
                igraph_integer_t edge = VECTOR(*inc_edges)[j];
                igraph_real_t weight = VECTOR(*lengths)[edge];

                /* Optimization: do not follow infinite-length edges. */
                if (weight == IGRAPH_INFINITY) {
                    continue;
                }

                igraph_integer_t to = IGRAPH_OTHER(graph, edge, vid);
                igraph_real_t altdist = mindist + weight;

                if (!igraph_2wheap_has_elem(&q, to)) {
                    /* This is the first non-infinite distance */
                    IGRAPH_CHECK(igraph_2wheap_push_with_index(&q, to, -altdist));
                } else if (igraph_2wheap_has_active(&q, to)) {
                    igraph_real_t curdist = -igraph_2wheap_get(&q, to);
                    if (altdist < curdist) {
                        /* This is a shorter path */
                        igraph_2wheap_modify(&q, to, -altdist);
                    }
                }
            }
        }

        /* All vertices have been excluded, no need to search further. */
        if (excluded_count == no_of_nodes) break;
    }

    if (rmax) {
        *rmax = radius_max;
    }

    igraph_2wheap_destroy(&q);
    igraph_inclist_destroy(&il);
    igraph_vector_bool_destroy(&excluded);
    igraph_vector_int_destroy(&ord);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}


/* Find the smallest and largest reasonable values of r to consider for the purpose
 * of choosing generator points. */
static igraph_error_t estimate_minmax_r(
        const igraph_t *graph,
        const igraph_vector_t *local_rel_dens,
        const igraph_vector_t *lengths,
        igraph_neimode_t mode,
        igraph_real_t *minr, igraph_real_t *maxr) {

    igraph_vector_int_t generators;

    /* As minimum distance, we use the shortest edge length. This may be shorter than the shortest
     * incident edge of a generator point, but underestimating the minimum distance does not affect
     * the radius optimization negatively. */
    *minr = igraph_vector_min(lengths);

    /* To determine the maximum distance, we run a generator selection with r=INFINITY,
     * and record the longest actual distance encountered in the process. */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&generators, 0);
    IGRAPH_CHECK(choose_generators(graph, &generators, maxr, local_rel_dens, lengths, mode, IGRAPH_INFINITY));
    igraph_vector_int_destroy(&generators);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


typedef igraph_error_t optfun_t(double x, double *res, void *extra);

/* This is the coefficient of the second order part when fitting a quadratic
 * polynomial to the points given in the argument. */
static igraph_real_t coeff2(
        igraph_real_t x1, igraph_real_t x2, igraph_real_t x3,
        igraph_real_t f1, igraph_real_t f2, igraph_real_t f3) {
    igraph_real_t num = x1*(f3 - f2) + x2*(f1 - f3) + x3*(f2 - f1);
    igraph_real_t denom = (x1 - x2)*(x1 - x3)*(x2 - x3);
    return num / denom;
}

/* Given the stationary point of the quadratic fit to the given points */
static igraph_real_t peakx(
        igraph_real_t x1, igraph_real_t x2, igraph_real_t x3,
        igraph_real_t f1, igraph_real_t f2, igraph_real_t f3) {
    igraph_real_t x1s = x1*x1, x2s = x2*x2, x3s = x3*x3;
    igraph_real_t num   = f3 * (x1s - x2s) + f1 * (x2s - x3s) + f2 * (x3s - x1s);
    igraph_real_t denom = f3 * (x1 - x2)   + f1 * (x2 - x3)   + f2 * (x3 - x1);
    return 0.5 * num / denom;
}

/**
 * Simple Brent's method optimizer, with some specializations for the use
 * case at hand (see code comments). It must be called with x2 > x1.
 * The optimal argument is the last one for which f() is invoked.
 * f() is expected to record this in 'extra'.
 */
static igraph_error_t brent_opt(optfun_t *f, igraph_real_t x1, igraph_real_t x2, void *extra) {
    igraph_real_t lo = x1, hi = x2;

    IGRAPH_ASSERT(isfinite(lo));
    IGRAPH_ASSERT(isfinite(hi));

    /* We choose the initial x3 point to be closer to x1 than x2.
     * This is so that if f1 == f2, the next computed point (newx)
     * would not coincide with x3. */
    igraph_real_t x3 = 0.6*x1 + 0.4*x2;
    igraph_real_t f1, f2, f3;

    IGRAPH_CHECK(f(x1, &f1, extra));

    /* Catch special case that would wreak havoc in the optimizer. */
    if (x1 == x2) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_CHECK(f(x2, &f2, extra));
    IGRAPH_CHECK(f(x3, &f3, extra));

    /* We expect that the middle point, f3, is greater than the boundary points. */

    /* Currently, we do not handle the case when f3 < f1. */
    if (f1 > f3) {
        IGRAPH_ERROR("Optimizer did not converge while maximizing modularity for Voronoi communities.",
                     IGRAPH_DIVERGED);
    }

    /* It sometimes happens in disconnected graphs that the maximum is reached at or near the
     * top of the radius range. If so, we bisect the (x3, x2) interval to search for a configuration
     * where f3 >= f2. */
    if (f2 > f3) {
        /* Limit iterations to 'maxiter'. */
        const int maxiter = 10;
        int i;
        for (i=0; i < maxiter; ++i) {
            x1 = x3; f1 = f3;
            x3 = 0.5 * (x1 + x2);
            IGRAPH_CHECK(f(x3, &f3, extra));

            if (f3 >= f2) break;
        }
        /* If no maximum was found in 'maxiter' bisections, just take the upper end of the range. */
        if (i == maxiter) {
            IGRAPH_CHECK(f(x2, &f2, extra));
            return IGRAPH_SUCCESS;
        }
    }

    /* Limit iterations to 20 */
    for (int i=0; i < 20; ++i) {
        igraph_real_t newx, newf;

        newx = peakx(x1, x2, x3, f1, f2, f3);
        IGRAPH_CHECK(f(newx, &newf, extra));

        /* We need to decide whether we drop (x1, f1) or (x2, f2) for the following iterations.
         * The sign of a1 (or a2) determines whether dropping x1 (or x2) yields a convex or concave
         * parabola in the next iteration. We need a negative sign = concave parabola,
         * as we are looking for a maximum. We always keep (x3, f3) as it was the last added point. */
        igraph_real_t a1 = coeff2(x2, x3, newx, f2, f3, newf);
        igraph_real_t a2 = coeff2(x1, x3, newx, f1, f3, newf);

        /* We cannot continue without the Brent optimizer switching to minimization.
         * Terminate search, accepting the current result. */
        if (a1 >= 0 && a2 >= 0) {
            break;
        }

        if (a1 <= a2) {
            x1 = x2;
            x2 = x3;
            x3 = newx;

            f1 = f2;
            f2 = f3;
            f3 = newf;
        } else {
            x2 = x1;
            x1 = x3;
            x3 = newx;

            f2 = f1;
            f1 = f3;
            f3 = newf;
        }

        /* Check if value goes out of initial interval. */
        if (x3 < lo || x3 > hi) {
            IGRAPH_ERROR("Optimizer did not converge while maximizing modularity for Voronoi communities.",
                         IGRAPH_DIVERGED);
        }

        /* We exploit the fact that we are optimizing a discrete valued function, and we can
         * detect convergence by checking that the function value stays exactly the same.
         *
         * As an optimization, we only check whether the two of the three f values are the same.
         * Almost always, when this is the case, another iteration would not yield a better
         * maximum, however, saving a call to f() improves performance noticeably.
         */
        const igraph_real_t eps = 1e-10;
        int c1 = igraph_cmp_epsilon(f1, f3, eps);
        int c2 = igraph_cmp_epsilon(f2, f3, eps);
        if (c1 == 0 || c2 == 0) {
            break;
        }
    }

    return IGRAPH_SUCCESS;
}


/* Work data for get_modularity() */
typedef struct {
    const igraph_t *graph;
    const igraph_vector_t *local_dens;
    const igraph_vector_t *lengths;
    const igraph_vector_t *weights;
    igraph_neimode_t mode;
    igraph_vector_int_t *generators;
    igraph_vector_int_t *membership;
    igraph_real_t modularity;
} get_modularity_work_t;


/* Objective function used with brent_opt(), it computes the modularity for a given radius. */
static igraph_error_t get_modularity(igraph_real_t r, igraph_real_t *modularity, void *extra) {
    get_modularity_work_t *gm = extra;

    IGRAPH_CHECK(choose_generators(gm->graph, gm->generators, NULL,
                                   gm->local_dens, gm->lengths, gm->mode,
                                   r));
    IGRAPH_CHECK(igraph_voronoi(gm->graph, gm->membership, NULL,
                                gm->generators, gm->lengths,
                                gm->mode, IGRAPH_VORONOI_RANDOM));
    IGRAPH_CHECK(igraph_modularity(gm->graph, gm->membership, gm->weights,
                                   1, gm->mode == IGRAPH_ALL ? IGRAPH_UNDIRECTED : IGRAPH_DIRECTED,
                                   &gm->modularity));
    *modularity = gm->modularity;

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_community_voronoi
 * \brief Finds communities using Voronoi partitioning.
 *
 * \experimental
 *
 * This function finds communities using a Voronoi partitioning of vertices based
 * on the given edge lengths divided by the edge clustering coefficient
 * (\ref igraph_ecc()). The generator vertices are chosen to be those with the
 * largest local relative density within a radius \p r, with the local relative
 * density of a vertex defined as
 * <code>s m / (m + k)</code>, where \c s is the strength of the vertex,
 * \c m is the number of edges within the vertex's first order neighborhood,
 * while \c k is the number of edges with only one endpoint within this
 * neighborhood.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Deritei et al., Community detection by graph Voronoi diagrams,
 * New Journal of Physics 16, 063007 (2014)
 * https://doi.org/10.1088/1367-2630/16/6/063007
 *
 * </para><para>
 * Molnár et al., Community Detection in Directed Weighted Networks using Voronoi Partitioning,
 * Scientific Reports 14, 8124 (2024)
 * https://doi.org/10.1038/s41598-024-58624-4
 *
 * \param graph The input graph. It must be simple.
 * \param membership If not \c NULL, the membership of each vertex is returned here.
 * \param generators If not \c NULL, the generator points used for Voronoi partitioning are returned here.
 * \param modularity If not \c NULL, the modularity score of the partitioning is returned here.
 * \param lengths Edge lengths, or \c NULL to consider all edges as having unit length.
 *   Voronoi partitioning will use edge lengths equal to lengths / ECC where ECC is the edge
 *   clustering coefficient.
 * \param weights Edge weights, or \c NULL to consider all edges as having unit weight.
 *   Weights are used when selecting generator points, as well as for computing modularity.
 * \param mode If \c IGRAPH_OUT, distances from generator points to all other nodes are considered.
 *   If \c IGRAPH_IN, the reverse distances are used. If \c IGRAPH_ALL, edge directions are ignored.
 *   This parameter is ignored for undirected graphs.
 * \param r The radius/resolution to use when selecting generator points. The larger this value, the
 *   fewer partitions there will be. Pass in a negative value to automatically select the radius
 *   that maximizes modularity.
 * \return Error code.
 *
 * \sa \ref igraph_voronoi(), \ref igraph_ecc().
 *
 * Time complexity: TODO.
 */
igraph_error_t igraph_community_voronoi(
        const igraph_t *graph,
        igraph_vector_int_t *membership, igraph_vector_int_t *generators,
        igraph_real_t *modularity,
        const igraph_vector_t *lengths, const igraph_vector_t *weights,
        igraph_neimode_t mode, igraph_real_t r) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_t local_rel_dens;
    igraph_vector_t lengths2; /* lengths2 = lengths / ecc */
    igraph_vector_int_t imembership, igenerators;
    igraph_vector_int_t *pmembership, *pgenerators;
    igraph_bool_t simple;

    if (! igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (lengths && igraph_vector_size(lengths) != no_of_edges) {
        IGRAPH_ERROR("Edge length vector size does not match edge count.", IGRAPH_EINVAL);
    }

    if (weights && igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERROR("Edge length vector size does not match edge count.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (! simple) {
        IGRAPH_ERROR("The graph must be simple for Voronoi communities.", IGRAPH_EINVAL);
    }
    if (igraph_is_directed(graph) && mode == IGRAPH_ALL) {
        igraph_bool_t has_mutual;
        /* When the graph is directed but edge directions are ignored,
         * mutual edges are effectively multi-edges. */
        IGRAPH_CHECK(igraph_has_mutual(graph, &has_mutual, false));
        if (has_mutual) {
            IGRAPH_ERROR("The graph must be simple for Voronoi communities. "
                         "Mutual directed edges are effectively multi-edges when ignoring edge directions.",
                         IGRAPH_EINVAL);
        }
    }

    if (no_of_edges == 0) {
        /* Also handles no_of_nodes <= 1 */
        if (membership) {
            IGRAPH_CHECK(igraph_vector_int_range(membership, 0, no_of_nodes));
        }
        if (generators) {
            IGRAPH_CHECK(igraph_vector_int_range(generators, 0, no_of_nodes));
        }
        return IGRAPH_SUCCESS;
    }

    if (! generators) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&igenerators, no_of_nodes);
        pgenerators = &igenerators;
    } else {
        pgenerators = generators;
    }

    if (! membership) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&imembership, no_of_nodes);
        pmembership = &imembership;
    } else {
        pmembership = membership;
    }

    if (lengths) {
        igraph_real_t m = igraph_vector_min(lengths);
        if (isnan(m)) {
            IGRAPH_ERROR("Edge lengths must not be NaN.", IGRAPH_EINVAL);
        }
        if (m < 0) {
            IGRAPH_ERROR("Edge lengths must be non-negative.", IGRAPH_EINVAL);
        }
    }

    if (weights) {
        igraph_real_t m = igraph_vector_min(weights);
        if (isnan(m)) {
            IGRAPH_ERROR("Edge weights must not be NaN.", IGRAPH_EINVAL);
        }
        if (m <= 0) {
            IGRAPH_ERROR("Edge weights must be positive.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&local_rel_dens, 0);

    IGRAPH_CHECK(weighted_local_density(graph, &local_rel_dens, weights));

    IGRAPH_VECTOR_INIT_FINALLY(&lengths2, 0);
    IGRAPH_CHECK(igraph_ecc(graph, &lengths2, igraph_ess_all(IGRAPH_EDGEORDER_ID), 3, true, true));

    /* Note: ECC is never NaN but it may be Inf */
    for (igraph_integer_t i=0; i < no_of_edges; i++) {
        VECTOR(lengths2)[i] = 1 / (VECTOR(lengths2)[i]);
    }
    if (lengths) {
        igraph_vector_mul(&lengths2, lengths);
    }

    if (r < 0) {
        igraph_real_t minr, maxr;

        IGRAPH_CHECK(estimate_minmax_r(graph, &local_rel_dens, &lengths2, mode, &minr, &maxr));

        get_modularity_work_t gm = {
                graph,
                &local_rel_dens,
                &lengths2,
                weights,
                mode,
                pgenerators,
                pmembership,
                /* modularity */ IGRAPH_NAN
        };
        IGRAPH_CHECK(brent_opt(get_modularity, minr, maxr, &gm));
        if (modularity) {
            *modularity = gm.modularity;
        }
    } else {
        IGRAPH_CHECK(choose_generators(graph, pgenerators, NULL, &local_rel_dens, &lengths2, mode, r));
        IGRAPH_CHECK(igraph_voronoi(graph, membership, NULL, pgenerators, &lengths2, mode, IGRAPH_VORONOI_RANDOM));
        if (modularity) {
            IGRAPH_CHECK(igraph_modularity(graph, membership, weights, 1,
                                           mode == IGRAPH_ALL ? IGRAPH_UNDIRECTED : IGRAPH_DIRECTED, modularity));
        }
    }

    if (! generators) {
        igraph_vector_int_destroy(&igenerators);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (! membership) {
        igraph_vector_int_destroy(&imembership);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&local_rel_dens);
    igraph_vector_destroy(&lengths2);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}
