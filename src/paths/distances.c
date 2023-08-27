/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "igraph_adjlist.h"
#include "igraph_datatype.h"
#include "igraph_dqueue.h"
#include "igraph_iterators.h"
#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_random.h"
#include "igraph_vector.h"

#include "core/interruption.h"
#include "core/indheap.h"

/* When vid_ecc is not NULL, only one vertex ID should be passed in vids.
 * vid_ecc will then return the id of the vertex farthest from the one in
 * vids. If unconn == false and not all other vertices were reachable from
 * the single given vertex, -1 is returned in vid_ecc. */
static igraph_error_t igraph_i_eccentricity(const igraph_t *graph,
                                 igraph_vector_t *res,
                                 igraph_vs_t vids,
                                 igraph_lazy_adjlist_t *adjlist,
                                 igraph_integer_t *vid_ecc,
                                 igraph_bool_t unconn) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_vector_int_t counted;
    igraph_integer_t i, mark = 1;
    igraph_integer_t min_degree = 0;

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&counted, no_of_nodes);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_fill(res, -1);

    for (i = 0, IGRAPH_VIT_RESET(vit);
         !IGRAPH_VIT_END(vit);
         IGRAPH_VIT_NEXT(vit), mark++, i++) {

        igraph_integer_t source;
        igraph_integer_t nodes_reached = 1;
        source = IGRAPH_VIT_GET(vit);
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, source));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        VECTOR(counted)[source] = mark;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t act = igraph_dqueue_int_pop(&q);
            igraph_integer_t dist = igraph_dqueue_int_pop(&q);
            igraph_vector_int_t *neis = igraph_lazy_adjlist_get(adjlist, act);
            igraph_integer_t j, n;

            IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");

            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                igraph_integer_t nei = VECTOR(*neis)[j];
                if (VECTOR(counted)[nei] != mark) {
                    VECTOR(counted)[nei] = mark;
                    nodes_reached++;
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, dist + 1));
                }
            }
            if (vid_ecc) {
                /* Return the vertex ID of the vertex which has the lowest
                 * degree of the vertices most distant from the starting
                 * vertex. Assumes there is only 1 vid in vids. Used for
                 * pseudo_diameter calculations. */
                if (dist > VECTOR(*res)[i] || (dist == VECTOR(*res)[i] && n < min_degree)) {
                    VECTOR(*res)[i] = dist;
                    *vid_ecc = act;
                    min_degree = n;
                }
            } else if (dist > VECTOR(*res)[i]) {
                VECTOR(*res)[i] = dist;
            }
        } /* while !igraph_dqueue_int_empty(dqueue) */

        if (nodes_reached != no_of_nodes && !unconn && vid_ecc) {
            *vid_ecc = -1;
            break;
        }
    } /* for IGRAPH_VIT_NEXT(vit) */

    igraph_vector_int_destroy(&counted);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * This function finds the weighted eccentricity and returns it via \p ecc.
 * It's used for igraph_pseudo_diameter_dijkstra() and igraph_eccentricity_dijkstra().
 * \p vid_ecc returns the vertex id of the ecc with the greatest
 * distance from \p vid_start. If two vertices have the same greatest distance,
 * the one with the lowest degree is chosen.
 * When the graph is not (strongly) connected and \p unconn is false, then \p ecc
 * wil be set to infinity, and \p vid_ecc to -1;
 */
static igraph_error_t igraph_i_eccentricity_dijkstra(
        const igraph_t *graph,
        const igraph_vector_t *weights,
        igraph_real_t *ecc,
        igraph_integer_t vid_start,
        igraph_integer_t *vid_ecc,
        igraph_bool_t unconn,
        igraph_lazy_inclist_t *inclist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_2wheap_t Q;
    igraph_vector_t vec_dist;
    igraph_integer_t i;
    igraph_real_t degree_ecc, dist;
    igraph_integer_t degree_i;
    igraph_vector_int_t *neis;

    IGRAPH_VECTOR_INIT_FINALLY(&vec_dist, no_of_nodes);
    igraph_vector_fill(&vec_dist, IGRAPH_INFINITY);
    IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
    IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);

    igraph_2wheap_clear(&Q);
    igraph_2wheap_push_with_index(&Q, vid_start, -1.0);

    while (!igraph_2wheap_empty(&Q)) {
        igraph_integer_t minnei = igraph_2wheap_max_index(&Q);
        igraph_real_t mindist = -igraph_2wheap_deactivate_max(&Q);
        igraph_integer_t nlen;

        VECTOR(vec_dist)[minnei] = mindist - 1.0;

        /* Now check all neighbors of 'minnei' for a shorter path */
        neis = igraph_lazy_inclist_get(inclist, minnei);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");
        nlen = igraph_vector_int_size(neis);
        for (i = 0; i < nlen; i++) {
            igraph_integer_t edge = VECTOR(*neis)[i];
            igraph_integer_t tto = IGRAPH_OTHER(graph, edge, minnei);
            igraph_real_t altdist = mindist + VECTOR(*weights)[edge];
            igraph_bool_t active = igraph_2wheap_has_active(&Q, tto);
            igraph_bool_t has = igraph_2wheap_has_elem(&Q, tto);
            igraph_real_t curdist = active ? -igraph_2wheap_get(&Q, tto) : 0.0;

            if (altdist == IGRAPH_INFINITY) {
                /* Ignore edges with positive infinite weights */
            } else if (!has) {
                /* This is the first non-infinite distance */
                IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
            } else if (altdist < curdist) {
                /* This is a shorter path */
                igraph_2wheap_modify(&Q, tto, -altdist);
            }
        }
    }

    *ecc = 0;
    *vid_ecc = vid_start;
    degree_ecc = 0;

    for (i = 0; i < no_of_nodes; i++) {
        if (i == vid_start) {
            continue;
        }
        dist = VECTOR(vec_dist)[i];

        /* inclist is used to ignore multiple edges when finding the degree */
        neis = igraph_lazy_inclist_get(inclist, i);
        IGRAPH_CHECK_OOM(neis, "Failed to query incident edges.");

        degree_i  = igraph_vector_int_size(neis);

        if (dist > *ecc) {
            if (!isfinite(dist)) {
                if (!unconn) {
                    *ecc = IGRAPH_INFINITY;
                    *vid_ecc = -1;
                    break;
                }
            } else {
                *ecc = dist;
                *vid_ecc = i;
                degree_ecc = degree_i;
            }
        } else if (dist == *ecc) {
            if (degree_i < degree_ecc) {
                degree_ecc = degree_i;
                *vid_ecc = i;
            }
        }
    }
    igraph_2wheap_destroy(&Q);
    igraph_vector_destroy(&vec_dist);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_eccentricity
 * \brief Eccentricity of some vertices.
 *
 * The eccentricity of a vertex is calculated by measuring the shortest
 * distance from (or to) the vertex, to (or from) all vertices in the
 * graph, and taking the maximum.
 *
 * </para><para>
 * This implementation ignores vertex pairs that are in different
 * components. Isolated vertices have eccentricity zero.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param res Pointer to an initialized vector, the result is stored
 *    here.
 * \param vids The vertices for which the eccentricity is calculated.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(v*(|V|+|E|)), where |V| is the number of
 * vertices, |E| is the number of edges and v is the number of
 * vertices for which eccentricity is calculated.
 *
 * \sa \ref igraph_radius().
 *
 * \example examples/simple/igraph_eccentricity.c
 */

igraph_error_t igraph_eccentricity(const igraph_t *graph,
                        igraph_vector_t *res,
                        igraph_vs_t vids,
                        igraph_neimode_t mode) {
    igraph_lazy_adjlist_t adjlist;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode,
                                          IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_i_eccentricity(graph, res, vids, &adjlist,
                                       /*vid_ecc*/ NULL, /*unconn*/ 1));
    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_eccentricity_dijkstra
 * \brief Eccentricity of some vertices, using weighted edges.
 *
 * The eccentricity of a vertex is calculated by measuring the shortest
 * distance from (or to) the vertex, to (or from) all vertices in the
 * graph, and taking the maximum.
 *
 * </para><para>
 * This implementation ignores vertex pairs that are in different
 * components. Isolated vertices have eccentricity zero.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_eccentricity() is called. Edges with positive
 *    infinite weights are ignored.
 * \param res Pointer to an initialized vector, the result is stored
 *    here.
 * \param vids The vertices for which the eccentricity is calculated.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V| |E| log|V| + |V|), where |V| is the number of
 * vertices, |E| the number of edges.
 *
 */

igraph_error_t igraph_eccentricity_dijkstra(const igraph_t *graph,
                        const igraph_vector_t *weights,
                        igraph_vector_t *res,
                        igraph_vs_t vids,
                        igraph_neimode_t mode) {
    igraph_lazy_inclist_t inclist;
    igraph_vit_t vit;
    igraph_integer_t dump;
    igraph_real_t ecc;
    igraph_integer_t no_of_edges = igraph_ecount(graph);

    if (weights == NULL) {
        return igraph_eccentricity(graph, res, vids, mode);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }

    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        } else if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode,
                                          IGRAPH_NO_LOOPS));
    IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_resize(res, 0));
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));

    for (IGRAPH_VIT_RESET(vit);
            !IGRAPH_VIT_END(vit);
            IGRAPH_VIT_NEXT(vit)) {
        IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc, IGRAPH_VIT_GET(vit), /*vid_ecc*/ &dump, /*unconn*/ 1, &inclist));
        IGRAPH_CHECK(igraph_vector_push_back(res, ecc));
    }
    igraph_lazy_inclist_destroy(&inclist);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_radius
 * \brief Radius of a graph.
 *
 * The radius of a graph is the defined as the minimum eccentricity of
 * its vertices, see \ref igraph_eccentricity().
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param radius Pointer to a real variable, the result is stored
 *   here.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V|(|V|+|E|)), where |V| is the number of
 * vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_radius_dijkstra() for the weighted version,
 * \ref igraph_diameter() for the maximum eccentricity,
 * \ref igraph_eccentricity() for the eccentricities of all vertices.
 *
 * \example examples/simple/igraph_radius.c
 */

igraph_error_t igraph_radius(const igraph_t *graph, igraph_real_t *radius,
                             igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (no_of_nodes == 0) {
        *radius = IGRAPH_NAN;
    } else {
        igraph_vector_t ecc;
        IGRAPH_VECTOR_INIT_FINALLY(&ecc, igraph_vcount(graph));
        IGRAPH_CHECK(igraph_eccentricity(graph, &ecc, igraph_vss_all(),
                                         mode));
        *radius = igraph_vector_min(&ecc);
        igraph_vector_destroy(&ecc);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_radius_dijkstra
 * \brief Radius of a graph, using weighted edges.
 *
 * \experimental
 *
 * The radius of a graph is the defined as the minimum eccentricity of
 * its vertices, see \ref igraph_eccentricity().
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_radius() is called. Edges with positive
 *    infinite weights are ignored.
 * \param radius Pointer to a real variable, the result is stored
 *   here.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V| |E| log|V| + |V|), where |V| is the number of
 * vertices, |E| the number of edges.
 *
 * \sa \ref igraph_radius() for the unweighted version,
 * \ref igraph_diameter_dijkstra() for the maximum weighted eccentricity,
 * \ref igraph_eccentricity_dijkstra() for weighted eccentricities of
 * all vertices.
 */

igraph_error_t igraph_radius_dijkstra(const igraph_t *graph, const igraph_vector_t *weights,
                                      igraph_real_t *radius, igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    if (weights == NULL) {
        return igraph_radius(graph, radius, mode);
    }

    if (no_of_nodes == 0) {
        *radius = IGRAPH_NAN;
    } else {
        igraph_vector_t ecc;
        IGRAPH_VECTOR_INIT_FINALLY(&ecc, igraph_vcount(graph));
        IGRAPH_CHECK(igraph_eccentricity_dijkstra(graph, weights, &ecc, igraph_vss_all(), mode));
        *radius = igraph_vector_min(&ecc);
        igraph_vector_destroy(&ecc);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_pseudo_diameter
 * \brief Approximation and lower bound of diameter.
 *
 * This algorithm finds a pseudo-peripheral vertex and returns its
 * eccentricity. This value can be used as an approximation
 * and lower bound of the diameter of a graph.
 *
 * </para><para>
 * A pseudo-peripheral vertex is a vertex v, such that for every
 * vertex u which is as far away from v as possible, v is also as
 * far away from u as possible. The process of finding one depends
 * on where the search starts, and for a disconnected graph the
 * maximum diameter found will be that of the component \p vid_start
 * is in.
 *
 * \param graph The input graph, if it is directed, its edge directions
 *        are ignored.
 * \param diameter Pointer to a real variable, the result is stored
 *        here.
 * \param vid_start Id of the starting vertex. If this is negative, a
 *        random starting vertex is chosen.
 * \param from Pointer to an integer, if not \c NULL it will be set to the
 *        source vertex of the diameter path. If \p unconn is \c false, and
 *        a disconnected graph is detected, this is set to -1.
 * \param to Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If \p unconn is \c false, and
 *        a disconnected graph is detected, this is set to -1.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c true the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is returned.
 * \return Error code.
 *
 * Time complexity: O(|V||E|)), where |V| is the number of
 * vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_eccentricity(), \ref igraph_diameter().
 *
 */
igraph_error_t igraph_pseudo_diameter(const igraph_t *graph,
                           igraph_real_t *diameter,
                           igraph_integer_t vid_start,
                           igraph_integer_t *from,
                           igraph_integer_t *to,
                           igraph_bool_t directed,
                           igraph_bool_t unconn) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_real_t ecc_v;
    igraph_real_t ecc_u;
    igraph_integer_t vid_ecc;
    igraph_integer_t ito, ifrom;
    igraph_bool_t inf = false;

    if (vid_start >= no_of_nodes) {
        IGRAPH_ERROR("Starting vertex ID for pseudo-diameter out of range.", IGRAPH_EINVAL);
    }

    /* We will reach here when vid_start < 0 and the graph has no vertices. */
    if (no_of_nodes == 0) {
        if (diameter) {
            *diameter = IGRAPH_NAN;
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (vid_start < 0) {
        RNG_BEGIN();
        vid_start = RNG_INTEGER(0, no_of_nodes - 1);
        RNG_END();
    }

    if (!igraph_is_directed(graph) || !directed) {
        igraph_lazy_adjlist_t adjlist;
        igraph_vector_t ecc_vec;

        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL,
                                              IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
        IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
        ifrom = vid_start;
        IGRAPH_VECTOR_INIT_FINALLY(&ecc_vec, no_of_nodes);

        IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_vec, igraph_vss_1(vid_start),
                                           &adjlist, &vid_ecc, unconn));
        ecc_u = VECTOR(ecc_vec)[0];

        if (!unconn && vid_ecc == -1) {
            inf = true;
        } else {
            while (true) {
                IGRAPH_ALLOW_INTERRUPTION();

                ito = vid_ecc;

                IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_vec, igraph_vss_1(vid_ecc),
                                                   &adjlist, &vid_ecc, 1));

                ecc_v = VECTOR(ecc_vec)[0];

                if (ecc_u < ecc_v) {
                    ecc_u = ecc_v;
                    ifrom = ito;
                } else {
                    break;
                }
            }
        }
        igraph_vector_destroy(&ecc_vec);
        igraph_lazy_adjlist_destroy(&adjlist);
        IGRAPH_FINALLY_CLEAN(2);
    } else {
        igraph_vector_t ecc_out;
        igraph_vector_t ecc_in;
        igraph_integer_t vid_ecc_in;
        igraph_integer_t vid_ecc_out;
        igraph_integer_t vid_end;
        igraph_bool_t direction;
        igraph_lazy_adjlist_t adjlist_in;
        igraph_lazy_adjlist_t adjlist_out;

        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist_in, IGRAPH_IN,
                                              IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
        IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist_in);
        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist_out, IGRAPH_OUT,
                                              IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
        IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist_out);

        IGRAPH_VECTOR_INIT_FINALLY(&ecc_in, igraph_vcount(graph));
        IGRAPH_VECTOR_INIT_FINALLY(&ecc_out, igraph_vcount(graph));

        IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_out, igraph_vss_1(vid_start),
                                           &adjlist_out, &vid_ecc_out, unconn));
        IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_in, igraph_vss_1(vid_start),
                                           &adjlist_in, &vid_ecc_in, unconn));

        /* A directed graph is strongly connected iff all vertices are reachable
         * from vid_start both when moving along or moving opposite the edge directions. */
        if (!unconn && (vid_ecc_out == -1 || vid_ecc_in == -1)) {
            inf = true;
        } else {
            if (VECTOR(ecc_out)[0] > VECTOR(ecc_in)[0]) {
                vid_ecc = vid_ecc_out;
                ecc_u = VECTOR(ecc_out)[0];
            } else {
                vid_ecc = vid_ecc_in;
                ecc_u = VECTOR(ecc_in)[0];
            }

            while (1) {
                IGRAPH_ALLOW_INTERRUPTION();

                vid_end = vid_ecc;

                /* TODO: In the undirected case, we break ties between vertices at the
                 * same distance based on their degree. In te directed case, should we
                 * use in-, out- or total degree? */
                IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_out, igraph_vss_1(vid_ecc),
                                                   &adjlist_out, &vid_ecc_out, 1));
                IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_in, igraph_vss_1(vid_ecc),
                                                   &adjlist_in, &vid_ecc_in, 1));

                if (VECTOR(ecc_out)[0] > VECTOR(ecc_in)[0]) {
                    vid_ecc = vid_ecc_out;
                    ecc_v = VECTOR(ecc_out)[0];
                    direction = 1;
                } else {
                    vid_ecc = vid_ecc_in;
                    ecc_v = VECTOR(ecc_in)[0];
                    direction = 0;
                }

                if (ecc_u < ecc_v) {
                    ecc_u = ecc_v;
                    vid_start = vid_end;
                } else {
                    break;
                }
            }

            if (direction) {
                ifrom = vid_end;
                ito   = vid_start;
            } else {
                ifrom = vid_start;
                ito   = vid_end;
            }

        }
        igraph_vector_destroy(&ecc_out);
        igraph_vector_destroy(&ecc_in);
        igraph_lazy_adjlist_destroy(&adjlist_in);
        igraph_lazy_adjlist_destroy(&adjlist_out);
        IGRAPH_FINALLY_CLEAN(4);
    }

    if (inf) {
        if (diameter) {
            *diameter = IGRAPH_INFINITY;
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
    } else {
        if (diameter) {
            *diameter = ecc_u;
        }
        if (from) {
            *from = ifrom;
        }
        if (to) {
            *to = ito;
        }
    }

    return IGRAPH_SUCCESS;
}


/**
 * \function igraph_pseudo_diameter_dijkstra
 * \brief Approximation and lower bound of the diameter of a weighted graph.
 *
 * This algorithm finds a pseudo-peripheral vertex and returns its
 * weighted eccentricity. This value can be used as an approximation
 * and lower bound of the diameter of a graph.
 *
 * </para><para>
 * A pseudo-peripheral vertex is a vertex v, such that for every
 * vertex u which is as far away from v as possible, v is also as
 * far away from u as possible. The process of finding one depends
 * on where the search starts, and for a disconnected graph the
 * maximum diameter found will be that of the component \p vid_start
 * is in.
 *
 * </para><para>
 * If the graph has no vertices, \c IGRAPH_NAN is returned.
 *
 * \param graph The input graph, can be directed or undirected.
 * \param weights The edge weights of the graph. Can be \c NULL for an
 *        unweighted graph. All weights should be non-negative. Edges with
 *        positive infinite weights are ignored.
 * \param diameter This will contain the weighted pseudo-diameter.
 * \param vid_start Id of the starting vertex. If this is negative, a
 *        random starting vertex is chosen.
 * \param from If not \c NULL this will be set to the
 *        source vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param to If not \c NULL this will be set to the
 *        target vertex of the diameter path. If the graph has no diameter path,
 *        it will be set to -1.
 * \param directed Boolean, whether to consider directed
 *        paths. Ignored for undirected graphs.
 * \param unconn What to do if the graph is not connected. If
 *        \c true the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is
 *        returned.
 * \return Error code.
 *
 * Time complexity: O(|V||E|*log|E|), |V| is the number of vertices,
 * |E| is the number of edges.
 *
 * \sa \ref igraph_diameter_dijkstra()
 */
igraph_error_t igraph_pseudo_diameter_dijkstra(const igraph_t *graph,
                                    const igraph_vector_t *weights,
                                    igraph_real_t *diameter,
                                    igraph_integer_t vid_start,
                                    igraph_integer_t *from,
                                    igraph_integer_t *to,
                                    igraph_bool_t directed,
                                    igraph_bool_t unconn) {


    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_real_t ecc_v;
    igraph_real_t ecc_u;
    igraph_integer_t vid_ecc;
    igraph_integer_t ito, ifrom;
    igraph_bool_t inf = false;

    if (vid_start >= no_of_nodes) {
        IGRAPH_ERROR("Starting vertex ID for pseudo-diameter out of range.", IGRAPH_EINVAL);
    }

    if (!weights) {
        return igraph_pseudo_diameter(graph, diameter, vid_start, from, to, directed, unconn);
    }

    if (igraph_vector_size(weights) != no_of_edges) {
        IGRAPH_ERRORF("Weight vector length (%" IGRAPH_PRId ") does not match number of edges (%" IGRAPH_PRId ").",
                      IGRAPH_EINVAL,
                      igraph_vector_size(weights), no_of_edges);
    }
    if (no_of_edges > 0) {
        igraph_real_t min = igraph_vector_min(weights);
        if (min < 0) {
            IGRAPH_ERRORF("Weight vector must be non-negative, got %g.", IGRAPH_EINVAL, min);
        }
        else if (isnan(min)) {
            IGRAPH_ERROR("Weight vector must not contain NaN values.", IGRAPH_EINVAL);
        }
    }

    /* We will reach here when vid_start < 0 and the graph has no vertices. */
    if (no_of_nodes == 0) {
        if (diameter) {
            *diameter = IGRAPH_NAN;
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
        return IGRAPH_SUCCESS;
    }

    if (vid_start < 0) {
        RNG_BEGIN();
        vid_start = RNG_INTEGER(0, no_of_nodes - 1);
        RNG_END();
    }

    if (!igraph_is_directed(graph) || !directed) {
        igraph_lazy_inclist_t inclist;
        IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

        ifrom = vid_start;

        IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_u, vid_start, &vid_ecc, unconn, &inclist));

        inf = !isfinite(ecc_u);

        if (!inf) {
            while (1) {
                IGRAPH_ALLOW_INTERRUPTION();

                ito = vid_ecc;
                IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_v, vid_ecc, &vid_ecc, unconn, &inclist));

                if (ecc_u < ecc_v) {
                    ecc_u = ecc_v;
                    ifrom = ito;
                } else {
                    break;
                }
            }
        }
        igraph_lazy_inclist_destroy(&inclist);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_real_t ecc_out;
        igraph_real_t ecc_in;
        igraph_integer_t vid_ecc_in;
        igraph_integer_t vid_ecc_out;
        igraph_integer_t vid_end;
        igraph_bool_t direction;
        igraph_lazy_inclist_t inclist_out;
        igraph_lazy_inclist_t inclist_in;

        IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist_out, IGRAPH_OUT, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist_out);
        IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist_in, IGRAPH_IN, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist_in);


        IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_out, vid_start, &vid_ecc_out, unconn, &inclist_out));
        IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_in, vid_start, &vid_ecc_in, unconn, &inclist_in));

        /* A directed graph is strongly connected iff all vertices are reachable
         * from vid_start both when moving along or moving opposite the edge directions. */
        if (!unconn && (vid_ecc_out == -1 || vid_ecc_in == -1)) {
            inf = true;
        } else {
            if (ecc_out > ecc_in) {
                vid_ecc = vid_ecc_out;
                ecc_u = ecc_out;
            } else {
                vid_ecc = vid_ecc_in;
                ecc_u = ecc_in;
            }

            while (1) {
                IGRAPH_ALLOW_INTERRUPTION();

                vid_end = vid_ecc;

                /* TODO: In the undirected case, we break ties between vertices at the
                 * same distance based on their degree. In te directed case, should we
                 * use in-, out- or total degree? */
                IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_out, vid_ecc, &vid_ecc_out, unconn, &inclist_out));
                IGRAPH_CHECK(igraph_i_eccentricity_dijkstra(graph, weights, &ecc_in, vid_ecc, &vid_ecc_in, unconn, &inclist_in));

                if (ecc_out > ecc_in) {
                    vid_ecc = vid_ecc_out;
                    ecc_v = ecc_out;
                    direction = 1;
                } else {
                    vid_ecc = vid_ecc_in;
                    ecc_v = ecc_in;
                    direction = 0;
                }

                if (ecc_u < ecc_v) {
                    ecc_u = ecc_v;
                    vid_start = vid_end;
                } else {
                    break;
                }
            }

            if (direction) {
                ifrom = vid_end;
                ito   = vid_start;
            } else {
                ifrom = vid_start;
                ito   = vid_end;
            }
        }
        igraph_lazy_inclist_destroy(&inclist_out);
        igraph_lazy_inclist_destroy(&inclist_in);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (inf) {
        if (diameter) {
            *diameter = IGRAPH_INFINITY;
        }
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
    } else {
        if (diameter) {
            *diameter = ecc_u;
        }
        if (from) {
            *from = ifrom;
        }
        if (to) {
            *to = ito;
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_graph_center
 * \brief Central vertices of a graph.
 *
 * The central vertices of a graph are calculated by finding the vertices
 * with the minimum eccentricity. This concept is typically applied to
 * (strongly) connected graphs. In disconnected graphs, the smallest
 * eccentricity is taken across all components.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param res Pointer to an initialized vector, the result is stored
 *    here.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V| (|V|+|E|)), where |V| is the number of
 * vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_graph_center_dijkstra(),
 * \ref igraph_eccentricity(), \ref igraph_radius()
 *
 */
igraph_error_t igraph_graph_center(
    const igraph_t *graph, igraph_vector_int_t *res, igraph_neimode_t mode
) {

    igraph_vector_t ecc;

    igraph_vector_int_clear(res);
    if (igraph_vcount(graph) == 0) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&ecc, 0);
    IGRAPH_CHECK(igraph_eccentricity(graph, &ecc, igraph_vss_all(), mode));

    /* igraph_eccentricity() does not return infinity or NaN, and the null graph
     * case was handled above, therefore calling vector_min() is safe. */
    igraph_real_t min_eccentricity = igraph_vector_min(&ecc);
    igraph_integer_t n = igraph_vector_size(&ecc);
    for (igraph_integer_t i = 0; i < n; i++) {
        if (VECTOR(ecc)[i] == min_eccentricity) {
            IGRAPH_CHECK(igraph_vector_int_push_back(res, i));
        }
    }

    igraph_vector_destroy(&ecc);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_graph_center_dijkstra
 * \brief Central vertices of a graph, using weighted edges.
 *
 * \experimental
 *
 * The central vertices of a graph are calculated by finding the vertices
 * with the minimum eccentricity. This function takes edge weights into
 * account and uses Dijkstra's algorithm for the shortest path calculation.
 * The concept of the graph center is typically applied to
 * (strongly) connected graphs. In disconnected graphs, the smallest
 * eccentricity is taken across all components.
 *
 * \param graph The input graph, it can be directed or undirected.
 * \param weights The edge weights. All edge weights must be
 *    non-negative for Dijkstra's algorithm to work. Additionally, no
 *    edge weight may be NaN. If either case does not hold, an error
 *    is returned. If this is a null pointer, then the unweighted
 *    version, \ref igraph_graph_center() is called. Edges with positive
 *    infinite weights are ignored.
 * \param res Pointer to an initialized vector, the result is stored
 *    here.
 * \param mode What kind of paths to consider for the calculation:
 *    \c IGRAPH_OUT, paths that follow edge directions;
 *    \c IGRAPH_IN, paths that follow the opposite directions; and
 *    \c IGRAPH_ALL, paths that ignore edge directions. This argument
 *    is ignored for undirected graphs.
 * \return Error code.
 *
 * Time complexity: O(|V| |E| log|V| + |V|), where |V| is the number of
 * vertices, |E| the number of edges.
 *
 * \sa \ref igraph_graph_center(),
 * \ref igraph_eccentricity_dijkstra(), \ref igraph_radius_dijkstra()
 *
 */
igraph_error_t igraph_graph_center_dijkstra(
    const igraph_t *graph, const igraph_vector_t *weights, igraph_vector_int_t *res, igraph_neimode_t mode
) {

    igraph_vector_t ecc;
    const igraph_real_t eps = IGRAPH_SHORTEST_PATH_EPSILON;

    if (weights == NULL) {
        return igraph_graph_center(graph, res, mode);
    }

    igraph_vector_int_clear(res);
    if (igraph_vcount(graph) == 0) {
        return IGRAPH_SUCCESS;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&ecc, 0);
    IGRAPH_CHECK(igraph_eccentricity_dijkstra(graph, weights, &ecc, igraph_vss_all(), mode));

    /* igraph_eccentricity_dijkstra() does not return infinity or NaN, and the null graph
     * case was handled above, therefore calling vector_min() is safe. */
    igraph_real_t min_eccentricity = igraph_vector_min(&ecc);
    igraph_integer_t n = igraph_vector_size(&ecc);
    for (igraph_integer_t i = 0; i < n; i++) {
        if (igraph_cmp_epsilon(VECTOR(ecc)[i], min_eccentricity, eps) == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(res, i));
        }
    }

    igraph_vector_destroy(&ecc);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
