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

#include "igraph_datatype.h"
#include "igraph_dqueue.h"
#include "igraph_iterators.h"
#include "igraph_vector.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"
#include "igraph_random.h"

#include "core/interruption.h"

/* When vid_ecc is not NULL, only one vertex id should be passed in vids.
 * vid_ecc will then return the id of the vertex farthest from the one in
 * vids. If unconn == FALSE and not all other vertices were reachable from
 * the single given vertex, -1 is returned n vid_ecc. */
static int igraph_i_eccentricity(const igraph_t *graph,
                                 igraph_vector_t *res,
                                 igraph_vs_t vids,
                                 igraph_lazy_adjlist_t *adjlist,
                                 igraph_integer_t *vid_ecc,
                                 igraph_bool_t unconn) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_long_t q;
    igraph_vit_t vit;
    igraph_vector_int_t counted;
    long int i, mark = 1;
    igraph_integer_t min_degree = 0;

    IGRAPH_CHECK(igraph_dqueue_long_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    IGRAPH_CHECK(igraph_vector_int_init(&counted, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &counted);

    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_fill(res, -1);

    for (i = 0, IGRAPH_VIT_RESET(vit);
         !IGRAPH_VIT_END(vit);
         IGRAPH_VIT_NEXT(vit), mark++, i++) {

        long int source;
        long int nodes_reached = 1;
        source = IGRAPH_VIT_GET(vit);
        IGRAPH_CHECK(igraph_dqueue_long_push(&q, source));
        IGRAPH_CHECK(igraph_dqueue_long_push(&q, 0));
        VECTOR(counted)[source] = mark;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_long_empty(&q)) {
            long int act = igraph_dqueue_long_pop(&q);
            long int dist = igraph_dqueue_long_pop(&q);
            igraph_vector_int_t *neis = igraph_lazy_adjlist_get(adjlist, act);
            long int j, n;

            n = igraph_vector_int_size(neis);
            for (j = 0; j < n; j++) {
                long int nei = VECTOR(*neis)[j];
                if (VECTOR(counted)[nei] != mark) {
                    VECTOR(counted)[nei] = mark;
                    nodes_reached++;
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, nei));
                    IGRAPH_CHECK(igraph_dqueue_long_push(&q, dist + 1));
                }
            }
            if (vid_ecc) {
                /* Return the vertex id of the vertex which has the lowest
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
        } /* while !igraph_dqueue_long_empty(dqueue) */

        if (nodes_reached != no_of_nodes && !unconn && vid_ecc) {
            *vid_ecc = -1;
            break;
        }
    } /* for IGRAPH_VIT_NEXT(vit) */

    igraph_vector_int_destroy(&counted);
    igraph_vit_destroy(&vit);
    igraph_dqueue_long_destroy(&q);
    IGRAPH_FINALLY_CLEAN(3);

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

int igraph_eccentricity(const igraph_t *graph,
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
 * \sa \ref igraph_eccentricity().
 *
 * \example examples/simple/igraph_radius.c
 */

int igraph_radius(const igraph_t *graph, igraph_real_t *radius,
                  igraph_neimode_t mode) {

    int no_of_nodes = igraph_vcount(graph);

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

    return 0;
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
 *        source vertex of the diameter path. If \p unconn is FALSE, and
 *        a disconnected graph is detected, this is set to -1.
 * \param to Pointer to an integer, if not \c NULL it will be set to the
 *        target vertex of the diameter path. If \p unconn is FALSE, and
 *        a disconnected graph is detected, this is set to -1.
 * \param unconn What to do if the graph is not connected. If
 *        \c TRUE the longest geodesic within a component
 *        will be returned, otherwise \c IGRAPH_INFINITY is returned.
 * \return Error code.
 *
 * Time complexity: O(|V||E|)), where |V| is the number of
 * vertices and |E| is the number of edges.
 *
 * \sa \ref igraph_eccentricity(), \ref igraph_diameter().
 *
 */
int igraph_pseudo_diameter(const igraph_t *graph,
                           igraph_real_t *diameter,
                           igraph_integer_t vid_start,
                           igraph_integer_t *from,
                           igraph_integer_t *to,
                           igraph_bool_t unconn) {

    int no_of_nodes = igraph_vcount(graph);
    igraph_real_t ecc_v;
    igraph_real_t ecc_u;
    igraph_integer_t vid_ecc;
    igraph_bool_t inf = 0;

    if (vid_start >= no_of_nodes) {
        IGRAPH_ERROR("Starting vertex id for pseudo-diameter out of range.", IGRAPH_EINVAL);
    }


    if (vid_start < 0) {
        RNG_BEGIN();
        vid_start = RNG_INTEGER(0, no_of_nodes - 1);
        RNG_END();
    }

    if (!igraph_is_directed(graph)) {
        igraph_lazy_adjlist_t adjlist;
        igraph_vector_t ecc_vec;

        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL,
                                              IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
        IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
        if (from) {
            *from = vid_start;
        }
        IGRAPH_VECTOR_INIT_FINALLY(&ecc_vec, no_of_nodes);

        IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_vec, igraph_vss_1(vid_start),
                                           &adjlist, &vid_ecc, unconn));
        ecc_u = VECTOR(ecc_vec)[0];

        if (!unconn && vid_ecc == -1) {
            inf = 1;
        } else {
            while (1) {
                if (to) {
                    *to = vid_ecc;
                }

                IGRAPH_CHECK(igraph_i_eccentricity(graph, &ecc_vec, igraph_vss_1(vid_ecc),
                                                   &adjlist, &vid_ecc, 1));

                ecc_v = VECTOR(ecc_vec)[0];

                if (ecc_u < ecc_v) {
                    ecc_u = ecc_v;
                    if (from) {
                        *from = *to;
                    }
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

        if (!unconn && (vid_ecc_out == -1 || vid_ecc_in == -1)) {
            inf = 1;
        } else {
            if (VECTOR(ecc_out)[0] > VECTOR(ecc_in)[0]) {
                vid_ecc = vid_ecc_out;
                ecc_u = VECTOR(ecc_out)[0];
            } else {
                vid_ecc = vid_ecc_in;
                ecc_u = VECTOR(ecc_in)[0];
            }

            while (1) {
                vid_end = vid_ecc;

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

            if (from) {
                if (direction) {
                    *from = vid_end;
                } else {
                    *from = vid_start;
                }
            }
            if (to) {
                if (direction) {
                    *to = vid_start;
                } else {
                    *to = vid_end;
                }
            }
        }
        igraph_vector_destroy(&ecc_out);
        igraph_vector_destroy(&ecc_in);
        igraph_lazy_adjlist_destroy(&adjlist_in);
        igraph_lazy_adjlist_destroy(&adjlist_out);
        IGRAPH_FINALLY_CLEAN(4);
    }

    if (inf) {
        *diameter = IGRAPH_INFINITY;
        if (from) {
            *from = -1;
        }
        if (to) {
            *to = -1;
        }
    } else {
        *diameter = ecc_u;
    }

    return IGRAPH_SUCCESS;
}
