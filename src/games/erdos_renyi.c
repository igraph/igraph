/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"
#include "igraph_nongraph.h"
#include "igraph_random.h"

/**
 * \section about_games
 *
 * <para>Games are randomized graph generators. Randomization means that
 * they generate a different graph every time you call them. </para>
 */

int igraph_erdos_renyi_game_gnp(
    igraph_t *graph, igraph_integer_t n, igraph_real_t p,
    igraph_bool_t directed, igraph_bool_t loops
) {

    long int no_of_nodes = n;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int retval = 0;
    long int vsize;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (p < 0.0 || p > 1.0) {
        IGRAPH_ERROR("Invalid probability given", IGRAPH_EINVAL);
    }

    if (p == 0.0 || no_of_nodes <= 1) {
        IGRAPH_CHECK(retval = igraph_empty(graph, n, directed));
    } else if (p == 1.0) {
        IGRAPH_CHECK(retval = igraph_full(graph, n, directed, loops));
    } else {

        long int i;
        double maxedges = n, last;
        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&s, (long int) (maxedges * p * 1.1)));

        RNG_BEGIN();

        last = RNG_GEOM(p);
        while (last < maxedges) {
            IGRAPH_CHECK(igraph_vector_push_back(&s, last));
            last += RNG_GEOM(p);
            last += 1;
        }

        RNG_END();

        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s) * 2));

        vsize = igraph_vector_size(&s);
        if (directed && loops) {
            for (i = 0; i < vsize; i++) {
                long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else if (directed && !loops) {
            for (i = 0; i < vsize; i++) {
                long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                if (from == to) {
                    to = no_of_nodes - 1;
                }
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else if (!directed && loops) {
            for (i = 0; i < vsize; i++) {
                long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        } else { /* !directed && !loops */
            for (i = 0; i < vsize; i++) {
                long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
                igraph_vector_push_back(&edges, from);
                igraph_vector_push_back(&edges, to);
            }
        }

        igraph_vector_destroy(&s);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_CHECK(retval = igraph_create(graph, &edges, n, directed));
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return retval;
}

int igraph_erdos_renyi_game_gnm(
    igraph_t *graph, igraph_integer_t n, igraph_real_t m,
    igraph_bool_t directed, igraph_bool_t loops
) {

    igraph_integer_t no_of_nodes = n;
    igraph_integer_t no_of_edges = (igraph_integer_t) m;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    igraph_vector_t s = IGRAPH_VECTOR_NULL;
    int retval = 0;

    if (n < 0) {
        IGRAPH_ERROR("Invalid number of vertices", IGRAPH_EINVAL);
    }
    if (m < 0) {
        IGRAPH_ERROR("Invalid number of edges", IGRAPH_EINVAL);
    }

    if (m == 0.0 || no_of_nodes <= 1) {
        IGRAPH_CHECK(retval = igraph_empty(graph, n, directed));
    } else {

        long int i;
        double maxedges = n;
        if (directed && loops) {
            maxedges *= n;
        } else if (directed && !loops) {
            maxedges *= (n - 1);
        } else if (!directed && loops) {
            maxedges *= (n + 1) / 2.0;
        } else {
            maxedges *= (n - 1) / 2.0;
        }

        if (no_of_edges > maxedges) {
            IGRAPH_ERROR("Invalid number (too large) of edges", IGRAPH_EINVAL);
        }

        if (maxedges == no_of_edges) {
            retval = igraph_full(graph, n, directed, loops);
        } else {

            long int slen;

            IGRAPH_VECTOR_INIT_FINALLY(&s, 0);
            IGRAPH_CHECK(igraph_random_sample(&s, 0, maxedges - 1,
                                              (igraph_integer_t) no_of_edges));

            IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
            IGRAPH_CHECK(igraph_vector_reserve(&edges, igraph_vector_size(&s) * 2));

            slen = igraph_vector_size(&s);
            if (directed && loops) {
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor(VECTOR(s)[i] / no_of_nodes);
                    long int from = (long int) (VECTOR(s)[i] - ((igraph_real_t)to) * no_of_nodes);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else if (directed && !loops) {
                for (i = 0; i < slen; i++) {
                    long int from = (long int) floor(VECTOR(s)[i] / (no_of_nodes - 1));
                    long int to = (long int) (VECTOR(s)[i] - ((igraph_real_t)from) * (no_of_nodes - 1));
                    if (from == to) {
                        to = no_of_nodes - 1;
                    }
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else if (!directed && loops) {
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) - 1) / 2);
                    long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to + 1)) / 2);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            } else { /* !directed && !loops */
                for (i = 0; i < slen; i++) {
                    long int to = (long int) floor((sqrt(8 * VECTOR(s)[i] + 1) + 1) / 2);
                    long int from = (long int) (VECTOR(s)[i] - (((igraph_real_t)to) * (to - 1)) / 2);
                    igraph_vector_push_back(&edges, from);
                    igraph_vector_push_back(&edges, to);
                }
            }

            igraph_vector_destroy(&s);
            IGRAPH_FINALLY_CLEAN(1);
            retval = igraph_create(graph, &edges, n, directed);
            igraph_vector_destroy(&edges);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return retval;
}

/**
 * \ingroup generators
 * \function igraph_erdos_renyi_game
 * \brief Generates a random (Erdős-Rényi) graph.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param type The type of the random graph, possible values:
 *        \clist
 *        \cli IGRAPH_ERDOS_RENYI_GNM
 *          G(n,m) graph,
 *          m edges are
 *          selected uniformly randomly in a graph with
 *          n vertices.
 *        \cli IGRAPH_ERDOS_RENYI_GNP
 *          G(n,p) graph,
 *          every possible edge is included in the graph with
 *          probability p.
 *        \endclist
 * \param n The number of vertices in the graph.
 * \param p_or_m This is the p parameter for
 *        G(n,p) graphs and the
 *        m
 *        parameter for G(n,m) graphs.
 * \param directed Logical, whether to generate a directed graph.
 * \param loops Logical, whether to generate loops (self) edges.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid
 *         \p type, \p n,
 *         \p p or \p m
 *          parameter.
 *         \c IGRAPH_ENOMEM: there is not enough
 *         memory for the operation.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges in the graph.
 *
 * \sa \ref igraph_barabasi_game(), \ref igraph_growing_random_game()
 *
 * \example examples/simple/igraph_erdos_renyi_game.c
 */
int igraph_erdos_renyi_game(igraph_t *graph, igraph_erdos_renyi_t type,
                            igraph_integer_t n, igraph_real_t p_or_m,
                            igraph_bool_t directed, igraph_bool_t loops) {
    int retval = 0;

    if (type == IGRAPH_ERDOS_RENYI_GNP) {
        retval = igraph_erdos_renyi_game_gnp(graph, n, p_or_m, directed, loops);
    } else if (type == IGRAPH_ERDOS_RENYI_GNM) {
        retval = igraph_erdos_renyi_game_gnm(graph, n, p_or_m, directed, loops);
    } else {
        IGRAPH_ERROR("Invalid type", IGRAPH_EINVAL);
    }

    return retval;
}
