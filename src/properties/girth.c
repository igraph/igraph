/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include "igraph_structural.h"

#include "igraph_adjlist.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"

#include "core/interruption.h"

#include <limits.h>

/**
 * \function igraph_girth
 * \brief The girth of a graph is the length of the shortest cycle in it.
 *
 * </para><para>
 * The current implementation works for undirected graphs only,
 * directed graphs are treated as undirected graphs. Self-loops and
 * multiple edges are ignored.
 *
 * </para><para>
 * For graphs that contain no cycles, and only for such graphs,
 * zero is returned. Note that in some applications, it is customary
 * to define the girth of acyclic graphs to be infinity. However, infinity
 * is not representable as an \c igraph_integer_t, therefore zero is used
 * for this case.
 *
 * </para><para>
 * This implementation is based on Alon Itai and Michael Rodeh:
 * Finding a minimum circuit in a graph
 * \emb Proceedings of the ninth annual ACM symposium on Theory of
 * computing \eme, 1-10, 1977. The first implementation of this
 * function was done by Keith Briggs, thanks Keith.
 * \param graph The input graph.
 * \param girth Pointer to an integer, if not \c NULL then the result
 *     will be stored here.
 * \param circle Pointer to an initialized vector, the vertex ids in
 *     the shortest circle will be stored here. If \c NULL then it is
 *     ignored.
 * \return Error code.
 *
 * Time complexity: O((|V|+|E|)^2), |V| is the number of vertices, |E|
 * is the number of edges in the general case. If the graph has no
 * cycles at all then the function needs O(|V|+|E|) time to realize
 * this and then it stops.
 *
 * \example examples/simple/igraph_girth.c
 */
int igraph_girth(const igraph_t *graph, igraph_integer_t *girth,
                 igraph_vector_t *circle) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_t q;
    igraph_lazy_adjlist_t adjlist;
    long int mincirc = LONG_MAX, minvertex = 0;
    long int node;
    igraph_bool_t triangle = 0;
    igraph_vector_int_t *neis;
    igraph_vector_long_t level;
    long int stoplevel = no_of_nodes + 1;
    igraph_bool_t anycircle = 0;
    long int t1 = 0, t2 = 0;

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vector_long_init(&level, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &level);

    for (node = 0; !triangle && node < no_of_nodes; node++) {

        /* Are there circles in this graph at all? */
        if (node == 1 && anycircle == 0) {
            igraph_bool_t conn;
            IGRAPH_CHECK(igraph_is_connected(graph, &conn, IGRAPH_WEAK));
            if (conn) {
                /* No, there are none */
                break;
            }
        }

        anycircle = 0;
        igraph_dqueue_clear(&q);
        igraph_vector_long_null(&level);
        IGRAPH_CHECK(igraph_dqueue_push(&q, node));
        VECTOR(level)[node] = 1;

        IGRAPH_ALLOW_INTERRUPTION();

        while (!igraph_dqueue_empty(&q)) {
            long int actnode = (long int) igraph_dqueue_pop(&q);
            long int actlevel = VECTOR(level)[actnode];
            long int i, n;

            if (actlevel >= stoplevel) {
                break;
            }

            neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
            n = igraph_vector_int_size(neis);
            for (i = 0; i < n; i++) {
                long int nei = (long int) VECTOR(*neis)[i];
                long int neilevel = VECTOR(level)[nei];
                if (neilevel != 0) {
                    if (neilevel == actlevel - 1) {
                        continue;
                    } else {
                        /* found circle */
                        stoplevel = neilevel;
                        anycircle = 1;
                        if (actlevel < mincirc) {
                            /* Is it a minimum circle? */
                            mincirc = actlevel + neilevel - 1;
                            minvertex = node;
                            t1 = actnode; t2 = nei;
                            if (neilevel == 2) {
                                /* Is it a triangle? */
                                triangle = 1;
                            }
                        }
                        if (neilevel == actlevel) {
                            break;
                        }
                    }
                } else {
                    igraph_dqueue_push(&q, nei);
                    VECTOR(level)[nei] = actlevel + 1;
                }
            }

        } /* while q !empty */
    } /* node */

    if (girth) {
        if (mincirc == LONG_MAX) {
            *girth = mincirc = 0;
        } else {
            *girth = (igraph_integer_t) mincirc;
        }
    }

    /* Store the actual circle, if needed */
    if (circle) {
        IGRAPH_CHECK(igraph_vector_resize(circle, mincirc));
        if (mincirc != 0) {
            long int i, n, idx = 0;
            igraph_dqueue_clear(&q);
            igraph_vector_long_null(&level); /* used for father pointers */
#define FATHER(x) (VECTOR(level)[(x)])
            IGRAPH_CHECK(igraph_dqueue_push(&q, minvertex));
            FATHER(minvertex) = minvertex;
            while (FATHER(t1) == 0 || FATHER(t2) == 0) {
                long int actnode = (long int) igraph_dqueue_pop(&q);
                neis = igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actnode);
                n = igraph_vector_int_size(neis);
                for (i = 0; i < n; i++) {
                    long int nei = (long int) VECTOR(*neis)[i];
                    if (FATHER(nei) == 0) {
                        FATHER(nei) = actnode + 1;
                        igraph_dqueue_push(&q, nei);
                    }
                }
            }  /* while q !empty */
            /* Ok, now use FATHER to create the path */
            while (t1 != minvertex) {
                VECTOR(*circle)[idx++] = t1;
                t1 = FATHER(t1) - 1;
            }
            VECTOR(*circle)[idx] = minvertex;
            idx = mincirc - 1;
            while (t2 != minvertex) {
                VECTOR(*circle)[idx--] = t2;
                t2 = FATHER(t2) - 1;
            }
        } /* anycircle */
    } /* circle */
#undef FATHER

    igraph_vector_long_destroy(&level);
    igraph_dqueue_destroy(&q);
    igraph_lazy_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}
