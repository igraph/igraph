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

#include "igraph_paths.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/interruption.h"

#include <string.h>  /* memset */

static void igraph_i_gasp_paths_destroy(igraph_vector_ptr_t *v) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(v); i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_destroy(VECTOR(*v)[i]);
            IGRAPH_FREE(VECTOR(*v)[i]);
        }
    }
    igraph_vector_ptr_destroy(v);
}

/**
 * \function igraph_get_all_shortest_paths
 * \brief All shortest paths (geodesics) from a vertex.
 *
 * When there is more than one shortest path between two vertices,
 * all of them will be returned.
 *
 * \param graph The graph object.
 * \param res Pointer to an initialized pointer vector, the result
 *   will be stored here in \ref igraph_vector_t objects. Each vector
 *   object contains the vertices along a shortest path from \p from
 *   to another vertex. The vectors are ordered according to their
 *   target vertex: first the shortest paths to vertex 0, then to
 *   vertex 1, etc. No data is included for unreachable vertices.
 * \param nrgeo Pointer to an initialized \ref igraph_vector_t object or
 *   \c NULL. If not \c NULL the number of shortest paths from \p from are
 *   stored here for every vertex in the graph. Note that the values
 *   will be accurate only for those vertices that are in the target
 *   vertex sequence (see \p to), since the search terminates as soon
 *   as all the target vertices have been found.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the ids of the vertices to/from which the
 *        shortest paths will be calculated. A vertex might be given multiple
 *        times.
 * \param mode The type of shortest paths to be use for the
 *        calculation in directed graphs. Possible values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the lengths of the outgoing paths are calculated.
 *        \cli IGRAPH_IN
 *          the lengths of the incoming paths are calculated.
 *        \cli IGRAPH_ALL
 *          the directed graph is considered as an
 *          undirected one for the computation.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *           not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p from is invalid vertex id.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|+|E|) for most graphs, O(|V|^2) in the worst
 * case.
 */
int igraph_get_all_shortest_paths(const igraph_t *graph,
                                  igraph_vector_ptr_t *res,
                                  igraph_vector_t *nrgeo,
                                  igraph_integer_t from, const igraph_vs_t to,
                                  igraph_neimode_t mode) {

    long int no_of_nodes = igraph_vcount(graph);
    long int *geodist;
    igraph_vector_ptr_t paths;
    igraph_dqueue_t q;
    igraph_vector_t *vptr;
    igraph_vector_t neis;
    igraph_vector_t ptrlist;
    igraph_vector_t ptrhead;
    long int n, j, i;
    long int to_reach, reached = 0, maxdist = 0;

    igraph_vit_t vit;

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("cannot get shortest paths", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    /* paths will store the shortest paths during the search */
    IGRAPH_CHECK(igraph_vector_ptr_init(&paths, 0));
    IGRAPH_FINALLY(igraph_i_gasp_paths_destroy, &paths);
    /* neis is a temporary vector holding the neighbors of the
     * node being examined */
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    /* ptrlist stores indices into the paths vector, in the order
     * of how they were found. ptrhead is a second-level index that
     * will be used to find paths that terminate in a given vertex */
    IGRAPH_VECTOR_INIT_FINALLY(&ptrlist, 0);
    /* ptrhead contains indices into ptrlist.
     * ptrhead[i] = j means that element #j-1 in ptrlist contains
     * the shortest path from the root to node i. ptrhead[i] = 0
     * means that node i was not reached so far */
    IGRAPH_VECTOR_INIT_FINALLY(&ptrhead, no_of_nodes);
    /* geodist[i] == 0 if i was not reached yet and it is not in the
     * target vertex sequence, or -1 if i was not reached yet and it
     * is in the target vertex sequence. Otherwise it is
     * one larger than the length of the shortest path from the
     * source */
    geodist = IGRAPH_CALLOC(no_of_nodes, long int);
    if (geodist == 0) {
        IGRAPH_ERROR("Cannot calculate shortest paths", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, geodist);
    /* dequeue to store the BFS queue -- odd elements are the vertex indices,
     * even elements are the distances from the root */
    IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

    if (nrgeo) {
        IGRAPH_CHECK(igraph_vector_resize(nrgeo, no_of_nodes));
        igraph_vector_null(nrgeo);
    }

    /* use geodist to count how many vertices we have to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (geodist[ (long int) IGRAPH_VIT_GET(vit) ] == 0) {
            geodist[ (long int) IGRAPH_VIT_GET(vit) ] = -1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    if (geodist[ (long int) from ] < 0) {
        reached++;
    }

    /* from -> from */
    vptr = IGRAPH_CALLOC(1, igraph_vector_t); /* TODO: dirty */
    IGRAPH_CHECK(igraph_vector_ptr_push_back(&paths, vptr));
    IGRAPH_CHECK(igraph_vector_init(vptr, 1));
    VECTOR(*vptr)[0] = from;
    geodist[(long int)from] = 1;
    VECTOR(ptrhead)[(long int)from] = 1;
    IGRAPH_CHECK(igraph_vector_push_back(&ptrlist, 0));
    if (nrgeo) {
        VECTOR(*nrgeo)[(long int)from] = 1;
    }

    /* Init queue */
    IGRAPH_CHECK(igraph_dqueue_push(&q, from));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0.0));
    while (!igraph_dqueue_empty(&q)) {
        long int actnode = (long int) igraph_dqueue_pop(&q);
        long int actdist = (long int) igraph_dqueue_pop(&q);

        IGRAPH_ALLOW_INTERRUPTION();

        if (reached >= to_reach) {
            /* all nodes were reached. Since we need all the shortest paths
             * to all these nodes, we can stop the search only if the distance
             * of the current node to the root is larger than the distance of
             * any of the nodes we wanted to reach */
            if (actdist > maxdist) {
                /* safety check, maxdist should have been set when we reached the last node */
                if (maxdist < 0) {
                    IGRAPH_ERROR("possible bug in igraph_get_all_shortest_paths, "
                                 "maxdist is negative", IGRAPH_EINVAL);
                }
                break;
            }
        }

        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) actnode,
                                      mode));
        n = igraph_vector_size(&neis);
        for (j = 0; j < n; j++) {
            long int neighbor = (long int) VECTOR(neis)[j];
            long int fatherptr;

            if (geodist[neighbor] > 0 &&
                geodist[neighbor] - 1 < actdist + 1) {
                /* this node was reached via a shorter path before */
                continue;
            }

            /* yay, found another shortest path to neighbor */

            if (nrgeo) {
                /* the number of geodesics leading to neighbor must be
                 * increased by the number of geodesics leading to actnode */
                VECTOR(*nrgeo)[neighbor] += VECTOR(*nrgeo)[actnode];
            }
            if (geodist[neighbor] <= 0) {
                /* this node was not reached yet, push it into the queue */
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                if (geodist[neighbor] < 0) {
                    reached++;
                }
                if (reached == to_reach) {
                    maxdist = actdist;
                }
            }
            geodist[neighbor] = actdist + 2;

            /* copy all existing paths to the parent */
            fatherptr = (long int) VECTOR(ptrhead)[actnode];
            while (fatherptr != 0) {
                /* allocate a new igraph_vector_t at the end of paths */
                vptr = IGRAPH_CALLOC(1, igraph_vector_t);
                IGRAPH_CHECK(igraph_vector_ptr_push_back(&paths, vptr));
                IGRAPH_CHECK(igraph_vector_copy(vptr, VECTOR(paths)[fatherptr - 1]));
                IGRAPH_CHECK(igraph_vector_reserve(vptr, actdist + 2));
                IGRAPH_CHECK(igraph_vector_push_back(vptr, neighbor));

                IGRAPH_CHECK(igraph_vector_push_back(&ptrlist,
                                                     VECTOR(ptrhead)[neighbor]));
                VECTOR(ptrhead)[neighbor] = igraph_vector_size(&ptrlist);

                fatherptr = (long int) VECTOR(ptrlist)[fatherptr - 1];
            }
        }
    }

    igraph_dqueue_destroy(&q);
    IGRAPH_FINALLY_CLEAN(1);

    /* mark the nodes for which we need the result */
    memset(geodist, 0, sizeof(long int) * (size_t) no_of_nodes);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        geodist[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
    }

    /* count the number of paths in the result */
    n = 0;
    for (i = 0; i < no_of_nodes; i++) {
        long int fatherptr = (long int) VECTOR(ptrhead)[i];
        if (geodist[i] > 0) {
            while (fatherptr != 0) {
                n++;
                fatherptr = (long int) VECTOR(ptrlist)[fatherptr - 1];
            }
        }
    }

    IGRAPH_CHECK(igraph_vector_ptr_resize(res, n));
    j = 0;
    for (i = 0; i < no_of_nodes; i++) {
        long int fatherptr = (long int) VECTOR(ptrhead)[i];

        IGRAPH_ALLOW_INTERRUPTION();

        /* do we need the paths leading to vertex i? */
        if (geodist[i] > 0) {
            /* yes, copy them to the result vector */
            while (fatherptr != 0) {
                VECTOR(*res)[j++] = VECTOR(paths)[fatherptr - 1];
                fatherptr = (long int) VECTOR(ptrlist)[fatherptr - 1];
            }
        } else {
            /* no, free them */
            while (fatherptr != 0) {
                igraph_vector_destroy(VECTOR(paths)[fatherptr - 1]);
                IGRAPH_FREE(VECTOR(paths)[fatherptr - 1]);
                fatherptr = (long int) VECTOR(ptrlist)[fatherptr - 1];
            }
        }
    }

    IGRAPH_FREE(geodist);
    igraph_vector_destroy(&ptrlist);
    igraph_vector_destroy(&ptrhead);
    igraph_vector_destroy(&neis);
    igraph_vector_ptr_destroy(&paths);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(6);

    return 0;
}
