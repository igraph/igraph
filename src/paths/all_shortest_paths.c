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

#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"

#include "core/interruption.h"

#include <string.h>  /* memset */

/**
 * \function igraph_get_all_shortest_paths
 * \brief All shortest paths (geodesics) from a vertex.
 *
 * When there is more than one shortest path between two vertices,
 * all of them will be returned. Every edge is considered separately,
 * therefore in graphs with multi-edges, this function may produce
 * a very large number of results.
 *
 * \param graph The graph object.
 * \param vertices The result, the IDs of the vertices along the paths.
 *   This is a list of integer vectors where each element is an
 *   \ref igraph_vector_int_t object. Each vector object contains the vertices
 *   along a shortest path from \p from to another vertex. The vectors are
 *   ordered according to their target vertex: first the shortest paths to
 *   vertex 0, then to vertex 1, etc. No data is included for unreachable
 *   vertices. The list will be resized as needed. Supply a null pointer here
 *   if you don't need these vectors.
 * \param edges The result, the IDs of the edges along the paths.
 *   This is a list of integer vectors where each element is an
 *   \ref igraph_vector_int_t object. Each vector object contains the edges
 *   along a shortest path from \p from to another vertex. The vectors are
 *   ordered according to their target vertex: first the shortest paths to
 *   vertex 0, then to vertex 1, etc. No data is included for unreachable
 *   vertices. The list will be resized as needed. Supply a null pointer here
 *   if you don't need these vectors.
 * \param nrgeo Pointer to an initialized \ref igraph_vector_int_t object or
 *   \c NULL. If not \c NULL the number of shortest paths from \p from are
 *   stored here for every vertex in the graph. Note that the values
 *   will be accurate only for those vertices that are in the target
 *   vertex sequence (see \p to), since the search terminates as soon
 *   as all the target vertices have been found.
 * \param from The id of the vertex from/to which the geodesics are
 *        calculated.
 * \param to Vertex sequence with the IDs of the vertices to/from which the
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
 *           \p from is invalid vertex ID.
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument.
 *        \endclist
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|+|E|) for most graphs, O(|V|^2) in the worst
 * case.
 */

igraph_error_t igraph_get_all_shortest_paths(const igraph_t *graph,
                                  igraph_vector_int_list_t *vertices,
                                  igraph_vector_int_list_t *edges,
                                  igraph_vector_int_t *nrgeo,
                                  igraph_integer_t from, const igraph_vs_t to,
                                  igraph_neimode_t mode) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t *geodist;
    igraph_vector_int_list_t paths;
    igraph_vector_int_list_t path_edge;
    igraph_dqueue_int_t q;
    igraph_vector_int_t *vptr;
    igraph_vector_int_t *vptr_e;
    igraph_vector_int_t neis;
    igraph_vector_int_t ptrlist;
    igraph_vector_int_t ptrhead;
    igraph_integer_t n;
    igraph_integer_t to_reach, reached = 0, maxdist = 0;

    igraph_vit_t vit;

    if (from < 0 || from >= no_of_nodes) {
        IGRAPH_ERROR("Index of source vertex is out of range.", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    /* paths will store the shortest paths during the search */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&paths, 0);
    /* path_edge will store the shortest paths during the search */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&path_edge, 0);
    /* neis is a temporary vector holding the neighbors of the
     * node being examined */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    /* ptrlist stores indices into the paths vector, in the order
     * of how they were found. ptrhead is a second-level index that
     * will be used to find paths that terminate in a given vertex */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptrlist, 0);
    /* ptrhead contains indices into ptrlist.
     * ptrhead[i] = j means that element #j-1 in ptrlist contains
     * the shortest path from the root to node i. ptrhead[i] = 0
     * means that node i was not reached so far */
    IGRAPH_VECTOR_INT_INIT_FINALLY(&ptrhead, no_of_nodes);
    /* geodist[i] == 0 if i was not reached yet and it is not in the
     * target vertex sequence, or -1 if i was not reached yet and it
     * is in the target vertex sequence. Otherwise it is
     * one larger than the length of the shortest path from the
     * source */
    geodist = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    IGRAPH_CHECK_OOM(geodist, "Insufficient memory for calculating shortest paths.");
    IGRAPH_FINALLY(igraph_free, geodist);
    /* dequeue to store the BFS queue -- odd elements are the vertex indices,
     * even elements are the distances from the root */
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    if (nrgeo) {
        IGRAPH_CHECK(igraph_vector_int_resize(nrgeo, no_of_nodes));
        igraph_vector_int_null(nrgeo);
    }

    /* use geodist to count how many vertices we have to reach */
    to_reach = IGRAPH_VIT_SIZE(vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        if (geodist[ IGRAPH_VIT_GET(vit) ] == 0) {
            geodist[ IGRAPH_VIT_GET(vit) ] = -1;
        } else {
            to_reach--;       /* this node was given multiple times */
        }
    }

    if (geodist[ from ] < 0) {
        reached++;
    }

    /* from -> from */
    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(&paths, &vptr));
    IGRAPH_CHECK(igraph_vector_int_push_back(vptr, from));
    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(&path_edge, &vptr_e));

    geodist[from] = 1;
    VECTOR(ptrhead)[from] = 1;
    IGRAPH_CHECK(igraph_vector_int_push_back(&ptrlist, 0));
    if (nrgeo) {
        VECTOR(*nrgeo)[from] = 1;
    }

    /* Init queue */
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, from));
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
    while (!igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
        igraph_integer_t actdist = igraph_dqueue_int_pop(&q);

        IGRAPH_ALLOW_INTERRUPTION();

        if (reached >= to_reach) {
            /* all nodes were reached. Since we need all the shortest paths
             * to all these nodes, we can stop the search only if the distance
             * of the current node to the root is larger than the distance of
             * any of the nodes we wanted to reach */
            if (actdist > maxdist) {
                /* safety check, maxdist should have been set when we reached the last node */
                IGRAPH_ASSERT(maxdist >= 0);
                break;
            }
        }

        /* If we need the edge-paths, we need to use igraph_incident() followed by an
         * IGRAPH_OTHER() macro in the main loop. This is going to be slower than
         * using igraph_neighbors() due to branch mispredictions in IGRAPH_OTHER(), so we
         * use igraph_incident() only if the user needs the edge-paths */
        if (edges) {
            IGRAPH_CHECK(igraph_incident(graph, &neis, actnode, mode));
        } else {
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
        }

        n = igraph_vector_int_size(&neis);
        for (igraph_integer_t j = 0; j < n; j++) {
            igraph_integer_t neighbor;
            igraph_integer_t parentptr;

            if (edges) {
                /* user needs the edge-paths, so 'neis' contains edge IDs, we need to resolve
                 * the next edge ID into a vertex ID */
                neighbor = IGRAPH_OTHER(graph, VECTOR(neis)[j], actnode);
            } else {
                /* user does not need the edge-paths, so 'neis' contains vertex IDs */
                neighbor = VECTOR(neis)[j];
            }

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
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                if (geodist[neighbor] < 0) {
                    reached++;
                }
                if (reached == to_reach) {
                    maxdist = actdist;
                }
            }
            geodist[neighbor] = actdist + 2;

            /* copy all existing paths to the parent */
            parentptr = VECTOR(ptrhead)[actnode];
            while (parentptr != 0) {
                /* allocate a new igraph_vector_int_t at the end of paths */
                IGRAPH_CHECK(igraph_vector_int_list_push_back_new(&paths, &vptr));
                IGRAPH_CHECK(igraph_vector_int_update(vptr, igraph_vector_int_list_get_ptr(&paths, parentptr - 1)));
                IGRAPH_CHECK(igraph_vector_int_push_back(vptr, neighbor));

                IGRAPH_CHECK(igraph_vector_int_list_push_back_new(&path_edge, &vptr_e));
                if (actnode != from) {
                    /* If the previous vertex was the source then there is no edge to add*/
                    IGRAPH_CHECK(igraph_vector_int_update(vptr_e, igraph_vector_int_list_get_ptr(&path_edge, parentptr - 1)));
                }
                IGRAPH_CHECK(igraph_vector_int_push_back(vptr_e, VECTOR(neis)[j]));

                IGRAPH_CHECK(igraph_vector_int_push_back(&ptrlist, VECTOR(ptrhead)[neighbor]));
                VECTOR(ptrhead)[neighbor] = igraph_vector_int_size(&ptrlist);

                parentptr = VECTOR(ptrlist)[parentptr - 1];
            }
        }
    }

    igraph_dqueue_int_destroy(&q);
    IGRAPH_FINALLY_CLEAN(1);

    /* mark the nodes for which we need the result */
    memset(geodist, 0, sizeof(geodist[0]) * (size_t) no_of_nodes);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        geodist[ IGRAPH_VIT_GET(vit) ] = 1;
    }

    if (vertices) {
        igraph_vector_int_list_clear(vertices);
    }
    if (edges) {
        igraph_vector_int_list_clear(edges);
    }

    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t parentptr = VECTOR(ptrhead)[i];

        IGRAPH_ALLOW_INTERRUPTION();

        /* do we need the paths leading to vertex i? */
        if (geodist[i] > 0) {
            /* yes, transfer them to the result vector */
            while (parentptr != 0) {
                /* Given two vector lists, list1 and list2, an efficient way to transfer
                 * a vector from list1 to the end of list2 is to extend list2 with an
                 * empty vector, then swap that empty vector with the given element of
                 * list1. This approach avoids creating a full copy of the vector. */
                if (vertices) {
                    igraph_vector_int_t *p;
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(vertices, &p));
                    igraph_vector_int_swap(p, igraph_vector_int_list_get_ptr(&paths, parentptr - 1));
                }
                if (edges) {
                    igraph_vector_int_t *p;
                    IGRAPH_CHECK(igraph_vector_int_list_push_back_new(edges, &p));
                    igraph_vector_int_swap(p, igraph_vector_int_list_get_ptr(&path_edge, parentptr - 1));
                }
                parentptr = VECTOR(ptrlist)[parentptr - 1];
            }
        }
    }

    IGRAPH_FREE(geodist);
    igraph_vector_int_destroy(&ptrlist);
    igraph_vector_int_destroy(&ptrhead);
    igraph_vector_int_destroy(&neis);
    igraph_vector_int_list_destroy(&paths);
    igraph_vector_int_list_destroy(&path_edge);
    igraph_vit_destroy(&vit);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}
