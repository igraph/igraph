/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "igraph_components.h"

#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_operators.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_structural.h"
#include "igraph_vector.h"

#include "core/interruption.h"
#include "operators/subgraph.h"

#include <limits.h>

static igraph_error_t igraph_i_clusters_weak(const igraph_t *graph, igraph_vector_int_t *membership,
                                  igraph_vector_int_t *csize, igraph_integer_t *no);

static igraph_error_t igraph_i_clusters_strong(const igraph_t *graph, igraph_vector_int_t *membership,
                                    igraph_vector_int_t *csize, igraph_integer_t *no);

/**
 * \ingroup structural
 * \function igraph_clusters
 * \brief Calculates the (weakly or strongly) connected components in a graph.
 *
 * \param graph The graph object to analyze.
 * \param membership First half of the result will be stored here. For
 *        every vertex the id of its component is given. The vector
 *        has to be preinitialized and will be resized. Alternatively
 *        this argument can be \c NULL, in which case it is ignored.
 * \param csize The second half of the result. For every component it
 *        gives its size, the order is defined by the component ids.
 *        The vector has to be preinitialized and will be resized.
 *        Alternatively this argument can be \c NULL, in which
 *        case it is ignored.
 * \param no Pointer to an integer, if not \c NULL then the number of
 *        clusters will be stored here.
 * \param mode For directed graph this specifies whether to calculate
 *        weakly or strongly connected components. Possible values:
 *        \c IGRAPH_WEAK,
 *        \c IGRAPH_STRONG. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *         \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 */

igraph_error_t igraph_clusters(const igraph_t *graph, igraph_vector_int_t *membership,
                    igraph_vector_int_t *csize, igraph_integer_t *no,
                    igraph_connectedness_t mode) {
    if (mode == IGRAPH_WEAK || !igraph_is_directed(graph)) {
        return igraph_i_clusters_weak(graph, membership, csize, no);
    } else if (mode == IGRAPH_STRONG) {
        return igraph_i_clusters_strong(graph, membership, csize, no);
    }

    IGRAPH_ERROR("Cannot calculate clusters", IGRAPH_EINVAL);
}

static igraph_error_t igraph_i_clusters_weak(const igraph_t *graph, igraph_vector_int_t *membership,
                                  igraph_vector_int_t *csize, igraph_integer_t *no) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    char *already_added;
    igraph_integer_t first_node, act_cluster_size = 0, no_of_clusters = 1;

    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

    igraph_integer_t i;
    igraph_vector_int_t neis = IGRAPH_VECTOR_NULL;

    already_added = IGRAPH_CALLOC(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("Cannot calculate clusters", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, no_of_nodes > 100000 ? 10000 : no_of_nodes / 10);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* Memory for result, csize is dynamically allocated */
    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
    }
    if (csize) {
        igraph_vector_int_clear(csize);
    }

    /* The algorithm */

    for (first_node = 0; first_node < no_of_nodes; ++first_node) {
        if (already_added[first_node] == 1) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        already_added[first_node] = 1;
        act_cluster_size = 1;
        if (membership) {
            VECTOR(*membership)[first_node] = no_of_clusters - 1;
        }
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, first_node));

        while ( !igraph_dqueue_int_empty(&q) ) {
            igraph_integer_t act_node = igraph_dqueue_int_pop(&q);
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, act_node, IGRAPH_ALL));
            for (i = 0; i < igraph_vector_int_size(&neis); i++) {
                igraph_integer_t neighbor = VECTOR(neis)[i];
                if (already_added[neighbor] == 1) {
                    continue;
                }
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                already_added[neighbor] = 1;
                act_cluster_size++;
                if (membership) {
                    VECTOR(*membership)[neighbor] = no_of_clusters - 1;
                }
            }
        }
        no_of_clusters++;
        if (csize) {
            IGRAPH_CHECK(igraph_vector_int_push_back(csize, act_cluster_size));
        }
    }

    /* Cleaning up */

    if (no) {
        *no = no_of_clusters - 1;
    }

    IGRAPH_FREE(already_added);
    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_clusters_strong(const igraph_t *graph, igraph_vector_int_t *membership,
                                    igraph_vector_int_t *csize, igraph_integer_t *no) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t next_nei = IGRAPH_VECTOR_NULL;

    igraph_integer_t i, n, num_seen;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

    igraph_integer_t no_of_clusters = 1;
    igraph_integer_t act_cluster_size;

    igraph_vector_int_t out = IGRAPH_VECTOR_NULL;
    const igraph_vector_int_t* tmp;

    igraph_adjlist_t adjlist;

    /* The result */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&next_nei, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&out, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    if (membership) {
        IGRAPH_CHECK(igraph_vector_int_resize(membership, no_of_nodes));
    }
    IGRAPH_CHECK(igraph_vector_int_reserve(&out, no_of_nodes));

    igraph_vector_int_null(&out);
    if (csize) {
        igraph_vector_int_clear(csize);
    }

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    num_seen = 0;
    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();

        tmp = igraph_adjlist_get(&adjlist, i);
        if (VECTOR(next_nei)[i] > igraph_vector_int_size(tmp)) {
            continue;
        }

        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));
        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t act_node = igraph_dqueue_int_back(&q);
            tmp = igraph_adjlist_get(&adjlist, act_node);
            if (VECTOR(next_nei)[act_node] == 0) {
                /* this is the first time we've met this vertex */
                VECTOR(next_nei)[act_node]++;
            } else if (VECTOR(next_nei)[act_node] <= igraph_vector_int_size(tmp)) {
                /* we've already met this vertex but it has more children */
                igraph_integer_t neighbor = VECTOR(*tmp)[VECTOR(next_nei)[act_node] - 1];
                if (VECTOR(next_nei)[neighbor] == 0) {
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                }
                VECTOR(next_nei)[act_node]++;
            } else {
                /* we've met this vertex and it has no more children */
                IGRAPH_CHECK(igraph_vector_int_push_back(&out, act_node));
                igraph_dqueue_int_pop_back(&q);
                num_seen++;

                if (num_seen % 10000 == 0) {
                    /* time to report progress and allow the user to interrupt */
                    IGRAPH_PROGRESS("Strongly connected components: ",
                                    num_seen * 50.0 / no_of_nodes, NULL);
                    IGRAPH_ALLOW_INTERRUPTION();
                }
            }
        } /* while q */
    }  /* for */

    IGRAPH_PROGRESS("Strongly connected components: ", 50.0, NULL);

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* OK, we've the 'out' values for the nodes, let's use them in
       decreasing order with the help of a heap */

    igraph_vector_int_null(&next_nei);             /* mark already added vertices */
    num_seen = 0;

    while (!igraph_vector_int_empty(&out)) {
        igraph_integer_t grandfather = igraph_vector_int_pop_back(&out);

        if (VECTOR(next_nei)[grandfather] != 0) {
            continue;
        }
        VECTOR(next_nei)[grandfather] = 1;
        act_cluster_size = 1;
        if (membership) {
            VECTOR(*membership)[grandfather] = no_of_clusters - 1;
        }
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, grandfather));

        num_seen++;
        if (num_seen % 10000 == 0) {
            /* time to report progress and allow the user to interrupt */
            IGRAPH_PROGRESS("Strongly connected components: ",
                            50.0 + num_seen * 50.0 / no_of_nodes, NULL);
            IGRAPH_ALLOW_INTERRUPTION();
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t act_node = igraph_dqueue_int_pop_back(&q);
            tmp = igraph_adjlist_get(&adjlist, act_node);
            n = igraph_vector_int_size(tmp);
            for (i = 0; i < n; i++) {
                igraph_integer_t neighbor = VECTOR(*tmp)[i];
                if (VECTOR(next_nei)[neighbor] != 0) {
                    continue;
                }
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                VECTOR(next_nei)[neighbor] = 1;
                act_cluster_size++;
                if (membership) {
                    VECTOR(*membership)[neighbor] = no_of_clusters - 1;
                }

                num_seen++;
                if (num_seen % 10000 == 0) {
                    /* time to report progress and allow the user to interrupt */
                    IGRAPH_PROGRESS("Strongly connected components: ",
                                    50.0 + num_seen * 50.0 / no_of_nodes, NULL);
                    IGRAPH_ALLOW_INTERRUPTION();
                }
            }
        }

        no_of_clusters++;
        if (csize) {
            IGRAPH_CHECK(igraph_vector_int_push_back(csize, act_cluster_size));
        }
    }

    IGRAPH_PROGRESS("Strongly connected components: ", 100.0, NULL);

    if (no) {
        *no = no_of_clusters - 1;
    }

    /* Clean up, return */

    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&out);
    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&next_nei);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_is_connected_weak(const igraph_t *graph, igraph_bool_t *res);

/**
 * \ingroup structural
 * \function igraph_is_connected
 * \brief Decides whether the graph is (weakly or strongly) connected.
 *
 * A graph with zero vertices (i.e. the null graph) is \em not connected by
 * definition. This behaviour changed in igraph 0.9; earlier versions assumed
 * that the null graph is connected. See the following issue on Github for the
 * argument that led us to change the definition:
 * https://github.com/igraph/igraph/issues/1538
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable, the result will be stored
 *        here.
 * \param mode For a directed graph this specifies whether to calculate
 *        weak or strong connectedness. Possible values:
 *        \c IGRAPH_WEAK,
 *        \c IGRAPH_STRONG. This argument is
 *        ignored for undirected graphs.
 * \return Error code:
 *        \c IGRAPH_EINVAL: invalid mode argument.
 *
 * Time complexity: O(|V|+|E|), the
 * number of vertices
 * plus the number of edges in the graph.
 */

igraph_error_t igraph_is_connected(const igraph_t *graph, igraph_bool_t *res,
                        igraph_connectedness_t mode) {
    if (igraph_vcount(graph) == 0) {
        /* Changed in igraph 0.9; see https://github.com/igraph/igraph/issues/1538
         * for the reasoning behind the change */
        *res = 0;
        return IGRAPH_SUCCESS;
    }

    if (mode == IGRAPH_WEAK || !igraph_is_directed(graph)) {
        return igraph_is_connected_weak(graph, res);
    } else if (mode == IGRAPH_STRONG) {
        igraph_error_t retval;
        igraph_integer_t no;
        retval = igraph_i_clusters_strong(graph, 0, 0, &no);
        *res = (no == 1);
        return retval;
    }

    IGRAPH_ERROR("mode argument", IGRAPH_EINVAL);
}

/**
 * \ingroup structural
 * \function igraph_is_connected_weak
 * \brief Query whether the graph is weakly connected.
 *
 * A graph with zero vertices (i.e. the null graph) is weakly connected by
 * definition. A directed graph is weakly connected if its undirected version
 * is connected. In the case of undirected graphs, weakly connected and
 * connected are equivalent.
 *
 * \param graph The graph object to analyze.
 * \param res Pointer to a logical variable; the result will be stored here.
 * \return Error code:
 *        \c IGRAPH_ENOMEM: unable to allocate requested memory.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of
 * edges in the graph.
 */

igraph_error_t igraph_is_connected_weak(const igraph_t *graph, igraph_bool_t *res) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    char *already_added;
    igraph_vector_int_t neis = IGRAPH_VECTOR_NULL;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

    igraph_integer_t i, j;

    if (no_of_nodes == 0) {
        *res = 1;
        return IGRAPH_SUCCESS;
    }

    already_added = IGRAPH_CALLOC(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("is connected (weak) failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 10);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    /* Try to find at least two clusters */
    already_added[0] = 1;
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));

    j = 1;
    while ( !igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, IGRAPH_ALL));
        for (i = 0; i < igraph_vector_int_size(&neis); i++) {
            igraph_integer_t neighbor = VECTOR(neis)[i];
            if (already_added[neighbor] != 0) {
                continue;
            }
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
            j++;
            already_added[neighbor]++;
        }
    }

    /* Connected? */
    *res = (j == no_of_nodes);

    IGRAPH_FREE(already_added);
    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&neis);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_decompose_destroy
 * \brief Free the memory allocated by \ref igraph_decompose().
 *
 * \param complist The list of graph components, as returned by
 *        \ref igraph_decompose().
 *
 * Time complexity: O(c), c is the number of components.
 */

void igraph_decompose_destroy(igraph_vector_ptr_t *complist) {
    igraph_integer_t i, n;

    n = igraph_vector_ptr_size(complist);
    for (i = 0; i < n; i++) {
        if (VECTOR(*complist)[i] != 0) {
            igraph_destroy(VECTOR(*complist)[i]);
            igraph_free(VECTOR(*complist)[i]);
        }
    }
}

static igraph_error_t igraph_i_decompose_weak(const igraph_t *graph,
                                   igraph_vector_ptr_t *components,
                                   igraph_integer_t maxcompno, igraph_integer_t minelements);

static igraph_error_t igraph_i_decompose_strong(const igraph_t *graph,
                                     igraph_vector_ptr_t *components,
                                     igraph_integer_t maxcompno, igraph_integer_t minelements);

/**
 * \function igraph_decompose
 * \brief Decompose a graph into connected components.
 *
 * Create separate graph for each component of a graph. Note that the
 * vertex IDs in the new graphs will be different than in the original
 * graph. (Except if there is only one component in the original graph.)
 *
 * \param graph The original graph.
 * \param components This pointer vector will contain pointers to the
 *   subcomponent graphs. It should be initialized before calling this
 *   function and will be resized to hold the graphs. Don't forget to
 *   call \ref igraph_destroy() and \ref igraph_free() on the elements of
 *   this pointer vector to free unneeded memory. Alternatively, you can
 *   simply call \ref igraph_decompose_destroy() that does this for you.
 * \param mode Either \c IGRAPH_WEAK or \c IGRAPH_STRONG for weakly
 *    and strongly connected components respectively.
 * \param maxcompno The maximum number of components to return. The
 *    first \p maxcompno components will be returned (which hold at
 *    least \p minelements vertices, see the next parameter), the
 *    others will be ignored. Supply -1 here if you don't want to limit
 *    the number of components.
 * \param minelements The minimum number of vertices a component
 *    should contain in order to place it in the \p components
 *    vector. Eg. supply 2 here to ignore isolated vertices.
 * \return Error code, \c IGRAPH_ENOMEM if there is not enough memory
 *   to perform the operation.
 *
 * Added in version 0.2.</para><para>
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges.
 *
 * \example examples/simple/igraph_decompose.c
 */

igraph_error_t igraph_decompose(const igraph_t *graph, igraph_vector_ptr_t *components,
                     igraph_connectedness_t mode,
                     igraph_integer_t maxcompno, igraph_integer_t minelements) {
    if (mode == IGRAPH_WEAK || !igraph_is_directed(graph)) {
        return igraph_i_decompose_weak(graph, components, maxcompno, minelements);
    } else if (mode == IGRAPH_STRONG) {
        return igraph_i_decompose_strong(graph, components, maxcompno, minelements);
    }

    IGRAPH_ERROR("Cannot decompose graph", IGRAPH_EINVAL);
}

static igraph_error_t igraph_i_decompose_weak(const igraph_t *graph,
                                   igraph_vector_ptr_t *components,
                                   igraph_integer_t maxcompno, igraph_integer_t minelements) {

    igraph_integer_t actstart;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t resco = 0;   /* number of graphs created so far */
    char *already_added;
    igraph_dqueue_int_t q;
    igraph_vector_int_t verts;
    igraph_vector_int_t neis;
    igraph_vector_int_t vids_old2new;
    igraph_integer_t i;
    igraph_t *newg;


    if (maxcompno < 0) {
        maxcompno = IGRAPH_INTEGER_MAX;
    }

    igraph_vector_ptr_clear(components);
    IGRAPH_FINALLY(igraph_decompose_destroy, components);

    /* already_added keeps track of what nodes made it into a graph already */
    already_added = IGRAPH_CALLOC(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("Cannot decompose graph", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&verts, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids_old2new, no_of_nodes);

    /* vids_old2new would have been created internally in igraph_induced_subgraph(),
       but it is slow if the graph is large and consists of many small components,
       so we create it once here and then re-use it */

    /* add a node and its neighbors at once, recursively
       then switch to next node that has not been added already */
    for (actstart = 0; resco < maxcompno && actstart < no_of_nodes; actstart++) {

        if (already_added[actstart]) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        igraph_vector_int_clear(&verts);

        /* add the node itself */
        already_added[actstart] = 1;
        IGRAPH_CHECK(igraph_vector_int_push_back(&verts, actstart));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actstart));

        /* add the neighbors, recursively */
        while (!igraph_dqueue_int_empty(&q) ) {
            /* pop from the queue of this component */
            igraph_integer_t actvert = igraph_dqueue_int_pop(&q);
            IGRAPH_CHECK(igraph_neighbors(graph, &neis, actvert, IGRAPH_ALL));
            /* iterate over the neighbors */
            for (i = 0; i < igraph_vector_int_size(&neis); i++) {
                igraph_integer_t neighbor = VECTOR(neis)[i];
                if (already_added[neighbor] == 1) {
                    continue;
                }
                /* add neighbor */
                already_added[neighbor] = 1;

                /* recursion: append neighbor to the queues */
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_vector_int_push_back(&verts, neighbor));
            }
        }

        /* ok, we have a component */
        if (igraph_vector_int_size(&verts) < minelements) {
            continue;
        }

        newg = IGRAPH_CALLOC(1, igraph_t);
        if (newg == 0) {
            IGRAPH_ERROR("Cannot decompose graph", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_ptr_push_back(components, newg));
        IGRAPH_CHECK(igraph_i_induced_subgraph_map(
            graph, newg, igraph_vss_vector(&verts),
            IGRAPH_SUBGRAPH_AUTO, &vids_old2new,
            /* invmap = */ 0, /* map_is_prepared = */ 1
        ));
        resco++;

        /* vids_old2new does not have to be cleaned up here; since we are doing
         * weak decomposition, each vertex will appear in only one of the
         * connected components so we won't ever touch an item in vids_old2new
         * if it was already set to a non-zero value in a previous component */

    } /* for actstart++ */

    igraph_vector_int_destroy(&vids_old2new);
    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&verts);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(already_added);
    IGRAPH_FINALLY_CLEAN(6);  /* + components */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_decompose_strong(const igraph_t *graph,
                                     igraph_vector_ptr_t *components,
                                     igraph_integer_t maxcompno, igraph_integer_t minelements) {


    igraph_integer_t no_of_nodes = igraph_vcount(graph);

    /* this is a heap used twice for checking what nodes have
     * been counted already */
    igraph_vector_int_t next_nei = IGRAPH_VECTOR_NULL;

    igraph_integer_t i, n, num_seen;
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;

    igraph_integer_t no_of_clusters = 0;
    igraph_integer_t act_cluster_size;

    igraph_vector_int_t out = IGRAPH_VECTOR_NULL;
    const igraph_vector_int_t* tmp;

    igraph_adjlist_t adjlist;
    igraph_vector_int_t verts;
    igraph_vector_int_t vids_old2new;
    igraph_t *newg;

    if (maxcompno < 0) {
        maxcompno = IGRAPH_INTEGER_MAX;
    }

    igraph_vector_ptr_clear(components);
    IGRAPH_FINALLY(igraph_decompose_destroy, components);

    /* The result */

    IGRAPH_VECTOR_INT_INIT_FINALLY(&vids_old2new, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&verts, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&next_nei, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&out, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_vector_int_reserve(&out, no_of_nodes));

    igraph_vector_int_null(&out);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* vids_old2new would have been created internally in igraph_induced_subgraph(),
       but it is slow if the graph is large and consists of many small components,
       so we create it once here and then re-use it */

    /* number of components seen */
    num_seen = 0;
    /* populate the 'out' vector by browsing a node and following up
       all its neighbors recursively, then switching to the next
       unassigned node */
    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* get all the 'out' neighbors of this node
         * NOTE: next_nei is initialized [0, 0, ...] */
        tmp = igraph_adjlist_get(&adjlist, i);
        if (VECTOR(next_nei)[i] > igraph_vector_int_size(tmp)) {
            continue;
        }

        /* add this node to the queue for this component */
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, i));

        /* consume the tree from this node ("root") recursively
         * until there is no more */
        while (!igraph_dqueue_int_empty(&q)) {
            /* this looks up but does NOT consume the queue */
            igraph_integer_t act_node = igraph_dqueue_int_back(&q);

            /* get all neighbors of this node */
            tmp = igraph_adjlist_get(&adjlist, act_node);
            if (VECTOR(next_nei)[act_node] == 0) {
                /* this is the first time we've met this vertex,
                     * because next_nei is initialized [0, 0, ...] */
                VECTOR(next_nei)[act_node]++;
                /* back to the queue, same vertex is up again */

            } else if (VECTOR(next_nei)[act_node] <= igraph_vector_int_size(tmp)) {
                /* we've already met this vertex but it has more children */
                igraph_integer_t neighbor = VECTOR(*tmp)[VECTOR(next_nei)[act_node] - 1];
                if (VECTOR(next_nei)[neighbor] == 0) {
                    /* add the root of the other children to the queue */
                    IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                }
                VECTOR(next_nei)[act_node]++;
            } else {
                /* we've met this vertex and it has no more children */
                IGRAPH_CHECK(igraph_vector_int_push_back(&out, act_node));
                /* this consumes the queue, since there's nowhere to go */
                igraph_dqueue_int_pop_back(&q);
                num_seen++;

                if (num_seen % 10000 == 0) {
                    /* time to report progress and allow the user to interrupt */
                    IGRAPH_PROGRESS("Strongly connected components: ",
                                    num_seen * 50.0 / no_of_nodes, NULL);
                    IGRAPH_ALLOW_INTERRUPTION();
                }
            }
        } /* while q */
    }  /* for */

    IGRAPH_PROGRESS("Strongly connected components: ", 50.0, NULL);

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    /* OK, we've the 'out' values for the nodes, let's use them in
     * decreasing order with the help of the next_nei heap */

    igraph_vector_int_null(&next_nei);             /* mark already added vertices */

    /* number of components built */
    num_seen = 0;
    while (!igraph_vector_int_empty(&out) && no_of_clusters < maxcompno) {
        /* consume the vector from the last element */
        igraph_integer_t grandfather = igraph_vector_int_pop_back(&out);

        /* been here, done that
         * NOTE: next_nei is initialized as [0, 0, ...] */
        if (VECTOR(next_nei)[grandfather] != 0) {
            continue;
        }

        /* collect all the members of this component */
        igraph_vector_int_clear(&verts);

        /* this node is gone for any future components */
        VECTOR(next_nei)[grandfather] = 1;
        act_cluster_size = 1;

        /* add to component */
        IGRAPH_CHECK(igraph_vector_int_push_back(&verts, grandfather));
        IGRAPH_CHECK(igraph_dqueue_int_push(&q, grandfather));

        num_seen++;
        if (num_seen % 10000 == 0) {
            /* time to report progress and allow the user to interrupt */
            IGRAPH_PROGRESS("Strongly connected components: ",
                            50.0 + num_seen * 50.0 / no_of_nodes, NULL);
            IGRAPH_ALLOW_INTERRUPTION();
        }

        while (!igraph_dqueue_int_empty(&q)) {
            /* consume the queue from this node */
            igraph_integer_t act_node = igraph_dqueue_int_pop_back(&q);
            tmp = igraph_adjlist_get(&adjlist, act_node);
            n = igraph_vector_int_size(tmp);
            for (i = 0; i < n; i++) {
                igraph_integer_t neighbor = VECTOR(*tmp)[i];
                if (VECTOR(next_nei)[neighbor] != 0) {
                    continue;
                }
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                VECTOR(next_nei)[neighbor] = 1;
                act_cluster_size++;

                /* add to component */
                IGRAPH_CHECK(igraph_vector_int_push_back(&verts, neighbor));

                num_seen++;
                if (num_seen % 10000 == 0) {
                    /* time to report progress and allow the user to interrupt */
                    IGRAPH_PROGRESS("Strongly connected components: ",
                                    50.0 + num_seen * 50.0 / no_of_nodes, NULL);
                    IGRAPH_ALLOW_INTERRUPTION();
                }
            }
        }

        /* ok, we have a component */
        if (igraph_vector_int_size(&verts) < minelements) {
            continue;
        }

        newg = IGRAPH_CALLOC(1, igraph_t);
        if (newg == 0) {
            IGRAPH_ERROR("Cannot decompose graph", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_ptr_push_back(components, newg));
        IGRAPH_CHECK(igraph_i_induced_subgraph_map(
            graph, newg, igraph_vss_vector(&verts),
            IGRAPH_SUBGRAPH_AUTO, &vids_old2new,
            /* invmap = */ 0, /* map_is_prepared = */ 1
        ));

        /* vids_old2new has to be cleaned up here because a vertex may appear
         * in multiple strongly connected components. Simply calling
         * igraph_vector_int_fill() would be an O(n) operation where n is the number
         * of vertices in the large graph so we cannot do that; we have to
         * iterate over 'verts' instead */
        n = igraph_vector_int_size(&verts);
        for (i = 0; i < n; i++) {
            VECTOR(vids_old2new)[VECTOR(verts)[i]] = 0;
        }

        no_of_clusters++;
    }

    IGRAPH_PROGRESS("Strongly connected components: ", 100.0, NULL);

    /* Clean up, return */

    igraph_vector_int_destroy(&vids_old2new);
    igraph_vector_int_destroy(&verts);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&out);
    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&next_nei);
    IGRAPH_FINALLY_CLEAN(7);  /* + components */

    return IGRAPH_SUCCESS;

}

/**
 * \function igraph_articulation_points
 * Find the articulation points in a graph.
 *
 * A vertex is an articulation point if its removal increases
 * the number of connected components in the graph.
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the
 *    articulation points will be stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and edges.
 *
 * \sa \ref igraph_biconnected_components(), \ref igraph_clusters(), \ref igraph_bridges()
 */

igraph_error_t igraph_articulation_points(const igraph_t *graph, igraph_vector_int_t *res) {

    return igraph_biconnected_components(graph, 0, 0, 0, 0, res);
}

static void igraph_i_free_vector_int_list(igraph_vector_ptr_t *list) {
    igraph_integer_t i, n = igraph_vector_ptr_size(list);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = VECTOR(*list)[i];
        if (v) {
            igraph_vector_int_destroy(v);
            IGRAPH_FREE(v);
        }
    }
    igraph_vector_ptr_destroy(list);
}

/**
 * \function igraph_biconnected_components
 * Calculate biconnected components
 *
 * A graph is biconnected if the removal of any single vertex (and
 * its incident edges) does not disconnect it.
 *
 * </para><para>
 * A biconnected component of a graph is a maximal biconnected
 * subgraph of it. The biconnected components of a graph can be given
 * by the partition of its edges: every edge is a member of exactly
 * one biconnected component. Note that this is not true for
 * vertices: the same vertex can be part of many biconnected
 * components.
 *
 * </para><para>
 * Somewhat arbitrarily, igraph does not consider components containing
 * a single vertex only as being biconnected. Isolated vertices will
 * not be part of any of the biconnected components.
 *
 * \param graph The input graph.
 * \param no If not a NULL pointer, the number of biconnected components will
 *     be stored here.
 * \param tree_edges If not a NULL pointer, then the found components
 *     are stored here, in a list of vectors. Every vector in the list
 *     is a biconnected component, represented by its edges. More precisely,
 *     a spanning tree of the biconnected component is returned.
 *     Note you'll have to
 *     destroy each vector first by calling \ref igraph_vector_destroy()
 *     and then \ref igraph_free() on it, plus you need to call
 *     \ref igraph_vector_ptr_destroy() on the list to regain all
 *     allocated memory.
 * \param component_edges If not a NULL pointer, then the edges of the
 *     biconnected components are stored here, in the same form as for
 *     \c tree_edges.
 * \param components If not a NULL pointer, then the vertices of the
 *     biconnected components are stored here, in the same format as
 *     for the previous two arguments.
 * \param articulation_points If not a NULL pointer, then the
 *     articulation points of the graph are stored in this vector.
 *     A vertex is an articulation point if its removal increases the
 *     number of (weakly) connected components in the graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges, but only if you do not calculate \c components and
 * \c component_edges. If you calculate \c components, then it is
 * quadratic in the number of vertices. If you calculate \c
 * component_edges as well, then it is cubic in the number of
 * vertices.
 *
 * \sa \ref igraph_articulation_points(), \ref igraph_clusters().
 *
 * \example examples/simple/igraph_biconnected_components.c
 */

igraph_error_t igraph_biconnected_components(const igraph_t *graph,
                                  igraph_integer_t *no,
                                  igraph_vector_ptr_t *tree_edges,
                                  igraph_vector_ptr_t *component_edges,
                                  igraph_vector_ptr_t *components,
                                  igraph_vector_int_t *articulation_points) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t nextptr;
    igraph_vector_int_t num, low;
    igraph_vector_bool_t found;
    igraph_vector_int_t *adjedges;
    igraph_stack_int_t path;
    igraph_stack_int_t edgestack;
    igraph_inclist_t inclist;
    igraph_integer_t i, counter, rootdfs = 0;
    igraph_vector_int_t vertex_added;
    igraph_integer_t comps = 0;
    igraph_vector_ptr_t *mycomponents = components, vcomponents;

    IGRAPH_CHECK(igraph_vector_int_init(&nextptr, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &nextptr);
    IGRAPH_CHECK(igraph_vector_int_init(&num, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &num);
    IGRAPH_CHECK(igraph_vector_int_init(&low, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &low);
    IGRAPH_CHECK(igraph_vector_bool_init(&found, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &found);

    IGRAPH_CHECK(igraph_stack_int_init(&path, 100));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &path);
    IGRAPH_CHECK(igraph_stack_int_init(&edgestack, 100));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &edgestack);

    IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

    IGRAPH_CHECK(igraph_vector_int_init(&vertex_added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &vertex_added);

    if (no) {
        *no = 0;
    }
    if (tree_edges) {
        igraph_vector_ptr_clear(tree_edges);
    }
    if (components) {
        igraph_vector_ptr_clear(components);
    }
    if (component_edges) {
        igraph_vector_ptr_clear(component_edges);
    }
    if (articulation_points) {
        igraph_vector_int_clear(articulation_points);
    }
    if (component_edges && !components) {
        mycomponents = &vcomponents;
        IGRAPH_CHECK(igraph_vector_ptr_init(mycomponents, 0));
        IGRAPH_FINALLY(igraph_i_free_vector_int_list, mycomponents);
    }

    for (i = 0; i < no_of_nodes; i++) {

        if (VECTOR(low)[i] != 0) {
            continue;    /* already visited */
        }

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_stack_int_push(&path, i));
        counter = 1;
        rootdfs = 0;
        VECTOR(low)[i] = VECTOR(num)[i] = counter++;
        while (!igraph_stack_int_empty(&path)) {
            igraph_integer_t n;
            igraph_integer_t act = igraph_stack_int_top(&path);
            igraph_integer_t actnext = VECTOR(nextptr)[act];

            adjedges = igraph_inclist_get(&inclist, act);
            n = igraph_vector_int_size(adjedges);
            if (actnext < n) {
                /* Step down (maybe) */
                igraph_integer_t edge = VECTOR(*adjedges)[actnext];
                igraph_integer_t nei = IGRAPH_OTHER(graph, edge, act);
                if (VECTOR(low)[nei] == 0) {
                    if (act == i) {
                        rootdfs++;
                    }
                    IGRAPH_CHECK(igraph_stack_int_push(&edgestack, edge));
                    IGRAPH_CHECK(igraph_stack_int_push(&path, nei));
                    VECTOR(low)[nei] = VECTOR(num)[nei] = counter++;
                } else {
                    /* Update low value if needed */
                    if (VECTOR(num)[nei] < VECTOR(low)[act]) {
                        VECTOR(low)[act] = VECTOR(num)[nei];
                    }
                }
                VECTOR(nextptr)[act] += 1;
            } else {
                /* Step up */
                igraph_stack_int_pop(&path);
                if (!igraph_stack_int_empty(&path)) {
                    igraph_integer_t prev = igraph_stack_int_top(&path);
                    /* Update LOW value if needed */
                    if (VECTOR(low)[act] < VECTOR(low)[prev]) {
                        VECTOR(low)[prev] = VECTOR(low)[act];
                    }
                    /* Check for articulation point */
                    if (VECTOR(low)[act] >= VECTOR(num)[prev]) {
                        if (articulation_points && !VECTOR(found)[prev]
                            && prev != i /* the root */) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(articulation_points, prev));
                            VECTOR(found)[prev] = 1;
                        }
                        if (no) {
                            *no += 1;
                        }

                        /*------------------------------------*/
                        /* Record the biconnected component just found */
                        if (tree_edges || mycomponents) {
                            igraph_vector_int_t *v = 0, *v2 = 0;
                            comps++;
                            if (tree_edges) {
                                v = IGRAPH_CALLOC(1, igraph_vector_int_t);
                                if (!v) {
                                    IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
                                }
                                IGRAPH_CHECK(igraph_vector_int_init(v, 0));
                                IGRAPH_FINALLY(igraph_vector_int_destroy, v);
                            }
                            if (mycomponents) {
                                v2 = IGRAPH_CALLOC(1, igraph_vector_int_t);
                                if (!v2) {
                                    IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
                                }
                                IGRAPH_CHECK(igraph_vector_int_init(v2, 0));
                                IGRAPH_FINALLY(igraph_vector_int_destroy, v2);
                            }

                            while (!igraph_stack_int_empty(&edgestack)) {
                                igraph_integer_t e = igraph_stack_int_pop(&edgestack);
                                igraph_integer_t from = IGRAPH_FROM(graph, e);
                                igraph_integer_t to = IGRAPH_TO(graph, e);
                                if (tree_edges) {
                                    IGRAPH_CHECK(igraph_vector_int_push_back(v, e));
                                }
                                if (mycomponents) {
                                    if (VECTOR(vertex_added)[from] != comps) {
                                        VECTOR(vertex_added)[from] = comps;
                                        IGRAPH_CHECK(igraph_vector_int_push_back(v2, from));
                                    }
                                    if (VECTOR(vertex_added)[to] != comps) {
                                        VECTOR(vertex_added)[to] = comps;
                                        IGRAPH_CHECK(igraph_vector_int_push_back(v2, to));
                                    }
                                }
                                if (from == prev || to == prev) {
                                    break;
                                }
                            }

                            if (mycomponents) {
                                IGRAPH_CHECK(igraph_vector_ptr_push_back(mycomponents, v2));
                                IGRAPH_FINALLY_CLEAN(1);
                            }
                            if (tree_edges) {
                                IGRAPH_CHECK(igraph_vector_ptr_push_back(tree_edges, v));
                                IGRAPH_FINALLY_CLEAN(1);
                            }
                            if (component_edges) {
                                igraph_vector_int_t *nodes = VECTOR(*mycomponents)[comps - 1];
                                igraph_vector_int_t *vv = IGRAPH_CALLOC(1, igraph_vector_int_t);
                                igraph_integer_t ii, no_vert = igraph_vector_int_size(nodes);
                                if (!vv) {
                                    IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
                                }
                                IGRAPH_CHECK(igraph_vector_int_init(vv, 0));
                                IGRAPH_FINALLY(igraph_vector_int_destroy, vv);
                                for (ii = 0; ii < no_vert; ii++) {
                                    igraph_integer_t vert = VECTOR(*nodes)[ii];
                                    igraph_vector_int_t *edges = igraph_inclist_get(&inclist,
                                                                 vert);
                                    igraph_integer_t j, nn = igraph_vector_int_size(edges);
                                    for (j = 0; j < nn; j++) {
                                        igraph_integer_t e = VECTOR(*edges)[j];
                                        igraph_integer_t nei = IGRAPH_OTHER(graph, e, vert);
                                        if (VECTOR(vertex_added)[nei] == comps && nei < vert) {
                                            IGRAPH_CHECK(igraph_vector_int_push_back(vv, e));
                                        }
                                    }
                                }
                                IGRAPH_CHECK(igraph_vector_ptr_push_back(component_edges, vv));
                                IGRAPH_FINALLY_CLEAN(1);
                            }
                        } /* record component if requested */
                        /*------------------------------------*/

                    }
                } /* !igraph_stack_int_empty(&path) */
            }

        } /* !igraph_stack_int_empty(&path) */

        if (articulation_points && rootdfs >= 2) {
            IGRAPH_CHECK(igraph_vector_int_push_back(articulation_points, i));
        }

    } /* i < no_of_nodes */

    if (mycomponents != components) {
        igraph_i_free_vector_int_list(mycomponents);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&vertex_added);
    igraph_inclist_destroy(&inclist);
    igraph_stack_int_destroy(&edgestack);
    igraph_stack_int_destroy(&path);
    igraph_vector_bool_destroy(&found);
    igraph_vector_int_destroy(&low);
    igraph_vector_int_destroy(&num);
    igraph_vector_int_destroy(&nextptr);
    IGRAPH_FINALLY_CLEAN(8);

    return IGRAPH_SUCCESS;
}


/* igraph_bridges -- find all bridges in the graph */
/* The algorithm is based on https://www.geeksforgeeks.org/bridge-in-a-graph/
   but instead of keeping track of the parent of each vertex in the DFS tree
   we keep track of its incoming edge. This is necessary to support multigraphs. */

static igraph_error_t igraph_i_bridges_rec(
        const igraph_t *graph, const igraph_inclist_t *il, igraph_integer_t u,
        igraph_integer_t *time, igraph_vector_int_t *bridges, igraph_vector_bool_t *visited,
        igraph_vector_int_t *disc, igraph_vector_int_t *low, igraph_vector_int_t *incoming_edge)
{
    igraph_vector_int_t *incedges;
    igraph_integer_t nc; /* neighbour count */
    igraph_integer_t i;

    VECTOR(*visited)[u] = 1;

    *time += 1;

    VECTOR(*disc)[u] = *time;
    VECTOR(*low)[u] = *time;

    incedges = igraph_inclist_get(il, u);
    nc = igraph_vector_int_size(incedges);
    for (i = 0; i < nc; ++i) {
        igraph_integer_t edge = VECTOR(*incedges)[i];
        igraph_integer_t v = IGRAPH_TO(graph, edge) == u ? IGRAPH_FROM(graph, edge) : IGRAPH_TO(graph, edge);

        if (! VECTOR(*visited)[v]) {
            VECTOR(*incoming_edge)[v] = edge;
            IGRAPH_CHECK(igraph_i_bridges_rec(graph, il, v, time, bridges, visited, disc, low, incoming_edge));

            VECTOR(*low)[u] = VECTOR(*low)[u] < VECTOR(*low)[v] ? VECTOR(*low)[u] : VECTOR(*low)[v];

            if (VECTOR(*low)[v] > VECTOR(*disc)[u]) {
                IGRAPH_CHECK(igraph_vector_int_push_back(bridges, edge));
            }
        } else if (edge != VECTOR(*incoming_edge)[u]) {
            VECTOR(*low)[u] = VECTOR(*low)[u] < VECTOR(*disc)[v] ? VECTOR(*low)[u] : VECTOR(*disc)[v];
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bridges
 * Find all bridges in a graph.
 *
 * An edge is a bridge if its removal increases the number of (weakly)
 * connected components in the graph.
 *
 * \param graph The input graph.
 * \param res Pointer to an initialized vector, the
 *    bridges will be stored here as edge indices.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and edges.
 *
 * \sa \ref igraph_articulation_points(), \ref igraph_biconnected_components(), \ref igraph_clusters()
 */

igraph_error_t igraph_bridges(const igraph_t *graph, igraph_vector_int_t *bridges) {
    igraph_inclist_t il;
    igraph_vector_bool_t visited;
    igraph_vector_int_t disc, low;
    igraph_vector_int_t incoming_edge;
    igraph_integer_t i, n;
    igraph_integer_t time;

    n = igraph_vcount(graph);

    IGRAPH_CHECK(igraph_inclist_init(graph, &il, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &il);

    IGRAPH_CHECK(igraph_vector_bool_init(&visited, n));
    IGRAPH_FINALLY(igraph_vector_bool_destroy, &visited);

    IGRAPH_CHECK(igraph_vector_int_init(&disc, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &disc);

    IGRAPH_CHECK(igraph_vector_int_init(&low, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &low);

    IGRAPH_CHECK(igraph_vector_int_init(&incoming_edge, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &incoming_edge);
    for (i = 0; i < n; ++i) {
        VECTOR(incoming_edge)[i] = -1;
    }

    igraph_vector_int_clear(bridges);

    time = 0;
    for (i = 0; i < n; ++i)
        if (! VECTOR(visited)[i]) {
            IGRAPH_CHECK(igraph_i_bridges_rec(graph, &il, i, &time, bridges, &visited, &disc, &low, &incoming_edge));
        }

    igraph_vector_int_destroy(&incoming_edge);
    igraph_vector_int_destroy(&low);
    igraph_vector_int_destroy(&disc);
    igraph_vector_bool_destroy(&visited);
    igraph_inclist_destroy(&il);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup structural
 * \function igraph_subcomponent
 * \brief The vertices in the same component as a given vertex.
 *
 * \param graph The graph object.
 * \param res The result, vector with the IDs of the vertices in the
 *        same component.
 * \param vertex The id of the vertex of which the component is
 *        searched.
 * \param mode Type of the component for directed graphs, possible
 *        values:
 *        \clist
 *        \cli IGRAPH_OUT
 *          the set of vertices reachable \em from the
 *          \p vertex,
 *        \cli IGRAPH_IN
 *          the set of vertices from which the
 *          \p vertex is reachable.
 *        \cli IGRAPH_ALL
 *          the graph is considered as an
 *          undirected graph. Note that this is \em not the same
 *          as the union of the previous two.
 *        \endclist
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_ENOMEM
 *          not enough memory for temporary data.
 *        \cli IGRAPH_EINVVID
 *           \p vertex is an invalid vertex ID
 *        \cli IGRAPH_EINVMODE
 *           invalid mode argument passed.
 *        \endclist
 *
 * Time complexity: O(|V|+|E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 *
 * \sa \ref igraph_induced_subgraph() if you want a graph object consisting only
 * a given set of vertices and the edges between them.
 */
igraph_error_t igraph_subcomponent(
    const igraph_t *graph, igraph_vector_int_t *res, igraph_integer_t vertex,
    igraph_neimode_t mode
) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q = IGRAPH_DQUEUE_NULL;
    char *already_added;
    igraph_integer_t i, vsize;
    igraph_vector_int_t tmp = IGRAPH_VECTOR_NULL;

    if (vertex < 0 || vertex >= no_of_nodes) {
        IGRAPH_ERROR("Vertex id out of range.", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    already_added = IGRAPH_CALLOC(no_of_nodes, char);
    if (already_added == 0) {
        IGRAPH_ERROR("Subcomponent failed.", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, already_added);

    igraph_vector_int_clear(res);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);

    IGRAPH_CHECK(igraph_dqueue_int_push(&q, vertex));
    IGRAPH_CHECK(igraph_vector_int_push_back(res, vertex));
    already_added[vertex] = 1;

    while (!igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actnode = igraph_dqueue_int_pop(&q);

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_neighbors(graph, &tmp, actnode, mode));
        vsize = igraph_vector_int_size(&tmp);
        for (i = 0; i < vsize; i++) {
            igraph_integer_t neighbor = VECTOR(tmp)[i];

            if (already_added[neighbor]) {
                continue;
            }
            already_added[neighbor] = 1;
            IGRAPH_CHECK(igraph_vector_int_push_back(res, neighbor));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
        }
    }

    igraph_dqueue_int_destroy(&q);
    igraph_vector_int_destroy(&tmp);
    IGRAPH_FREE(already_added);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
