/*
   IGraph library.
   Copyright (C) 2006-2023  The igraph development team <igraph@igraph.org>

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

#include "igraph_visitor.h"
#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_dqueue.h"
#include "igraph_stack.h"

/**
 * \function igraph_bfs
 * \brief Breadth-first search.
 *
 * A simple breadth-first search, with a lot of different results and
 * the possibility to call a callback whenever a vertex is visited.
 * It is allowed to supply null pointers as the output arguments the
 * user is not interested in, in this case they will be ignored.
 *
 * </para><para>
 * If not all vertices can be reached from the supplied root vertex,
 * then additional root vertices will be used, in the order of their
 * vertex IDs.
 *
 * </para><para>
 * Consider using \ref igraph_bfs_simple instead if you set most of the output
 * arguments provided by this function to a null pointer.
 *
 * \param graph The input graph.
 * \param root The id of the root vertex. It is ignored if the \c
 *        roots argument is not a null pointer.
 * \param roots Pointer to an initialized vector, or a null
 *        pointer. If not a null pointer, then it is a vector
 *        containing root vertices to start the BFS from. The vertices
 *        are considered in the order they appear. If a root vertex
 *        was already found while searching from another one, then no
 *        search is conducted from it.
 * \param mode For directed graphs, it defines which edges to follow.
 *        \c IGRAPH_OUT means following the direction of the edges,
 *        \c IGRAPH_IN means the opposite, and
 *        \c IGRAPH_ALL ignores the direction of the edges.
 *        This parameter is ignored for undirected graphs.
 * \param unreachable Logical scalar, whether the search should visit
 *        the vertices that are unreachable from the given root
 *        node(s). If true, then additional searches are performed
 *        until all vertices are visited.
 * \param restricted If not a null pointer, then it must be a pointer
 *        to a vector containing vertex IDs. The BFS is carried out
 *        only on these vertices.
 * \param order If not null pointer, then the vertex IDs of the graph are
 *        stored here, in the same order as they were visited.
 * \param rank If not a null pointer, then the rank of each vertex is
 *        stored here.
 * \param parents If not a null pointer, then the id of the parent of
 *        each vertex is stored here. When a vertex was not visited
 *        during the traversal, -2 will be stored as the ID of its parent.
 *        When a vertex was visited during the traversal and it was one of
 *        the roots of the search trees, -1 will be stored as the ID of
 *        its parent.
 * \param pred If not a null pointer, then the id of vertex that was
 *        visited before the current one is stored here. If there is
 *        no such vertex (the current vertex is the root of a search
 *        tree), then -1 is stored as the predecessor of the vertex.
 *        If the vertex was not visited at all, then -2 is stored for
 *        the predecessor of the vertex.
 * \param succ If not a null pointer, then the id of the vertex that
 *        was visited after the current one is stored here. If there
 *        is no such vertex (the current one is the last in a search
 *        tree), then -1 is stored as the successor of the vertex.
 *        If the vertex was not visited at all, then -2 is stored for
 *        the successor of the vertex.
 * \param dist If not a null pointer, then the distance from the root of
 *        the current search tree is stored here for each vertex. If a
 *        vertex was not reached during the traversal, its distance will
 *        be -1 in this vector.
 * \param callback If not null, then it should be a pointer to a
 *        function of type \ref igraph_bfshandler_t. This function
 *        will be called, whenever a new vertex is visited.
 * \param extra Extra argument to pass to the callback function.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_bfs.c
 * \example examples/simple/igraph_bfs_callback.c
 */
igraph_error_t igraph_bfs(const igraph_t *graph,
               igraph_integer_t root, const igraph_vector_int_t *roots,
               igraph_neimode_t mode, igraph_bool_t unreachable,
               const igraph_vector_int_t *restricted,
               igraph_vector_int_t *order, igraph_vector_int_t *rank,
               igraph_vector_int_t *parents,
               igraph_vector_int_t *pred, igraph_vector_int_t *succ,
               igraph_vector_int_t *dist, igraph_bfshandler_t *callback,
               void *extra) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);

    igraph_error_t ret;

    igraph_dqueue_int_t Q;
    igraph_integer_t actroot = 0;
    igraph_vector_char_t added;

    igraph_lazy_adjlist_t adjlist;

    igraph_integer_t act_rank = 0;
    igraph_integer_t pred_vec = -1;

    igraph_integer_t rootpos = 0;
    igraph_integer_t noroots = roots ? igraph_vector_int_size(roots) : 1;

    if (!roots && (root < 0 || root >= no_of_nodes)) {
        IGRAPH_ERROR("Invalid root vertex in BFS.", IGRAPH_EINVVID);
    }

    if (roots && !igraph_vector_int_isininterval(roots, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("Invalid root vertex in BFS.", IGRAPH_EINVVID);
    }

    if (restricted && !igraph_vector_int_isininterval(restricted, 0, no_of_nodes-1)) {
        IGRAPH_ERROR("Invalid vertex ID in restricted set.", IGRAPH_EINVVID);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&added, no_of_nodes);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&Q, 100);

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    /* Mark the vertices that are not in the restricted set, as already
       found. Special care must be taken for vertices that are not in
       the restricted set, but are to be used as 'root' vertices. */
    if (restricted) {
        igraph_integer_t i, n = igraph_vector_int_size(restricted);
        igraph_vector_char_fill(&added, true);
        for (i = 0; i < n; i++) {
            igraph_integer_t v = VECTOR(*restricted)[i];
            VECTOR(added)[v] = false;
        }
    }

    /* Resize result vectors, and fill them with the initial value. */

# define VINIT(v, initial) \
    if (v) { \
        IGRAPH_CHECK(igraph_vector_int_resize((v), no_of_nodes)); \
        igraph_vector_int_fill((v), initial); \
    }

    VINIT(order, -1);
    VINIT(rank, -1);
    VINIT(parents, -2);
    VINIT(pred, -2);
    VINIT(succ, -2);
    VINIT(dist, -1);
# undef VINIT

    while (1) {

        /* Get the next root vertex, if any */

        if (roots && rootpos < noroots) {
            /* We are still going through the 'roots' vector */
            actroot = VECTOR(*roots)[rootpos++];
        } else if (!roots && rootpos == 0) {
            /* We have a single root vertex given, and start now */
            actroot = root;
            rootpos++;
        } else if (rootpos == noroots && unreachable) {
            /* We finished the given root(s), but other vertices are also
            tried as root */
            actroot = 0;
            rootpos++;
        } else if (unreachable && actroot + 1 < no_of_nodes) {
            /* We are already doing the other vertices, take the next one */
            actroot++;
        } else {
            /* No more root nodes to do */
            break;
        }

        /* OK, we have a new root, start BFS */
        if (VECTOR(added)[actroot]) {
            continue;
        }
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, actroot));
        IGRAPH_CHECK(igraph_dqueue_int_push(&Q, 0));
        VECTOR(added)[actroot] = true;
        if (parents) {
            VECTOR(*parents)[actroot] = -1;
        }

        pred_vec = -1;

        while (!igraph_dqueue_int_empty(&Q)) {
            igraph_integer_t actvect = igraph_dqueue_int_pop(&Q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&Q);
            igraph_integer_t succ_vec;
            igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, actvect);

            IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");
            const igraph_integer_t n = igraph_vector_int_size(neis);

            if (pred) {
                VECTOR(*pred)[actvect] = pred_vec;
            }
            if (rank) {
                VECTOR(*rank)[actvect] = act_rank;
            }
            if (order) {
                VECTOR(*order)[act_rank++] = actvect;
            }
            if (dist) {
                VECTOR(*dist)[actvect] = actdist;
            }

            for (igraph_integer_t i = 0; i < n; i++) {
                igraph_integer_t nei = VECTOR(*neis)[i];
                if (! VECTOR(added)[nei]) {
                    VECTOR(added)[nei] = true;
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, nei));
                    IGRAPH_CHECK(igraph_dqueue_int_push(&Q, actdist + 1));
                    if (parents) {
                        VECTOR(*parents)[nei] = actvect;
                    }
                }
            }

            succ_vec = igraph_dqueue_int_empty(&Q)
                           ? -1
                           : igraph_dqueue_int_head(&Q);
            if (callback) {
                IGRAPH_CHECK_CALLBACK(
                    callback(graph, actvect, pred_vec, succ_vec, act_rank - 1, actdist, extra),
                    &ret
                );

                if (ret == IGRAPH_STOP) {
                    goto cleanup;
                }
            }

            if (succ) {
                VECTOR(*succ)[actvect] = succ_vec;
            }
            pred_vec = actvect;

        } /* while Q !empty */

    } /* for actroot < no_of_nodes */

cleanup:

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_dqueue_int_destroy(&Q);
    igraph_vector_char_destroy(&added);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_bfs_simple
 * Breadth-first search, single-source version
 *
 * An alternative breadth-first search implementation to cater for the
 * simpler use-cases when only a single breadth-first search has to be conducted
 * from a source node and most of the output arguments from \ref igraph_bfs
 * are not needed. It is allowed to supply null pointers as
 * the output arguments the user is not interested in, in this case they will
 * be ignored.
 *
 * \param graph The input graph.
 * \param root The id of the root vertex.
 * \param mode For directed graphs, it defines which edges to follow.
 *        \c IGRAPH_OUT means following the direction of the edges,
 *        \c IGRAPH_IN means the opposite, and
 *        \c IGRAPH_ALL ignores the direction of the edges.
 *        This parameter is ignored for undirected graphs.
 * \param order If not a null pointer, then an initialized vector must be passed
 *        here. The IDs of the vertices visited during the traversal will be
 *        stored here, in the same order as they were visited.
 * \param layers If not a null pointer, then an initialized vector must be
 *        passed here. The i-th element of the vector will contain the index
 *        into \c order where the vertices that are at distance i from the root
 *        are stored. In other words, if you are interested in the vertices that
 *        are at distance i from the root, you need to look in the \c order
 *        vector from \c layers[i] to \c layers[i+1].
 * \param parents If not a null pointer, then an initialized vector must be
 *        passed here. The vector will be resized so its length is equal to the
 *        number of nodes, and it will contain the index of the parent node for
 *        each \em visited node. The values in the vector are set to -2 for
 *        vertices that were \em not visited, and -1 for the root vertex.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_bfs_simple.c
 */
igraph_error_t igraph_bfs_simple(
    const igraph_t *graph, igraph_integer_t root, igraph_neimode_t mode,
    igraph_vector_int_t *order, igraph_vector_int_t *layers,
    igraph_vector_int_t *parents
) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_integer_t num_visited = 0;
    igraph_vector_int_t neis;
    bool *added;
    igraph_integer_t lastlayer = -1;

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    /* temporary storage */

    added = IGRAPH_CALLOC(no_of_nodes, bool);
    IGRAPH_CHECK_OOM(added, "Insufficient memory for BFS.");
    IGRAPH_FINALLY(igraph_free, added);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_int_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_int_destroy, &q);

    /* results */
    if (order) {
        igraph_vector_int_clear(order);
    }
    if (layers) {
        igraph_vector_int_clear(layers);
    }
    if (parents) {
        IGRAPH_CHECK(igraph_vector_int_resize(parents, no_of_nodes));
        igraph_vector_int_fill(parents, -2);
    }

    /* ok start with root */
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, root));
    IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
    if (layers) {
        IGRAPH_CHECK(igraph_vector_int_push_back(layers, num_visited));
    }
    if (order) {
        IGRAPH_CHECK(igraph_vector_int_push_back(order, root));
    }
    if (parents) {
        VECTOR(*parents)[root] = -1;
    }
    num_visited++;
    added[root] = true;

    while (!igraph_dqueue_int_empty(&q)) {
        igraph_integer_t actvect = igraph_dqueue_int_pop(&q);
        igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, actvect,
                                      mode));
        igraph_integer_t nei_count = igraph_vector_int_size(&neis);
        for (igraph_integer_t i = 0; i < nei_count; i++) {
            const igraph_integer_t neighbor = VECTOR(neis)[i];
            if (! added[neighbor]) {
                added[neighbor] = true;
                if (parents) {
                    VECTOR(*parents)[neighbor] = actvect;
                }
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                if (layers && lastlayer != actdist + 1) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(layers, num_visited));
                }
                if (order) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(order, neighbor));
                }
                num_visited++;
                lastlayer = actdist + 1;
            }
        } /* for i in neis */
    } /* while ! dqueue_int_empty */

    if (layers) {
        IGRAPH_CHECK(igraph_vector_int_push_back(layers, num_visited));
    }

    igraph_vector_int_destroy(&neis);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_dfs
 * \brief Depth-first search.
 *
 * A simple depth-first search, with
 * the possibility to call a callback whenever a vertex is discovered
 * and/or whenever a subtree is finished.
 * It is allowed to supply null pointers as the output arguments the
 * user is not interested in, in this case they will be ignored.
 *
 * </para><para>
 * If not all vertices can be reached from the supplied root vertex,
 * then additional root vertices will be used, in the order of their
 * vertex IDs.
 *
 * \param graph The input graph.
 * \param root The id of the root vertex.
 * \param mode For directed graphs, it defines which edges to follow.
 *        \c IGRAPH_OUT means following the direction of the edges,
 *        \c IGRAPH_IN means the opposite, and
 *        \c IGRAPH_ALL ignores the direction of the edges.
 *        This parameter is ignored for undirected graphs.
 * \param unreachable Logical scalar, whether the search should visit
 *        the vertices that are unreachable from the given root
 *        node(s). If true, then additional searches are performed
 *        until all vertices are visited.
 * \param order If not null pointer, then the vertex IDs of the graph are
 *        stored here, in the same order as they were discovered. The tail of
 *        the vector will be padded with -1 to ensure that the length of the
 *        vector is the same as the number of vertices, even if some vertices
 *        were not visited during the traversal.
 * \param order_out If not a null pointer, then the vertex IDs of the
 *        graphs are stored here, in the order of the completion of
 *        their subtree. The tail of the vector will be padded with -1 to ensure
 *        that the length of the vector is the same as the number of vertices,
 *        even if some vertices were not visited during the traversal.
 * \param parents If not a null pointer, then the id of the parent of
 *        each vertex is stored here. -1 will be stored for the root of the
 *        search tree; -2 will be stored for vertices that were not visited.
 * \param dist If not a null pointer, then the distance from the root of
 *        the current search tree is stored here. -1 will be stored for vertices
 *        that were not visited.
 * \param in_callback If not null, then it should be a pointer to a
 *        function of type \ref igraph_dfshandler_t. This function
 *        will be called, whenever a new vertex is discovered.
 * \param out_callback If not null, then it should be a pointer to a
 *        function of type \ref igraph_dfshandler_t. This function
 *        will be called, whenever the subtree of a vertex is completed.
 * \param extra Extra argument to pass to the callback function(s).
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 */

igraph_error_t igraph_dfs(const igraph_t *graph, igraph_integer_t root,
               igraph_neimode_t mode, igraph_bool_t unreachable,
               igraph_vector_int_t *order,
               igraph_vector_int_t *order_out, igraph_vector_int_t *parents,
               igraph_vector_int_t *dist, igraph_dfshandler_t *in_callback,
               igraph_dfshandler_t *out_callback,
               void *extra) {

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_lazy_adjlist_t adjlist;
    igraph_stack_int_t stack;
    igraph_vector_char_t added;
    igraph_vector_int_t nptr;
    igraph_error_t ret;
    igraph_integer_t act_rank = 0;
    igraph_integer_t rank_out = 0;
    igraph_integer_t act_dist = 0;

    if (root < 0 || root >= no_of_nodes) {
        IGRAPH_ERROR("Invalid root vertex for DFS.", IGRAPH_EINVAL);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument.", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&added, no_of_nodes);
    IGRAPH_STACK_INT_INIT_FINALLY(&stack, 100);

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&nptr, no_of_nodes);

# define FREE_ALL() do { \
        igraph_vector_int_destroy(&nptr); \
        igraph_lazy_adjlist_destroy(&adjlist); \
        igraph_stack_int_destroy(&stack); \
        igraph_vector_char_destroy(&added); \
        IGRAPH_FINALLY_CLEAN(4); } while (0)

    /* Resize result vectors and fill them with the initial value */

# define VINIT(v, initial) if (v) { \
        IGRAPH_CHECK(igraph_vector_int_resize(v, no_of_nodes)); \
        igraph_vector_int_fill(v, initial); }

    VINIT(order, -1);
    VINIT(order_out, -1);
    VINIT(parents, -2);
    VINIT(dist, -1);

# undef VINIT

    IGRAPH_CHECK(igraph_stack_int_push(&stack, root));
    VECTOR(added)[root] = true;
    if (parents) {
        VECTOR(*parents)[root] = -1;
    }
    if (order) {
        VECTOR(*order)[act_rank++] = root;
    }
    if (dist) {
        VECTOR(*dist)[root] = 0;
    }
    if (in_callback) {
        IGRAPH_CHECK_CALLBACK(in_callback(graph, root, 0, extra), &ret);
        if (ret == IGRAPH_STOP) {
            FREE_ALL();
            return IGRAPH_SUCCESS;
        }
    }

    for (igraph_integer_t actroot = 0; actroot < no_of_nodes; ) {

        /* 'root' first, then all other vertices */
        if (igraph_stack_int_empty(&stack)) {
            if (!unreachable) {
                break;
            }
            if (VECTOR(added)[actroot]) {
                actroot++;
                continue;
            }
            IGRAPH_CHECK(igraph_stack_int_push(&stack, actroot));
            VECTOR(added)[actroot] = true;
            if (parents) {
                VECTOR(*parents)[actroot] = -1;
            }
            if (order) {
                VECTOR(*order)[act_rank++] = actroot;
            }
            if (dist) {
                VECTOR(*dist)[actroot] = 0;
            }

            if (in_callback) {
                IGRAPH_CHECK_CALLBACK(in_callback(graph, actroot, 0, extra), &ret);
                if (ret == IGRAPH_STOP) {
                    FREE_ALL();
                    return IGRAPH_SUCCESS;
                }
            }

            actroot++;
        }

        while (!igraph_stack_int_empty(&stack)) {
            igraph_integer_t actvect = igraph_stack_int_top(&stack);
            igraph_integer_t *ptr = igraph_vector_int_get_ptr(&nptr, actvect);

            igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist, actvect);
            IGRAPH_CHECK_OOM(neis, "Failed to query neighbors.");

            const igraph_integer_t n = igraph_vector_int_size(neis);

            /* Search for a neighbor that was not yet visited */
            igraph_bool_t any = false;
            igraph_integer_t nei = 0;
            while (!any && (*ptr) < n) {
                nei = VECTOR(*neis)[(*ptr)];
                any = !VECTOR(added)[nei];
                (*ptr) ++;
            }
            if (any) {
                /* There is such a neighbor, add it */
                IGRAPH_CHECK(igraph_stack_int_push(&stack, nei));
                VECTOR(added)[nei] = true;
                if (parents) {
                    VECTOR(*parents)[ nei ] = actvect;
                }
                if (order) {
                    VECTOR(*order)[act_rank++] = nei;
                }
                act_dist++;
                if (dist) {
                    VECTOR(*dist)[nei] = act_dist;
                }

                if (in_callback) {
                    IGRAPH_CHECK_CALLBACK(
                        in_callback(graph, nei, act_dist, extra),
                        &ret
                    );
                    if (ret == IGRAPH_STOP) {
                        FREE_ALL();
                        return IGRAPH_SUCCESS;
                    }
                }

            } else {
                /* There is no such neighbor, finished with the subtree */
                igraph_stack_int_pop(&stack);
                if (order_out) {
                    VECTOR(*order_out)[rank_out++] = actvect;
                }
                act_dist--;

                if (out_callback) {
                    IGRAPH_CHECK_CALLBACK(
                        out_callback(graph, actvect, act_dist, extra),
                        &ret
                    );

                    if (ret == IGRAPH_STOP) {
                        FREE_ALL();
                        return IGRAPH_SUCCESS;
                    }
                }
            }
        }
    }

    FREE_ALL();
# undef FREE_ALL

    return IGRAPH_SUCCESS;
}
