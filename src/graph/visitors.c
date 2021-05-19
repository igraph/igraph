/* -*- mode: C -*-  */
/*
   IGraph R package.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_visitor.h"
#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_dqueue.h"
#include "igraph_stack.h"

/**
 * \function igraph_bfs
 * Breadth-first search
 *
 * A simple breadth-first search, with a lot of different results and
 * the possibility to call a callback whenever a vertex is visited.
 * It is allowed to supply null pointers as the output arguments the
 * user is not interested in, in this case they will be ignored.
 *
 * </para><para>
 * If not all vertices can be reached from the supplied root vertex,
 * then additional root vertices will be used, in the order of their
 * vertex ids.
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
 *        to a vector containing vertex ids. The BFS is carried out
 *        only on these vertices.
 * \param order If not null pointer, then the vertex ids of the graph are
 *        stored here, in the same order as they were visited.
 * \param rank If not a null pointer, then the rank of each vertex is
 *        stored here.
 * \param father If not a null pointer, then the id of the father of
 *        each vertex is stored here.
 * \param pred If not a null pointer, then the id of vertex that was
 *        visited before the current one is stored here. If there is
 *        no such vertex (the current vertex is the root of a search
 *        tree), then -1 is stored.
 * \param succ If not a null pointer, then the id of the vertex that
 *        was visited after the current one is stored here. If there
 *        is no such vertex (the current one is the last in a search
 *        tree), then -1 is stored.
 * \param dist If not a null pointer, then the distance from the root of
 *        the current search tree is stored here.
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
int igraph_bfs(const igraph_t *graph,
               igraph_integer_t root, const igraph_vector_t *roots,
               igraph_neimode_t mode, igraph_bool_t unreachable,
               const igraph_vector_t *restricted,
               igraph_vector_t *order, igraph_vector_t *rank,
               igraph_vector_t *father,
               igraph_vector_t *pred, igraph_vector_t *succ,
               igraph_vector_t *dist, igraph_bfshandler_t *callback,
               void *extra) {

    igraph_dqueue_t Q;
    long int no_of_nodes = igraph_vcount(graph);
    long int actroot = 0;
    igraph_vector_char_t added;

    igraph_lazy_adjlist_t adjlist;

    long int act_rank = 0;
    long int pred_vec = -1;

    long int rootpos = 0;
    long int noroots = roots ? igraph_vector_size(roots) : 1;

    if (!roots && (root < 0 || root >= no_of_nodes)) {
        IGRAPH_ERROR("Invalid root vertex in BFS", IGRAPH_EINVAL);
    }

    if (roots) {
        igraph_real_t min, max;
        igraph_vector_minmax(roots, &min, &max);
        if (min < 0 || max >= no_of_nodes) {
            IGRAPH_ERROR("Invalid root vertex in BFS", IGRAPH_EINVAL);
        }
    }

    if (restricted) {
        igraph_real_t min, max;
        igraph_vector_minmax(restricted, &min, &max);
        if (min < 0 || max >= no_of_nodes) {
            IGRAPH_ERROR("Invalid vertex id in restricted set", IGRAPH_EINVAL);
        }
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
    IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

    /* Mark the vertices that are not in the restricted set, as already
       found. Special care must be taken for vertices that are not in
       the restricted set, but are to be used as 'root' vertices. */
    if (restricted) {
        long int i, n = igraph_vector_size(restricted);
        igraph_vector_char_fill(&added, 1);
        for (i = 0; i < n; i++) {
            long int v = (long int) VECTOR(*restricted)[i];
            VECTOR(added)[v] = 0;
        }
    }

    /* Resize result vectors, and fill them with IGRAPH_NAN */

# define VINIT(v) if (v) {                      \
        igraph_vector_resize((v), no_of_nodes);   \
        igraph_vector_fill((v), IGRAPH_NAN); }

    VINIT(order);
    VINIT(rank);
    VINIT(father);
    VINIT(pred);
    VINIT(succ);
    VINIT(dist);
# undef VINIT

    while (1) {

        /* Get the next root vertex, if any */

        if (roots && rootpos < noroots) {
            /* We are still going through the 'roots' vector */
            actroot = (long int) VECTOR(*roots)[rootpos++];
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
        IGRAPH_CHECK(igraph_dqueue_push(&Q, actroot));
        IGRAPH_CHECK(igraph_dqueue_push(&Q, 0));
        VECTOR(added)[actroot] = 1;
        if (father) {
            VECTOR(*father)[actroot] = -1;
        }

        pred_vec = -1;

        while (!igraph_dqueue_empty(&Q)) {
            long int actvect = (long int) igraph_dqueue_pop(&Q);
            long int actdist = (long int) igraph_dqueue_pop(&Q);
            long int succ_vec;
            igraph_vector_int_t *neis = igraph_lazy_adjlist_get(&adjlist,
                                    (igraph_integer_t) actvect);
            long int i, n = igraph_vector_int_size(neis);

            if (pred) {
                VECTOR(*pred)[actvect] = pred_vec;
            }
            if (rank) {
                VECTOR(*rank) [actvect] = act_rank;
            }
            if (order) {
                VECTOR(*order)[act_rank++] = actvect;
            }
            if (dist) {
                VECTOR(*dist)[actvect] = actdist;
            }

            for (i = 0; i < n; i++) {
                long int nei = (long int) VECTOR(*neis)[i];
                if (! VECTOR(added)[nei]) {
                    VECTOR(added)[nei] = 1;
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
                    IGRAPH_CHECK(igraph_dqueue_push(&Q, actdist + 1));
                    if (father) {
                        VECTOR(*father)[nei] = actvect;
                    }
                }
            }

            succ_vec = igraph_dqueue_empty(&Q) ? -1L :
                       (long int) igraph_dqueue_head(&Q);
            if (callback) {
                igraph_bool_t terminate =
                    callback(graph, (igraph_integer_t) actvect, (igraph_integer_t)
                             pred_vec, (igraph_integer_t) succ_vec,
                             (igraph_integer_t) act_rank - 1, (igraph_integer_t) actdist,
                             extra);
                if (terminate) {
                    igraph_lazy_adjlist_destroy(&adjlist);
                    igraph_dqueue_destroy(&Q);
                    igraph_vector_char_destroy(&added);
                    IGRAPH_FINALLY_CLEAN(3);
                    return 0;
                }
            }

            if (succ) {
                VECTOR(*succ)[actvect] = succ_vec;
            }
            pred_vec = actvect;

        } /* while Q !empty */

    } /* for actroot < no_of_nodes */

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_dqueue_destroy(&Q);
    igraph_vector_char_destroy(&added);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
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
 * \param vid The id of the root vertex.
 * \param mode For directed graphs, it defines which edges to follow.
 *        \c IGRAPH_OUT means following the direction of the edges,
 *        \c IGRAPH_IN means the opposite, and
 *        \c IGRAPH_ALL ignores the direction of the edges.
 *        This parameter is ignored for undirected graphs.
 * \param vids If not a null pointer, then an initialized vector must be passed
 *        here. The ids of the vertices visited during the traversal will be
 *        stored here, in the same order as they were visited.
 * \param layers If not a null pointer, then an initialized vector must be
 *        passed here. The i-th element of the vector will contain the index
 *        into \c vids where the vertices that are at distance i from the root
 *        are stored. In other words, if you are interested in the vertices that
 *        are at distance i from the root, you need to look in the \c vids
 *        vector from \c layers[i] to \c layers[i+1].
 * \param parents If not a null pointer, then an initialized vector must be
 *        passed here. The vector will be resized so its length is equal to the
 *        number of nodes, and it will contain the index of the parent node for
 *        each \em visited node. The values in the vector are undefined for
 *        vertices that were \em not visited.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number of vertices and
 * edges.
 *
 * \example examples/simple/igraph_bfs_simple.c
 */
int igraph_bfs_simple(igraph_t *graph, igraph_integer_t vid, igraph_neimode_t mode,
                      igraph_vector_t *vids, igraph_vector_t *layers,
                      igraph_vector_t *parents) {

    igraph_dqueue_t q;
    long int num_visited = 0;
    igraph_vector_t neis;
    long int no_of_nodes = igraph_vcount(graph);
    long int i;
    char *added;
    long int lastlayer = -1;

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    /* temporary storage */
    added = IGRAPH_CALLOC(no_of_nodes, char);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate BFS", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
    IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
    IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

    /* results */
    if (vids) {
        igraph_vector_clear(vids);
    }
    if (layers) {
        igraph_vector_clear(layers);
    }
    if (parents) {
        IGRAPH_CHECK(igraph_vector_resize(parents, no_of_nodes));
    }

    /* ok start with vid */
    IGRAPH_CHECK(igraph_dqueue_push(&q, vid));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    if (layers) {
        IGRAPH_CHECK(igraph_vector_push_back(layers, num_visited));
    }
    if (vids) {
        IGRAPH_CHECK(igraph_vector_push_back(vids, vid));
    }
    if (parents) {
        VECTOR(*parents)[(long int)vid] = vid;
    }
    num_visited++;
    added[(long int)vid] = 1;

    while (!igraph_dqueue_empty(&q)) {
        long int actvect = (long int) igraph_dqueue_pop(&q);
        long int actdist = (long int) igraph_dqueue_pop(&q);
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) actvect,
                                      mode));
        for (i = 0; i < igraph_vector_size(&neis); i++) {
            long int neighbor = (long int) VECTOR(neis)[i];
            if (added[neighbor] == 0) {
                added[neighbor] = 1;
                if (parents) {
                    VECTOR(*parents)[neighbor] = actvect;
                }
                IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
                IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                if (layers && lastlayer != actdist + 1) {
                    IGRAPH_CHECK(igraph_vector_push_back(layers, num_visited));
                }
                if (vids) {
                    IGRAPH_CHECK(igraph_vector_push_back(vids, neighbor));
                }
                num_visited++;
                lastlayer = actdist + 1;
            }
        } /* for i in neis */
    } /* while ! dqueue_empty */

    if (layers) {
        IGRAPH_CHECK(igraph_vector_push_back(layers, num_visited));
    }

    igraph_vector_destroy(&neis);
    igraph_dqueue_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_dfs
 * Depth-first search
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
 * vertex ids.
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
 * \param order If not null pointer, then the vertex ids of the graph are
 *        stored here, in the same order as they were discovered.
 * \param order_out If not a null pointer, then the vertex ids of the
 *        graphs are stored here, in the order of the completion of
 *        their subtree.
 * \param father If not a null pointer, then the id of the father of
 *        each vertex is stored here.
 * \param dist If not a null pointer, then the distance from the root of
 *        the current search tree is stored here.
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

int igraph_dfs(const igraph_t *graph, igraph_integer_t root,
               igraph_neimode_t mode, igraph_bool_t unreachable,
               igraph_vector_t *order,
               igraph_vector_t *order_out, igraph_vector_t *father,
               igraph_vector_t *dist, igraph_dfshandler_t *in_callback,
               igraph_dfshandler_t *out_callback,
               void *extra) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_lazy_adjlist_t adjlist;
    igraph_stack_t stack;
    igraph_vector_char_t added;
    igraph_vector_long_t nptr;
    long int actroot;
    long int act_rank = 0;
    long int rank_out = 0;
    long int act_dist = 0;

    if (root < 0 || root >= no_of_nodes) {
        IGRAPH_ERROR("Invalid root vertex for DFS", IGRAPH_EINVAL);
    }

    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
    IGRAPH_CHECK(igraph_stack_init(&stack, 100));
    IGRAPH_FINALLY(igraph_stack_destroy, &stack);
    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    IGRAPH_CHECK(igraph_vector_long_init(&nptr, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &nptr);

# define FREE_ALL() do {            \
        igraph_vector_long_destroy(&nptr);            \
        igraph_lazy_adjlist_destroy(&adjlist);        \
        igraph_stack_destroy(&stack);                 \
        igraph_vector_char_destroy(&added);           \
        IGRAPH_FINALLY_CLEAN(4); } while (0)

    /* Resize result vectors and fill them with IGRAPH_NAN */

# define VINIT(v) if (v) {                      \
        igraph_vector_resize(v, no_of_nodes);       \
        igraph_vector_fill(v, IGRAPH_NAN); }

    VINIT(order);
    VINIT(order_out);
    VINIT(father);
    VINIT(dist);

# undef VINIT

    IGRAPH_CHECK(igraph_stack_push(&stack, root));
    VECTOR(added)[(long int)root] = 1;
    if (father) {
        VECTOR(*father)[(long int)root] = -1;
    }
    if (order) {
        VECTOR(*order)[act_rank++] = root;
    }
    if (dist) {
        VECTOR(*dist)[(long int)root] = 0;
    }
    if (in_callback) {
        igraph_bool_t terminate = in_callback(graph, root, 0, extra);
        if (terminate) {
            FREE_ALL();
            return 0;
        }
    }

    for (actroot = 0; actroot < no_of_nodes; ) {

        /* 'root' first, then all other vertices */
        if (igraph_stack_empty(&stack)) {
            if (!unreachable) {
                break;
            }
            if (VECTOR(added)[actroot]) {
                actroot++;
                continue;
            }
            IGRAPH_CHECK(igraph_stack_push(&stack, actroot));
            VECTOR(added)[actroot] = 1;
            if (father) {
                VECTOR(*father)[actroot] = -1;
            }
            if (order) {
                VECTOR(*order)[act_rank++] = actroot;
            }
            if (dist) {
                VECTOR(*dist)[actroot] = 0;
            }

            if (in_callback) {
                igraph_bool_t terminate = in_callback(graph, (igraph_integer_t) actroot,
                                                      0, extra);
                if (terminate) {
                    FREE_ALL();
                    return 0;
                }
            }
            actroot++;
        }

        while (!igraph_stack_empty(&stack)) {
            long int actvect = (long int) igraph_stack_top(&stack);
            igraph_vector_int_t *neis =
                igraph_lazy_adjlist_get(&adjlist, (igraph_integer_t) actvect);
            long int n = igraph_vector_int_size(neis);
            long int *ptr = igraph_vector_long_e_ptr(&nptr, actvect);

            /* Search for a neighbor that was not yet visited */
            igraph_bool_t any = 0;
            long int nei = 0;
            while (!any && (*ptr) < n) {
                nei = (long int) VECTOR(*neis)[(*ptr)];
                any = !VECTOR(added)[nei];
                (*ptr) ++;
            }
            if (any) {
                /* There is such a neighbor, add it */
                IGRAPH_CHECK(igraph_stack_push(&stack, nei));
                VECTOR(added)[nei] = 1;
                if (father) {
                    VECTOR(*father)[ nei ] = actvect;
                }
                if (order) {
                    VECTOR(*order)[act_rank++] = nei;
                }
                act_dist++;
                if (dist) {
                    VECTOR(*dist)[nei] = act_dist;
                }

                if (in_callback) {
                    igraph_bool_t terminate = in_callback(graph, (igraph_integer_t) nei,
                                                          (igraph_integer_t) act_dist,
                                                          extra);
                    if (terminate) {
                        FREE_ALL();
                        return 0;
                    }
                }

            } else {
                /* There is no such neighbor, finished with the subtree */
                igraph_stack_pop(&stack);
                if (order_out) {
                    VECTOR(*order_out)[rank_out++] = actvect;
                }
                act_dist--;

                if (out_callback) {
                    igraph_bool_t terminate = out_callback(graph, (igraph_integer_t)
                                                           actvect, (igraph_integer_t)
                                                           act_dist, extra);
                    if (terminate) {
                        FREE_ALL();
                        return 0;
                    }
                }
            }
        }
    }

    FREE_ALL();
# undef FREE_ALL

    return 0;
}
