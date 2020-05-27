/* -*- mode: C -*-  */
/*
   IGraph library.
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

#include "igraph_operators.h"
#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_interrupt_internal.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_adjlist.h"
#include "igraph_attributes.h"
#include "igraph_conversion.h"
#include "igraph_qsort.h"
#include "config.h"
#include <limits.h>

/**
 * \function igraph_disjoint_union
 * \brief Creates the union of two disjoint graphs
 *
 * </para><para>
 * First the vertices of the second graph will be relabeled with new
 * vertex ids to have two disjoint sets of vertex ids, then the union
 * of the two graphs will be formed.
 * If the two graphs have |V1| and |V2| vertices and |E1| and |E2|
 * edges respectively then the new graph will have |V1|+|V2| vertices
 * and |E1|+|E2| edges.
 *
 * </para><para>
 * Both graphs need to have the same directedness, ie. either both
 * directed or both undirected.
 *
 * </para><para>
 * The current version of this function cannot handle graph, vertex
 * and edge attributes, they will be lost.
 *
 * \param res  Pointer to an uninitialized graph object, the result
 *        will stored here.
 * \param left The first graph.
 * \param right The second graph.
 * \return Error code.
 * \sa \ref igraph_disjoint_union_many() for creating the disjoint union
 * of more than two graphs, \ref igraph_union() for non-disjoint
 * union.
 *
 * Time complexity: O(|V1|+|V2|+|E1|+|E2|).
 *
 * \example examples/simple/igraph_disjoint_union.c
 */

int igraph_disjoint_union(igraph_t *res, const igraph_t *left,
                          const igraph_t *right) {

    long int no_of_nodes_left = igraph_vcount(left);
    long int no_of_nodes_right = igraph_vcount(right);
    long int no_of_edges_left = igraph_ecount(left);
    long int no_of_edges_right = igraph_ecount(right);
    igraph_vector_t edges;
    igraph_bool_t directed_left = igraph_is_directed(left);
    igraph_integer_t from, to;
    long int i;

    if (directed_left != igraph_is_directed(right)) {
        IGRAPH_ERROR("Cannot union directed and undirected graphs",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges,
                                       2 * (no_of_edges_left + no_of_edges_right)));
    for (i = 0; i < no_of_edges_left; i++) {
        igraph_edge(left, (igraph_integer_t) i, &from, &to);
        igraph_vector_push_back(&edges, from);
        igraph_vector_push_back(&edges, to);
    }
    for (i = 0; i < no_of_edges_right; i++) {
        igraph_edge(right, (igraph_integer_t) i, &from, &to);
        igraph_vector_push_back(&edges, from + no_of_nodes_left);
        igraph_vector_push_back(&edges, to + no_of_nodes_left);
    }

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t)
                               (no_of_nodes_left + no_of_nodes_right),
                               directed_left));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_disjoint_union_many
 * \brief The disjint union of many graphs.
 *
 * </para><para>
 * First the vertices in the graphs will be relabeled with new vertex
 * ids to have pairwise disjoint vertex id sets and then the union of
 * the graphs is formed.
 * The number of vertices and edges in the result is the total number
 * of vertices and edges in the graphs.
 *
 * </para><para>
 * Both graphs need to have the same directedness, ie. either both
 * directed or both undirected.
 *
 * </para><para>
 * The current version of this function cannot handle graph, vertex
 * and edge attributes, they will be lost.
 *
 * \param res Pointer to an uninitialized graph object, the result of
 *        the operation will be stored here.
 * \param graphs Pointer vector, contains pointers to initialized
 *        graph objects.
 * \return Error code.
 * \sa \ref igraph_disjoint_union() for an easier syntax if you have
 * only two graphs, \ref igraph_union_many() for non-disjoint union.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number
 * of edges in the result.
 */

int igraph_disjoint_union_many(igraph_t *res,
                               const igraph_vector_ptr_t *graphs) {
    long int no_of_graphs = igraph_vector_ptr_size(graphs);
    igraph_bool_t directed = 1;
    igraph_vector_t edges;
    long int no_of_edges = 0;
    long int shift = 0;
    igraph_t *graph;
    long int i, j;
    igraph_integer_t from, to;

    if (no_of_graphs != 0) {
        graph = VECTOR(*graphs)[0];
        directed = igraph_is_directed(graph);
        for (i = 0; i < no_of_graphs; i++) {
            graph = VECTOR(*graphs)[i];
            no_of_edges += igraph_ecount(graph);
            if (directed != igraph_is_directed(graph)) {
                IGRAPH_ERROR("Cannot union directed and undirected graphs",
                             IGRAPH_EINVAL);
            }
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, 2 * no_of_edges));

    for (i = 0; i < no_of_graphs; i++) {
        long int ec;
        graph = VECTOR(*graphs)[i];
        ec = igraph_ecount(graph);
        for (j = 0; j < ec; j++) {
            igraph_edge(graph, (igraph_integer_t) j, &from, &to);
            igraph_vector_push_back(&edges, from + shift);
            igraph_vector_push_back(&edges, to + shift);
        }
        shift += igraph_vcount(graph);
    }

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) shift, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

static int igraph_i_order_edgelist_cmp(void *edges, const void *e1, const void *e2) {
    igraph_vector_t *edgelist = edges;
    long int edge1 = (*(const long int*) e1) * 2;
    long int edge2 = (*(const long int*) e2) * 2;
    long int from1 = VECTOR(*edgelist)[edge1];
    long int from2 = VECTOR(*edgelist)[edge2];
    if (from1 < from2) {
        return -1;
    } else if (from1 > from2) {
        return 1;
    } else {
        long int to1 = VECTOR(*edgelist)[edge1 + 1];
        long int to2 = VECTOR(*edgelist)[edge2 + 1];
        if (to1 < to2) {
            return -1;
        } else if (to1 > to2) {
            return 1;
        } else {
            return 0;
        }
    }
}

#define IGRAPH_MODE_UNION        1
#define IGRAPH_MODE_INTERSECTION 2

static int igraph_i_merge(igraph_t *res, int mode,
                          const igraph_t *left, const igraph_t *right,
                          igraph_vector_t *edge_map1, igraph_vector_t *edge_map2) {

    long int no_of_nodes_left = igraph_vcount(left);
    long int no_of_nodes_right = igraph_vcount(right);
    long int no_of_nodes;
    long int no_edges_left = igraph_ecount(left);
    long int no_edges_right = igraph_ecount(right);
    igraph_bool_t directed = igraph_is_directed(left);
    igraph_vector_t edges;
    igraph_vector_t edges1, edges2;
    igraph_vector_long_t order1, order2;
    long int i, j, eptr = 0;
    long int idx1, idx2, edge1 = -1, edge2 = -1, from1 = -1, from2 = -1, to1 = -1, to2 = -1;
    igraph_bool_t l;

    if (directed != igraph_is_directed(right)) {
        IGRAPH_ERROR("Cannot make union or intersection of directed "
                     "and undirected graph", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges1, no_edges_left * 2);
    IGRAPH_VECTOR_INIT_FINALLY(&edges2, no_edges_right * 2);
    IGRAPH_CHECK(igraph_vector_long_init(&order1, no_edges_left));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &order1);
    IGRAPH_CHECK(igraph_vector_long_init(&order2, no_edges_right));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &order2);

    if (edge_map1) {
        switch (mode) {
        case IGRAPH_MODE_UNION:
            IGRAPH_CHECK(igraph_vector_resize(edge_map1, no_edges_left));
            break;
        case IGRAPH_MODE_INTERSECTION:
            igraph_vector_clear(edge_map1);
            break;
        }
    }
    if (edge_map2) {
        switch (mode) {
        case IGRAPH_MODE_UNION:
            IGRAPH_CHECK(igraph_vector_resize(edge_map2, no_edges_right));
            break;
        case IGRAPH_MODE_INTERSECTION:
            igraph_vector_clear(edge_map2);
            break;
        }
    }

    no_of_nodes = no_of_nodes_left > no_of_nodes_right ?
                  no_of_nodes_left : no_of_nodes_right;

    /* We merge the two edge lists. We need to sort them first.
       For undirected graphs, we also need to make sure that
       for every edge, that larger (non-smaller) vertex id is in the
       second column. */

    IGRAPH_CHECK(igraph_get_edgelist(left, &edges1, /*bycol=*/ 0));
    IGRAPH_CHECK(igraph_get_edgelist(right, &edges2, /*bycol=*/ 0));
    if (!directed) {
        for (i = 0, j = 0; i < no_edges_left; i++, j += 2) {
            if (VECTOR(edges1)[j] > VECTOR(edges1)[j + 1]) {
                long int tmp = VECTOR(edges1)[j];
                VECTOR(edges1)[j] = VECTOR(edges1)[j + 1];
                VECTOR(edges1)[j + 1] = tmp;
            }
        }
        for (i = 0, j = 0; i < no_edges_right; i++, j += 2) {
            if (VECTOR(edges2)[j] > VECTOR(edges2)[j + 1]) {
                long int tmp = VECTOR(edges2)[j];
                VECTOR(edges2)[j] = VECTOR(edges2)[j + 1];
                VECTOR(edges2)[j + 1] = tmp;
            }
        }
    }

    for (i = 0; i < no_edges_left; i++) {
        VECTOR(order1)[i] = i;
    }
    for (i = 0; i < no_edges_right; i++) {
        VECTOR(order2)[i] = i;
    }

    igraph_qsort_r(VECTOR(order1), no_edges_left, sizeof(VECTOR(order1)[0]),
                   &edges1, igraph_i_order_edgelist_cmp);
    igraph_qsort_r(VECTOR(order2), no_edges_right, sizeof(VECTOR(order2)[0]),
                   &edges2, igraph_i_order_edgelist_cmp);

#define INC1() if ( (++idx1) < no_edges_left) {          \
        edge1 = VECTOR(order1)[idx1];                \
        from1 = VECTOR(edges1)[2*edge1];                 \
        to1 = VECTOR(edges1)[2*edge1+1];                 \
    }
#define INC2() if ( (++idx2) < no_edges_right) {         \
        edge2 = VECTOR(order2)[idx2];                \
        from2 = VECTOR(edges2)[2*edge2];                 \
        to2 = VECTOR(edges2)[2*edge2+1];                 \
    }

    idx1 = idx2 = -1;
    INC1();
    INC2();

#define CONT() switch (mode) {              \
    case IGRAPH_MODE_UNION:                \
        l = idx1 < no_edges_left || idx2 < no_edges_right;   \
        break;                       \
    case IGRAPH_MODE_INTERSECTION:             \
        l = idx1 < no_edges_left && idx2 < no_edges_right;   \
        break;                       \
    }

    CONT();
    while (l) {
        if (idx2 >= no_edges_right ||
            (idx1 < no_edges_left && from1 < from2) ||
            (idx1 < no_edges_left && from1 == from2 && to1 < to2)) {
            /* Edge from first graph */
            if (mode == IGRAPH_MODE_UNION) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, from1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to1));
                if (edge_map1) {
                    VECTOR(*edge_map1)[edge1] = eptr;
                }
                eptr++;
            }
            INC1();
        } else if (idx1 >= no_edges_left ||
                   (idx2 < no_edges_right && from2 < from1) ||
                   (idx2 < no_edges_right && from1 == from2 && to2 < to1)) {
            /* Edge from second graph */
            if (mode == IGRAPH_MODE_UNION) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, from2));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, to2));
                if (edge_map2) {
                    VECTOR(*edge_map2)[edge2] = eptr;
                }
                eptr++;
            }
            INC2();
        } else {
            /* Edge from both */
            IGRAPH_CHECK(igraph_vector_push_back(&edges, from1));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, to1));
            if (mode == IGRAPH_MODE_UNION) {
                if (edge_map1) {
                    VECTOR(*edge_map1)[edge1] = eptr;
                }
                if (edge_map2) {
                    VECTOR(*edge_map2)[edge2] = eptr;
                }
            } else if (mode == IGRAPH_MODE_INTERSECTION) {
                if (edge_map1) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map1, edge1));
                }
                if (edge_map2) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map2, edge2));
                }
            }
            eptr++;
            INC1();
            INC2();
        }
        CONT();
    }

#undef INC1
#undef INC2

    igraph_vector_long_destroy(&order2);
    igraph_vector_long_destroy(&order1);
    igraph_vector_destroy(&edges2);
    igraph_vector_destroy(&edges1);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \function igraph_intersection
 * \brief Collect the common edges from two graphs.
 *
 * </para><para>
 * The result graph contains only edges present both in the first and
 * the second graph. The number of vertices in the result graph is the
 * same as the larger from the two arguments.
 *
 * \param res Pointer to an uninitialized graph object. This will
 * contain the result of the operation.
 * \param left The first operand, a graph object.
 * \param right The second operand, a graph object.
 * \param edge_map1 Null pointer, or an initialized \type igraph_vector_t.
 *    If the latter, then a mapping from the edges of the result graph, to
 *    the edges of the \p left input graph is stored here.
 * \param edge_map2 Null pointer, or an \type igraph_vector_t. The same
 *    as \p edge_map1, but for the \p right input graph.
 * \return Error code.
 * \sa \ref igraph_intersection_many() to calculate the intersection
 * of many graphs at once, \ref igraph_union(), \ref
 * igraph_difference() for other operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of nodes, |E|
 * is the number of edges in the smaller graph of the two. (The one
 * containing less vertices is considered smaller.)
 *
 * \example examples/simple/igraph_intersection.c
 */

int igraph_intersection(igraph_t *res,
                        const igraph_t *left, const igraph_t *right,
                        igraph_vector_t *edge_map1,
                        igraph_vector_t *edge_map2) {
    return igraph_i_merge(res, IGRAPH_MODE_INTERSECTION, left, right,
                          edge_map1, edge_map2);
}

static void igraph_i_union_many_free(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_destroy(VECTOR(*v)[i]);
            igraph_Free(VECTOR(*v)[i]);
        }
    }
    igraph_vector_ptr_destroy(v);
}

static void igraph_i_union_many_free2(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_long_destroy(VECTOR(*v)[i]);
            igraph_Free(VECTOR(*v)[i]);
        }
    }
    igraph_vector_ptr_destroy(v);
}

static void igraph_i_union_many_free3(igraph_vector_ptr_t *v) {
    long int i, n = igraph_vector_ptr_size(v);
    for (i = 0; i < n; i++) {
        if (VECTOR(*v)[i] != 0) {
            igraph_vector_destroy(VECTOR(*v)[i]);
            igraph_Free(VECTOR(*v)[i]);
        }
    }
}

/**
 * \function igraph_intersection_many
 * \brief The intersection of more than two graphs.
 *
 * </para><para>
 * This function calculates the intersection of the graphs stored in
 * the \c graphs argument. Only those edges will be included in the
 * result graph which are part of every graph in \c graphs.
 *
 * </para><para>
 * The number of vertices in the result graph will be the maximum
 * number of vertices in the argument graphs.
 *
 * \param res Pointer to an uninitialized graph object, the result of
 *        the operation will be stored here.
 * \param graphs Pointer vector, contains pointers to graphs objects,
 *        the operands of the intersection operator.
 * \param edgemaps If not a null pointer, then it must be an initialized
 *        pointer vector and the mappings of edges from the graphs to the
 *        result graph will be stored here, in the same order as
 *        \p graphs. Each mapping is stored in a separate
 *        \type igraph_vector_t object. For the edges that are not in
 *        the intersection, -1 is stored.
 * \return Error code.
 * \sa \ref igraph_intersection() for the intersection of two graphs,
 * \ref igraph_union_many(), \ref igraph_union() and \ref
 * igraph_difference() for other operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges in the smallest graph (ie. the graph having
 * the less vertices).
 */

int igraph_intersection_many(igraph_t *res,
                             const igraph_vector_ptr_t *graphs,
                             igraph_vector_ptr_t *edgemaps) {

    long int no_of_graphs = igraph_vector_ptr_size(graphs);
    long int no_of_nodes = 0;
    igraph_bool_t directed = 1;
    igraph_vector_t edges;
    igraph_vector_ptr_t edge_vects, order_vects;
    long int i, j, tailfrom = no_of_graphs > 0 ? 0 : -1, tailto = -1;
    igraph_vector_long_t no_edges;
    igraph_bool_t allne = no_of_graphs == 0 ? 0 : 1, allsame = 0;
    long int idx = 0;

    /* Check directedness */
    if (no_of_graphs != 0) {
        directed = igraph_is_directed(VECTOR(*graphs)[0]);
    }
    for (i = 1; i < no_of_graphs; i++) {
        if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
            IGRAPH_ERROR("Cannot intersect directed and undirected graphs",
                         IGRAPH_EINVAL);
        }
    }

    if (edgemaps) {
        IGRAPH_CHECK(igraph_vector_ptr_resize(edgemaps, no_of_graphs));
        igraph_vector_ptr_null(edgemaps);
        IGRAPH_FINALLY(igraph_i_union_many_free3, edgemaps);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_long_init(&no_edges, no_of_graphs));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &no_edges);

    /* Calculate number of nodes, query number of edges */
    for (i = 0; i < no_of_graphs; i++) {
        long int n = igraph_vcount(VECTOR(*graphs)[i]);
        if (n > no_of_nodes) {
            no_of_nodes = n;
        }
        VECTOR(no_edges)[i] = igraph_ecount(VECTOR(*graphs)[i]);
        allne = allne && VECTOR(no_edges)[i] > 0;
    }

    if (edgemaps) {
        for (i = 0; i < no_of_graphs; i++) {
            VECTOR(*edgemaps)[i] = igraph_Calloc(1, igraph_vector_t);
            if (!VECTOR(*edgemaps)[i]) {
                IGRAPH_ERROR("Cannot intersect graphs", IGRAPH_ENOMEM);
            }
            IGRAPH_CHECK(igraph_vector_init(VECTOR(*edgemaps)[i],
                                            VECTOR(no_edges)[i]));
            igraph_vector_fill(VECTOR(*edgemaps)[i], -1);
        }
    }

    /* Allocate memory for the edge lists and their index vectors */
    if (no_of_graphs != 0) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&edge_vects, no_of_graphs));
        IGRAPH_FINALLY(igraph_i_union_many_free, &edge_vects);
        IGRAPH_CHECK(igraph_vector_ptr_init(&order_vects, no_of_graphs));
        IGRAPH_FINALLY(igraph_i_union_many_free2, &order_vects);
    }
    for (i = 0; i < no_of_graphs; i++) {
        VECTOR(edge_vects)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(order_vects)[i] = igraph_Calloc(1, igraph_vector_long_t);
        if (! VECTOR(edge_vects)[i] || ! VECTOR(order_vects)[i]) {
            IGRAPH_ERROR("Cannot intersect graphs", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(edge_vects)[i],
                                        2 * VECTOR(no_edges)[i]));
        IGRAPH_CHECK(igraph_vector_long_init(VECTOR(order_vects)[i],
                                             VECTOR(no_edges)[i]));
    }

    /* Query and sort the edge lists */
    for (i = 0; i < no_of_graphs; i++) {
        long int k, j, n = VECTOR(no_edges)[i];
        igraph_vector_t *edges = VECTOR(edge_vects)[i];
        igraph_vector_long_t *order = VECTOR(order_vects)[i];
        IGRAPH_CHECK(igraph_get_edgelist(VECTOR(*graphs)[i], edges, /*bycol=*/0));
        if (!directed) {
            for (k = 0, j = 0; k < n; k++, j += 2) {
                if (VECTOR(*edges)[j] > VECTOR(*edges)[j + 1]) {
                    long int tmp = VECTOR(*edges)[j];
                    VECTOR(*edges)[j] = VECTOR(*edges)[j + 1];
                    VECTOR(*edges)[j + 1] = tmp;
                }
            }
        }
        for (k = 0; k < n; k++) {
            VECTOR(*order)[k] = k;
        }
        igraph_qsort_r(VECTOR(*order), n, sizeof(VECTOR(*order)[0]), edges,
                       igraph_i_order_edgelist_cmp);
    }

    /* Do the merge. We work from the end of the edge lists,
       because then we don't have to keep track of where we are right
       now in the edge and order lists. We find the "largest" edge,
       and if it is present in all graphs, then we copy it to the
       result. We remove all instances of this edge.  */

    while (allne) {

        /* Look for the smallest tail element */
        for (j = 0, tailfrom = LONG_MAX, tailto = LONG_MAX; j < no_of_graphs; j++) {
            long int edge = igraph_vector_long_tail(VECTOR(order_vects)[j]);
            igraph_vector_t *ev = VECTOR(edge_vects)[j];
            long int from = VECTOR(*ev)[2 * edge];
            long int to = VECTOR(*ev)[2 * edge + 1];
            if (from < tailfrom || (from == tailfrom && to < tailto)) {
                tailfrom = from; tailto = to;
            }
        }

        /* OK, now remove all elements from the tail(s) that are bigger
           than the smallest tail element. */
        for (j = 0, allsame = 1; j < no_of_graphs; j++) {
            long int from = -1, to = -1;
            while (1) {
                long int edge = igraph_vector_long_tail(VECTOR(order_vects)[j]);
                igraph_vector_t *ev = VECTOR(edge_vects)[j];
                from = VECTOR(*ev)[2 * edge];
                to = VECTOR(*ev)[2 * edge + 1];
                if (from > tailfrom || (from == tailfrom && to > tailto)) {
                    igraph_vector_long_pop_back(VECTOR(order_vects)[j]);
                    if (igraph_vector_long_empty(VECTOR(order_vects)[j])) {
                        allne = 0;
                        break;
                    }
                } else {
                    break;
                }
            }
            if (from != tailfrom || to != tailto) {
                allsame = 0;
            }
        }

        /* Add the edge, if the smallest tail element was present
           in all graphs. */
        if (allsame) {
            IGRAPH_CHECK(igraph_vector_push_back(&edges, tailfrom));
            IGRAPH_CHECK(igraph_vector_push_back(&edges, tailto));
        }

        /* Drop edges matching the smalles tail elements
           from the order vectors, build edge maps */
        if (allne) {
            for (j = 0; j < no_of_graphs; j++) {
                long int edge = igraph_vector_long_tail(VECTOR(order_vects)[j]);
                igraph_vector_t *ev = VECTOR(edge_vects)[j];
                long int from = VECTOR(*ev)[2 * edge];
                long int to = VECTOR(*ev)[2 * edge + 1];
                if (from == tailfrom && to == tailto) {
                    igraph_vector_long_pop_back(VECTOR(order_vects)[j]);
                    if (igraph_vector_long_empty(VECTOR(order_vects)[j])) {
                        allne = 0;
                    }
                    if (edgemaps && allsame) {
                        igraph_vector_t *map = VECTOR(*edgemaps)[j];
                        VECTOR(*map)[edge] = idx;
                    }
                }
            }
            if (allsame) {
                idx++;
            }
        }

    } /* while allne */

    if (no_of_graphs > 0) {
        igraph_i_union_many_free2(&order_vects);
        igraph_i_union_many_free(&edge_vects);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_long_destroy(&no_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    if (edgemaps) {
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_union
 * \brief Calculates the union of two graphs.
 *
 * </para><para>
 * The number of vertices in the result is that of the larger graph
 * from the two arguments. The result graph contains edges which are
 * present in at least one of the operand graphs.
 *
 * \param res Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param left The first graph.
 * \param right The second graph.
 * \param edge_map1 Pointer to an initialized vector or a null pointer.
 *     If not a null pointer, it will contain a mapping from the edges
 *     of the first argument graph (\p left) to the edges of the
 *     result graph.
 * \param edge_map2 The same as \p edge_map1, but for the second
 *     graph, \p right.
 * \return Error code.
 * \sa \ref igraph_union_many() for the union of many graphs,
 * \ref igraph_intersection() and \ref igraph_difference() for other
 * operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of
 * vertices, |E| the number of edges in the result graph.
 *
 * \example examples/simple/igraph_union.c
 */

int igraph_union(igraph_t *res,
                 const igraph_t *left, const igraph_t *right,
                 igraph_vector_t *edge_map1, igraph_vector_t *edge_map2) {
    return igraph_i_merge(res, IGRAPH_MODE_UNION, left, right,
                          edge_map1, edge_map2);
}

/**
 * \function igraph_union_many
 * \brief Creates the union of many graphs.
 *
 * </para><para>
 * The result graph will contain as many vertices as the largest graph
 * among the arguments does, and an edge will be included in it if it
 * is part of at least one operand graph.
 *
 * </para><para>
 * The directedness of the operand graphs must be the same.
 *
 * \param res Pointer to an uninitialized graph object, this will
 *        contain the result.
 * \param graphs Pointer vector, contains pointers to the operands of
 *        the union operator, graph objects of course.
 * \param edgemaps If not a null pointer, then it must be an initialized
 *        pointer vector and the mappings of edges from the graphs to the
 *        result graph will be stored here, in the same order as
 *        \p graphs. Each mapping is stored in a separate
 *        \type igraph_vector_t object.
 * \return Error code.
 * \sa \ref igraph_union() for the union of two graphs, \ref
 * igraph_intersection_many(), \ref igraph_intersection() and \ref
 * igraph_difference for other operators.
 *
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices
 * in largest graph and |E| is the number of edges in the result graph.
 *
 * \example examples/simple/igraph_union.c
 */

int igraph_union_many(igraph_t *res, const igraph_vector_ptr_t *graphs,
                      igraph_vector_ptr_t *edgemaps) {

    long int no_of_graphs = igraph_vector_ptr_size(graphs);
    long int no_of_nodes = 0;
    igraph_bool_t directed = 1;
    igraph_vector_t edges;
    igraph_vector_ptr_t edge_vects, order_vects;
    igraph_vector_long_t no_edges;
    long int i, j, tailfrom = no_of_graphs > 0 ? 0 : -1, tailto = -1;
    long int idx = 0;

    /* Check directedness */
    if (no_of_graphs != 0) {
        directed = igraph_is_directed(VECTOR(*graphs)[0]);
        no_of_nodes = igraph_vcount(VECTOR(*graphs)[0]);
    }
    for (i = 1; i < no_of_graphs; i++) {
        if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
            IGRAPH_ERROR("Cannot union directed and undirected graphs",
                         IGRAPH_EINVAL);
        }
    }

    if (edgemaps) {
        IGRAPH_CHECK(igraph_vector_ptr_resize(edgemaps, no_of_graphs));
        igraph_vector_ptr_null(edgemaps);
        IGRAPH_FINALLY(igraph_i_union_many_free3, edgemaps);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_long_init(&no_edges, no_of_graphs));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &no_edges);

    /* Calculate number of nodes, query number of edges */
    for (i = 0; i < no_of_graphs; i++) {
        long int n = igraph_vcount(VECTOR(*graphs)[i]);
        if (n > no_of_nodes) {
            no_of_nodes = n;
        }
        VECTOR(no_edges)[i] = igraph_ecount(VECTOR(*graphs)[i]);
    }

    if (edgemaps) {
        for (i = 0; i < no_of_graphs; i++) {
            VECTOR(*edgemaps)[i] = igraph_Calloc(1, igraph_vector_t);
            if (!VECTOR(*edgemaps)[i]) {
                IGRAPH_ERROR("Cannot union graphs", IGRAPH_ENOMEM);
            }
            IGRAPH_CHECK(igraph_vector_init(VECTOR(*edgemaps)[i],
                                            VECTOR(no_edges)[i]));
        }
    }

    /* Allocate memory for the edge lists and their index vectors */
    if (no_of_graphs != 0) {
        IGRAPH_CHECK(igraph_vector_ptr_init(&edge_vects, no_of_graphs));
        IGRAPH_FINALLY(igraph_i_union_many_free, &edge_vects);
        IGRAPH_CHECK(igraph_vector_ptr_init(&order_vects, no_of_graphs));
        IGRAPH_FINALLY(igraph_i_union_many_free2, &order_vects);
    }
    for (i = 0; i < no_of_graphs; i++) {
        VECTOR(edge_vects)[i] = igraph_Calloc(1, igraph_vector_t);
        VECTOR(order_vects)[i] = igraph_Calloc(1, igraph_vector_long_t);
        if (! VECTOR(edge_vects)[i] || ! VECTOR(order_vects)[i]) {
            IGRAPH_ERROR("Cannot union graphs", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_vector_init(VECTOR(edge_vects)[i],
                                        2 * VECTOR(no_edges)[i]));
        IGRAPH_CHECK(igraph_vector_long_init(VECTOR(order_vects)[i],
                                             VECTOR(no_edges)[i]));
    }

    /* Query and sort the edge lists */
    for (i = 0; i < no_of_graphs; i++) {
        long int k, j, n = VECTOR(no_edges)[i];
        igraph_vector_t *edges = VECTOR(edge_vects)[i];
        igraph_vector_long_t *order = VECTOR(order_vects)[i];
        IGRAPH_CHECK(igraph_get_edgelist(VECTOR(*graphs)[i], edges, /*bycol=*/0));
        if (!directed) {
            for (k = 0, j = 0; k < n; k++, j += 2) {
                if (VECTOR(*edges)[j] > VECTOR(*edges)[j + 1]) {
                    long int tmp = VECTOR(*edges)[j];
                    VECTOR(*edges)[j] = VECTOR(*edges)[j + 1];
                    VECTOR(*edges)[j + 1] = tmp;
                }
            }
        }
        for (k = 0; k < n; k++) {
            VECTOR(*order)[k] = k;
        }
        igraph_qsort_r(VECTOR(*order), n, sizeof(VECTOR(*order)[0]), edges,
                       igraph_i_order_edgelist_cmp);
    }

    while (tailfrom >= 0) {

        /* Get the largest tail element */
        tailfrom = tailto = -1;
        for (j = 0; j < no_of_graphs; j++) {
            if (!igraph_vector_long_empty(VECTOR(order_vects)[j])) {
                long int edge = igraph_vector_long_tail(VECTOR(order_vects)[j]);
                igraph_vector_t *ev = VECTOR(edge_vects)[j];
                long int from = VECTOR(*ev)[2 * edge];
                long int to = VECTOR(*ev)[2 * edge + 1];
                if (from > tailfrom || (from == tailfrom && to > tailto)) {
                    tailfrom = from; tailto = to;
                }
            }
        }
        if (tailfrom < 0) {
            continue;
        }

        /* add the edge */
        IGRAPH_CHECK(igraph_vector_push_back(&edges, tailfrom));
        IGRAPH_CHECK(igraph_vector_push_back(&edges, tailto));

        /* update edge lists, we just modify the 'order' vectors */
        for (j = 0; j < no_of_graphs; j++) {
            if (!igraph_vector_long_empty(VECTOR(order_vects)[j])) {
                long int edge = igraph_vector_long_tail(VECTOR(order_vects)[j]);
                igraph_vector_t *ev = VECTOR(edge_vects)[j];
                long int from = VECTOR(*ev)[2 * edge];
                long int to = VECTOR(*ev)[2 * edge + 1];
                if (from == tailfrom && to == tailto) {
                    igraph_vector_long_pop_back(VECTOR(order_vects)[j]);
                    if (edgemaps) {
                        igraph_vector_t *map = VECTOR(*edgemaps)[j];
                        VECTOR(*map)[edge] = idx;
                    }
                }
            }
        }
        idx++;

    }

    if (no_of_graphs > 0) {
        igraph_i_union_many_free2(&order_vects);
        igraph_i_union_many_free(&edge_vects);
        IGRAPH_FINALLY_CLEAN(2);
    }

    igraph_vector_long_destroy(&no_edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    if (edgemaps) {
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_difference
 * \brief Calculate the difference of two graphs
 *
 * </para><para>
 * The number of vertices in the result is the number of vertices in
 * the original graph, ie. the left, first operand. In the results
 * graph only edges will be included from \c orig which are not
 * present in \c sub.
 *
 * \param res Pointer to an uninitialized graph object, the result
 * will be stored here.
 * \param orig The left operand of the operator, a graph object.
 * \param sub The right operand of the operator, a graph object.
 * \return Error code.
 * \sa \ref igraph_intersection() and \ref igraph_union() for other
 * operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number vertices in
 * the smaller graph, |E| is the
 * number of edges in the result graph.
 *
 * \example examples/simple/igraph_difference.c
 */

int igraph_difference(igraph_t *res,
                      const igraph_t *orig, const igraph_t *sub) {

    /* Quite nasty, but we will use that an edge adjacency list
       contains the vertices according to the order of the
       vertex ids at the "other" end of the edge. */

    long int no_of_nodes_orig = igraph_vcount(orig);
    long int no_of_nodes_sub = igraph_vcount(sub);
    long int no_of_nodes = no_of_nodes_orig;
    long int smaller_nodes;
    igraph_bool_t directed = igraph_is_directed(orig);
    igraph_vector_t edges;
    igraph_vector_t edge_ids;
    igraph_vector_int_t *nei1, *nei2;
    igraph_inclist_t inc_orig, inc_sub;
    long int i;
    igraph_integer_t v1, v2;

    if (directed != igraph_is_directed(sub)) {
        IGRAPH_ERROR("Cannot subtract directed and undirected graphs",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edge_ids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_inclist_init(orig, &inc_orig, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inc_orig);
    IGRAPH_CHECK(igraph_inclist_init(sub, &inc_sub, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inc_sub);

    smaller_nodes = no_of_nodes_orig > no_of_nodes_sub ?
                    no_of_nodes_sub : no_of_nodes_orig;

    for (i = 0; i < smaller_nodes; i++) {
        long int n1, n2, e1, e2;
        IGRAPH_ALLOW_INTERRUPTION();
        nei1 = igraph_inclist_get(&inc_orig, i);
        nei2 = igraph_inclist_get(&inc_sub, i);
        n1 = igraph_vector_int_size(nei1) - 1;
        n2 = igraph_vector_int_size(nei2) - 1;
        while (n1 >= 0 && n2 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            e2 = (long int) VECTOR(*nei2)[n2];
            v1 = IGRAPH_OTHER(orig, e1, i);
            v2 = IGRAPH_OTHER(sub, e2, i);

            if (!directed && v1 < i) {
                n1--;
            } else if (!directed && v2 < i) {
                n2--;
            } else if (v1 > v2) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
                n1--;
            } else if (v2 > v1) {
                n2--;
            } else {
                n1--;
                n2--;
            }
        }

        /* Copy remaining edges */
        while (n1 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            v1 = IGRAPH_OTHER(orig, e1, i);
            if (directed || v1 >= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
            }
            n1--;
        }
    }

    /* copy remaining edges, use the previous value of 'i' */
    for (; i < no_of_nodes_orig; i++) {
        long int n1, e1;
        nei1 = igraph_inclist_get(&inc_orig, i);
        n1 = igraph_vector_int_size(nei1) - 1;
        while (n1 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            v1 = IGRAPH_OTHER(orig, e1, i);
            if (directed || v1 >= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
            }
            n1--;
        }
    }

    igraph_inclist_destroy(&inc_sub);
    igraph_inclist_destroy(&inc_orig);
    IGRAPH_FINALLY_CLEAN(2);
    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    /* Attributes */
    if (orig->attr) {
        IGRAPH_I_ATTRIBUTE_DESTROY(res);
        IGRAPH_I_ATTRIBUTE_COPY(res, orig, /*graph=*/1, /*vertex=*/1, /*edge=*/0);
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(orig, res, &edge_ids));
    }

    igraph_vector_destroy(&edge_ids);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_complementer
 * \brief Create the complementer of a graph
 *
 * </para><para>The complementer graph means that all edges which are
 * not part of the original graph will be included in the result.
 *
 * \param res Pointer to an uninitialized graph object.
 * \param graph The original graph.
 * \param loops Whether to add loop edges to the complementer graph.
 * \return Error code.
 * \sa \ref igraph_union(), \ref igraph_intersection() and \ref
 * igraph_difference().
 *
 * Time complexity: O(|V|+|E1|+|E2|), |V| is the number of
 * vertices in the graph, |E1| is the number of edges in the original
 * and |E2| in the complementer graph.
 *
 * \example examples/simple/igraph_complementer.c
 */

int igraph_complementer(igraph_t *res, const igraph_t *graph,
                        igraph_bool_t loops) {

    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t edges;
    igraph_vector_t neis;
    long int i, j;
    long int zero = 0, *limit;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

    if (igraph_is_directed(graph)) {
        limit = &zero;
    } else {
        limit = &i;
    }

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
                                      IGRAPH_OUT));
        if (loops) {
            for (j = no_of_nodes - 1; j >= *limit; j--) {
                if (igraph_vector_empty(&neis) || j > igraph_vector_tail(&neis)) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                } else {
                    igraph_vector_pop_back(&neis);
                }
            }
        } else {
            for (j = no_of_nodes - 1; j >= *limit; j--) {
                if (igraph_vector_empty(&neis) || j > igraph_vector_tail(&neis)) {
                    if (i != j) {
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, j));
                    }
                } else {
                    igraph_vector_pop_back(&neis);
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&neis);
    IGRAPH_I_ATTRIBUTE_DESTROY(res);
    IGRAPH_I_ATTRIBUTE_COPY(res, graph, /*graph=*/1, /*vertex=*/1, /*edge=*/0);
    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \function igraph_compose
 * \brief Calculates the composition of two graphs
 *
 * The composition of graphs contains the same number of vertices as
 * the bigger graph of the two operands. It contains an (i,j) edge if
 * and only if there is a k vertex, such that the first graphs
 * contains an (i,k) edge and the second graph a (k,j) edge.
 *
 * </para><para>This is of course exactly the composition of two
 * binary relations.
 *
 * </para><para>Two two graphs must have the same directedness,
 * otherwise the function returns with an error message.
 * Note that for undirected graphs the two relations are by definition
 * symmetric.
 *
 * \param res Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param g1 The firs operand, a graph object.
 * \param g2 The second operand, another graph object.
 * \param edge_map1 If not a null pointer, then it must be a pointer
 *        to an initialized vector, and a mapping from the edges of
 *        the result graph to the edges of the first graph is stored
 *        here.
 * \param edge_map1 If not a null pointer, then it must be a pointer
 *        to an initialized vector, and a mapping from the edges of
 *        the result graph to the edges of the second graph is stored
 *        here.
 * \return Error code.
 *
 * Time complexity: O(|V|*d1*d2), |V| is the number of vertices in the
 * first graph, d1 and d2 the average degree in the first and second
 * graphs.
 *
 * \example examples/simple/igraph_compose.c
 */

int igraph_compose(igraph_t *res, const igraph_t *g1, const igraph_t *g2,
                   igraph_vector_t *edge_map1, igraph_vector_t *edge_map2) {

    long int no_of_nodes_left = igraph_vcount(g1);
    long int no_of_nodes_right = igraph_vcount(g2);
    long int no_of_nodes;
    igraph_bool_t directed = igraph_is_directed(g1);
    igraph_vector_t edges;
    igraph_vector_t neis1, neis2;
    long int i;

    if (directed != igraph_is_directed(g2)) {
        IGRAPH_ERROR("Cannot compose directed and undirected graph",
                     IGRAPH_EINVAL);
    }

    no_of_nodes = no_of_nodes_left > no_of_nodes_right ?
                  no_of_nodes_left : no_of_nodes_right;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis2, 0);

    if (edge_map1) {
        igraph_vector_clear(edge_map1);
    }
    if (edge_map2) {
        igraph_vector_clear(edge_map2);
    }

    for (i = 0; i < no_of_nodes_left; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_incident(g1, &neis1, (igraph_integer_t) i,
                                     IGRAPH_OUT));
        while (!igraph_vector_empty(&neis1)) {
            long int con = (long int) igraph_vector_pop_back(&neis1);
            long int v1 = IGRAPH_OTHER(g1, con, i);
            if (v1 < no_of_nodes_right) {
                IGRAPH_CHECK(igraph_incident(g2, &neis2, (igraph_integer_t) v1,
                                             IGRAPH_OUT));
            } else {
                continue;
            }
            while (!igraph_vector_empty(&neis2)) {
                long int con2 = igraph_vector_pop_back(&neis2);
                long int v2 = IGRAPH_OTHER(g2, con2, v1);
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v2));
                if (edge_map1) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map1, con));
                }
                if (edge_map2) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map2, con2));
                }
            }
        }
    }

    igraph_vector_destroy(&neis1);
    igraph_vector_destroy(&neis2);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
