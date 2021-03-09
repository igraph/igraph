/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020 The igraph development team

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

#include <limits.h>

#include "igraph_operators.h"

#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"

#include "operators/misc_internal.h"

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
    return igraph_i_merge(res, IGRAPH_MERGE_MODE_UNION, left, right,
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
 * If the graph list has length zero, the result will be a \em directed
 * graph with no vertices.
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
        IGRAPH_FINALLY(igraph_i_union_intersection_destroy_vectors, edgemaps);
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
            VECTOR(*edgemaps)[i] = IGRAPH_CALLOC(1, igraph_vector_t);
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
        IGRAPH_FINALLY(igraph_i_union_intersection_destroy_vectors, &edge_vects);
        IGRAPH_CHECK(igraph_vector_ptr_init(&order_vects, no_of_graphs));
        IGRAPH_FINALLY(igraph_i_union_intersection_destroy_vector_longs, &order_vects);
    }
    for (i = 0; i < no_of_graphs; i++) {
        VECTOR(edge_vects)[i] = IGRAPH_CALLOC(1, igraph_vector_t);
        VECTOR(order_vects)[i] = IGRAPH_CALLOC(1, igraph_vector_long_t);
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
        igraph_i_union_intersection_destroy_vector_longs(&order_vects);
        igraph_i_union_intersection_destroy_vectors(&edge_vects);
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
