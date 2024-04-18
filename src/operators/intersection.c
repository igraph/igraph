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

#include "igraph_operators.h"

#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_interface.h"
#include "igraph_qsort.h"
#include "igraph_vector_list.h"

#include "operators/misc_internal.h"

/**
 * \function igraph_intersection
 * \brief Collect the common edges from two graphs.
 *
 * The result graph contains only edges present both in the first and
 * the second graph. The number of vertices in the result graph is the
 * same as the larger from the two arguments.
 *
 * </para><para>
 * The directedness of the operand graphs must be the same.
 *
 * </para><para>
 * Edge multiplicities are handled by taking the \em smaller of the two
 * multiplicities in the input graphs. In other words, if the first graph
 * has N edges between a vertex pair (u, v) and the second graph has M edges,
 * the result graph will have min(N, M) edges between them.
 *
 * \param res Pointer to an uninitialized graph object. This will
 * contain the result of the operation.
 * \param left The first operand, a graph object.
 * \param right The second operand, a graph object.
 * \param edge_map1 Null pointer, or an initialized vector.
 *    If the latter, then a mapping from the edges of the result graph, to
 *    the edges of the \p left input graph is stored here. For the edges that
 *    are not in the intersection, -1 is stored.
 * \param edge_map2 Null pointer, or an initialized vector. The same
 *    as \p edge_map1, but for the \p right input graph. For the edges that
 *    are not in the intersection, -1 is stored.
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
igraph_error_t igraph_intersection(
        igraph_t *res,
        const igraph_t *left, const igraph_t *right,
        igraph_vector_int_t *edge_map1, igraph_vector_int_t *edge_map2) {
    return igraph_i_merge(res, IGRAPH_MERGE_MODE_INTERSECTION, left, right,
                          edge_map1, edge_map2);
}

/**
 * \function igraph_intersection_many
 * \brief The intersection of more than two graphs.
 *
 * This function calculates the intersection of the graphs stored in
 * the \p graphs argument. Only those edges will be included in the
 * result graph which are part of every graph in \p graphs.
 *
 * </para><para>
 * The number of vertices in the result graph will be the maximum
 * number of vertices in the argument graphs.
 *
 * </para><para>
 * The directedness of the argument graphs must be the same.
 * If the graph list has length zero, the result will be a \em directed
 * graph with no vertices.
 *
 * </para><para>
 * Edge multiplicities are handled by taking the \em minimum multiplicity of the
 * all multiplicities for the same vertex pair (u, v) in the input graphs; this
 * will be the multiplicity of (u, v) in the result graph.
 *
 * \param res Pointer to an uninitialized graph object, the result of
 *        the operation will be stored here.
 * \param graphs Pointer vector, contains pointers to graphs objects,
 *        the operands of the intersection operator.
 * \param edgemaps If not a null pointer, then it must be an initialized
 *        list of integer vectors, and the mappings of edges from the graphs to
 *        the result graph will be stored here, in the same order as
 *        \p graphs. Each mapping is stored in a separate
 *        \type igraph_vector_int_t object. For the edges that are not in
 *        the intersection, -1 is stored.
 * \return Error code.
 * \sa \ref igraph_intersection() for the intersection of two graphs,
 * \ref igraph_union_many(), \ref igraph_union() and \ref
 * igraph_difference() for other operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number of vertices,
 * |E| is the number of edges in the smallest graph (i.e. the graph having
 * the less vertices).
 */
igraph_error_t igraph_intersection_many(
    igraph_t *res, const igraph_vector_ptr_t *graphs,
    igraph_vector_int_list_t *edgemaps
) {

    igraph_integer_t no_of_graphs = igraph_vector_ptr_size(graphs);
    igraph_integer_t no_of_nodes = 0;
    igraph_bool_t directed = true;
    igraph_vector_int_t edges;
    igraph_vector_int_list_t edge_vects, order_vects;
    igraph_integer_t i, j, tailfrom = no_of_graphs > 0 ? 0 : -1, tailto = -1;
    igraph_vector_int_t no_edges;
    igraph_bool_t allne = no_of_graphs > 0;
    igraph_bool_t allsame = false;
    igraph_integer_t idx = 0;

    /* Check directedness */
    if (no_of_graphs != 0) {
        directed = igraph_is_directed(VECTOR(*graphs)[0]);
    }
    for (i = 1; i < no_of_graphs; i++) {
        if (directed != igraph_is_directed(VECTOR(*graphs)[i])) {
            IGRAPH_ERROR("Cannot create intersection of directed and undirected graphs.",
                         IGRAPH_EINVAL);
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_init(&no_edges, no_of_graphs));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &no_edges);

    /* Calculate number of nodes, query number of edges */
    for (i = 0; i < no_of_graphs; i++) {
        igraph_integer_t n = igraph_vcount(VECTOR(*graphs)[i]);
        if (n > no_of_nodes) {
            no_of_nodes = n;
        }
        VECTOR(no_edges)[i] = igraph_ecount(VECTOR(*graphs)[i]);
        allne = allne && VECTOR(no_edges)[i] > 0;
    }

    if (edgemaps) {
        IGRAPH_CHECK(igraph_vector_int_list_resize(edgemaps, no_of_graphs));
        for (i = 0; i < no_of_graphs; i++) {
            igraph_vector_int_t* v = igraph_vector_int_list_get_ptr(edgemaps, i);
            IGRAPH_CHECK(igraph_vector_int_resize(v, VECTOR(no_edges)[i]));
            igraph_vector_int_fill(v, -1);
        }
    }

    /* Allocate memory for the edge lists and their index vectors */
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&edge_vects, no_of_graphs);
    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&order_vects, no_of_graphs);

    /* Query and sort the edge lists */
    for (i = 0; i < no_of_graphs; i++) {
        igraph_integer_t k, j, n = VECTOR(no_edges)[i];
        igraph_vector_int_t *ev = igraph_vector_int_list_get_ptr(&edge_vects, i);
        igraph_vector_int_t *order = igraph_vector_int_list_get_ptr(&order_vects, i);
        IGRAPH_CHECK(igraph_get_edgelist(VECTOR(*graphs)[i], ev, /*bycol=*/ false));
        if (!directed) {
            for (k = 0, j = 0; k < n; k++, j += 2) {
                if (VECTOR(*ev)[j] > VECTOR(*ev)[j + 1]) {
                    igraph_integer_t tmp = VECTOR(*ev)[j];
                    VECTOR(*ev)[j] = VECTOR(*ev)[j + 1];
                    VECTOR(*ev)[j + 1] = tmp;
                }
            }
        }
        IGRAPH_CHECK(igraph_vector_int_resize(order, n));
        for (k = 0; k < n; k++) {
            VECTOR(*order)[k] = k;
        }
        igraph_qsort_r(VECTOR(*order), n, sizeof(VECTOR(*order)[0]), ev,
                       igraph_i_order_edgelist_cmp);
    }

    /* Do the merge. We work from the end of the edge lists,
       because then we don't have to keep track of where we are right
       now in the edge and order lists. We find the "largest" edge,
       and if it is present in all graphs, then we copy it to the
       result. We remove all instances of this edge.  */

    while (allne) {

        /* Look for the smallest tail element */
        for (j = 0, tailfrom = IGRAPH_INTEGER_MAX, tailto = IGRAPH_INTEGER_MAX; j < no_of_graphs; j++) {
            igraph_vector_int_t *order = igraph_vector_int_list_get_ptr(&order_vects, j);
            igraph_vector_int_t *ev = igraph_vector_int_list_get_ptr(&edge_vects, j);
            igraph_integer_t edge = igraph_vector_int_tail(order);
            igraph_integer_t from = VECTOR(*ev)[2 * edge];
            igraph_integer_t to = VECTOR(*ev)[2 * edge + 1];
            if (from < tailfrom || (from == tailfrom && to < tailto)) {
                tailfrom = from; tailto = to;
            }
        }

        /* OK, now remove all elements from the tail(s) that are bigger
           than the smallest tail element. */
        for (j = 0, allsame = true; j < no_of_graphs; j++) {
            igraph_integer_t from = -1, to = -1;
            igraph_vector_int_t *order = igraph_vector_int_list_get_ptr(&order_vects, j);
            while (true) {
                igraph_integer_t edge = igraph_vector_int_tail(order);
                igraph_vector_int_t *ev = igraph_vector_int_list_get_ptr(&edge_vects, j);
                from = VECTOR(*ev)[2 * edge];
                to = VECTOR(*ev)[2 * edge + 1];
                if (from > tailfrom || (from == tailfrom && to > tailto)) {
                    igraph_vector_int_pop_back(order);
                    if (igraph_vector_int_empty(order)) {
                        allne = false;
                        break;
                    }
                } else {
                    break;
                }
            }
            if (from != tailfrom || to != tailto) {
                allsame = false;
            }
        }

        /* Add the edge, if the smallest tail element was present
           in all graphs. */
        if (allsame) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tailfrom));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, tailto));
        }

        /* Drop edges matching the smallest tail elements
           from the order vectors, build edge maps */
        if (allne) {
            for (j = 0; j < no_of_graphs; j++) {
                igraph_vector_int_t *order = igraph_vector_int_list_get_ptr(&order_vects, j);
                igraph_integer_t edge = igraph_vector_int_tail(order);
                igraph_vector_int_t *ev = igraph_vector_int_list_get_ptr(&edge_vects, j);
                igraph_integer_t from = VECTOR(*ev)[2 * edge];
                igraph_integer_t to = VECTOR(*ev)[2 * edge + 1];
                if (from == tailfrom && to == tailto) {
                    igraph_vector_int_pop_back(order);
                    if (igraph_vector_int_empty(order)) {
                        allne = false;
                    }
                    if (edgemaps && allsame) {
                        igraph_vector_int_t *map = igraph_vector_int_list_get_ptr(edgemaps, j);
                        VECTOR(*map)[edge] = idx;
                    }
                }
            }
            if (allsame) {
                idx++;
            }
        }

    } /* while allne */

    igraph_vector_int_list_destroy(&order_vects);
    igraph_vector_int_list_destroy(&edge_vects);
    igraph_vector_int_destroy(&no_edges);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, directed));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
