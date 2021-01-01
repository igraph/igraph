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
#include "igraph_interface.h"

#include "operators/misc_internal.h"

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
 * Both graphs need to have the same directedness, i.e. either both
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
 * All graphs need to have the same directedness, i.e. either all
 * directed or all undirected. If the graph list has length zero,
 * the result will be a \em directed graph with no vertices.
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
