/* -*- mode: C -*-  */
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

#include "igraph_constructors.h"

#include "igraph_interface.h"

#include "core/interruption.h"

/* Note to self: tried using adjacency lists instead of igraph_incident queries,
 * with minimal performance improvements on a graph with 70K vertices and 360K
 * edges. (1.09s instead of 1.10s). I think it's not worth the fuss. */
static igraph_error_t igraph_i_linegraph_undirected(const igraph_t *graph, igraph_t *linegraph) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t adjedges, adjedges2;
    igraph_vector_int_t edges;
    igraph_integer_t prev = -1;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjedges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjedges2, 0);

    for (igraph_integer_t e1 = 0; e1 < no_of_edges; e1++) {
        igraph_integer_t from = IGRAPH_FROM(graph, e1);
        igraph_integer_t to = IGRAPH_TO(graph, e1);
        igraph_integer_t n;

        IGRAPH_ALLOW_INTERRUPTION();

        if (from != prev) {
            IGRAPH_CHECK(igraph_incident(graph, &adjedges, from, IGRAPH_ALL));
        }
        n = igraph_vector_int_size(&adjedges);
        for (igraph_integer_t i = 0; i < n; i++) {
            igraph_integer_t e2 = VECTOR(adjedges)[i];
            if (e2 < e1) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e2));
            }
        }

        IGRAPH_CHECK(igraph_incident(graph, &adjedges2, to, IGRAPH_ALL));
        n = igraph_vector_int_size(&adjedges2);
        for (igraph_integer_t i = 0; i < n; i++) {
            igraph_integer_t e2 = VECTOR(adjedges2)[i];
            if (e2 < e1) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e1));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e2));
            }
        }

        /* Self-loops are considered self-adjacent. */
        if (from == to) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e1));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e1));
        }

        prev = from;
    }

    igraph_vector_int_destroy(&adjedges);
    igraph_vector_int_destroy(&adjedges2);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(linegraph, &edges, no_of_edges, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_linegraph_directed(const igraph_t *graph, igraph_t *linegraph) {
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t i, j, n;
    igraph_vector_int_t adjedges;
    igraph_vector_int_t edges;
    igraph_integer_t prev = -1;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&adjedges, 0);

    for (i = 0; i < no_of_edges; i++) {
        igraph_integer_t from = IGRAPH_FROM(graph, i);

        IGRAPH_ALLOW_INTERRUPTION();

        if (from != prev) {
            IGRAPH_CHECK(igraph_incident(graph, &adjedges, from, IGRAPH_IN));
        }
        n = igraph_vector_int_size(&adjedges);
        for (j = 0; j < n; j++) {
            igraph_integer_t e = VECTOR(adjedges)[j];
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, e));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
        }

        prev = from;
    }

    igraph_vector_int_destroy(&adjedges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_create(linegraph, &edges, no_of_edges, igraph_is_directed(graph)));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_linegraph
 * \brief Create the line graph of a graph.
 *
 * The line graph L(G) of a G undirected graph is defined as follows.
 * L(G) has one vertex for each edge in G and two different vertices in L(G)
 * are connected by an edge if their corresponding edges share an end point.
 * In a multigraph, if two end points are shared, two edges are created.
 * The single vertex of an undirected self-loop is counted as two end points.
 *
 * </para><para>
 * The line graph L(G) of a G directed graph is slightly different:
 * L(G) has one vertex for each edge in G and two vertices in L(G) are connected
 * by a directed edge if the target of the first vertex's corresponding edge
 * is the same as the source of the second vertex's corresponding edge.
 *
 * </para><para>
 * Self-loops are considered self-adjacent, thus their corresponding vertex
 * in the line graph will also a have a single self-loop, in both undirected
 * and directed graphs.
 *
 * </para><para>
 * Edge \em i  in the original graph will correspond to vertex \em i
 * in the line graph.
 *
 * </para><para>
 * The first version of this function was contributed by Vincent Matossian,
 * thanks.
 *
 * \param graph The input graph, may be directed or undirected.
 * \param linegraph Pointer to an uninitialized graph object, the
 *        result is stored here.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), the number of edges plus the number of vertices.
 */

igraph_error_t igraph_linegraph(const igraph_t *graph, igraph_t *linegraph) {

    if (igraph_is_directed(graph)) {
        return igraph_i_linegraph_directed(graph, linegraph);
    } else {
        return igraph_i_linegraph_undirected(graph, linegraph);
    }
}
