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

#include "graph/attributes.h"
#include "core/interruption.h"

/**
 * \function igraph_complementer
 * \brief Creates the complementer of a graph.
 *
 * The complementer graph means that all edges which are
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
igraph_error_t igraph_complementer(igraph_t *res, const igraph_t *graph,
                        igraph_bool_t loops) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_vector_int_t edges;
    igraph_vector_int_t neis;
    igraph_integer_t i, j;
    igraph_integer_t zero = 0, *limit;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);

    if (igraph_is_directed(graph)) {
        limit = &zero;
    } else {
        limit = &i;
    }

    for (i = 0; i < no_of_nodes; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
        if (loops) {
            for (j = no_of_nodes - 1; j >= *limit; j--) {
                if (igraph_vector_int_empty(&neis) || j > igraph_vector_int_tail(&neis)) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                    IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
                } else {
                    igraph_vector_int_pop_back(&neis);
                }
            }
        } else {
            for (j = no_of_nodes - 1; j >= *limit; j--) {
                if (igraph_vector_int_empty(&neis) || j > igraph_vector_int_tail(&neis)) {
                    if (i != j) {
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, j));
                    }
                } else {
                    igraph_vector_int_pop_back(&neis);
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, igraph_is_directed(graph)));
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&neis);
    IGRAPH_I_ATTRIBUTE_DESTROY(res);
    IGRAPH_I_ATTRIBUTE_COPY(res, graph, /*graph=*/true, /*vertex=*/true, /*edge=*/false);
    IGRAPH_FINALLY_CLEAN(2);
    return IGRAPH_SUCCESS;
}
