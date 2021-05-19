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

/**
 * \function igraph_permute_vertices
 * Permute the vertices
 *
 * This function creates a new graph from the input graph by permuting
 * its vertices according to the specified mapping. Call this function
 * with the output of \ref igraph_canonical_permutation() to create
 * the canonical form of a graph.
 * \param graph The input graph.
 * \param res Pointer to an uninitialized graph object. The new graph
 *    is created here.
 * \param permutation The permutation to apply. Vertex 0 is mapped to
 *    the first element of the vector, vertex 1 to the second,
 * etc. Note that it is not checked that the vector contains every
 *    element only once, and no range checking is performed either.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in terms of the number of
 * vertices and edges.
 */
int igraph_permute_vertices(const igraph_t *graph, igraph_t *res,
                            const igraph_vector_t *permutation) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_vector_t edges;
    long int i, p = 0;

    if (igraph_vector_size(permutation) != no_of_nodes) {
        IGRAPH_ERROR("Permute vertices: invalid permutation vector size", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    for (i = 0; i < no_of_edges; i++) {
        VECTOR(edges)[p++] = VECTOR(*permutation)[ (long int) IGRAPH_FROM(graph, i) ];
        VECTOR(edges)[p++] = VECTOR(*permutation)[ (long int) IGRAPH_TO(graph, i) ];
    }

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));

    /* Attributes */
    if (graph->attr) {
        igraph_vector_t index;
        igraph_vector_t vtypes;
        IGRAPH_I_ATTRIBUTE_DESTROY(res);
        IGRAPH_I_ATTRIBUTE_COPY(res, graph, /*graph=*/1, /*vertex=*/0, /*edge=*/1);
        IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
        IGRAPH_CHECK(igraph_i_attribute_get_info(graph, 0, 0, 0, &vtypes, 0, 0));
        if (igraph_vector_size(&vtypes) != 0) {
            IGRAPH_VECTOR_INIT_FINALLY(&index, no_of_nodes);
            for (i = 0; i < no_of_nodes; i++) {
                VECTOR(index)[ (long int) VECTOR(*permutation)[i] ] = i;
            }
            IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, res, &index));
            igraph_vector_destroy(&index);
            IGRAPH_FINALLY_CLEAN(1);
        }
        igraph_vector_destroy(&vtypes);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
