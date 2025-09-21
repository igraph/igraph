/*
   igraph library.
   Copyright (C) 2006-2025  The igraph development team <igraph@igraph.org>

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

#include "igraph_operators.h"
#include "igraph_isomorphism.h"

#include "igraph_constructors.h"
#include "igraph_interface.h"

#include "graph/attributes.h"

/**
 * \function igraph_invert_permutation
 * \brief Inverts a permutation.
 *
 * Produces the inverse of \p permutation into \p inverse and at the same time it checks
 * that the permutation vector is valid, i.e. all indices are within range and there are
 * no duplicate entries.
 *
 * \param permutation A permutation vector containing 0-based integer indices.
 * \param inverse An initialized vector. The inverse of \p permutation will be stored here.
 * \return Error code.
 */
igraph_error_t igraph_invert_permutation(const igraph_vector_int_t *permutation, igraph_vector_int_t *inverse) {
    const igraph_int_t n = igraph_vector_int_size(permutation);

    IGRAPH_CHECK(igraph_vector_int_resize(inverse, n));
    igraph_vector_int_fill(inverse, -1);

    for (igraph_int_t i = 0; i < n; i++) {
        igraph_int_t j = VECTOR(*permutation)[i];
        if (j < 0 || j >= n) {
            IGRAPH_ERROR("Invalid index in permutation vector.", IGRAPH_EINVAL);
        }
        if (VECTOR(*inverse)[j] != -1) {
            /* This element of 'inverse' has already been set, 'j' is a duplicate value. */
            IGRAPH_ERROR("Duplicate entry in permutation vector.", IGRAPH_EINVAL);
        }
        VECTOR(*inverse)[j] = i;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_permute_vertices
 * \brief Permute the vertices.
 *
 * This function creates a new graph from the input graph by permuting
 * its vertices according to the specified mapping. Call this function
 * with the output of \ref igraph_canonical_permutation() to create
 * the canonical form of a graph.
 *
 * \param graph The input graph.
 * \param res Pointer to an uninitialized graph object. The new graph
 *    is created here.
 * \param permutation The permutation to apply. The i-th element of the
 *    vector specifies the index of the vertex in the original graph that
 *    will become vertex i in the new graph.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in terms of the number of
 * vertices and edges.
 */
igraph_error_t igraph_permute_vertices(const igraph_t *graph, igraph_t *res,
                                       const igraph_vector_int_t *permutation) {

    igraph_int_t no_of_nodes = igraph_vcount(graph);
    igraph_int_t no_of_edges = igraph_ecount(graph);
    igraph_vector_int_t edges;
    igraph_vector_int_t index;
    igraph_int_t p;

    if (igraph_vector_int_size(permutation) != no_of_nodes) {
        IGRAPH_ERROR("Permute vertices: invalid permutation vector size.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&index, no_of_nodes);

    /* Also checks that 'permutation' is valid: */
    IGRAPH_CHECK(igraph_invert_permutation(permutation, &index));

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    p = 0;
    for (igraph_int_t i = 0; i < no_of_edges; i++) {
        VECTOR(edges)[p++] = VECTOR(index)[ IGRAPH_FROM(graph, i) ];
        VECTOR(edges)[p++] = VECTOR(index)[ IGRAPH_TO(graph, i) ];
    }

    IGRAPH_CHECK(igraph_create(res, &edges, no_of_nodes, igraph_is_directed(graph)));
    IGRAPH_FINALLY(igraph_destroy, res);

    /* Attributes */
    if (graph->attr) {
        igraph_vector_int_t vtypes;
        IGRAPH_CHECK(igraph_i_attribute_copy(res, graph, true, /* vertex= */ false, true));
        IGRAPH_VECTOR_INT_INIT_FINALLY(&vtypes, 0);
        IGRAPH_CHECK(igraph_i_attribute_get_info(graph, 0, 0, 0, &vtypes, 0, 0));
        if (igraph_vector_int_size(&vtypes) != 0) {
            IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph, res, permutation));
        }
        igraph_vector_int_destroy(&vtypes);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&index);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(3); /* +1 for res */

    return IGRAPH_SUCCESS;
}
