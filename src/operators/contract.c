/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2021 The igraph development team

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
 * \function igraph_contract_vertices
 * \brief Replace multiple vertices with a single one.
 *
 * This function modifies the graph by merging several vertices
 * into one. The vertices in the modified graph correspond
 * to groups of vertices in the input graph. No edges are removed,
 * thus the modified graph will typically have self-loops
 * (corresponding to in-group edges) and multi-edges
 * (corresponding to multiple connections between two groups).
 * Use \ref igraph_simplify() to eliminate self-loops and
 * merge multi-edges.
 *
 * \param graph The input graph. It will be modified in-place.
 * \param mapping A vector giving the mapping. For each
 *        vertex in the original graph, it should contain
 *        its desired ID in the result graph. In order not to create
 *        "orphan vertices" that have no corresponding vertices
 *        in the original graph, ensure that the IDs are consecutive
 *        integers starting from zero.
 * \param vertex_comb What to do with the vertex attributes.
 *        \c NULL means that vertex attributes are not kept
 *        after the contraction (not even for unaffected
 *        vertices). See the igraph manual section about attributes
 *        for details.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number
 * or vertices plus edges.
 */

igraph_error_t igraph_contract_vertices(igraph_t *graph,
                             const igraph_vector_int_t *mapping,
                             const igraph_attribute_combination_t *vertex_comb) {
    igraph_vector_int_t edges;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_bool_t vattr = vertex_comb && igraph_has_attribute_table();
    igraph_t res;
    igraph_integer_t last;
    igraph_integer_t no_new_vertices;

    if (igraph_vector_int_size(mapping) != no_of_nodes) {
        IGRAPH_ERRORF("Mapping vector length (%" IGRAPH_PRId ") "
                "not equal to number of nodes (%" IGRAPH_PRId ").",
                IGRAPH_EINVAL, igraph_vector_int_size(mapping), no_of_nodes);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));

    if (no_of_nodes > 0) {
        last = igraph_vector_int_max(mapping);
    } else {
        /* Ensure that no_new_vertices will be zero
         * when the input graph has no vertices. */
        last = -1;
    }

    for (igraph_integer_t edge = 0; edge < no_of_edges; edge++) {
        igraph_integer_t from = IGRAPH_FROM(graph, edge);
        igraph_integer_t to = IGRAPH_TO(graph, edge);

        igraph_integer_t nfrom = VECTOR(*mapping)[from];
        igraph_integer_t nto = VECTOR(*mapping)[to];

        igraph_vector_int_push_back(&edges, nfrom);
        igraph_vector_int_push_back(&edges, nto);

        if (nfrom > last) {
            last = nfrom;
        }
        if (nto   > last) {
            last = nto;
        }
    }

    no_new_vertices = last + 1;

    IGRAPH_CHECK(igraph_create(&res, &edges, no_new_vertices,
                               igraph_is_directed(graph)));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ true, /*vertex=*/ false, /*edge=*/ true);

    if (vattr) {
        igraph_vector_int_list_t merges;
        igraph_vector_int_t sizes;

        IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&merges, no_new_vertices);
        IGRAPH_VECTOR_INT_INIT_FINALLY(&sizes, no_new_vertices);

        for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
            igraph_integer_t to = VECTOR(*mapping)[i];
            igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(&merges, to);
            VECTOR(sizes)[to] += 1;
            IGRAPH_CHECK(igraph_vector_int_push_back(v, i));
        }

        IGRAPH_CHECK(igraph_i_attribute_combine_vertices(graph, &res,
                     &merges,
                     vertex_comb));

        igraph_vector_int_destroy(&sizes);
        igraph_vector_int_list_destroy(&merges);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    return IGRAPH_SUCCESS;
}
