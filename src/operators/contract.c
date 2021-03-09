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
#include "igraph_memory.h"

#include "graph/attributes.h"

static void igraph_i_simplify_free(igraph_vector_ptr_t *p) {
    long int i, n = igraph_vector_ptr_size(p);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*p)[i];
        if (v) {
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_destroy(p);
}

/**
 * \function igraph_contract_vertices
 * Replace multiple vertices with a single one.
 *
 * This function creates a new graph, by merging several
 * vertices into one. The vertices in the new graph correspond
 * to sets of vertices in the input graph.
 * \param graph The input graph, it can be directed or
 *        undirected.
 * \param mapping A vector giving the mapping. For each
 *        vertex in the original graph, it should contain
 *        its id in the new graph.
 * \param vertex_comb What to do with the vertex attributes.
 *        See the igraph manual section about attributes for
 *        details.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|), linear in the number
 * or vertices plus edges.
 */

int igraph_contract_vertices(igraph_t *graph,
                             const igraph_vector_t *mapping,
                             const igraph_attribute_combination_t *vertex_comb) {
    igraph_vector_t edges;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t vattr = vertex_comb && igraph_has_attribute_table();
    igraph_t res;
    long int e, last = -1;
    long int no_new_vertices;

    if (igraph_vector_size(mapping) != no_of_nodes) {
        IGRAPH_ERROR("Invalid mapping vector length",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    if (no_of_nodes > 0) {
        last = (long int) igraph_vector_max(mapping);
    }

    for (e = 0; e < no_of_edges; e++) {
        long int from = IGRAPH_FROM(graph, e);
        long int to = IGRAPH_TO(graph, e);

        long int nfrom = (long int) VECTOR(*mapping)[from];
        long int nto = (long int) VECTOR(*mapping)[to];

        igraph_vector_push_back(&edges, nfrom);
        igraph_vector_push_back(&edges, nto);

        if (nfrom > last) {
            last = nfrom;
        }
        if (nto   > last) {
            last = nto;
        }
    }

    no_new_vertices = last + 1;

    IGRAPH_CHECK(igraph_create(&res, &edges, (igraph_integer_t) no_new_vertices,
                               igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ 1,
                            /*vertex=*/ 0, /*edge=*/ 1);

    if (vattr) {
        long int i;
        igraph_vector_ptr_t merges;
        igraph_vector_t sizes;
        igraph_vector_t *vecs;

        vecs = IGRAPH_CALLOC(no_new_vertices, igraph_vector_t);
        if (!vecs) {
            IGRAPH_ERROR("Cannot combine attributes while contracting"
                         " vertices", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, vecs);
        IGRAPH_CHECK(igraph_vector_ptr_init(&merges, no_new_vertices));
        IGRAPH_FINALLY(igraph_i_simplify_free, &merges);
        IGRAPH_VECTOR_INIT_FINALLY(&sizes, no_new_vertices);

        for (i = 0; i < no_of_nodes; i++) {
            long int to = (long int) VECTOR(*mapping)[i];
            VECTOR(sizes)[to] += 1;
        }
        for (i = 0; i < no_new_vertices; i++) {
            igraph_vector_t *v = &vecs[i];
            IGRAPH_CHECK(igraph_vector_init(v, (long int) VECTOR(sizes)[i]));
            igraph_vector_clear(v);
            VECTOR(merges)[i] = v;
        }
        for (i = 0; i < no_of_nodes; i++) {
            long int to = (long int) VECTOR(*mapping)[i];
            igraph_vector_t *v = &vecs[to];
            igraph_vector_push_back(v, i);
        }

        IGRAPH_CHECK(igraph_i_attribute_combine_vertices(graph, &res,
                     &merges,
                     vertex_comb));

        igraph_vector_destroy(&sizes);
        igraph_i_simplify_free(&merges);
        igraph_free(vecs);
        IGRAPH_FINALLY_CLEAN(3);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    return 0;
}
