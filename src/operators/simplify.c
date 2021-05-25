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

#include "core/fixed_vectorlist.h"
#include "graph/attributes.h"

/**
 * \ingroup structural
 * \function igraph_simplify
 * \brief Removes loop and/or multiple edges from the graph.
 *
 * \param graph The graph object.
 * \param multiple Logical, if true, multiple edges will be removed.
 * \param loops Logical, if true, loops (self edges) will be removed.
 * \param edge_comb What to do with the edge attributes. See the igraph
 *        manual section about attributes for details.
 * \return Error code:
 *    \c IGRAPH_ENOMEM if we are out of memory.
 *
 * Time complexity: O(|V|+|E|).
 *
 * \example examples/simple/igraph_simplify.c
 */

int igraph_simplify(igraph_t *graph, igraph_bool_t multiple,
                    igraph_bool_t loops,
                    const igraph_attribute_combination_t *edge_comb) {

    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    long int edge;
    igraph_bool_t attr = edge_comb && igraph_has_attribute_table();
    long int from, to, pfrom = -1, pto = -2;
    igraph_t res;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_vector_t mergeinto;
    long int actedge;

    if (!multiple && !loops)
        /* nothing to do */
    {
        return IGRAPH_SUCCESS;
    }

    if (!multiple) {
        /* removing loop edges only, this is simple. No need to combine anything
         * and the whole process can be done in-place */
        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
        IGRAPH_FINALLY(igraph_es_destroy, &es);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        while (!IGRAPH_EIT_END(eit)) {
            edge = IGRAPH_EIT_GET(eit);
            from = IGRAPH_FROM(graph, edge);
            to = IGRAPH_TO(graph, edge);
            if (from == to) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, edge));
            }
            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
        IGRAPH_FINALLY_CLEAN(2);

        if (igraph_vector_size(&edges) > 0) {
            IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_vector(&edges)));
        }

        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);

        return IGRAPH_SUCCESS;
    }

    if (attr) {
        IGRAPH_VECTOR_INIT_FINALLY(&mergeinto, no_of_edges);
    }
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_FROM));
    IGRAPH_FINALLY(igraph_es_destroy, &es);
    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    for (actedge = -1; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        edge = IGRAPH_EIT_GET(eit);
        from = IGRAPH_FROM(graph, edge);
        to = IGRAPH_TO(graph, edge);

        if (loops && from == to) {
            /* Loop edge to be removed */
            if (attr) {
                VECTOR(mergeinto)[edge] = -1;
            }
        } else if (multiple && from == pfrom && to == pto) {
            /* Multiple edge to be contracted */
            if (attr) {
                VECTOR(mergeinto)[edge] = actedge;
            }
        } else {
            /* Edge to be kept */
            igraph_vector_push_back(&edges, from);
            igraph_vector_push_back(&edges, to);
            if (attr) {
                actedge++;
                VECTOR(mergeinto)[edge] = actedge;
            }
        }
        pfrom = from; pto = to;
    }

    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(&res, &edges, (igraph_integer_t) no_of_nodes,
                               igraph_is_directed(graph)));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ 1,
                            /*vertex=*/ 1, /*edge=*/ 0);

    if (attr) {
        igraph_fixed_vectorlist_t vl;
        IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto,
                     actedge + 1));
        IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

        IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &res, &vl.v,
                     edge_comb));

        igraph_fixed_vectorlist_destroy(&vl);
        igraph_vector_destroy(&mergeinto);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    return 0;
}
