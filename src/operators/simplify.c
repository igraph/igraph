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
 * This function merges parallel edges and removes self-loops, according
 * to the \p multiple and \p loops parameters. Note that this function
 * may change the edge order, even if the input was already a simple graph.
 *
 * \param graph The graph object.
 * \param multiple Logical, if true, multiple edges will be removed.
 * \param loops Logical, if true, loops (self edges) will be removed.
 * \param edge_comb What to do with the edge attributes. \c NULL means to
 *        discard the edge attributes after the operation, even for edges
 *        that were unaffected. See the igraph manual section about attributes
 *        for details.
 * \return Error code:
 *    \c IGRAPH_ENOMEM if we are out of memory.
 *
 * Time complexity: O(|V|+|E|).
 *
 * \example examples/simple/igraph_simplify.c
 */

igraph_error_t igraph_simplify(igraph_t *graph,
                               igraph_bool_t multiple, igraph_bool_t loops,
                               const igraph_attribute_combination_t *edge_comb) {

    igraph_vector_int_t edges;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_integer_t edge;
    igraph_bool_t attr = edge_comb && igraph_has_attribute_table();
    igraph_integer_t from, to, pfrom = -1, pto = -2;
    igraph_t res;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_vector_int_t mergeinto;
    igraph_integer_t actedge;

    /* if we already know there are no multi-edges, they don't need to be removed */
    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_MULTI) &&
        !igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_MULTI)) {
        multiple = false;
    }

    /* if we already know there are no loops, they don't need to be removed */
    if (igraph_i_property_cache_has(graph, IGRAPH_PROP_HAS_LOOP) &&
        !igraph_i_property_cache_get_bool(graph, IGRAPH_PROP_HAS_LOOP)) {
        loops = false;
    }

    if (!multiple && !loops)
        /* nothing to do */
    {
        return IGRAPH_SUCCESS;
    }

    if (!multiple) {
        igraph_vector_int_t edges_to_delete;

        /* removing loop edges only, this is simple. No need to combine anything
         * and the whole process can be done in-place */
        IGRAPH_VECTOR_INT_INIT_FINALLY(&edges_to_delete, 0);
        IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
        IGRAPH_FINALLY(igraph_es_destroy, &es);
        IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
        IGRAPH_FINALLY(igraph_eit_destroy, &eit);

        while (!IGRAPH_EIT_END(eit)) {
            edge = IGRAPH_EIT_GET(eit);
            from = IGRAPH_FROM(graph, edge);
            to = IGRAPH_TO(graph, edge);
            if (from == to) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges_to_delete, edge));
            }
            IGRAPH_EIT_NEXT(eit);
        }

        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
        IGRAPH_FINALLY_CLEAN(2);

        if (igraph_vector_int_size(&edges_to_delete) > 0) {
            IGRAPH_CHECK(igraph_delete_edges(graph, igraph_ess_vector(&edges_to_delete)));
        }

        igraph_vector_int_destroy(&edges_to_delete);
        IGRAPH_FINALLY_CLEAN(1);

        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, false);

        return IGRAPH_SUCCESS;
    }

    if (attr) {
        IGRAPH_VECTOR_INT_INIT_FINALLY(&mergeinto, no_of_edges);
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges * 2));

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
            igraph_vector_int_push_back(&edges, from);  /* reserved */
            igraph_vector_int_push_back(&edges, to);  /* reserved */
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

    IGRAPH_CHECK(igraph_create(&res, &edges, no_of_nodes, igraph_is_directed(graph)));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_FINALLY(igraph_destroy, &res);

    IGRAPH_I_ATTRIBUTE_DESTROY(&res);
    IGRAPH_I_ATTRIBUTE_COPY(&res, graph, /*graph=*/ true, /*vertex=*/ true, /*edge=*/ false);

    if (attr) {
        igraph_fixed_vectorlist_t vl;
        IGRAPH_CHECK(igraph_fixed_vectorlist_convert(&vl, &mergeinto, actedge + 1));
        IGRAPH_FINALLY(igraph_fixed_vectorlist_destroy, &vl);

        IGRAPH_CHECK(igraph_i_attribute_combine_edges(graph, &res, &vl.vecs, edge_comb));

        igraph_fixed_vectorlist_destroy(&vl);
        igraph_vector_int_destroy(&mergeinto);
        IGRAPH_FINALLY_CLEAN(2);
    }

    IGRAPH_FINALLY_CLEAN(1);
    igraph_destroy(graph);
    *graph = res;

    /* The cache must be set as the very last step, only after all functions that can
     * potentially return with an error have finished. */

    if (loops) {
        /* Loop edges were removed so we know for sure that there aren't any
         * loop edges now */
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_LOOP, false);
    }

    if (multiple) {
        /* Multi-edges were removed so we know for sure that there aren't any
         * multi-edges now */
        igraph_i_property_cache_set_bool_checked(graph, IGRAPH_PROP_HAS_MULTI, false);
    }

    return IGRAPH_SUCCESS;
}
