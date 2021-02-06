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

#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_interface.h"

#include "graph/attributes.h"
#include "core/interruption.h"

/**
 * \function igraph_difference
 * \brief Calculate the difference of two graphs
 *
 * </para><para>
 * The number of vertices in the result is the number of vertices in
 * the original graph, i.e. the left, first operand. In the results
 * graph only edges will be included from \p orig which are not
 * present in \p sub.
 *
 * \param res Pointer to an uninitialized graph object, the result
 * will be stored here.
 * \param orig The left operand of the operator, a graph object.
 * \param sub The right operand of the operator, a graph object.
 * \return Error code.
 * \sa \ref igraph_intersection() and \ref igraph_union() for other
 * operators.
 *
 * Time complexity: O(|V|+|E|), |V| is the number vertices in
 * the smaller graph, |E| is the
 * number of edges in the result graph.
 *
 * \example examples/simple/igraph_difference.c
 */
int igraph_difference(igraph_t *res,
                      const igraph_t *orig, const igraph_t *sub) {

    /* Quite nasty, but we will use that an edge adjacency list
       contains the vertices according to the order of the
       vertex ids at the "other" end of the edge. */

    long int no_of_nodes_orig = igraph_vcount(orig);
    long int no_of_nodes_sub = igraph_vcount(sub);
    long int no_of_nodes = no_of_nodes_orig;
    long int smaller_nodes;
    igraph_bool_t directed = igraph_is_directed(orig);
    igraph_vector_t edges;
    igraph_vector_t edge_ids;
    igraph_vector_int_t *nei1, *nei2;
    igraph_inclist_t inc_orig, inc_sub;
    long int i;
    igraph_integer_t v1, v2;

    if (directed != igraph_is_directed(sub)) {
        IGRAPH_ERROR("Cannot subtract directed and undirected graphs",
                     IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edge_ids, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_inclist_init(orig, &inc_orig, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inc_orig);
    IGRAPH_CHECK(igraph_inclist_init(sub, &inc_sub, IGRAPH_OUT, IGRAPH_LOOPS_ONCE));
    IGRAPH_FINALLY(igraph_inclist_destroy, &inc_sub);

    smaller_nodes = no_of_nodes_orig > no_of_nodes_sub ?
                    no_of_nodes_sub : no_of_nodes_orig;

    for (i = 0; i < smaller_nodes; i++) {
        long int n1, n2, e1, e2;
        IGRAPH_ALLOW_INTERRUPTION();
        nei1 = igraph_inclist_get(&inc_orig, i);
        nei2 = igraph_inclist_get(&inc_sub, i);
        n1 = igraph_vector_int_size(nei1) - 1;
        n2 = igraph_vector_int_size(nei2) - 1;
        while (n1 >= 0 && n2 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            e2 = (long int) VECTOR(*nei2)[n2];
            v1 = IGRAPH_OTHER(orig, e1, i);
            v2 = IGRAPH_OTHER(sub, e2, i);

            if (!directed && v1 < i) {
                n1--;
            } else if (!directed && v2 < i) {
                n2--;
            } else if (v1 > v2) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));
                n1--;
                /* handle loop edges properly in undirected graphs */
                if (!directed && i == v1) {
                    n1--;
                }
            } else if (v2 > v1) {
                n2--;
            } else {
                n1--;
                n2--;
            }
        }

        /* Copy remaining edges */
        while (n1 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            v1 = IGRAPH_OTHER(orig, e1, i);
            if (directed || v1 >= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));

                /* handle loop edges properly in undirected graphs */
                if (!directed && v1 == i) {
                    n1--;
                }
            }
            n1--;
        }
    }

    /* copy remaining edges, use the previous value of 'i' */
    for (; i < no_of_nodes_orig; i++) {
        long int n1, e1;
        nei1 = igraph_inclist_get(&inc_orig, i);
        n1 = igraph_vector_int_size(nei1) - 1;
        while (n1 >= 0) {
            e1 = (long int) VECTOR(*nei1)[n1];
            v1 = IGRAPH_OTHER(orig, e1, i);
            if (directed || v1 >= i) {
                IGRAPH_CHECK(igraph_vector_push_back(&edge_ids, e1));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v1));

                /* handle loop edges properly in undirected graphs */
                if (!directed && v1 == i) {
                    n1--;
                }
            }
            n1--;
        }
    }

    igraph_inclist_destroy(&inc_sub);
    igraph_inclist_destroy(&inc_orig);
    IGRAPH_FINALLY_CLEAN(2);
    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    /* Attributes */
    if (orig->attr) {
        IGRAPH_I_ATTRIBUTE_DESTROY(res);
        IGRAPH_I_ATTRIBUTE_COPY(res, orig, /*graph=*/1, /*vertex=*/1, /*edge=*/0);
        IGRAPH_CHECK(igraph_i_attribute_permute_edges(orig, res, &edge_ids));
    }

    igraph_vector_destroy(&edge_ids);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
