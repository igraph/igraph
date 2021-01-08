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

#include "core/interruption.h"

/**
 * \function igraph_compose
 * \brief Calculates the composition of two graphs
 *
 * The composition of graphs contains the same number of vertices as
 * the bigger graph of the two operands. It contains an (i,j) edge if
 * and only if there is a k vertex, such that the first graphs
 * contains an (i,k) edge and the second graph a (k,j) edge.
 *
 * </para><para>This is of course exactly the composition of two
 * binary relations.
 *
 * </para><para>Two two graphs must have the same directedness,
 * otherwise the function returns with an error message.
 * Note that for undirected graphs the two relations are by definition
 * symmetric.
 *
 * \param res Pointer to an uninitialized graph object, the result
 *        will be stored here.
 * \param g1 The firs operand, a graph object.
 * \param g2 The second operand, another graph object.
 * \param edge_map1 If not a null pointer, then it must be a pointer
 *        to an initialized vector, and a mapping from the edges of
 *        the result graph to the edges of the first graph is stored
 *        here.
 * \param edge_map1 If not a null pointer, then it must be a pointer
 *        to an initialized vector, and a mapping from the edges of
 *        the result graph to the edges of the second graph is stored
 *        here.
 * \return Error code.
 *
 * Time complexity: O(|V|*d1*d2), |V| is the number of vertices in the
 * first graph, d1 and d2 the average degree in the first and second
 * graphs.
 *
 * \example examples/simple/igraph_compose.c
 */
int igraph_compose(igraph_t *res, const igraph_t *g1, const igraph_t *g2,
                   igraph_vector_t *edge_map1, igraph_vector_t *edge_map2) {

    long int no_of_nodes_left = igraph_vcount(g1);
    long int no_of_nodes_right = igraph_vcount(g2);
    long int no_of_nodes;
    igraph_bool_t directed = igraph_is_directed(g1);
    igraph_vector_t edges;
    igraph_vector_t neis1, neis2;
    long int i;

    if (directed != igraph_is_directed(g2)) {
        IGRAPH_ERROR("Cannot compose directed and undirected graph",
                     IGRAPH_EINVAL);
    }

    no_of_nodes = no_of_nodes_left > no_of_nodes_right ?
                  no_of_nodes_left : no_of_nodes_right;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&neis2, 0);

    if (edge_map1) {
        igraph_vector_clear(edge_map1);
    }
    if (edge_map2) {
        igraph_vector_clear(edge_map2);
    }

    for (i = 0; i < no_of_nodes_left; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        IGRAPH_CHECK(igraph_incident(g1, &neis1, (igraph_integer_t) i,
                                     IGRAPH_OUT));
        while (!igraph_vector_empty(&neis1)) {
            long int con = (long int) igraph_vector_pop_back(&neis1);
            long int v1 = IGRAPH_OTHER(g1, con, i);
            if (v1 < no_of_nodes_right) {
                IGRAPH_CHECK(igraph_incident(g2, &neis2, (igraph_integer_t) v1,
                                             IGRAPH_OUT));
            } else {
                continue;
            }
            while (!igraph_vector_empty(&neis2)) {
                long int con2 = igraph_vector_pop_back(&neis2);
                long int v2 = IGRAPH_OTHER(g2, con2, v1);
                IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, v2));
                if (edge_map1) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map1, con));
                }
                if (edge_map2) {
                    IGRAPH_CHECK(igraph_vector_push_back(edge_map2, con2));
                }
            }
        }
    }

    igraph_vector_destroy(&neis1);
    igraph_vector_destroy(&neis2);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(res, &edges, (igraph_integer_t) no_of_nodes,
                               directed));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
