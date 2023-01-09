/*
   IGraph library.
   Copyright (C) 2015-2022  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    /* this is a directed graph with 10 vertices and 20 edges: */
    igraph_integer_t vc = 10, ec = 20;
    igraph_integer_t edge_data[] = {
        0, 2, 1, 4, 2, 5, 2, 7, 3, 7, 3, 8, 4, 2, 5, 8, 6, 0, 6, 1, 6, 2, 7,
        0, 8, 0, 8, 2, 8, 3, 8, 5, 9, 2, 9, 3, 9, 4, 9, 5
    };

    igraph_vector_int_t edges;
    igraph_vector_t tri;
    igraph_t graph;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    igraph_vector_int_view(&edges, edge_data, 2 * ec);

    igraph_create(&graph, &edges, vc, 1 /* directed=true */);

    igraph_vector_init(&tri, 0);

    igraph_triad_census(&graph, &tri);
    print_vector_round(&tri);

    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, NULL); /* convert to undirected */
    igraph_triad_census(&graph, &tri);
    print_vector_round(&tri);

    igraph_vector_destroy(&tri);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
