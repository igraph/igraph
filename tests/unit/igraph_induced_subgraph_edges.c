/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_vector_int_init(&edges, 0);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
              // 0    1    2    3    4    5    6    7    8    9    10
                 4,4, 0,3, 0,3, 2,3, 1,4, 2,4, 3,4, 2,3, 4,4, 2,3, 1,1,
                 -1);

    igraph_induced_subgraph_edges(&graph, igraph_vss_range(2,5), &edges);
    igraph_vector_int_sort(&edges); /* canonicalize */
    print_vector_int(&edges);

    igraph_induced_subgraph_edges(&graph, igraph_vss_1(0), &edges);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);

    igraph_induced_subgraph_edges(&graph, igraph_vss_1(1), &edges);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);

    igraph_induced_subgraph_edges(&graph, igraph_vss_all(), &edges);
    igraph_vector_int_sort(&edges);
    print_vector_int(&edges);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&edges);

    VERIFY_FINALLY_STACK();
    return 0;
}
