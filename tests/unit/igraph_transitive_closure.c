/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

void print_closure(const igraph_t *graph) {
    igraph_t closure;
    igraph_bool_t simple;

    igraph_transitive_closure(graph, &closure);
    print_graph_canon(&closure);
    IGRAPH_ASSERT(igraph_vcount(graph) == igraph_vcount(&closure));
    IGRAPH_ASSERT(igraph_is_directed(graph) == igraph_is_directed(&closure));
    igraph_is_simple(&closure, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&closure);
}

int main(void) {
    igraph_t graph;

    printf("Directed null graph\n");
    igraph_empty(&graph, 0, IGRAPH_DIRECTED);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nUndirected null graph\n");
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nDirected singleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_DIRECTED);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nUndirected singleton graph\n");
    igraph_empty(&graph, 1, IGRAPH_UNDIRECTED);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nEdgeless graph\n");
    igraph_empty(&graph, 3, IGRAPH_DIRECTED);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nSmall DAG\n");
    igraph_small(&graph, 9, IGRAPH_DIRECTED,
                 8, 7, 7, 6, 6, 3, 6, 0, 3, 2, 3, 1, 5, 0, 4, 1,
                 -1);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nSmall directed multigraph with cycles \n");
    igraph_small(&graph, 10, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,0, 2,0, 0,3, 3,4, 4,3, 3,3, 0,5, 5,6, 6,6, 8,7,
                 -1);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nSmall undirected graph\n");
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 3,4,
                 -1);
    print_closure(&graph);
    igraph_destroy(&graph);

    printf("\nGraph and vertex attributes are kept, edge attributes are dropped\n");
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_small(&graph, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,1, 3,3,
                 -1);
    SETGAS(&graph, "name", "test graph");
    SETVAN(&graph, "weight", 0, 10);
    SETVAN(&graph, "weight", 1, 11);
    SETVAN(&graph, "weight", 2, 12);
    SETVAN(&graph, "weight", 3, 13);
    SETVAS(&graph, "label", 2, "two");
    SETEAN(&graph, "capacity", 0, 5);
    {
        igraph_t closure;
        igraph_transitive_closure(&graph, &closure);
        print_graph_canon(&closure);
        print_attributes(&closure);
        IGRAPH_ASSERT(!igraph_cattribute_has_attr(&closure, IGRAPH_ATTRIBUTE_EDGE, "capacity"));
        igraph_destroy(&closure);
    }
    igraph_destroy(&graph);
    igraph_set_attribute_table(NULL);

    VERIFY_FINALLY_STACK();

    return 0;
}
