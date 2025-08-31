/*
   igraph library.
   Copyright (C) 2006-2024  The igraph development team <igraph@igraph.org>

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

/* Prints a vector of edge IDs as u--v vertex pairs. */
void print_edge_vector(const igraph_t *graph, const igraph_vector_int_t *edges) {
    const igraph_int_t n = igraph_vector_int_size(edges);
    for (igraph_int_t i=0; i < n; i++) {
        igraph_int_t edge = VECTOR(*edges)[i];
        printf("%" IGRAPH_PRId "--%" IGRAPH_PRId " ", IGRAPH_FROM(graph, edge), IGRAPH_TO(graph, edge));
    }
    printf("\n");
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t component_vertices, component_edges;
    igraph_int_t no;

    /* Initialize the library. */
    igraph_setup();

    /* Create an example graph. */
    igraph_small(&graph, 7, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,3, 3,0,
                 2,4, 4,5, 5,2,
                 0,6,
                 0,7,
                 -1);

    /* The data structures that the result will be written to must be initialized first. */
    igraph_vector_int_list_init(&component_vertices, 0);
    igraph_vector_int_list_init(&component_edges, 0);

    igraph_biconnected_components(&graph, &no, NULL, &component_edges, &component_vertices, NULL);

    printf("Number of components: %" IGRAPH_PRId "\n", no);
    for (igraph_int_t i=0; i < no; i++) {
        printf("\n");
        printf("Component %" IGRAPH_PRId ":\n", i);
        printf("Vertices: ");
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(&component_vertices, i));
        printf("Edges: ");
        print_edge_vector(&graph, igraph_vector_int_list_get_ptr(&component_edges, i));
    }

    /* Destroy data structures after we no longer need them. */

    igraph_vector_int_list_destroy(&component_edges);
    igraph_vector_int_list_destroy(&component_vertices);

    igraph_destroy(&graph);

    return 0;
}
