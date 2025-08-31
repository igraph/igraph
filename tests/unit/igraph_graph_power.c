/* igraph library.
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

void print_result(igraph_t *g, igraph_int_t order, igraph_bool_t directed) {
    igraph_t res;
    igraph_graph_power(g, &res, order, directed);
    print_graph_canon(&res);
    if ((order == 0 || order == 1) && igraph_has_attribute_table()) {
        print_attributes(&res);
    }
    igraph_destroy(&res);
}

int main(void) {
    igraph_t g;


    printf("Directed graph with no vertices:\n");
    igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
    print_result(&g, 10, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 10, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A->B->C, directed, order 2:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 2, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A->B->C, undirected, order 2:\n");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 1,2, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A->B<-C, directed, order 2:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 2,1, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 2, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A<-B<-C, directed, order 2:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, 1,0, 2,1, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 2, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A->B->C->A, directed, order 2:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 2, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nBasic example, A->B<->C->A, directed, order 2:\n");
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, 2,1, 2,0, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);
    printf("Ignore edge directions:\n");
    print_result(&g, 2, IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    printf("\nDirected graph with loops and multiple edges, order 0:\n");
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    print_result(&g, 0, IGRAPH_DIRECTED);
    igraph_destroy(&g);

    printf("\nDirected graph with loops and multiple edges, order 0, with attributes set:\n");
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    SETGAB(&g, "bool_graph_attr", 1);
    SETEAB(&g, "bool_edge_attr", 0, 1);
    SETVAB(&g, "bool_vertex_attr", 0, 1);

    print_result(&g, 0, IGRAPH_DIRECTED);

    printf("Same graph, order 1:\n");
    print_result(&g, 1, IGRAPH_DIRECTED);

    printf("Same starting graph, order 2:\n");
    print_result(&g, 2, IGRAPH_DIRECTED);

    printf("Same starting graph, order 3:\n");
    print_result(&g, 3, IGRAPH_DIRECTED);

    printf("Same starting graph, order 12:\n");
    print_result(&g, 12, IGRAPH_DIRECTED);

    printf("Same starting graph, but undirected, order 2:\n");
    igraph_destroy(&g);
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0,1, 0,2, 1,1, 1,3, 2,3, 3,4, 3,4, -1);
    print_result(&g, 2, IGRAPH_DIRECTED);

    VERIFY_FINALLY_STACK();

    printf("\nCheck negative order error.\n");
    CHECK_ERROR(igraph_graph_power(&g, NULL, -1, IGRAPH_DIRECTED), IGRAPH_EINVAL);
    igraph_destroy(&g);

    return 0;
}
