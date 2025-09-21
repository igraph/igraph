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

void print_and_destroy(igraph_t *g, igraph_int_t source, igraph_int_t target, igraph_vconn_nei_t neighbors) {
    igraph_int_t res;
    igraph_st_vertex_connectivity(g, &res, source, target, neighbors);
    printf("%" IGRAPH_PRId "\n", res);
    igraph_destroy(g);
}

int main(void) {
    igraph_t g;
    igraph_int_t res;

    printf("graph with two unconnected vertices: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, -1);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_ERROR);

    printf("graph with two connected vertices, IGRAPH_VCONN_NEI_NEGATIVE: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_NEGATIVE);

    printf("graph with two connected vertices, IGRAPH_VCONN_NEI_NUMBER_OF_NODES: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_NUMBER_OF_NODES);

    printf("graph with two connected vertices, IGRAPH_VCONN_NEI_IGNORE: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, 0,1, 0,1, -1);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_IGNORE);

    printf("directed graph with two connected vertices, IGRAPH_VCONN_NEI_IGNORE: ");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, 0,1, 0,1, 1,0, 1,0, -1);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_IGNORE);

    printf("line graph with 6 vertices: ");
    igraph_small(&g, 6, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,5, -1);
    print_and_destroy(&g, 0, 5, IGRAPH_VCONN_NEI_ERROR);

    printf("full graph with 6 vertices IGRAPH_VCONN_NEI_IGNORE: ");
    igraph_full(&g, 6, IGRAPH_UNDIRECTED, 0);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_IGNORE);

    printf("full graph with 6 vertices IGRAPH_VCONN_NEI_IGNORE, directed: ");
    igraph_full(&g, 6, IGRAPH_DIRECTED, 0);
    print_and_destroy(&g, 0, 1, IGRAPH_VCONN_NEI_IGNORE);

    printf("line graph with 3 vertices, 6 edges: ");
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 0,1, 1,2, 1,2, 1,2, 1,2, -1);
    print_and_destroy(&g, 0, 2, IGRAPH_VCONN_NEI_ERROR);

    VERIFY_FINALLY_STACK();

    printf("check error on graph with two connected vertices.\n");
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0,1, -1);
    CHECK_ERROR(igraph_st_vertex_connectivity(&g, &res, 0, 0, IGRAPH_VCONN_NEI_ERROR), IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("check error on graph with no vertices.\n");
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, -1);
    CHECK_ERROR(igraph_st_vertex_connectivity(&g, &res, 0, 0, IGRAPH_VCONN_NEI_ERROR), IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("check error on graph with one vertex.\n");
    igraph_small(&g, 1, IGRAPH_UNDIRECTED, -1);
    CHECK_ERROR(igraph_st_vertex_connectivity(&g, &res, 0, 0, IGRAPH_VCONN_NEI_ERROR), IGRAPH_EINVAL);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
}
