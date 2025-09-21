/*
   igraph library.
   Copyright (C) 2023  The igraph development team <igraph@igraph.org>

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

void print_and_destroy(igraph_t *g1, igraph_t *g2) {
    igraph_bool_t iso;
    igraph_isomorphic(g1, g2, &iso);
    printf("%d\n", iso);
    igraph_destroy(g1);
    igraph_destroy(g2);
}

int main(void) {
    igraph_t g1, g2;
    igraph_bool_t iso;
    printf("Two multigraphs: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,1, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,2, 1,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Only second graph is multigraph: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,2, 1,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Only first graph is multigraph: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,1, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Two multigraphs, first has loop: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,0, 0,1, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,2, 1,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Two multigraphs, both have loops: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,0, 0,1, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,1, 1,2, 1,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Two multigraphs, only loops: ");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,0, 0,0, 1,1, -1);
    igraph_small(&g2, 3, IGRAPH_UNDIRECTED, 1,1, 2,2, 2,2, -1);
    print_and_destroy(&g1, &g2);

    printf("Check error for two multigraphs, different directedness\n");
    igraph_small(&g1, 3, IGRAPH_UNDIRECTED, 0,1, 0,1, -1);
    igraph_small(&g2, 3, IGRAPH_DIRECTED, 1,2, 1,2, -1);
    CHECK_ERROR(igraph_isomorphic(&g1, &g2, &iso), IGRAPH_EINVAL);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    VERIFY_FINALLY_STACK();
    return 0;
}
