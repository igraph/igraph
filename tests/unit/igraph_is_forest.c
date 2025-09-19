/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

void check_output(const igraph_t *graph, igraph_bool_t *res, igraph_neimode_t mode){
    igraph_bool_t result;
    igraph_vector_int_t roots;

    igraph_vector_int_init(&roots, 0);
    igraph_is_forest(graph, &result, &roots, mode);
    if (result) {
        printf("Root nodes: ");
        igraph_vector_int_print(&roots);
    }
    else {
        printf("Not a forest.\n");
    }
    IGRAPH_ASSERT(*res == result);
    igraph_vector_int_destroy(&roots);
    printf("\n");
}

int main(void) {
    igraph_t graph;
    igraph_bool_t res;
    igraph_neimode_t mode;

    printf("Empty Graph\n");
    mode=IGRAPH_ALL;
    res=1;
    igraph_empty(&graph, 0, 0);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Graph with 0 edges\n");
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Graph with 1 vertex\n");
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph, 1, IGRAPH_UNDIRECTED, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Undirected Graph\n");
    mode=IGRAPH_ALL;
    res=1;
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 3, 4, 3, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Directed Graph out trees\n");
    mode=IGRAPH_OUT;
    res=1;
    igraph_small(&graph, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 3, 4, 3, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Directed Graph in trees\n");
    mode=IGRAPH_IN;
    res=0;
    igraph_small(&graph, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 3, 4, 3, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Undirected Graph with cycle\n");
    mode=IGRAPH_ALL;
    res=0;
    igraph_small(&graph, 7, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 3, 4, 3, 5, 4, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Undirected Graph with ecount>vcount-1\n");
    mode=IGRAPH_ALL;
    res=0;
    igraph_small(&graph, 6, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 3, 4, 3, 5, 4, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Undirected Graph with self-loop\n");
    mode=IGRAPH_ALL;
    res=0;
    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 4, 3, 3, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Directed Multigraph\n");
    mode=IGRAPH_OUT;
    res=0;
    igraph_small(&graph, 6, IGRAPH_DIRECTED, 0, 1, 1, 2, 3, 4, 3, 5, 3, 5, -1);
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    printf("Disjoint union of a cycle with two trees, directed\n");
    mode=IGRAPH_OUT;
    res=0;
    igraph_small(&graph, 10, IGRAPH_DIRECTED,
                 0, 1, 1, 2, 1, 3, 4, 5, 5, 6, 6, 7, 7, 4, 8, 9,
                 -1);
    check_output(&graph, &res, mode);

    printf("Disjoint union of a cycle with two trees, undirected\n");
    mode=IGRAPH_ALL;
    res=0;
    check_output(&graph, &res, mode);
    igraph_destroy(&graph);

    /* Cache testing */
    /* 1 <- 0 -> 2 <- 3 */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0,1, 0,2, 3,2, -1);
    /* This must not cache that the graph is not a forest,
     * as we are only checking the directed case: */
    igraph_is_forest(&graph, &res, NULL, IGRAPH_OUT);
    IGRAPH_ASSERT(!res);
    igraph_is_forest(&graph, &res, NULL, IGRAPH_ALL);
    IGRAPH_ASSERT(res);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
