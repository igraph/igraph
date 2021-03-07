/*
   IGraph library.
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
#include "test_utilities.inc"

void print_result(igraph_t *g, igraph_integer_t start_vid) {
    igraph_real_t result;
    igraph_integer_t from, to;
    IGRAPH_ASSERT(igraph_pseudo_diameter(g, &result, start_vid, &from, &to) == IGRAPH_SUCCESS);
    printf("result:");
    print_real(stdout, result, "%g");
    printf(", from %" IGRAPH_PRId " to %" IGRAPH_PRId "\n\n", from, to);
}

int main() {
    igraph_t g;
    igraph_real_t result;
    int i;

    igraph_set_error_handler(igraph_error_handler_ignore);

    printf("No vertices, no allowed starting vertex id:\n");
    igraph_small(&g, 0, 0, -1);
    IGRAPH_ASSERT(igraph_pseudo_diameter(&g, &result, 0, NULL, NULL) == IGRAPH_EINVAL);
    igraph_destroy(&g);

    printf("1 vertex:\n");
    igraph_small(&g, 1, 0, -1);
    print_result(&g, 0);
    igraph_destroy(&g);

    printf("2 vertices:\n");
    igraph_small(&g, 2, 0, -1);
    print_result(&g, 0);
    igraph_destroy(&g);

    printf("Disconnected graph with loops and multiple edges.\n");
    igraph_small(&g, 6, 1, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
    for (i = 0; i < 6; i ++) {
        print_result(&g, i);
    }
    igraph_destroy(&g);

    printf("Petersen graph.\n");
    igraph_famous(&g, "petersen");
    for (i = 0; i < 6; i ++) {
        print_result(&g, i);
    }
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();
    return 0;
}
