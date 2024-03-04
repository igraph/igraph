/*
   IGraph library.
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

void print(const char *str) {
    igraph_t g;
    igraph_from_nauty(&g, str);
    print_graph_canon(&g);
    igraph_destroy(&g);
}

int main(void) {
    //graph6
    const char *str1 = "Gr`HOk";

    //digraph6
    const char *singleton = "&@?";
    const char *one_vertex_with_loop = "&@_";

    printf("some graph6 example\n");
    print(str1);

    printf("\ndigraph6 singleton\n");
    print(singleton);

    printf("\ndigraph6 single vertex with loop\n");
    print(one_vertex_with_loop);

    VERIFY_FINALLY_STACK();
    return 0;
}
