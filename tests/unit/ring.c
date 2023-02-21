/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2022  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    int directed, circular, mutual;
    int n;

    for (n=0; n < 4; ++n) {
        for (directed=0; directed < 2; ++directed) {
            for (circular=0; circular < 2; ++circular) {
                for (mutual=0; mutual < 2; ++mutual) {
                    igraph_t graph;
                    igraph_ring(&graph, n, directed, mutual, circular);
                    printf("n=%d, %s, %s, %s\n", n,
                           directed ? "directed" : "undirected",
                           circular ? "cycle" : "path",
                           mutual ? "bidirectional" : "unidirectional"
                           );
                    print_graph_canon(&graph);
                    printf("\n");
                    igraph_destroy(&graph);
                }
            }
        }
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
