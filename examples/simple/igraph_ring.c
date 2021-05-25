/*
   IGraph library.
   Copyright (C) 2020  The igraph development team <igraph@igraph.org>

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
#include <stdio.h>

int main() {
    igraph_t graph;

    /* Create a directed path graph on 10 vertices. */
    igraph_ring(&graph, 10, IGRAPH_DIRECTED, /* mutual= */ 0, /* circular= */ 0);

    /* Output the edge list of the graph. */
    printf("10-path graph:\n");
    igraph_write_graph_edgelist(&graph, stdout);

    /* Destroy the graph. */
    igraph_destroy(&graph);

    /* Create a 4-cycle graph. */
    igraph_ring(&graph, 4, IGRAPH_UNDIRECTED, /* mutual= */ 0, /* circular= */ 1);

    /* Output the edge list of the graph. */
    printf("\n4-cycle graph:\n");
    igraph_write_graph_edgelist(&graph, stdout);

    /* Destroy the graph. */
    igraph_destroy(&graph);

    return 0;
}
