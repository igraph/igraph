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

    /* Create an undirected 6-star, with the 0th node as the centre. */
    igraph_star(&graph, 7, IGRAPH_STAR_UNDIRECTED, 0);

    /* Output the edge list of the graph. */
    igraph_write_graph_edgelist(&graph, stdout);

    /* Destroy the graph when we are done using it. */
    igraph_destroy(&graph);

    return 0;
}
