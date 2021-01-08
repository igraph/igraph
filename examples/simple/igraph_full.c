/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2020  The igraph development team <igraph@igraph.org>

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
    long int n_vertices = 10;

    /* Create an undirected complete graph. */
    /* Use IGRAPH_UNDIRECTED and IGRAPH_NO_LOOPS instead of 1/TRUE and 0/FALSE for better readability. */
    igraph_full(&graph, n_vertices, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    printf("The undirected complete graph on %ld vertices has %ld edges.\n",
           (long int) igraph_vcount(&graph), (long int) igraph_ecount(&graph));

    /* Remember to destroy the object at the end. */
    igraph_destroy(&graph);

    /* Create a directed complete graph. */
    igraph_full(&graph, n_vertices, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    printf("The directed complete graph on %ld vertices has %ld edges.\n",
           (long int) igraph_vcount(&graph), (long int) igraph_ecount(&graph));

    igraph_destroy(&graph);

    /* Create an undirected complete graph with self-loops. */
    igraph_full(&graph, n_vertices, IGRAPH_UNDIRECTED, IGRAPH_LOOPS);
    printf("The undirected complete graph on %ld vertices with self-loops has %ld edges.\n",
           (long int) igraph_vcount(&graph), (long int) igraph_ecount(&graph));

    igraph_destroy(&graph);

    /* Create a directed graph with self-loops. */
    igraph_full(&graph, n_vertices, IGRAPH_DIRECTED, IGRAPH_LOOPS);
    printf("The directed complete graph on %ld vertices with self-loops has %ld edges.\n",
           (long int) igraph_vcount(&graph), (long int) igraph_ecount(&graph));

    igraph_destroy(&graph);

    return 0;

}
