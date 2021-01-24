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
    igraph_matrix_t matrix;

    /* Create a small test graph. */
    igraph_small(&graph, 0, IGRAPH_DIRECTED,
                 0, 1, 2, 1, 2, 0, 3, 0,
                 -1);

    /* As usual with igraph functions, the data structure in which the result
       will be returned must be initialized in advance. */
    igraph_matrix_init(&matrix, 0, 0);
    igraph_bibcoupling(&graph, &matrix, igraph_vss_all());
    printf("Bibliographic coupling matrix:\n");
    igraph_matrix_print(&matrix);

    igraph_cocitation(&graph, &matrix, igraph_vss_all());
    printf("\nCocitation matrix:\n");
    igraph_matrix_print(&matrix);

    /* Destroy data structures when we are done with them. */
    igraph_matrix_destroy(&matrix);
    igraph_destroy(&graph);

    return 0;
}
