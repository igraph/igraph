/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

int main(void) {
    igraph_t left, right, joined;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&left, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, -1);
    igraph_small(&right, 5, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,4, -1);

    igraph_join(&joined, &left, &right);
    igraph_write_graph_edgelist(&joined, stdout);
    printf("\n");

    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&joined);


    igraph_small(&left, 2, IGRAPH_DIRECTED, 0,1, -1);
    igraph_small(&right, 3, IGRAPH_DIRECTED, 0,1, 2,1, -1);

    igraph_join(&joined, &left, &right);
    igraph_write_graph_edgelist(&joined, stdout);
    printf("\n");

    igraph_destroy(&left);
    igraph_destroy(&right);
    igraph_destroy(&joined);

    return 0;
}
