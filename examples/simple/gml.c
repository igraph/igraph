/*
   igraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
#include <stdio.h>

int main(void) {
    igraph_t graph;
    FILE *ifile;

    /* Initialize the library. */
    igraph_setup();

    ifile = fopen("karate.gml", "r");
    if (ifile == 0) {
        return 1;
    }

    igraph_read_graph_gml(&graph, ifile);
    fclose(ifile);

    printf("The graph is %s.\n", igraph_is_directed(&graph) ? "directed" : "undirected");

    /* Output as edge list */
    printf("\n-----------------\n");
    igraph_write_graph_edgelist(&graph, stdout);

    /* Output as GML */
    printf("\n-----------------\n");
    igraph_write_graph_gml(&graph, stdout,IGRAPH_WRITE_GML_DEFAULT_SW,  0, "");

    igraph_destroy(&graph);

    return 0;
}
