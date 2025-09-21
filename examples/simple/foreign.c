/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
    igraph_t g;
    FILE *ifile;

    /* Initialize the library. */
    igraph_setup();

    /* Turn on attribute handling. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Read a Pajek file. */
    ifile = fopen("links.net", "r");
    if (ifile == 0) {
        return 10;
    }
    igraph_read_graph_pajek(&g, ifile);
    fclose(ifile);

    /* Write it in edgelist format. */
    printf("The graph:\n");
    printf("Vertices: %" IGRAPH_PRId "\n", igraph_vcount(&g));
    printf("Edges: %" IGRAPH_PRId "\n", igraph_ecount(&g));
    printf("Directed: %i\n", igraph_is_directed(&g) ? 1 : 0);
    igraph_write_graph_edgelist(&g, stdout);
    igraph_destroy(&g);

    return 0;
}
