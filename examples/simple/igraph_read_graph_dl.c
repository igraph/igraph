/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <stdlib.h>

int main(void) {

    const char *files[] = { "fullmatrix1.dl", "fullmatrix2.dl",
                            "fullmatrix3.dl", "fullmatrix4.dl",
                            "edgelist1.dl", "edgelist2.dl", "edgelist3.dl",
                            "edgelist4.dl", "edgelist5.dl", "edgelist6.dl",
                            "edgelist7.dl", "nodelist1.dl", "nodelist2.dl" };
    igraph_t graph;
    FILE *infile;

    /* Turn on attribute handling. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    for (size_t i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
        printf("Doing %s\n", files[i]);
        infile = fopen(files[i], "r");
        if (!infile) {
            printf("Cannot open file: %s\n", files[i]);
            abort();
        }
        igraph_read_graph_dl(&graph, infile, IGRAPH_DIRECTED);
        fclose(infile);
        igraph_write_graph_edgelist(&graph, stdout);
        igraph_destroy(&graph);
    }

    if (IGRAPH_FINALLY_STACK_SIZE() != 0) {
        return 1;
    }

    return 0;
}
