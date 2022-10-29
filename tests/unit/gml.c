/*
   IGraph library.
   Copyright (C) 2022  The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#include <igraph.h>
#include <stdio.h>

#include "test_utilities.h"

void test_input(const char *filename) {
    igraph_t graph;
    FILE *ifile;

    ifile = fopen(filename, "r");
    if (! ifile) {
        fprintf(stderr, "Cannot open '%s'.\n", filename);
        abort();
    }

    printf("===== %s =====\n", filename);

    IGRAPH_ASSERT(igraph_read_graph_gml(&graph, ifile) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_write_graph_gml(&graph, stdout, IGRAPH_WRITE_GML_DEFAULT_SW, NULL, "igraph") == IGRAPH_SUCCESS);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    printf("======================\n\n");
}

int main(void) {

    /* Enable attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    test_input("graph1.gml");
    test_input("graph2.gml");
    test_input("graph3.gml");

    return 0;
}
