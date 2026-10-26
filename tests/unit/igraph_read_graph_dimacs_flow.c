/*
   igraph library.
   Copyright (C) 2026  The igraph development team

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
#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t label;
    FILE *ifile;

    igraph_vector_int_init(&label, 0);

    /* EDGE problem with node labels. Node indices are 1-based in the file. */
    ifile = fopen("dimacs_edge_labels.dimacs", "r");
    IGRAPH_ASSERT(ifile != NULL);
    IGRAPH_ASSERT(igraph_read_graph_dimacs_flow(&g, ifile, NULL, &label, NULL, NULL, NULL, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    fclose(ifile);
    print_graph(&g);
    print_vector_int(&label);
    igraph_destroy(&g);

    igraph_error_handler_t *handler = igraph_set_error_handler(&igraph_error_handler_printignore);

    /* Node line without a label. */
    ifile = fopen("dimacs_edge_nolabel.dimacs", "r");
    IGRAPH_ASSERT(ifile != NULL);
    IGRAPH_ASSERT(igraph_read_graph_dimacs_flow(&g, ifile, NULL, &label, NULL, NULL, NULL, IGRAPH_UNDIRECTED) == IGRAPH_PARSEERROR);
    fclose(ifile);

    /* Node index out of range. */
    ifile = fopen("dimacs_edge_badindex.dimacs", "r");
    IGRAPH_ASSERT(ifile != NULL);
    IGRAPH_ASSERT(igraph_read_graph_dimacs_flow(&g, ifile, NULL, &label, NULL, NULL, NULL, IGRAPH_UNDIRECTED) == IGRAPH_PARSEERROR);
    fclose(ifile);

    igraph_set_error_handler(handler);

    igraph_vector_int_destroy(&label);

    VERIFY_FINALLY_STACK();

    return 0;
}
