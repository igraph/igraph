/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

/* This test verifies that trying to read an empty file does not cause
 * a crash or fatal error in any of the file format readers. */

int main(void) {
    igraph_t graph;
    FILE *file;

    file = fopen("empty", "r");
    IGRAPH_ASSERT(file != NULL);

    /* Formats for which an emtpy file is valid */

    /* Edge list */
    IGRAPH_ASSERT(igraph_read_graph_edgelist(&graph, file, 0, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* NCOL */
    IGRAPH_ASSERT(igraph_read_graph_ncol(&graph, file, NULL, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* LGL */
    IGRAPH_ASSERT(igraph_read_graph_lgl(&graph, file, 1, IGRAPH_ADD_WEIGHTS_IF_PRESENT, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vcount(&graph) == 0);
    igraph_destroy(&graph);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* Formats for which an empty file is invalid */

    igraph_set_error_handler(&igraph_error_handler_printignore);

    /* GraphML */
    IGRAPH_ASSERT(igraph_read_graph_graphml(&graph, file, 0) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* Pajek */
    IGRAPH_ASSERT(igraph_read_graph_pajek(&graph, file) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* DL */
    IGRAPH_ASSERT(igraph_read_graph_dl(&graph, file, IGRAPH_UNDIRECTED) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* GML */
    IGRAPH_ASSERT(igraph_read_graph_gml(&graph, file) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* graphdb */
    IGRAPH_ASSERT(igraph_read_graph_graphdb(&graph, file, IGRAPH_UNDIRECTED) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    /* DIMACS flow problem */
    IGRAPH_ASSERT(igraph_read_graph_dimacs_flow(&graph, file, NULL, NULL, NULL, NULL, NULL, IGRAPH_DIRECTED) != IGRAPH_SUCCESS);
    rewind(file);
    VERIFY_FINALLY_STACK();

    return 0;
}
