/*
   igraph library.
   Copyright (C) 2006-2022  The igraph development team

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
#include <unistd.h> /* unlink */

int main(void) {
    igraph_t graph;
    const char *infilename  = "test.graphml";
    const char *outfilename = "test2.graphml";

    /* Initialize the library. */
    igraph_setup();

    /* Set up attribute handling, so graph attributes can be imported
     * from the GraphML file. */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Problems in the GraphML file may cause igraph to print warnings.
     * If this is not desired, set a silent warning handler: */
    igraph_set_warning_handler(&igraph_warning_handler_ignore);

    /* Read the contents of a GraphML file. */

    /* GraphML */
    FILE *infile = fopen("test.graphml", "r");
    if (! infile) {
        fprintf(stderr, "Could not open input file '%s'.", infilename);
        exit(1);
    }

    /* GraphML support is an optional feature in igraph. If igraph was compiled
     * without GraphML support, igraph_read_graph_graphml() returns IGRAPH_UNIMPLEMENTED.
     * We temporarily disable the default error handler so we can test for this condition. */
    igraph_error_handler_t *oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_error_t ret = igraph_read_graph_graphml(&graph, infile, 0);
    if (ret == IGRAPH_UNIMPLEMENTED) {
        fprintf(stderr, "igraph was compiled without GraphML support.");
        exit(77);
    }
    if (ret != IGRAPH_SUCCESS) {
        fprintf(stderr, "Unexpected error while reading GraphML.");
        exit(1);
    }
    igraph_set_error_handler(oldhandler);


    fclose(infile);

    /* Write it back into another file. */

    FILE *outfile = fopen(outfilename, "w");
    if (outfile) {
        igraph_write_graph_graphml(&graph, outfile, true);
        fclose(outfile);

        /* Clean up after ourselves */
        unlink(outfilename);
    } else {
        fprintf(stderr, "Could not write output file '%s'.", outfilename);
    }

    /* Destroy the graph */
    igraph_destroy(&graph);

    return 0;
}
