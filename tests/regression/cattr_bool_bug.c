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

#include "../unit/test_utilities.h"

void check_attr(igraph_t *graph) {
    IGRAPH_ASSERT(igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "name"));
    IGRAPH_ASSERT(igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "type"));
    IGRAPH_ASSERT(igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_GRAPH, "p"));
    IGRAPH_ASSERT(igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, "name"));
    IGRAPH_ASSERT(igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, "weight"));
}

int main(void) {

    igraph_t graph;
    igraph_error_handler_t* oldhandler;
    int result;
    FILE *ifile = fopen("cattr_bool_bug.graphml", "r");

    if (!ifile) {
        printf("Cannot open input file\n");
        return 1;
    }

    igraph_set_attribute_table(&igraph_cattribute_table);

    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if ((result = igraph_read_graph_graphml(&graph, ifile, 0))) {
        /* Maybe it is simply disabled at compile-time. If so, skip test. */
        if (result == IGRAPH_UNIMPLEMENTED) {
            return 77;
        }
        printf("Failed to read GraphML file\n");
        return 1;
    }
    igraph_set_error_handler(oldhandler);

    fclose(ifile);

    printf("Checkng attributes of original graph\n");
    check_attr(&graph);

    printf("Checkng attributes of directed version\n");
    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_ARBITRARY);
    check_attr(&graph);

    IGRAPH_ASSERT(! GAB(&graph, "loops"));

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
