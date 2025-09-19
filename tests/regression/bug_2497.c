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

#include "../unit/test_utilities.h"

int main(void) {
    FILE *file;
    igraph_t graph;

    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    file = fopen("bug_2497.gml", "r");
    IGRAPH_ASSERT(file != NULL);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);
    IGRAPH_ASSERT(igraph_read_graph_gml(&graph, file) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    fclose(file);

    igraph_write_graph_gml(&graph, stdout, IGRAPH_WRITE_GML_DEFAULT_SW, NULL, "");
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
