/*
azert hallod, az a resze nem kis mertekben baszott fel, hogy "majd ha tobbet leteszunk az asztalra"
*/

/*
   IGraph library.
   Copyright (C) 2024  The igraph development team

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

#include <stdio.h>
#include "igraph.h"

int main(int argc, char* argv[]) {
    igraph_t g;
    igraph_error_t result;
    FILE *ifile;

    igraph_set_error_handler(igraph_error_handler_ignore);

    igraph_set_attribute_table(&igraph_cattribute_table);

    ifile = fopen("bug_2506_1.graphml", "r");
    IGRAPH_ASSERT(ifile != NULL);

    result = igraph_read_graph_graphml(&g, ifile, 0);
    fclose(ifile);

    if (result == IGRAPH_UNIMPLEMENTED) {
        /* maybe it is simply disabled at compile-time */
        return 77;
    }

    IGRAPH_ASSERT(result == IGRAPH_PARSEERROR);
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    ifile = fopen("bug_2506_2.graphml", "r");
    IGRAPH_ASSERT(ifile != NULL);

    result = igraph_read_graph_graphml(&g, ifile, 0);
    fclose(ifile);

    IGRAPH_ASSERT(result == IGRAPH_PARSEERROR);
    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;
}
