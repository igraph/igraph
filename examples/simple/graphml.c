/* -*- mode: C -*-  */
/*
   IGraph library.
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
#include <unistd.h>     /* unlink */

int main() {
    igraph_t g;
    igraph_error_handler_t* oldhandler;
    int result;
    FILE *ifile, *ofile;

    igraph_set_attribute_table(&igraph_cattribute_table);

    /* GraphML */
    ifile = fopen("test.graphml", "r");
    if (ifile == 0) {
        return 10;
    }

    /* GraphML support is an optional feature so we disable the error handler
     * temporarily to handle the case when it is not implemented */
    oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
    if ((result = igraph_read_graph_graphml(&g, ifile, 0))) {
        /* maybe it is simply disabled at compile-time */
        if (result == IGRAPH_UNIMPLEMENTED) {
            return 77;
        }
        return 1;
    }
    igraph_set_error_handler(oldhandler);

    fclose(ifile);

    /* Write it back */
    ofile = fopen("test2.graphml", "w");
    if (ofile) {
        if ((result = igraph_write_graph_graphml(&g, ofile, /*prefixattr=*/ 1))) {
            printf("Received unexpected return code: %d\n", result);
            return 1;
        }
        fclose(ofile);

        /* Clean up after ourselves */
        unlink("test2.graphml");
    }

    /* Destroy the graph */
    igraph_destroy(&g);

    return 0;
}
