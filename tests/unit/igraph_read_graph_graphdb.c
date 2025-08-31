/*
   igraph library.
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/
#include <igraph.h>
#include <stdio.h>
#include "test_utilities.h"

int main(void) {
    igraph_t g;
    FILE *ifile;

    ifile = fopen("si2_b06m_s20.A98", "rb");
    IGRAPH_ASSERT(ifile != NULL);

    igraph_read_graph_graphdb(&g, ifile, IGRAPH_DIRECTED);
    fclose(ifile);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    print_graph(&g);
    igraph_destroy(&g);

    igraph_error_handler_t *handler = igraph_set_error_handler(&igraph_error_handler_printignore);

    /* Node count too low, extra bytes at end. */
    ifile = fopen("si2_b06m_s20-bad1.A98", "rb");
    IGRAPH_ASSERT(ifile != NULL);
    IGRAPH_ASSERT(igraph_read_graph_graphdb(&g, ifile, IGRAPH_DIRECTED) == IGRAPH_PARSEERROR);
    fclose(ifile);

    /* Truncated file. */
    ifile = fopen("si2_b06m_s20-bad2.A98", "rb");
    IGRAPH_ASSERT(ifile != NULL);
    IGRAPH_ASSERT(igraph_read_graph_graphdb(&g, ifile, IGRAPH_DIRECTED) == IGRAPH_PARSEERROR);
    fclose(ifile);

    igraph_set_error_handler(handler);

    VERIFY_FINALLY_STACK();

    return 0;
}
