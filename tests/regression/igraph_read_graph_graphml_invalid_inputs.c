/*
   IGraph library.
   Copyright (C) 2021  The igraph development team

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

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "igraph.h"
#include <stdio.h>

#include "../unit/test_utilities.h"

int test_file(const char* fname, igraph_bool_t should_parse) {
    FILE *ifile;
    igraph_t g;
    igraph_error_t retval;

    ifile = fopen(fname, "r");
    if (ifile == 0) {
        printf("Cannot open input file: %s\n", fname);
        return 1;
    }

    retval = igraph_read_graph_graphml(&g, ifile, 0);
    if (retval == IGRAPH_SUCCESS) {
        /* destroy the graph if we managed to parse it */
        igraph_destroy(&g);
    }

    fclose(ifile);

    if (!should_parse && retval == IGRAPH_SUCCESS) {
        /* input was accepted but it should not have been, this is a bug */
        return 2;
    } else {
        return 0;
    }
}

#undef RUN_TEST
#define RUN_TEST(fname, should_parse) {   \
    index++;                              \
    if (test_file(fname, should_parse)) { \
        return index;                     \
    }                                     \
    VERIFY_FINALLY_STACK();               \
}

int main(void) {
    int index = 0;

    /* We do not care about errors; all we care about is that the library
     * should not segfault, should not accept invalid input and should not
     * print anything to stdout or stderr while parsing, apart from warnings */
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);

    RUN_TEST("invalid1.graphml", /* should_parse = */ 0);
    RUN_TEST("invalid2.graphml", /* should_parse = */ 1);
    RUN_TEST("invalid3.graphml", /* should_parse = */ 0);
    RUN_TEST("invalid4.graphml", /* should_parse = */ 0);
    RUN_TEST("invalid5.graphml", /* should_parse = */ 0);

    return 0;
}
