/*
   igraph library.
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

#include <igraph.h>
#include <stdio.h>

int test_file(const char* fname) {
    FILE *ifile;
    igraph_t g;
    int retval;

    ifile = fopen(fname, "r");
    if (ifile == 0) {
        return 1;
    }

    retval = igraph_read_graph_pajek(&g, ifile);
    if (!retval) {
        /* input was accepted, this is a bug; attempt to clean up after
         * ourselves nevertheless */
        igraph_destroy(&g);
        fclose(ifile);
        return 2;
    }

    fclose(ifile);

    return 0;
}

#define RUN_TEST(fname) {   \
    index++;                \
    if (test_file(fname)) { \
        return index;       \
    }                       \
}

int main(void) {
    int index = 0;

    /* Turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* We do not care about errors; all we care about is that the library
     * should not segfault and should not accept invalid input either */
    igraph_set_error_handler(igraph_error_handler_ignore);

    RUN_TEST("invalid_pajek1.net");
    RUN_TEST("invalid_pajek2.net");
    RUN_TEST("invalid_pajek3.net");

    return 0;
}
