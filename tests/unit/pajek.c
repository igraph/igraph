/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    FILE *ifile;

    /* turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    ifile = fopen("pajek5.net", "r");
    IGRAPH_ASSERT(ifile != NULL);

    igraph_read_graph_pajek(&g, ifile);
    fclose(ifile);

    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 9);

    igraph_destroy(&g);

    ifile = fopen("pajek6.net", "r");
    IGRAPH_ASSERT(ifile != NULL);

    igraph_read_graph_pajek(&g, ifile);
    fclose(ifile);

    IGRAPH_ASSERT(igraph_vcount(&g) == 10);
    IGRAPH_ASSERT(igraph_ecount(&g) == 9);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
