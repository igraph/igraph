/*
   igraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

    /* empty directed graph, zero vertices */
    igraph_empty(&g, 0, 1);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /* empty undirected graph, zero vertices */
    igraph_empty(&g, 0, 0);
    IGRAPH_ASSERT(igraph_vcount(&g) == 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /* empty directed graph, 20 vertices */
    igraph_empty(&g, 20, 1);
    IGRAPH_ASSERT(igraph_vcount(&g) == 20);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /* empty undirected graph, 30 vertices */
    igraph_empty(&g, 30, 0);
    IGRAPH_ASSERT(igraph_vcount(&g) == 30);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /* error: negative number of vertices */
    CHECK_ERROR(igraph_empty(&g, -1, 0), IGRAPH_EINVAL);

    VERIFY_FINALLY_STACK();

    return 0;
}
