/* -*- mode: C -*-  */
/*
   IGraph library.
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

int main() {

    igraph_t g;
    int ret;

    /* empty directed graph, zero vertices */
    igraph_empty(&g, 0, 1);
    if (igraph_vcount(&g) != 0) {
        return 1;
    }
    if (igraph_ecount(&g) != 0) {
        return 2;
    }
    igraph_destroy(&g);

    /* empty undirected graph, zero vertices */
    igraph_empty(&g, 0, 0);
    if (igraph_vcount(&g) != 0) {
        return 3;
    }
    if (igraph_ecount(&g) != 0) {
        return 4;
    }
    igraph_destroy(&g);

    /* empty directed graph, 20 vertices */
    igraph_empty(&g, 20, 1);
    if (igraph_vcount(&g) != 20) {
        return 5;
    }
    if (igraph_ecount(&g) != 0) {
        return 6;
    }
    igraph_destroy(&g);

    /* empty undirected graph, 30 vertices */
    igraph_empty(&g, 30, 0);
    if (igraph_vcount(&g) != 30) {
        return 7;
    }
    if (igraph_ecount(&g) != 0) {
        return 8;
    }
    igraph_destroy(&g);

    /* error: negative number of vertices */
    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_empty(&g, -1, 0);
    if (ret != IGRAPH_EINVAL) {
        return 9;
    }

    return 0;
}
