/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_vector_t deg;
    igraph_bool_t is_simple;

    igraph_set_error_handler(&igraph_error_handler_ignore);

    igraph_vector_init(&deg, 0);

    /* k-regular undirected graph, even degrees, no multiple edges */
    igraph_k_regular_game(&g, 10, 4, 0, 0);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_ALL, 1);
    igraph_vector_print(&deg);
    igraph_is_simple(&g, &is_simple);
    if (!is_simple) {
        return 1;
    }
    if (igraph_is_directed(&g)) {
        return 1;
    }
    igraph_destroy(&g);

    /* k-regular undirected graph, odd degrees, even number of vertices, no multiple edges */
    igraph_k_regular_game(&g, 10, 3, 0, 0);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_ALL, 1);
    igraph_vector_print(&deg);
    igraph_is_simple(&g, &is_simple);
    if (!is_simple) {
        return 2;
    }
    if (igraph_is_directed(&g)) {
        return 2;
    }
    igraph_destroy(&g);

    /* k-regular undirected graph, odd degrees, odd number of vertices, no multiple edges */
    if (!igraph_k_regular_game(&g, 9, 3, 0, 0)) {
        return 3;
    }

    /* k-regular undirected graph, even degrees, multiple edges */
    igraph_k_regular_game(&g, 10, 4, 0, 1);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_ALL, 1);
    igraph_vector_print(&deg);
    if (igraph_is_directed(&g)) {
        return 14;
    }
    igraph_destroy(&g);

    /* k-regular undirected graph, odd degrees, even number of vertices, multiple edges */
    igraph_k_regular_game(&g, 10, 3, 0, 1);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_ALL, 1);
    igraph_vector_print(&deg);
    if (igraph_is_directed(&g)) {
        return 15;
    }
    igraph_destroy(&g);

    /* k-regular undirected graph, odd degrees, odd number of vertices, multiple edges */
    if (!igraph_k_regular_game(&g, 9, 3, 0, 1)) {
        return 4;
    }

    /* k-regular directed graph, even degrees, no multiple edges */
    igraph_k_regular_game(&g, 10, 4, 1, 0);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    igraph_is_simple(&g, &is_simple);
    if (!is_simple) {
        return 5;
    }
    if (!igraph_is_directed(&g)) {
        return 5;
    }
    igraph_destroy(&g);

    /* k-regular directed graph, odd degrees, even number of vertices, no multiple edges */
    igraph_k_regular_game(&g, 10, 3, 1, 0);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    igraph_is_simple(&g, &is_simple);
    if (!is_simple) {
        return 6;
    }
    if (!igraph_is_directed(&g)) {
        return 6;
    }
    igraph_destroy(&g);

    /* k-regular directed graph, odd degrees, odd number of vertices, no multiple edges */
    igraph_k_regular_game(&g, 9, 3, 1, 0);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    igraph_is_simple(&g, &is_simple);
    if (!is_simple) {
        return 7;
    }
    if (!igraph_is_directed(&g)) {
        return 7;
    }
    igraph_destroy(&g);

    /* k-regular directed graph, even degrees, multiple edges */
    igraph_k_regular_game(&g, 10, 4, 1, 1);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    if (!igraph_is_directed(&g)) {
        return 16;
    }
    igraph_destroy(&g);

    /* k-regular directed graph, odd degrees, even number of vertices, multiple edges */
    igraph_k_regular_game(&g, 10, 3, 1, 1);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    if (!igraph_is_directed(&g)) {
        return 17;
    }
    igraph_destroy(&g);

    /* k-regular directed graph, odd degrees, odd number of vertices, multiple edges */
    igraph_k_regular_game(&g, 9, 3, 1, 1);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_IN, 1);
    igraph_vector_print(&deg);
    igraph_degree(&g, &deg, igraph_vss_all(), IGRAPH_OUT, 1);
    igraph_vector_print(&deg);
    if (!igraph_is_directed(&g)) {
        return 18;
    }
    igraph_destroy(&g);

    /* k-regular undirected graph, too large degree, no multiple edges */
    if (!igraph_k_regular_game(&g, 10, 10, 0, 0)) {
        return 8;
    }

    /* k-regular directed graph, too large degree, no multiple edges */
    if (!igraph_k_regular_game(&g, 10, 10, 1, 0)) {
        return 9;
    }

    /* empty graph */
    if (igraph_k_regular_game(&g, 0, 0, 0, 0)) {
        return 10;
    }
    if (igraph_vcount(&g) != 0 || igraph_ecount(&g) != 0 || igraph_is_directed(&g)) {
        return 11;
    }
    igraph_destroy(&g);
    if (igraph_k_regular_game(&g, 0, 0, 1, 0)) {
        return 12;
    }
    if (igraph_vcount(&g) != 0 || igraph_ecount(&g) != 0 || !igraph_is_directed(&g)) {
        return 13;
    }
    igraph_destroy(&g);

    igraph_vector_destroy(&deg);

    VERIFY_FINALLY_STACK();

    return 0;
}

