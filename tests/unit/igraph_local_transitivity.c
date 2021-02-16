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

#include "test_utilities.inc"

int main() {

    igraph_t g;
    igraph_vs_t vertices;
    igraph_vector_t result1, result2;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_init(&result1, 0);
    igraph_vector_init(&result2, 0);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, 100, .1,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_vs_seq(&vertices, 0, 99);

    igraph_transitivity_local_undirected(&g, &result1, igraph_vss_all(),
                                         IGRAPH_TRANSITIVITY_NAN);
    igraph_transitivity_local_undirected(&g, &result2, vertices,
                                         IGRAPH_TRANSITIVITY_NAN);

    if (!igraph_vector_all_e(&result1, &result2)) {
        igraph_vector_print(&result1);
        igraph_vector_print(&result2);
        return 1;
    }

    igraph_vector_destroy(&result1);
    igraph_vector_destroy(&result2);
    igraph_vs_destroy(&vertices);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
