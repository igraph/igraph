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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {

    igraph_t graph;
    igraph_vector_int_t sep;
    igraph_bool_t result;

    /* Simple star graph, remove the center */
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_is_separator(&graph, igraph_vss_1(0), &result);
    IGRAPH_ASSERT(result);

    /* Same graph, but another vertex */
    igraph_is_separator(&graph, igraph_vss_1(6), &result);
    IGRAPH_ASSERT(! result);

    /* Same graph, all vertices but the center */
    igraph_is_separator(&graph, igraph_vss_range(1, 10), &result);
    IGRAPH_ASSERT(! result);

    /* Same graph, all vertices */
    igraph_is_separator(&graph, igraph_vss_range(0, 10), &result);
    IGRAPH_ASSERT(! result);
    igraph_destroy(&graph);

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 0,3, 4,5, 5,6,
                 -1);

    igraph_vector_int_init_int_end(&sep, -1,
                                   1, 0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   1, 3, 0, 3, 3, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(!result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   4, 1, 0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   4, 1, 0, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(! result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   4, 1, 0, 5, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(! result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   1, 0, 5, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(result);
    igraph_vector_int_destroy(&sep);

    igraph_vector_int_init_int_end(&sep, -1,
                                   5, 6, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(!result);
    igraph_vector_int_destroy(&sep);

    igraph_destroy(&graph);

    /* Karate club */
    igraph_famous(&graph, "zachary");
    igraph_vector_int_init(&sep, 0);
    igraph_vector_int_push_back(&sep, 32);
    igraph_vector_int_push_back(&sep, 33);
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(result);

    igraph_vector_int_resize(&sep, 5);
    VECTOR(sep)[0] = 8;
    VECTOR(sep)[1] = 9;
    VECTOR(sep)[2] = 19;
    VECTOR(sep)[3] = 30;
    VECTOR(sep)[4] = 31;
    igraph_is_separator(&graph, igraph_vss_vector(&sep), &result);
    IGRAPH_ASSERT(!result);

    igraph_destroy(&graph);
    igraph_vector_int_destroy(&sep);

    return 0;
}
