/*
   igraph library.
   Copyright (C) 2010-2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

/* Brute force test of minimality. */
igraph_bool_t is_minimal_sep(const igraph_t *graph, const igraph_vector_int_t *S) {
    igraph_bool_t res;
    igraph_is_separator(graph, igraph_vss_vector(S), &res);

    if (res) {
        for (igraph_int_t i=0; i < igraph_vector_int_size(S); i++) {
            igraph_bool_t sep;
            igraph_vector_int_t S2;
            igraph_vector_int_init_copy(&S2, S);
            igraph_vector_int_remove_fast(&S2, i);
            igraph_is_separator(graph, igraph_vss_vector(&S2), &sep);
            igraph_vector_int_destroy(&S2);
            if (sep) {
                res = false;
                break;
            }
        }
    }

    return res;
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_t S;
    igraph_bool_t is_sep, is_min;

    /* Simple star graph, remove the center */
    igraph_star(&graph, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_is_separator(&graph, igraph_vss_1(0), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_1(0), &is_min);
    IGRAPH_ASSERT(is_min);

    /* Same graph, but another vertex */
    igraph_is_separator(&graph, igraph_vss_1(6), &is_sep);
    IGRAPH_ASSERT(! is_sep);

    /* Same graph, all vertices but the center */
    igraph_is_separator(&graph, igraph_vss_range(1, 10), &is_sep);
    IGRAPH_ASSERT(! is_sep);

    /* Same graph, all vertices */
    igraph_is_separator(&graph, igraph_vss_range(0, 10), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_destroy(&graph);

    /* Same graph, no vertices */
    igraph_is_separator(&graph, igraph_vss_range(0, 0), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_destroy(&graph);

    /* Test graph 2 */

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 2,0, 0,3, 4,5, 5,6, 1,1, 5,6,
                 -1);

    igraph_vector_int_init_int_end(&S, -1,
                                   0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   5, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   1, 0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    /* Pass in vertices multiple times */
    igraph_vector_int_init_int_end(&S, -1,
                                   1, 3, 0, 1, 3, 3, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(!is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   4, 1, 0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   4, 1, 0, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   4, 1, 0, 5, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   1, 0, 5, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   5, 6, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   5, 0, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_is_separator(&graph, igraph_vss_range(0, 0), &is_sep);
    IGRAPH_ASSERT(!is_sep);

    igraph_destroy(&graph);

    /* Test graph 3 */

    igraph_small(&graph, 6, IGRAPH_UNDIRECTED,
                 0,1, 1,3, 0,2, 2,3, 2,4,
                 -1);

    igraph_vector_int_init_int_end(&S, -1,
                                   5, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(! is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   1, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);

    igraph_add_edge(&graph, 1, 4);

    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    /* Pass in vertices multiple times */
    igraph_vector_int_init_int_end(&S, -1,
                                   1, 2, 1, 2, 2, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_destroy(&graph);

    /* Test graph 4 */

    igraph_small(&graph, 0, IGRAPH_UNDIRECTED,
                 0,1, 0,2, 0,3, 1,2, 1,3, 2,3,
                 4,5, 4,6, 4,7, 5,6, 5,7, 6,7,
                 2,5, 3,4,
                 -1);

    igraph_vector_int_init_int_end(&S, -1,
                                   2, 3, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   2, 4, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   3, 5, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   2, 3, 5, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(! is_min);
    igraph_vector_int_destroy(&S);

    igraph_destroy(&graph);

    /* Karate club network */

    igraph_famous(&graph, "Zachary");

    igraph_vector_int_init_int_end(&S, -1, 32, 33, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_sep);
    IGRAPH_ASSERT(is_sep);
    igraph_is_minimal_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(is_min);
    igraph_vector_int_destroy(&S);

    igraph_vector_int_init_int_end(&S, -1,
                                   8, 9, 19, 30, 31, -1);
    igraph_is_separator(&graph, igraph_vss_vector(&S), &is_min);
    IGRAPH_ASSERT(!is_min);
    igraph_vector_int_destroy(&S);

    /* Caution:
     * igraph_all_minimal_st_separators() returns minimal (s,t) separators,
     * i.e. separators that are minimal with respect to separating a specific
     * (s,t) pair, but not necessarily minimal separators. See the documentation
     * for an example.
     */

    igraph_vector_int_list_t separators;
    igraph_vector_int_list_init(&separators, 0);
    igraph_all_minimal_st_separators(&graph, &separators);
    for (igraph_int_t i=0; i < igraph_vector_int_list_size(&separators); i++) {
        const igraph_vector_int_t *sep = igraph_vector_int_list_get_ptr(&separators, i);
        igraph_is_separator(&graph, igraph_vss_vector(sep), &is_sep);
        IGRAPH_ASSERT(is_sep);
        igraph_is_minimal_separator( &graph, igraph_vss_vector(sep), &is_min);

        igraph_bool_t is_min2 = is_minimal_sep(&graph, sep);
        IGRAPH_ASSERT((!is_min) == (! is_min2));
    }
    igraph_vector_int_list_destroy(&separators);

    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
