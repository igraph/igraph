/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/* This is the exact 11-vertex, 13-edge toy citation graph used by
 * fast-cdindex's own regression tests (tests/tests.py's ctimes/cedges,
 * also used by src/main.cpp), with vertex IDs 0..10 matching igraph's
 * natural numbering. The expected CD-index and mCD-index values below are
 * taken verbatim from fast-cdindex's own golden output (tests/py.ok /
 * tests/bin.ok) for time_window = 157680000. The I-index has no golden
 * value upstream (fast-cdindex's own tests never print iindex()); the
 * expected values here were derived by hand from the same timestamps and
 * window, and cross-checked against mCD-index / CD-index from py.ok
 * wherever CD-index is nonzero. */

int main(void) {
    igraph_t g;
    igraph_vector_int_t timestamps;
    igraph_vector_t cd, i_index, mcd;
    static const igraph_int_t ts[] = {
        694224000, 694224000, 725846400, 725846400, 788918400,
        852098400, 883612800, 915148800, 915148800, 883612800, 852076800
    };
    const igraph_int_t time_window = 157680000;

    igraph_small(&g, 11, IGRAPH_DIRECTED,
                 4, 2, 4, 0, 4, 1, 4, 3, 5, 2, 6, 2, 6, 4, 7, 4, 8, 4, 9, 4, 9, 1, 9, 3, 10, 4,
                 -1);
    igraph_vector_int_init_array(&timestamps, ts, sizeof(ts) / sizeof(ts[0]));

    igraph_vector_init(&cd, 0);
    igraph_vector_init(&i_index, 0);
    igraph_vector_init(&mcd, 0);

    /* All vertices at once. */
    igraph_cd_index(&g, &timestamps, &cd, &i_index, &mcd, igraph_vss_all(), time_window);
    printf("CD index:  "); print_vector(&cd);
    printf("I index:   "); print_vector(&i_index);
    printf("mCD index: "); print_vector(&mcd);

    /* Single vertex via a vs selector; must match the vertex-4 entries above. */
    igraph_cd_index(&g, &timestamps, &cd, &i_index, &mcd, igraph_vss_1(4), time_window);
    printf("\nVertex 4 only:\n");
    printf("CD index:  "); print_vector(&cd);
    printf("I index:   "); print_vector(&i_index);
    printf("mCD index: "); print_vector(&mcd);

    /* Error: undirected graph. */
    {
        igraph_t ug;
        igraph_copy(&ug, &g);
        igraph_to_undirected(&ug, IGRAPH_TO_UNDIRECTED_EACH, NULL);
        CHECK_ERROR(igraph_cd_index(&ug, &timestamps, &cd, NULL, NULL, igraph_vss_all(), time_window), IGRAPH_EINVAL);
        igraph_destroy(&ug);
    }

    /* Error: timestamp vector length mismatch. */
    {
        igraph_vector_int_t bad_ts;
        igraph_vector_int_init(&bad_ts, 3);
        CHECK_ERROR(igraph_cd_index(&g, &bad_ts, &cd, NULL, NULL, igraph_vss_all(), time_window), IGRAPH_EINVAL);
        igraph_vector_int_destroy(&bad_ts);
    }

    /* Error: self-loop present. */
    {
        igraph_t lg;
        igraph_copy(&lg, &g);
        igraph_add_edge(&lg, 0, 0);
        CHECK_ERROR(igraph_cd_index(&lg, &timestamps, &cd, NULL, NULL, igraph_vss_all(), time_window), IGRAPH_EINVAL);
        igraph_destroy(&lg);
    }

    igraph_vector_destroy(&mcd);
    igraph_vector_destroy(&i_index);
    igraph_vector_destroy(&cd);
    igraph_vector_int_destroy(&timestamps);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
