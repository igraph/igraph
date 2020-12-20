/* -*- mode: C -*-  */
/*
  Test suite for whether two vertices are connected by an edge.
  Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>

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
#include <stdio.h>

#include "test_utilities.inc"

#define R_INTEGER(a,b) (igraph_rng_get_integer(igraph_rng_default(), (a), (b)))

/* Crash the library function here. We expect error codes to be returned here.
 */
int error_test() {
    igraph_t g;
    igraph_bool_t connected;
    igraph_integer_t nvert, u, v;
    int ret;

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */
    igraph_small(&g, /*nvert*/ 0, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    nvert = igraph_vcount(&g);
    u = (igraph_integer_t)R_INTEGER(-100 * nvert, 100 * nvert);
    v = (igraph_integer_t)R_INTEGER(nvert, 100 * nvert);

    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_are_connected(&g, u, v, &connected);
    if (ret != IGRAPH_EINVVID) {
        printf("Error test failed.\n");
        return IGRAPH_FAILURE;
    }
    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}

/* Testing for two vertices being connected by an edge in various graphs.
 */
int connected_test() {
    igraph_t gcomplete, gempty;
    igraph_bool_t connected;
    igraph_integer_t nvert, u, v;

    igraph_rng_seed(igraph_rng_default(), 57); /* make tests deterministic */

    /* A complete graph on n vertices. Any two distinct vertices are connected */
    /* by an edge. Hence we expect the test to return true for any given pair */
    /* of distinct vertices. */
    nvert = (igraph_integer_t)R_INTEGER(2, 100);
    igraph_full(&gcomplete, nvert, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    u = (igraph_integer_t)R_INTEGER(0, nvert - 1);
    do {
        v = (igraph_integer_t)R_INTEGER(0, nvert - 1);
    } while (v == u);
    igraph_are_connected(&gcomplete, u, v, &connected);
    if (!connected) {
        printf("Expected connected = true, but received connected = false.\n");
        return IGRAPH_FAILURE;
    }
    igraph_destroy(&gcomplete);

    /* A graph with n vertices, but no edges. Any two distinct vertices are */
    /* not joined by an edge. Thus we expect the test to return false for any */
    /* given pair of distinct vertices. */
    nvert = (igraph_integer_t)R_INTEGER(2, 100);
    igraph_empty(&gempty, nvert, IGRAPH_DIRECTED);
    u = (igraph_integer_t)R_INTEGER(0, nvert - 1);
    do {
        v = (igraph_integer_t)R_INTEGER(0, nvert - 1);
    } while (v == u);
    igraph_are_connected(&gempty, u, v, &connected);
    if (connected) {
        printf("Expected connected = false, but received connected = true.\n");
        return IGRAPH_FAILURE;
    }
    igraph_destroy(&gempty);

    return IGRAPH_SUCCESS;
}

int main() {
    int ret;

    ret = error_test();
    if (ret) {
        return 1;
    }
    ret = connected_test();
    if (ret) {
        return 1;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
