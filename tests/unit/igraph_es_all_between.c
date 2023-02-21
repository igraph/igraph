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
#include <stdlib.h>

#include "test_utilities.h"

igraph_error_t test_directed(void) {
    igraph_t g;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_integer_t size;

    IGRAPH_CHECK(igraph_small(&g, 5, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 3, 1, 2, 3, 4, 1, 2, 4, 0, -1));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 3, 4));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));

    /* single edge between 0 and 1 */
    igraph_es_all_between(&es, 0, 1, IGRAPH_DIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 1);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 0);
        IGRAPH_ASSERT(to == 1);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* no edge between 3 and 0 */
    igraph_es_all_between(&es, 3, 0, IGRAPH_DIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 0);
    IGRAPH_ASSERT(IGRAPH_EIT_END(eit));
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* many edges between 1 and 2 */
    igraph_es_all_between(&es, 1, 2, IGRAPH_DIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 8);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 1);
        IGRAPH_ASSERT(to == 2);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* two edges between 3 and 4, using IGRAPH_UNDIRECTED */
    igraph_es_all_between(&es, 4, 3, IGRAPH_UNDIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 2);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 3);
        IGRAPH_ASSERT(to == 4);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return IGRAPH_SUCCESS;
}

igraph_error_t test_undirected(void) {
    igraph_t g;
    igraph_es_t es;
    igraph_eit_t eit;
    igraph_integer_t size;

    IGRAPH_CHECK(igraph_small(&g, 5, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, 1, 2, 3, 4, 1, 2, 4, 0, -1));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 2, 1));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));
    IGRAPH_CHECK(igraph_add_edge(&g, 3, 4));
    IGRAPH_CHECK(igraph_add_edge(&g, 2, 1));
    IGRAPH_CHECK(igraph_add_edge(&g, 1, 2));

    /* single edge between 0 and 1 */
    igraph_es_all_between(&es, 0, 1, IGRAPH_UNDIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 1);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 0);
        IGRAPH_ASSERT(to == 1);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* no edge between 3 and 0 */
    igraph_es_all_between(&es, 3, 0, IGRAPH_UNDIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 0);
    IGRAPH_ASSERT(IGRAPH_EIT_END(eit));
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* many edges between 1 and 2 */
    igraph_es_all_between(&es, 1, 2, IGRAPH_UNDIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 8);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 1);
        IGRAPH_ASSERT(to == 2);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    /* two edges between 3 and 4 */
    igraph_es_all_between(&es, 4, 3, IGRAPH_UNDIRECTED);
    IGRAPH_CHECK(igraph_eit_create(&g, es, &eit));
    IGRAPH_CHECK(igraph_es_size(&g, &es, &size));
    IGRAPH_ASSERT(size == 2);
    while (!IGRAPH_EIT_END(eit)) {
        igraph_integer_t edge = IGRAPH_EIT_GET(eit);
        igraph_integer_t from, to;
        IGRAPH_CHECK(igraph_edge(&g, edge, &from, &to));
        IGRAPH_ASSERT(from == 3);
        IGRAPH_ASSERT(to == 4);
        IGRAPH_EIT_NEXT(eit);
        size--;
    }
    IGRAPH_ASSERT(size == 0);
    igraph_eit_destroy(&eit);
    igraph_es_destroy(&es);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return IGRAPH_SUCCESS;
}

int main(void) {
    IGRAPH_ASSERT(test_undirected() == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(test_directed() == IGRAPH_SUCCESS);

    return 0;
}
