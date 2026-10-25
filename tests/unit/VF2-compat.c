/*
   igraph library.
   Copyright (C) 2012-2026  The igraph development team <igraph@igraph.org>

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

/* ----------------------------------------------------------- */

/* Vertices/edges with the same parity match */
igraph_bool_t compat_parity(const igraph_t *graph1,
                            const igraph_t *graph2,
                            const igraph_int_t g1_num,
                            const igraph_int_t g2_num,
                            void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    return (g1_num % 2) == (g2_num % 2);
}

/* Vertex/edge 0 in graph1 never matches */
igraph_bool_t compat_not0(const igraph_t *graph1,
                          const igraph_t *graph2,
                          const igraph_int_t g1_num,
                          const igraph_int_t g2_num,
                          void *arg) {
    IGRAPH_UNUSED(graph1);
    IGRAPH_UNUSED(graph2);
    IGRAPH_UNUSED(arg);
    IGRAPH_UNUSED(g2_num);
    return g1_num != 0;
}

void match_rings(void) {
    igraph_t r1, r2;
    igraph_bool_t iso;
    igraph_int_t count;
    igraph_ring(&r1, 10, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ true);
    igraph_ring(&r2, 10, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ true);

    igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                          &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                          /*node_compat_fn=*/ NULL, /*edge_compat_fn=*/ NULL,
                          /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                          &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                          compat_parity, /*edge_compat_fn=*/ NULL, /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                          &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                          compat_not0, /*edge_compat_fn=*/ NULL, /*arg=*/ NULL);
    IGRAPH_ASSERT(!iso);

    /* ------- */

    igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                          &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                          /*node_compat_fn=*/ NULL, compat_parity, /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_isomorphic_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                          &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                          /*node_compat_fn=*/ NULL, compat_not0, /*arg=*/ NULL);
    IGRAPH_ASSERT(!iso);

    /* ------- */

    igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                  &count, /*node_compat_fn=*/ NULL,
                                  /*edge_compat_fn=*/ NULL, /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 20);

    igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                  &count, compat_parity, /*edge_compat_fn=*/ NULL,
                                  /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 10);

    igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                  &count, compat_not0, /*edge_compat_fn=*/ NULL,
                                  /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 0);

    /* ------- */

    igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                  &count, /*node_compat_fn=*/ NULL, compat_parity,
                                  /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 10);

    igraph_count_isomorphisms_vf2(&r1, &r2, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                  &count, /*node_compat_fn=*/ NULL, compat_not0,
                                  /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 0);

    igraph_destroy(&r1);
    igraph_destroy(&r2);
}

void match_rings_open_closed(void) {
    igraph_t ro, rc;
    igraph_bool_t iso;
    igraph_int_t count;
    igraph_ring(&ro, 10, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ false);
    igraph_ring(&rc, 10, IGRAPH_UNDIRECTED, /*mutual=*/ false, /*circular=*/ true);

    igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                             &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                             /*node_compat_fn=*/ NULL, /*edge_compat_fn=*/ NULL,
                             /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                             &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                             compat_parity, /*edge_compat_fn=*/ NULL,
                             /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                             &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                             compat_not0, /*edge_compat_fn=*/ NULL,
                             /*arg=*/ NULL);
    IGRAPH_ASSERT(!iso);

    /* ------- */

    igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                             &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                             /*node_compat_fn=*/ NULL, compat_parity,
                             /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    igraph_subisomorphic_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                             &iso, /*map12=*/ NULL, /*map21=*/ NULL,
                             /*node_compat_fn=*/ NULL, compat_not0,
                             /*arg=*/ NULL);
    IGRAPH_ASSERT(iso);

    /* ------- */

    igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                     &count, /*node_compat_fn=*/ NULL,
                                     /*edge_compat_fn=*/ NULL, /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 20);

    igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                     &count, compat_parity,
                                     /*edge_compat_fn=*/ NULL, /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 10);

    igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                     &count, compat_not0, /*edge_compat_fn=*/ NULL,
                                     /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 0);

    /* ------- */

    igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                     &count, /*node_compat_fn=*/ NULL,
                                     compat_parity, /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 10);

    igraph_count_subisomorphisms_vf2(&rc, &ro, /*colors(4x)*/ NULL, NULL, NULL, NULL,
                                     &count, /*node_compat_fn=*/ NULL, compat_not0,
                                     /*arg=*/ NULL);

    IGRAPH_ASSERT(count == 2);

    igraph_destroy(&ro);
    igraph_destroy(&rc);
}

/* ----------------------------------------------------------- */

int main(void) {
    match_rings();
    match_rings_open_closed();

    VERIFY_FINALLY_STACK();

    return 0;
}
