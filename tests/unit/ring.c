/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   Ring test suite
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

#include "test_utilities.inc"

typedef struct {
    int n, m;
    igraph_bool_t directed, mutual, circular;
    igraph_real_t *edges;
} ring_test_t;

#define RING_TEST(id, n, m, di, mu, ci, ...) \
    igraph_real_t ring_ ## id ## _edges[] = { __VA_ARGS__ };      \
    ring_test_t ring_ ## id = { n, m, di, mu, ci, ring_ ## id ## _edges }

/*---------------n--m--di-mu-ci--edges-------------------------------------*/
RING_TEST(uc_6,  6, 6, 0, 0, 1,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0 );
RING_TEST(uc_0,  0, 0, 0, 0, 1,  -1 );
RING_TEST(uc_1,  1, 0, 0, 0, 1,  -1 );
RING_TEST(uc_2,  2, 1, 0, 0, 1,  0, 1 );

RING_TEST(u_6,   6, 5, 0, 0, 0,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5 );
RING_TEST(u_0,   0, 0, 0, 0, 0,  -1 );
RING_TEST(u_1,   1, 0, 0, 0, 0,  -1 );
RING_TEST(u_2,   2, 1, 0, 0, 0,  0, 1 );

RING_TEST(umc_6, 6, 6, 0, 1, 1,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0 );
RING_TEST(umc_0, 0, 0, 0, 1, 1,  -1 );
RING_TEST(umc_1, 1, 0, 0, 1, 1,  -1 );
RING_TEST(umc_2, 2, 1, 0, 1, 1,  0, 1 );

RING_TEST(um_6,  6, 5, 0, 1, 0,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5 );
RING_TEST(um_0,  0, 0, 0, 1, 0,  -1 );
RING_TEST(um_1,  1, 0, 0, 1, 0,  -1 );
RING_TEST(um_2,  2, 1, 0, 1, 0,  0, 1 );

RING_TEST(dc_6,  6, 6, 1, 0, 1,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0 );
RING_TEST(dc_0,  0, 0, 1, 0, 1,  -1 );
RING_TEST(dc_1,  1, 0, 1, 0, 1,  -1 );
RING_TEST(dc_2,  2, 2, 1, 0, 1,  0, 1, 1, 0 );

RING_TEST(d_6,   6, 5, 1, 0, 1,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5 );
RING_TEST(d_0,   0, 0, 1, 0, 1,  -1 );
RING_TEST(d_1,   1, 0, 1, 0, 1,  -1 );
RING_TEST(d_2,   2, 1, 1, 0, 1,  0, 1 );

RING_TEST(dmc_6,  6, 12, 1, 1, 1, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0,
          1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 0, 5 );
RING_TEST(dmc_0,  0, 0, 1, 1, 1, -1 );
RING_TEST(dmc_1,  1, 0, 1, 1, 1, -1 );
RING_TEST(dmc_2,  2, 2, 1, 1, 1, 0, 1, 1, 0 );

RING_TEST(dm_6,  6, 10, 1, 1, 0,  0, 1, 1, 2, 2, 3, 3, 4, 4, 5,
          1, 0, 2, 1, 3, 2, 4, 3, 5, 4 );
RING_TEST(dm_0,  0, 0, 1, 1, 0,  -1 );
RING_TEST(dm_1,  1, 0, 1, 1, 0,  -1 );
RING_TEST(dm_2,  2, 2, 1, 1, 0,  0, 1, 1, 0 );
/*---------------n--m--di-mu-ci--edges-------------------------------------*/

ring_test_t *all_checks[] = { /*  1 */ &ring_uc_6,   /*  2 */ &ring_uc_0,
                                       /*  3 */ &ring_uc_1,   /*  4 */ &ring_uc_2,
                                       /*  5 */ &ring_u_6,    /*  6 */ &ring_u_0,
                                       /*  7 */ &ring_u_1,    /*  8 */ &ring_u_2,
                                       /*  9 */ &ring_umc_6,  /* 10 */ &ring_umc_0,
                                       /* 11 */ &ring_umc_1,  /* 12 */ &ring_umc_2,
                                       /* 13 */ &ring_um_6,   /* 14 */ &ring_um_0,
                                       /* 15 */ &ring_um_1,   /* 16 */ &ring_um_2,
                                       /* 17 */ &ring_dc_6,   /* 18 */ &ring_dc_0,
                                       /* 19 */ &ring_dc_1,   /* 20 */ &ring_dc_2,
                                       /* 21 */ &ring_dmc_6,  /* 22 */ &ring_dmc_0,
                                       /* 23 */ &ring_dmc_1,  /* 24 */ &ring_dmc_2,
                                       /* 25 */ &ring_dm_6,   /* 26 */ &ring_dm_0,
                                       /* 27 */ &ring_dm_1,   /* 28 */ &ring_dm_2,
                                       0
                            };

int check_ring_properties(const igraph_t *ring, igraph_bool_t directed,
                          igraph_bool_t mutual, igraph_bool_t circular) {

    igraph_bool_t res;

    /* Connected */
    igraph_is_connected(ring, &res, IGRAPH_WEAK);
    if (!res && igraph_vcount(ring) > 0) {
        printf("Not connected\n");
        return 1;
    }

    /* Simple */
    igraph_is_simple(ring, &res);
    if (!res) {
        printf("Not simple\n");
        return 2;
    }

    /* Girth, for big enough circular graphs */
    if (circular && igraph_vcount(ring) > 2) {
        igraph_integer_t girth;
        igraph_girth(ring, &girth, NULL);
        if (girth != igraph_vcount(ring)) {
            printf("Wrong girth\n");
            return 3;
        }
    }

    return 0;
}

int check_ring(const ring_test_t *test) {
    igraph_t graph, othergraph;
    igraph_vector_t otheredges;
    igraph_bool_t iso;
    int ret;

    /* Create ring */
    igraph_ring(&graph, test->n, test->directed, test->mutual, test->circular);

    /* Check its properties */
    if ((ret = check_ring_properties(&graph, test->directed, test->mutual,
                                     test->circular))) {
        return ret;
    }

    /* Check that it is isomorphic to the stored graph */
    igraph_vector_view(&otheredges, test->edges, test->m * 2);
    igraph_create(&othergraph, &otheredges, test->n, test->directed);
    igraph_isomorphic(&graph, &othergraph, &iso);
    if (!iso) {
        return 50;
    }

    /* Clean up */
    igraph_destroy(&graph);
    igraph_destroy(&othergraph);

    return 0;
}

int main() {
    int i, ret;

    i = 0;
    while (all_checks[i]) {
        if ((ret = check_ring(all_checks[i]))) {
            printf("Check no #%d failed.\n", (int) (i + 1));
            return ret;
        }
        i++;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
