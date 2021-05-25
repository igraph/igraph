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

typedef struct {
    int dim;
    int m;
    int nei;
    igraph_bool_t directed, mutual, circular;
    igraph_real_t *dimedges;
} lat_test_t;

#define LAT_TEST(id, d, m, ne, di, mu, ci, ...) \
    igraph_real_t lat_ ## id ## _edges[] = { __VA_ARGS__ } ; \
    lat_test_t lat_ ## id = { d, m, ne, di, mu, ci, lat_ ## id ## _edges }

/*----------------d--m--ne-di-mu-ci-dimedges------------------------*/
LAT_TEST(u_0,     0, 0, 1, 0, 0, 0,  -1 );
LAT_TEST(u_01,    1, 0, 1, 0, 0, 0,  0 );
LAT_TEST(u_02,    2, 0, 1, 0, 0, 0,  0, 1 );
LAT_TEST(u_03,    2, 0, 1, 0, 0, 0,  1, 0 );

LAT_TEST(d_0,     0, 0, 1, 1, 0, 0,  -1 );
LAT_TEST(d_01,    1, 0, 1, 1, 0, 0,  0 );
LAT_TEST(d_02,    2, 0, 1, 1, 0, 0,  0, 1 );
LAT_TEST(d_03,    2, 0, 1, 1, 0, 0,  1, 0 );

LAT_TEST(u_1,     1, 0, 1, 0, 0, 0,  1 );
LAT_TEST(u_1x1,   2, 0, 1, 0, 0, 0,  1, 1 );
LAT_TEST(u_2,     1, 1, 1, 0, 0, 0,  2, 0, 1 );
LAT_TEST(u_2x1,   2, 1, 1, 0, 0, 0,  2, 1, 0, 1 );
LAT_TEST(u_2x2,   2, 4, 1, 0, 0, 0,  2, 2, 0, 1, 0, 2, 1, 3, 2, 3 );

LAT_TEST(uc_1,    1, 0, 1, 0, 0, 1,  1 );
LAT_TEST(uc_1x1,  2, 0, 1, 0, 0, 1,  1, 1 );
LAT_TEST(uc_2,    1, 1, 1, 0, 0, 1,  2, 0, 1 );
LAT_TEST(uc_2x1,  2, 1, 1, 0, 0, 1,  2, 1, 0, 1 );
LAT_TEST(uc_2x2,  2, 4, 1, 0, 0, 1,  2, 2, 0, 1, 0, 2, 1, 3, 2, 3 );

LAT_TEST(dc_1,    1, 0, 1, 1, 0, 1,  1 );
LAT_TEST(dc_1x1,  2, 0, 1, 1, 0, 1,  1, 1 );
LAT_TEST(dc_2,    1, 2, 1, 1, 0, 1,  2, 0, 1, 1, 0 );
LAT_TEST(dc_2x1,  2, 2, 1, 1, 0, 1,  2, 1, 0, 1, 1, 0 );
LAT_TEST(dc_2x2,  2, 8, 1, 1, 0, 1,  2, 2, 0, 1, 0, 2, 1, 3, 2, 3,
         1, 0, 2, 0, 3, 1, 3, 2, );

LAT_TEST(d_1,     1, 0, 1, 1, 0, 0,  1 );
LAT_TEST(d_1x1,   2, 0, 1, 1, 0, 0,  1, 1 );
LAT_TEST(d_2,     1, 1, 1, 1, 0, 0,  2, 0, 1 );
LAT_TEST(d_2x1,   2, 1, 1, 1, 0, 0,  2, 1, 0, 1 );
LAT_TEST(d_2x2,   2, 4, 1, 1, 0, 0,  2, 2, 0, 1, 0, 2, 1, 3, 2, 3 );

LAT_TEST(dmc_1,   1, 0, 1, 1, 0, 1,  1 );
LAT_TEST(dmc_1x1, 2, 0, 1, 1, 0, 1,  1, 1 );
LAT_TEST(dmc_2,   1, 2, 1, 1, 0, 1,  2, 0, 1, 1, 0 );
LAT_TEST(dmc_2x1, 2, 2, 1, 1, 0, 1,  2, 1, 0, 1, 1, 0 );
LAT_TEST(dmc_2x2, 2, 4, 1, 1, 0, 1,  2, 2, 0, 1, 0, 2, 1, 3, 2, 3,
         1, 0, 3, 2, );
/*----------------d--m--ne-di-mu-ci-dimedges------------------------*/

/* TODO: add more */

lat_test_t *all_checks[] = { /*  1 */ &lat_u_0,   /*  2 */ &lat_u_01,
                                      /*  3 */ &lat_u_02,  /*  4 */ &lat_u_03,
                                      /*  5 */ &lat_d_0,   /*  6 */ &lat_d_01,
                                      /*  7 */ &lat_d_02,  /*  8 */ &lat_d_03,
                                      /*  9 */ &lat_u_1,   /* 10 */ &lat_u_1x1,
                                      /* 11 */ &lat_u_2,   /* 12 */ &lat_u_2x1,
                                      /* 13 */ &lat_u_2x2, /* 14 */ &lat_u_1,
                                      /* 15 */ &lat_u_1x1, /* 16 */ &lat_u_2,
                                      /* 17 */ &lat_u_2x1, /* 18 */ &lat_uc_2x2,
                                      /* 19 */ &lat_dc_1,  /* 20 */ &lat_dc_1x1,
                                      /* 21 */ &lat_dc_2,  /* 22 */ &lat_dc_2x1,
                                      /* 23 */ &lat_dc_2x2,/* 24 */ &lat_d_1,
                                      /* 25 */ &lat_d_1x1, /* 26 */ &lat_d_2,
                                      /* 27 */ &lat_d_2x1, /* 28 */ &lat_d_2x2,
                                      /* 29 */ &lat_dc_2x2,/* 30 */ &lat_d_1,
                                      /* 31 */ &lat_d_1x1, /* 32 */ &lat_d_2,
                                      /* 33 */ &lat_d_2x1, /* 34 */ &lat_d_2x2,
                                      0
                           };

int check_lattice_properties(const igraph_t *lattice,
                             const igraph_vector_t *dim,
                             igraph_bool_t directed,
                             igraph_bool_t mutual,
                             igraph_bool_t circular) {
    igraph_bool_t res;

    /* Connected */
    igraph_is_connected(lattice, &res, IGRAPH_WEAK);
    if (!res && igraph_vcount(lattice) > 0) {
        printf("Not connected\n");
        return 1;
    }

    /* Simple */
    igraph_is_simple(lattice, &res);
    if (!res) {
        printf("Not simple\n");
        return 2;
    }

    return 0;
}

int check_lattice(const lat_test_t *test) {
    igraph_t graph, othergraph;
    igraph_vector_t otheredges;
    igraph_vector_t dimvector;
    igraph_bool_t iso;
    int ret;

    /* Create lattice */
    igraph_vector_view(&dimvector, test->dimedges, test->dim);
    igraph_lattice(&graph, &dimvector, test->nei, test->directed,
                   test->mutual, test->circular);

    /* Check its properties */
    if ((ret = check_lattice_properties(&graph, &dimvector, test->directed,
                                        test->mutual, test->circular))) {
        igraph_destroy(&graph);
        printf("Lattice properties are not satisfied\n");
        return ret;
    }

    /* Check that it is isomorphic to the stored graph */
    igraph_vector_view(&otheredges, test->dimedges + test->dim, test->m * 2);
    igraph_create(&othergraph, &otheredges, igraph_vector_prod(&dimvector),
                  test->directed);
    igraph_isomorphic(&graph, &othergraph, &iso);
    if (!iso) {
        printf("--\n");
        igraph_write_graph_edgelist(&graph, stdout);
        printf("--\n");
        igraph_write_graph_edgelist(&othergraph, stdout);
        igraph_destroy(&graph);
        igraph_destroy(&othergraph);
        return 50;
    }

    igraph_destroy(&graph);
    igraph_destroy(&othergraph);
    return 0;
}

int main() {
    int i, ret;

    i = 0;
    while (all_checks[i]) {
        if ((ret = check_lattice(all_checks[i]))) {
            printf("Check no #%d failed.\n", (int) (i + 1));
            return ret;
        }
        i++;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
