/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

void print_vector_int(igraph_vector_int_t *v, FILE *f) {
    igraph_integer_t i;
    for (i = 0; i < igraph_vector_int_size(v); i++) {
        fprintf(f, " %" IGRAPH_PRId, VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

int check_simple() {

    igraph_t g;
    igraph_integer_t nodes = 100;
    igraph_integer_t edges = 1000;
    igraph_real_t p = 3.0 / nodes;
    igraph_integer_t runs = 10;
    igraph_integer_t r, e, ecount;
    igraph_vector_int_t eids, pairs, path;

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */

    igraph_vector_int_init(&pairs, edges * 2);
    igraph_vector_int_init(&path, 0);
    igraph_vector_int_init(&eids, 0);

    for (r = 0; r < runs; r++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, nodes, p,
                                /*directed=*/ 0, /*loops=*/ 0);
        ecount = igraph_ecount(&g);
        for (e = 0; e < edges; e++) {
            igraph_integer_t edge = RNG_INTEGER(0, ecount - 1);
            VECTOR(pairs)[2 * e] = IGRAPH_FROM(&g, edge);
            VECTOR(pairs)[2 * e + 1] = IGRAPH_TO(&g, edge);
        }
        igraph_get_eids(&g, &eids, &pairs, /*path=*/ 0, 0, /*error=*/ 1);
        for (e = 0; e < edges; e++) {
            igraph_integer_t edge = VECTOR(eids)[e];
            igraph_integer_t from1 = VECTOR(pairs)[2 * e];
            igraph_integer_t to1 = VECTOR(pairs)[2 * e + 1];
            igraph_integer_t from2 = IGRAPH_FROM(&g, edge);
            igraph_integer_t to2 = IGRAPH_TO(&g, edge);
            igraph_integer_t min1 = from1 < to1 ? from1 : to1;
            igraph_integer_t max1 = from1 < to1 ? to1 : from1;
            igraph_integer_t min2 = from2 < to2 ? from2 : to2;
            igraph_integer_t max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 11;
            }
        }

        igraph_diameter(&g, /*res=*/ 0, /*from=*/ 0, /*to=*/ 0, &path, NULL,
                        IGRAPH_UNDIRECTED, /*unconn=*/ 1);
        igraph_get_eids(&g, &eids, /*pairs=*/ 0, &path, 0, /*error=*/ 1);
        for (e = 0; e < igraph_vector_int_size(&path) - 1; e++) {
            igraph_integer_t edge = VECTOR(eids)[e];
            igraph_integer_t from1 = VECTOR(path)[e];
            igraph_integer_t to1 = VECTOR(path)[e + 1];
            igraph_integer_t from2 = IGRAPH_FROM(&g, edge);
            igraph_integer_t to2 = IGRAPH_TO(&g, edge);
            igraph_integer_t min1 = from1 < to1 ? from1 : to1;
            igraph_integer_t max1 = from1 < to1 ? to1 : from1;
            igraph_integer_t min2 = from2 < to2 ? from2 : to2;
            igraph_integer_t max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 12;
            }
        }

        igraph_destroy(&g);
    }

    igraph_vector_int_destroy(&path);
    igraph_vector_int_destroy(&pairs);
    igraph_vector_int_destroy(&eids);

    return 0;
}

int check_multi() {

    igraph_t g;
    igraph_vector_int_t vec;
    igraph_vector_int_t eids;
    igraph_vector_int_t eids2;
    igraph_error_t ret;
    igraph_integer_t i;

    igraph_integer_t q1[] = { 0, 1, 0, 1 };
    igraph_integer_t q2[] = { 0, 1, 0, 1, 0, 1 };
    igraph_integer_t q3[] = { 1, 0, 3, 4, 1, 0, 0, 1, 3, 4, 0, 1 };

    igraph_vector_int_init(&eids, 0);

    /*********************************/
    igraph_small(&g, /*n=*/ 10, /*directed=*/ 1,
                 0, 1, 0, 1, 1, 0, 1, 2, 3, 4, 3, 4, 3, 4, 3, 5, 3, 7,
                 9, 8,
                 -1);

    igraph_vector_int_view(&vec, q1, sizeof(q1) / sizeof(q1[0]));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 1, /*error=*/ 1);
    igraph_vector_int_sort(&eids);
    print_vector_int(&eids, stdout);

    igraph_vector_int_view(&vec, q2, sizeof(q2) / sizeof(q2[0]));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 0, /*error=*/ 1);
    igraph_vector_int_sort(&eids);
    print_vector_int(&eids, stdout);

    igraph_vector_int_view(&vec, q2, sizeof(q2) / sizeof(q2[0]));
    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 1, /*error=*/1);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }
    igraph_set_error_handler(igraph_error_handler_abort);

    igraph_destroy(&g);
    /*********************************/

    /*********************************/
    igraph_small(&g, /*n=*/10, /*directed=*/0,
                 0, 1, 1, 0, 0, 1, 3, 4, 3, 4, 5, 4, 9, 8,
                 -1);

    igraph_vector_int_view(&vec, q1, sizeof(q1) / sizeof(q1[0]));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/1, /*error=*/ 1);
    igraph_vector_int_sort(&eids);
    print_vector_int(&eids, stdout);

    igraph_vector_int_view(&vec, q3, sizeof(q3) / sizeof(q3[0]));
    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/0, /*error=*/ 1);
    if (ret != IGRAPH_EINVAL) {
        return 2;
    }
    igraph_set_error_handler(igraph_error_handler_abort);

    igraph_destroy(&g);

    /*********************************/

    igraph_vector_int_destroy(&eids);

    /*********************************/
    /* Speed tests */

#define NODES 10000
    igraph_barabasi_game(&g, /*n=*/ NODES, /*power=*/ 1.0, /*m=*/ 3,
                         /*outseq=*/ 0, /*outpref=*/ 0, /*A=*/ 1,
                         /*directed=*/ 1, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

    igraph_vector_int_init(&eids2, NODES / 2);
    igraph_random_sample(&eids2, 0, igraph_ecount(&g) - 1, NODES / 2);
    igraph_vector_int_init(&vec, NODES);
    for (i = 0; i < NODES / 2; i++) {
        VECTOR(vec)[2 * i]   = IGRAPH_FROM(&g, VECTOR(eids2)[i]);
        VECTOR(vec)[2 * i + 1] = IGRAPH_TO(&g, VECTOR(eids2)[i]);
    }
    igraph_vector_int_init(&eids, 0);
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 1, /*error=*/ 1);
    for (i = 0; i < NODES / 2; i++) {
        if (VECTOR(eids)[i] != VECTOR(eids2)[i]) {
            return 3;
        }
    }

    /**/

    for (i = 0; i < NODES / 2; i++) {
        VECTOR(vec)[2 * i]   = IGRAPH_TO(&g, VECTOR(eids)[i]);
        VECTOR(vec)[2 * i + 1] = IGRAPH_FROM(&g, VECTOR(eids)[i]);
    }
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 0, /*error=*/ 1);
    for (i = 0; i < NODES / 2; i++) {
        if (VECTOR(eids)[i] != VECTOR(eids2)[i]) {
            return 4;
        }
    }

    igraph_vector_int_destroy(&eids);
    igraph_vector_int_destroy(&eids2);
    igraph_vector_int_destroy(&vec);
    igraph_destroy(&g);

    /*********************************/

    return 0;
}

int main() {
    int ret;

    if ( (ret = check_simple()) != 0) {
        return ret;
    }
    if ( (ret = check_multi()) != 0)  {
        return ret;
    }

    return 0;
}
