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

void print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

int check_simple() {

    igraph_t g;
    long int nodes = 100;
    long int edges = 1000;
    igraph_real_t p = 3.0 / nodes;
    long int runs = 10;
    long int r, e, ecount;
    igraph_vector_t eids, pairs, path;

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */

    igraph_vector_init(&pairs, edges * 2);
    igraph_vector_init(&path, 0);
    igraph_vector_init(&eids, 0);

    for (r = 0; r < runs; r++) {
        igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNP, nodes, p,
                                /*directed=*/ 0, /*loops=*/ 0);
        ecount = igraph_ecount(&g);
        for (e = 0; e < edges; e++) {
            long int edge = RNG_INTEGER(0, ecount - 1);
            VECTOR(pairs)[2 * e] = IGRAPH_FROM(&g, edge);
            VECTOR(pairs)[2 * e + 1] = IGRAPH_TO(&g, edge);
        }
        igraph_get_eids(&g, &eids, &pairs, /*path=*/ 0, 0, /*error=*/ 1);
        for (e = 0; e < edges; e++) {
            long int edge = VECTOR(eids)[e];
            long int from1 = VECTOR(pairs)[2 * e];
            long int to1 = VECTOR(pairs)[2 * e + 1];
            long int from2 = IGRAPH_FROM(&g, edge);
            long int to2 = IGRAPH_TO(&g, edge);
            long int min1 = from1 < to1 ? from1 : to1;
            long int max1 = from1 < to1 ? to1 : from1;
            long int min2 = from2 < to2 ? from2 : to2;
            long int max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 11;
            }
        }

        igraph_diameter(&g, /*res=*/ 0, /*from=*/ 0, /*to=*/ 0, &path,
                        IGRAPH_UNDIRECTED, /*unconn=*/ 1);
        igraph_get_eids(&g, &eids, /*pairs=*/ 0, &path, 0, /*error=*/ 1);
        for (e = 0; e < igraph_vector_size(&path) - 1; e++) {
            long int edge = VECTOR(eids)[e];
            long int from1 = VECTOR(path)[e];
            long int to1 = VECTOR(path)[e + 1];
            long int from2 = IGRAPH_FROM(&g, edge);
            long int to2 = IGRAPH_TO(&g, edge);
            long int min1 = from1 < to1 ? from1 : to1;
            long int max1 = from1 < to1 ? to1 : from1;
            long int min2 = from2 < to2 ? from2 : to2;
            long int max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 12;
            }
        }

        igraph_destroy(&g);
    }

    igraph_vector_destroy(&path);
    igraph_vector_destroy(&pairs);
    igraph_vector_destroy(&eids);

    return 0;
}

int check_multi() {

    igraph_t g;
    igraph_vector_t vec;
    igraph_vector_t eids, eids2;
    int ret;
    long int i;

    igraph_real_t q1[] = { 0, 1, 0, 1 };
    igraph_real_t q2[] = { 0, 1, 0, 1, 0, 1 };
    igraph_real_t q3[] = { 1, 0, 3, 4, 1, 0, 0, 1, 3, 4, 0, 1 };

    igraph_vector_init(&eids, 0);

    /*********************************/
    igraph_small(&g, /*n=*/ 10, /*directed=*/ 1,
                 0, 1, 0, 1, 1, 0, 1, 2, 3, 4, 3, 4, 3, 4, 3, 5, 3, 7,
                 9, 8,
                 -1);

    igraph_vector_view(&vec, q1, sizeof(q1) / sizeof(igraph_real_t));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 1, /*error=*/ 1);
    igraph_vector_sort(&eids);
    print_vector(&eids, stdout);

    igraph_vector_view(&vec, q2, sizeof(q2) / sizeof(igraph_real_t));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/ 0, /*error=*/ 1);
    igraph_vector_sort(&eids);
    print_vector(&eids, stdout);

    igraph_vector_view(&vec, q2, sizeof(q2) / sizeof(igraph_real_t));
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

    igraph_vector_view(&vec, q1, sizeof(q1) / sizeof(igraph_real_t));
    igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/1, /*error=*/ 1);
    igraph_vector_sort(&eids);
    print_vector(&eids, stdout);

    igraph_vector_view(&vec, q3, sizeof(q3) / sizeof(igraph_real_t));
    igraph_set_error_handler(igraph_error_handler_ignore);
    ret = igraph_get_eids_multi(&g, &eids, &vec, 0, /*directed=*/0, /*error=*/ 1);
    if (ret != IGRAPH_EINVAL) {
        return 2;
    }
    igraph_set_error_handler(igraph_error_handler_abort);

    igraph_destroy(&g);

    /*********************************/

    igraph_vector_destroy(&eids);

    /*********************************/
    /* Speed tests */

#define NODES 10000
    igraph_barabasi_game(&g, /*n=*/ NODES, /*power=*/ 1.0, /*m=*/ 3,
                         /*outseq=*/ 0, /*outpref=*/ 0, /*A=*/ 1,
                         /*directed=*/ 1, IGRAPH_BARABASI_BAG,
                         /*start_from=*/ 0);
    igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 0, /*edge_comb=*/ 0);

    igraph_vector_init(&eids, NODES / 2);
    igraph_random_sample(&eids, 0, igraph_ecount(&g) - 1, NODES / 2);
    igraph_vector_init(&vec, NODES);
    for (i = 0; i < NODES / 2; i++) {
        VECTOR(vec)[2 * i]   = IGRAPH_FROM(&g, VECTOR(eids)[i]);
        VECTOR(vec)[2 * i + 1] = IGRAPH_TO(&g, VECTOR(eids)[i]);
    }
    igraph_vector_init(&eids2, 0);
    igraph_get_eids_multi(&g, &eids2, &vec, 0, /*directed=*/ 1, /*error=*/ 1);
    if (!igraph_vector_all_e(&eids, &eids2)) {
        return 3;
    }

    /**/

    for (i = 0; i < NODES / 2; i++) {
        VECTOR(vec)[2 * i]   = IGRAPH_TO(&g, VECTOR(eids)[i]);
        VECTOR(vec)[2 * i + 1] = IGRAPH_FROM(&g, VECTOR(eids)[i]);
    }
    igraph_get_eids_multi(&g, &eids2, &vec, 0, /*directed=*/ 0, /*error=*/ 1);
    if (!igraph_vector_all_e(&eids, &eids2)) {
        return 4;
    }

    igraph_vector_destroy(&eids);
    igraph_vector_destroy(&eids2);
    igraph_vector_destroy(&vec);
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
