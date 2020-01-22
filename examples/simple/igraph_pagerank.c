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
#include <unistd.h>
#include <libgen.h>

void warning_handler_stdout (const char *reason, const char *file,
                             int line, int igraph_errno) {
    IGRAPH_UNUSED(igraph_errno);
    printf("Warning: %s\n", reason);
}

void print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %4.2f", VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

igraph_warning_handler_t *oldwarn;

int main() {

    igraph_t g;
    igraph_vector_t v, res, reset, weights;
    igraph_arpack_options_t arpack_options;
    igraph_real_t value;
    int ret;
    igraph_pagerank_power_options_t power_options;

    /* Test graphs taken from http://www.iprcom.com/papers/pagerank/ */
    igraph_vector_init(&v, 10);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 2;
    VECTOR(v)[4] = 2;
    VECTOR(v)[5] = 0;
    VECTOR(v)[6] = 3;
    VECTOR(v)[7] = 2;
    VECTOR(v)[8] = 0;
    VECTOR(v)[9] = 2;
    igraph_create(&g, &v, 0, 1);

    igraph_vector_init(&res, 0);
    oldwarn = igraph_set_warning_handler(warning_handler_stdout);
    igraph_pagerank_old(&g, &res, igraph_vss_all(), 1, 1000, 0.001, 0.85, 0);
    print_vector(&res, stdout);
    igraph_vector_destroy(&res);
    igraph_vector_destroy(&v);

    igraph_destroy(&g);

    igraph_vector_init(&v, 28);
    VECTOR(v)[ 0] = 0;
    VECTOR(v)[ 1] = 1;
    VECTOR(v)[ 2] = 0;
    VECTOR(v)[ 3] = 2;
    VECTOR(v)[ 4] = 0;
    VECTOR(v)[ 5] = 3;
    VECTOR(v)[ 6] = 1;
    VECTOR(v)[ 7] = 0;
    VECTOR(v)[ 8] = 2;
    VECTOR(v)[ 9] = 0;
    VECTOR(v)[10] = 3;
    VECTOR(v)[11] = 0;
    VECTOR(v)[12] = 3;
    VECTOR(v)[13] = 4;
    VECTOR(v)[14] = 3;
    VECTOR(v)[15] = 5;
    VECTOR(v)[16] = 3;
    VECTOR(v)[17] = 6;
    VECTOR(v)[18] = 3;
    VECTOR(v)[19] = 7;
    VECTOR(v)[20] = 4;
    VECTOR(v)[21] = 0;
    VECTOR(v)[22] = 5;
    VECTOR(v)[23] = 0;
    VECTOR(v)[24] = 6;
    VECTOR(v)[25] = 0;
    VECTOR(v)[26] = 7;
    VECTOR(v)[27] = 0;
    igraph_create(&g, &v, 0, 1);

    igraph_vector_init(&res, 0);
    igraph_pagerank_old(&g, &res, igraph_vss_all(), 1, 10000, 0.0001, 0.85, 0);
    print_vector(&res, stdout);
    igraph_vector_destroy(&res);
    igraph_vector_destroy(&v);
    igraph_destroy(&g);

    igraph_set_warning_handler(oldwarn);

    /* New PageRank */
    igraph_star(&g, 11, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_vector_init(&res, 0);
    igraph_arpack_options_init(&arpack_options);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    print_vector(&res, stdout);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    print_vector(&res, stdout);
    /* Check twice more for consistency */
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    print_vector(&res, stdout);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    print_vector(&res, stdout);

    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, &arpack_options);
    print_vector(&res, stdout);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                    igraph_vss_all(), 0, 0.85, 0, 0);
    print_vector(&res, stdout);

    /* Check personalized PageRank */
    igraph_personalized_pagerank_vs(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                                    igraph_vss_all(), 0, 0.5,
                                    igraph_vss_1(1), 0, &arpack_options);
    print_vector(&res, stdout);
    igraph_personalized_pagerank_vs(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                                    igraph_vss_all(), 0, 0.5,
                                    igraph_vss_1(1), 0, 0);
    print_vector(&res, stdout);

    /* Errors */
    power_options.niter = -1;
    power_options.eps = 0.0001;
    igraph_set_error_handler(igraph_error_handler_ignore);
    igraph_set_warning_handler(igraph_warning_handler_ignore);
    ret = igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_POWER, &res,
                          /*value=*/ 0, igraph_vss_all(), 1, 0.85,
                          /*weights=*/ 0, &power_options);
    if (ret != IGRAPH_EINVAL) {
        return 1;
    }

    power_options.niter = 10000;
    power_options.eps = -1;
    ret = igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_POWER, &res,
                          /*value=*/ 0, igraph_vss_all(), 1, 0.85,
                          /*weights=*/ 0, &power_options);
    if (ret != IGRAPH_EINVAL) {
        return 2;
    }

    power_options.niter = 10000;
    power_options.eps = 0.0001;
    ret = igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_POWER, &res,
                          /*value=*/ 0, igraph_vss_all(), 1, 1.2,
                          /*weights=*/ 0, &power_options);
    if (ret != IGRAPH_EINVAL) {
        return 3;
    }

    igraph_vector_init(&reset, 2);
    ret = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                                       igraph_vss_all(), 0, 0.85, &reset, 0,
                                       &arpack_options);
    if (ret != IGRAPH_EINVAL) {
        return 4;
    }
    ret = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                                       igraph_vss_all(), 0, 0.85, &reset, 0, 0);
    if (ret != IGRAPH_EINVAL) {
        return 4;
    }
    igraph_vector_resize(&reset, 10);
    igraph_vector_fill(&reset, 0);
    ret = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK,
                                       &res, 0, igraph_vss_all(), 0, 0.85,
                                       &reset, 0, &arpack_options);
    if (ret != IGRAPH_EINVAL) {
        return 5;
    }
    ret = igraph_personalized_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK,
                                       &res, 0, igraph_vss_all(), 0, 0.85,
                                       &reset, 0, 0);
    if (ret != IGRAPH_EINVAL) {
        return 5;
    }
    igraph_vector_destroy(&reset);
    igraph_destroy(&g);
    igraph_set_error_handler(igraph_error_handler_abort);

    /* Special cases: check for empty graph */
    igraph_empty(&g, 10, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, 0, &arpack_options);
    if (value != 1.0) {
        return 6;
    }
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, 0, 0);
    if (value != 1.0) {
        return 6;
    }
    print_vector(&res, stdout);
    igraph_destroy(&g);

    /* Special cases: check for full graph, zero weights */
    igraph_full(&g, 10, 0, 0);
    igraph_vector_init(&v, 45);
    igraph_vector_fill(&v, 0);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &v, &arpack_options);
    if (value != 1.0) {
        return 7;
    }
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, &value,
                    igraph_vss_all(), 1, 0.85, &v, 0);
    if (value != 1.0) {
        return 7;
    }
    igraph_vector_destroy(&v);
    print_vector(&res, stdout);
    igraph_destroy(&g);

    /* Another test case for PageRank (bug #792352) */
    igraph_small(&g, 9, 1, 0, 5, 1, 5, 2, 0, 3, 1, 5, 4, 5, 7, 6, 0, 8, 0, 8, 1, -1);
    igraph_vector_init(&weights, 9);
    VECTOR(weights)[0] = 4;
    VECTOR(weights)[1] = 5;
    VECTOR(weights)[2] = 5;
    VECTOR(weights)[3] = 4;
    VECTOR(weights)[4] = 4;
    VECTOR(weights)[5] = 4;
    VECTOR(weights)[6] = 3;
    VECTOR(weights)[7] = 4;
    VECTOR(weights)[8] = 4;
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_ARPACK, &res, 0,
                    igraph_vss_all(), 1, 0.85, &weights, &arpack_options);
    print_vector(&res, stdout);
    igraph_pagerank(&g, IGRAPH_PAGERANK_ALGO_PRPACK, &res, 0,
                    igraph_vss_all(), 1, 0.85, &weights, 0);
    print_vector(&res, stdout);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g);

    igraph_vector_destroy(&res);
    return 0;
}
