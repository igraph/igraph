/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

void print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i = 0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

int main() {

    igraph_t g;
    igraph_vector_t v, seq;
    int ret;
    igraph_integer_t mdeg, nedges;
    long int i;
    long int ndeg;

    /* Create graph */
    igraph_vector_init(&v, 8);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 2;
    VECTOR(v)[4] = 2;
    VECTOR(v)[5] = 3;
    VECTOR(v)[6] = 2;
    VECTOR(v)[7] = 2;
    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_set_error_handler(igraph_error_handler_ignore);

    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    ndeg = 0;
    nedges = igraph_ecount(&g);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        ndeg += (long int) VECTOR(v)[i];
    }
    if (ndeg != 2 * nedges) {
        return 1;
    }

    igraph_destroy(&g);

    igraph_vector_resize(&v, 8);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 2;
    VECTOR(v)[4] = 2;
    VECTOR(v)[5] = 3;
    VECTOR(v)[6] = 2;
    VECTOR(v)[7] = 2;
    igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_vector(&v, stdout);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    ndeg = 0;
    nedges = igraph_ecount(&g);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        ndeg += (long int) VECTOR(v)[i];
    }
    if (ndeg != 2 * nedges) {
        return 2;
    }

    /* Degree of the same vertex multiple times */

    igraph_vector_init(&seq, 3);
    VECTOR(seq)[0] = 2;
    VECTOR(seq)[1] = 0;
    VECTOR(seq)[2] = 2;
    igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    /* Errors */
    ret = igraph_degree(&g, &v, igraph_vss_vector(&seq), (igraph_neimode_t)0,
                        IGRAPH_LOOPS);
    if (ret != IGRAPH_EINVMODE) {
        return 3;
    }

    VECTOR(seq)[0] = 4;
    ret = igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS);
    if (ret != IGRAPH_EINVVID) {
        return 4;
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&seq);

    /* Maximum degree */

    igraph_ring(&g, 10, 0 /*undirected*/, 0 /*undirected*/, 0/*uncircular*/);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 2) {
        return 5;
    }
    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    ndeg = 0;
    nedges = igraph_ecount(&g);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        ndeg += (long int) VECTOR(v)[i];
    }
    if (ndeg != 2 * nedges) {
        return 6;
    }
    igraph_destroy(&g);

    igraph_full(&g, 10, 0 /*undirected*/, 0/*no loops*/);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 9) {
        return 7;
    }
    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    ndeg = 0;
    nedges = igraph_ecount(&g);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        ndeg += (long int) VECTOR(v)[i];
    }
    if (ndeg != 2 * nedges) {
        return 8;
    }
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    if (mdeg != 9) {
        return 9;
    }
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    if (mdeg != 1) {
        return 10;
    }
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 9) {
        return 11;
    }
    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    ndeg = 0;
    nedges = igraph_ecount(&g);
    for (i = 0; i < igraph_vector_size(&v); i++) {
        ndeg += (long int) VECTOR(v)[i];
    }
    if (ndeg != 2 * nedges) {
        return 12;
    }
    igraph_destroy(&g);

    igraph_vector_destroy(&v);

    return 0;
}
