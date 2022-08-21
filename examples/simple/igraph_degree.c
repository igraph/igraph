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

void print_vector(igraph_vector_int_t *v, FILE *f) {
    igraph_integer_t i;
    for (i = 0; i < igraph_vector_int_size(v); i++) {
        fprintf(f, " %" IGRAPH_PRId "", VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

void handshaking_lemma(igraph_t *g, igraph_vector_int_t *v) {
    igraph_integer_t i, ndeg, nedges;
    /* Consistency check of the handshaking lemma. */
    /* If d is the sum of all vertex degrees, then d = 2|E|. */
    ndeg = 0;
    nedges = igraph_ecount(g);
    for (i = 0; i < igraph_vector_int_size(v); i++) {
        ndeg += VECTOR(*v)[i];
    }
    IGRAPH_ASSERT(ndeg == 2 * nedges);
}

int main() {

    igraph_t g;
    igraph_vector_int_t v;
    igraph_vector_int_t seq;
    igraph_integer_t mdeg;

    /* Create graph */
    igraph_vector_int_init(&v, 8);
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 2,2, -1);

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

    handshaking_lemma(&g, &v);

    igraph_destroy(&g);

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 2,2, -1);

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

    handshaking_lemma(&g, &v);

    /* Degree of the same vertex multiple times */

    igraph_vector_int_init(&seq, 3);
    VECTOR(seq)[0] = 2;
    VECTOR(seq)[1] = 0;
    VECTOR(seq)[2] = 2;
    igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS);
    print_vector(&v, stdout);

    igraph_destroy(&g);
    igraph_vector_int_destroy(&seq);

    /* Maximum degree */

    igraph_ring(&g, 10, 0 /*undirected*/, 0 /*undirected*/, 0/*uncircular*/);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 2) {
        return 5;
    }

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    handshaking_lemma(&g, &v);
    igraph_destroy(&g);

    igraph_full(&g, 10, 0 /*undirected*/, 0/*no loops*/);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 9) {
        return 7;
    }

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    handshaking_lemma(&g, &v);
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

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    handshaking_lemma(&g, &v);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&v);

    return 0;
}
