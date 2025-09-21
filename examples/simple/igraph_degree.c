/*
   igraph library.
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

igraph_bool_t handshaking_lemma(igraph_t *g, igraph_vector_int_t *v) {
    /* Consistency check of the handshaking lemma:
     * If d is the sum of all vertex degrees, then d = 2|E|. */
    return igraph_vector_int_sum(v) == 2 * igraph_ecount(g);
}

int main(void) {
    igraph_t g;
    igraph_vector_int_t v;
    igraph_vector_int_t seq;
    igraph_int_t mdeg;

    /* Initialize the library. */
    igraph_setup();

    /* Create graph */
    igraph_vector_int_init(&v, 8);
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 2,2, -1);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    if (!handshaking_lemma(&g, &v)) {
        exit(3);
    }

    igraph_destroy(&g);

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 2,2, -1);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_int_print(&v);

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    if (!handshaking_lemma(&g, &v)) {
        exit(4);
    }

    /* Degree of the same vertex multiple times */

    igraph_vector_int_init(&seq, 3);
    VECTOR(seq)[0] = 2;
    VECTOR(seq)[1] = 0;
    VECTOR(seq)[2] = 2;
    igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_vector_int_print(&v);

    igraph_destroy(&g);
    igraph_vector_int_destroy(&seq);

    /* Maximum degree */

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, /* mutual */ false, /* circular */ false);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 2) {
        exit(5);
    }

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (! handshaking_lemma(&g, &v)) {
        exit(6);
    }
    igraph_destroy(&g);

    igraph_full(&g, 10, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 9) {
        exit(7);
    }

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (! handshaking_lemma(&g, &v)) {
        exit(8);
    }
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, /* center */ 0);
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS);
    if (mdeg != 9) {
        exit(9);
    }
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS);
    if (mdeg != 1) {
        exit(10);
    }
    igraph_maxdegree(&g, &mdeg, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (mdeg != 9) {
        exit(11);
    }

    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    if (! handshaking_lemma(&g, &v)) {
        exit(12);
    }
    igraph_destroy(&g);

    igraph_vector_int_destroy(&v);

    return 0;
}
