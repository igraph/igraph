/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

#include "bench.h"

#define BENCH_BLOCK(G, N, Q, TEXT1, CODE1, TEXT2, CODE2) \
    do { \
        char msg[128]; \
        const igraph_int_t vcount = pow(Q, N); \
        const igraph_int_t ecount = (vcount * (Q - 1) * N) / 2; \
        snprintf(msg, sizeof(msg) / sizeof(msg[0]), \
             " 1 %s: vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId "", \
             TEXT1, vcount, ecount); \
        BENCH(msg, CODE1); \
        IGRAPH_ASSERT(igraph_vcount(&G) == vcount); \
        IGRAPH_ASSERT(igraph_ecount(&G) == ecount); \
        igraph_destroy(&G); \
        snprintf(msg, sizeof(msg) / sizeof(msg[0]), \
             " 2 %s: vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId "", \
             TEXT2, vcount, ecount); \
        BENCH(msg, CODE2); \
        IGRAPH_ASSERT(igraph_vcount(&G) == vcount); \
        IGRAPH_ASSERT(igraph_ecount(&G) == ecount); \
        igraph_destroy(&G); \
        printf("\n"); \
    } while (0)

/*
 * This is an alternative Hamming graph generator based on the fact that
 * the Hamming graph H(n,q) is, equivalently, the cartesian product of n
 * copies of complete graphs K(q).
 */
igraph_error_t cartesian_hamming(igraph_t *graph, igraph_int_t n, igraph_int_t q,
                              igraph_bool_t directed) {
    igraph_t full;

    if (n <= 0 || q <= 0) {
        IGRAPH_ERROR("n and q must be greater than zero.", IGRAPH_EINVAL);
    }

    igraph_full(&full, q, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_copy(graph, &full);

    for (igraph_int_t i = 1; i < n; i++) {
        igraph_t temp;
        igraph_product(&temp, graph, &full, IGRAPH_PRODUCT_CARTESIAN);
        igraph_destroy(graph);
        igraph_copy(graph, &temp);
        igraph_destroy(&temp);
    }

    igraph_destroy(&full);

    return IGRAPH_SUCCESS;
}

/*
 * This is a rook graph generator as a special case of cartesian_hamming for n = 2.
 */
igraph_error_t lattice_rook(igraph_t *graph, igraph_int_t q,
                              igraph_bool_t directed) {
    igraph_t full;

    if (q <= 0) {
        IGRAPH_ERROR("q must be greater than zero.", IGRAPH_EINVAL);
    }

    igraph_full(&full, q, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_product(graph, &full, &full, IGRAPH_PRODUCT_CARTESIAN);

    igraph_destroy(&full);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t g;
    igraph_int_t d, q;

    BENCH_INIT();

    /* Compare H_1_q to the complete(full) graph K_q. */
    d = 1, q = 10000; // d must be 1
    printf("Compare H(1,q) and K(q) for q=%"IGRAPH_PRId"\n", q);
    BENCH_BLOCK(g, d, q,
              "HAMMING   ", igraph_hamming(&g, d, q, IGRAPH_UNDIRECTED),
              "FULL      ", igraph_full(&g, q, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS)
             );

    /* Compare H_d_2 to the hypercube graph Q_d. */
    d = 22, q = 2; // q must be 2
    printf("Compare H(n,2) and Q(n) for n=%"IGRAPH_PRId"\n", d);
    BENCH_BLOCK(g, d, q,
              "HAMMING   ", igraph_hamming(&g, d, q, IGRAPH_UNDIRECTED),
              "HYPERCUBE ", igraph_hypercube(&g, d, IGRAPH_UNDIRECTED)
             );

    /* Compare H_2_q to the lattice(rook) graph L_q_q. */
    d = 2, q = 300; // n must be 2
    printf("Compare H(2,q) and L(q,q) for q=%"IGRAPH_PRId"\n", q);
    BENCH_BLOCK(g, d, q,
              "HAMMING  ", igraph_hamming(&g, d, q, IGRAPH_UNDIRECTED),
              "ROOK     ", lattice_rook(&g, q, IGRAPH_UNDIRECTED)
             );

    /* Compare H_n_q to d cartesian products of K_q. */
    d = 9, q = 5;
    printf("Compare H(n,q) and K(q)^n for n=%"IGRAPH_PRId" and q=%"IGRAPH_PRId"\n", d, q);
    BENCH_BLOCK(g, d, q,
              "HAMMING   ", igraph_hamming(&g, d, q, IGRAPH_UNDIRECTED),
              "CARTESIAN ", cartesian_hamming(&g, d, q, IGRAPH_UNDIRECTED)
             );

    return 0;
}
