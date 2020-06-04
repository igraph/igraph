/* -*- mode: C -*-  */
/*
  Test suite for the Fisher-Yates shuffle.
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

#define R_INTEGER(a,b) (igraph_rng_get_integer(igraph_rng_default(), (a), (b)))
#define R_UNIF(a,b) (igraph_rng_get_unif(igraph_rng_default(), (a), (b)))

int main() {
    igraph_real_t d;
    igraph_vector_t u, v;
    int ret;
    long int i, k, n;

    /********************************
     * Example usage
     ********************************/

    /* Sequences with one element. Such sequences are trivially permuted.
     * The result of any Fisher-Yates shuffle on a sequence with one element
     * must be the original sequence itself.
     */
    n = 1;
    igraph_vector_init(&v, n);
    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */
    k = R_INTEGER(-1000, 1000);
    VECTOR(v)[0] = k;
    igraph_vector_shuffle(&v);
    if (VECTOR(v)[0] != k) {
        return 1;
    }
    d = R_UNIF(-1000.0, 1000.0);

    VECTOR(v)[0] = d;
    igraph_vector_shuffle(&v);
    if (VECTOR(v)[0] != d) {
        return 2;
    }
    igraph_vector_destroy(&v);

    /* Sequences with multiple elements. A Fisher-Yates shuffle of a sequence S
     * is a random permutation \pi(S) of S. Thus \pi(S) must have the same
     * length and elements as the original sequence S. A major difference between
     * S and its random permutation \pi(S) is that the order in which elements
     * appear in \pi(S) is probably different from how elements are ordered in S.
     * If S has length n = 1, then both \pi(S) and S are equivalent sequences in
     * that \pi(S) is merely S and no permutation has taken place. If S has
     * length n > 1, then there are n! possible permutations of S. Assume that
     * each such permutation is equally likely to appear as a result of the
     * Fisher-Yates shuffle. As n increases, the probability that S is different
     * from \pi(S) also increases. We have a probability of 1 / n! that S and
     * \pi(S) are equivalent sequences.
     */
    n = 100;
    igraph_vector_init(&u, n);
    igraph_vector_init(&v, n);

    for (i = 0; i < n; i++) {
        k = R_INTEGER(-1000, 1000);
        VECTOR(u)[i] = k;
        VECTOR(v)[i] = k;
    }

    igraph_vector_shuffle(&v);
    /* must have same length */
    if (igraph_vector_size(&v) != n) {
        return 3;
    }
    if (igraph_vector_size(&u) != igraph_vector_size(&v)) {
        return 4;
    }
    /* must have same elements */
    igraph_vector_sort(&u);
    igraph_vector_sort(&v);
    if (!igraph_vector_all_e(&u, &v)) {
        return 5;
    }
    igraph_vector_destroy(&u);
    igraph_vector_destroy(&v);

    /* empty sequence */
    igraph_vector_init(&v, 0);
    ret = igraph_vector_shuffle(&v);
    igraph_vector_destroy(&v);

    return ret == 0 ? 0 : 6;
}
