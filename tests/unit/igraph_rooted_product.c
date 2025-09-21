/*
    igraph library.
    Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

// P2 X P4 with root 0 is P8
void test_p2_p4(void) {
    igraph_t p2, p4, p8, product;
    igraph_bool_t is_iso;

    igraph_ring(&p2, 2, IGRAPH_UNDIRECTED, false, false);
    igraph_ring(&p4, 4, IGRAPH_UNDIRECTED, false, false);
    igraph_ring(&p8, 8, IGRAPH_UNDIRECTED, false, false);

    igraph_rooted_product(&product, &p2, &p4, 0);

    igraph_isomorphic(&product, &p8, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&product);
    igraph_destroy(&p2);
    igraph_destroy(&p4);
    igraph_destroy(&p8);
}

void test_multigraph(void) {
    igraph_t g1, g2, exp_prod, product;
    igraph_bool_t is_iso;

    igraph_small(&g1, 2, IGRAPH_UNDIRECTED, 0, 1,    0, 1, -1);
    igraph_small(&g2, 2, IGRAPH_UNDIRECTED, 0, 1,    0, 1,    0, 1, -1);
    igraph_small(&exp_prod, 4, IGRAPH_UNDIRECTED, 0, 1,    0, 1,    0,1,
                                            2, 3,   2, 3,   2,3,
                                            0,2,    0, 2, -1, IGRAPH_UNDIRECTED);

    igraph_rooted_product(&product, &g1, &g2, 0);
    igraph_isomorphic(&product, &exp_prod, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&product);
    igraph_destroy(&g1);
    igraph_destroy(&g2);
    igraph_destroy(&exp_prod);
}

void test_directed(void) {
    igraph_t p3, c3, exp_prod, product;
    igraph_bool_t is_iso;

    igraph_ring(&p3, 3, IGRAPH_DIRECTED, false, false);
    igraph_ring(&c3, 3, IGRAPH_DIRECTED, false, true);

    // calculated by hand
    igraph_small(&exp_prod, 9, IGRAPH_DIRECTED, 0,1,  1,2,
                                          3,0,  0,4,  4,3,
                                          5,1,  1,6,  6,5,
                                          7,2,  2,8,  8,7, -1);

    igraph_rooted_product(&product, &p3, &c3, 0);
    igraph_isomorphic(&exp_prod, &product, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&product);
    igraph_destroy(&p3);
    igraph_destroy(&c3);
    igraph_destroy(&exp_prod);
}

void test_null_and_singleton(void) {
    igraph_t null_g, singleton_g, product;

    igraph_empty(&null_g, 0, IGRAPH_UNDIRECTED);
    igraph_empty(&singleton_g, 1, IGRAPH_UNDIRECTED);

    igraph_rooted_product(&product, &null_g, &singleton_g, 0);
    IGRAPH_ASSERT(igraph_vcount(&product) == 0);
    igraph_destroy(&product);

    igraph_rooted_product(&product, &singleton_g, &singleton_g, 0);
    IGRAPH_ASSERT(igraph_vcount(&product) == 1);
    IGRAPH_ASSERT(igraph_ecount(&product) == 0);
    igraph_destroy(&product);

    CHECK_ERROR(igraph_rooted_product(&product, &singleton_g, &null_g, 0), IGRAPH_EINVVID);

    igraph_destroy(&singleton_g);
    igraph_destroy(&null_g);
}

int main(void) {
    test_p2_p4();
    test_multigraph();
    test_directed();
    test_null_and_singleton();

    VERIFY_FINALLY_STACK();

    return 0;
}
