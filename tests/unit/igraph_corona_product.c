/*
    IGraph library.
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

// P1 X P2 = C3
void test_p1_p2(void) {
    igraph_t p1, p2, c3, res;
    igraph_bool_t is_iso;

    igraph_ring(&p1, 1, IGRAPH_UNDIRECTED, false, false);
    igraph_ring(&p2, 2, IGRAPH_UNDIRECTED, false, false);
    igraph_ring(&c3, 3, IGRAPH_UNDIRECTED, false, true);

    igraph_corona_product(&res, &p1, &p2, IGRAPH_OUT);

    igraph_isomorphic(&c3, &res, &is_iso);
    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&c3);
    igraph_destroy(&res);
}

void test_c6_p4(void) {
    igraph_t k6, c4, res;

    igraph_full(&k6, 6, IGRAPH_UNDIRECTED, false);
    igraph_ring(&c4, 4, IGRAPH_UNDIRECTED, false, true);

    // 0, 5, 15, 20, 25 are the vertices of the centre graph
    // 0 is connected to 1, 2, 3, 4; 5 is connected to 6, 7, 8, 9; and so on
    // 1, 2, 3, 4 are the 1st copy of c4, 6, 7, 8, 9 are the 2nd copy of c4, and so on
    igraph_corona_product(&res, &k6, &c4, IGRAPH_IN);

    print_graph(&res);

    igraph_destroy(&k6);
    igraph_destroy(&c4);
    igraph_destroy(&res);
}

int main(void) {
    test_p1_p2();
    test_c6_p4();

    VERIFY_FINALLY_STACK();
    return 0;
}
