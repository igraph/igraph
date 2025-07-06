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

int main(void) {
    test_p2_p4();

    VERIFY_FINALLY_STACK();

    return 0;
}
