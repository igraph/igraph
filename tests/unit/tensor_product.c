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

// K2 X G(5,2) = G(10,3)
void test_petersen(void) {
    igraph_t k2, g_5_2, g_10_3, product;
    igraph_bool_t is_iso;

    igraph_full(&k2, 2, IGRAPH_UNDIRECTED, false);
    igraph_generalized_petersen(&g_5_2, 5, 2);
    igraph_generalized_petersen(&g_10_3, 10, 3);

    igraph_product(&product, &k2, &g_5_2, IGRAPH_PRODUCT_TENSOR);

    igraph_isomorphic(&product, &g_10_3, &is_iso);

    IGRAPH_ASSERT(is_iso);

    igraph_destroy(&k2);
    igraph_destroy(&g_5_2);
    igraph_destroy(&g_10_3);
    igraph_destroy(&product);
}

int main(void) {
    test_petersen();

    VERIFY_FINALLY_STACK();

    return 0;
}