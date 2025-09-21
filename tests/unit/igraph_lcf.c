/*
   igraph library.
   Copyright (C) 2007-2025  The igraph development team <igraph@igraph.org>

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

int main(void) {

    igraph_t g, g2;
    igraph_bool_t iso;

    // Franklin graph: [5, -5]^6x
    igraph_lcf_small(&g, 12, 5, -5, 6, 0);
    igraph_famous(&g2, "franklin");

    igraph_isomorphic(&g, &g2, &iso);
    IGRAPH_ASSERT(iso);

    igraph_destroy(&g);
    igraph_destroy(&g2);

    // [3, -2]^4, n=8
    igraph_lcf_small(&g, 8, 3, -2, 4, 0);

    IGRAPH_ASSERT(igraph_ecount(&g) == 16);

    igraph_destroy(&g);

    // [2, -2]^2, n=2
    igraph_lcf_small(&g, 2, 2, -2, 2, 0);

    IGRAPH_ASSERT(igraph_ecount(&g) == 1);

    igraph_destroy(&g);

    // [2]^2, n=2
    igraph_lcf_small(&g, 2, 2, 2, 0);

    IGRAPH_ASSERT(igraph_ecount(&g) == 1);

    igraph_destroy(&g);

    // Regression test for bug #996
    igraph_lcf_small(&g, 0, 0);

    IGRAPH_ASSERT(igraph_ecount(&g) == 0 && igraph_vcount(&g) == 0);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
