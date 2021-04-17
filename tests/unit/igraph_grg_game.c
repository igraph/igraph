/*
   IGraph library.
   Copyright (C) 2006-2021  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.inc"

int main() {

    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Empty graph */
    igraph_grg_game(&g, 100, 0, 0, 0, 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == 0);
    igraph_destroy(&g);

    /* Full graph */
    igraph_grg_game(&g, 10, sqrt(2.0) / 2, 1, 0, 0);
    IGRAPH_ASSERT(igraph_ecount(&g) == igraph_vcount(&g) * (igraph_vcount(&g) - 1) / 2);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
