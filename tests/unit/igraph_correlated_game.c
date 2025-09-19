/*
   igraph library.
   Copyright (C) 2003-2024  The igraph development team <igraph@igraph.org>

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
    igraph_t g1, g2;

    igraph_rng_seed(igraph_rng_default(), 9275);

    igraph_erdos_renyi_game_gnp(&g1, 10, 0.3, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);

    igraph_correlated_game(&g2, &g1, 0.9, 0.3, /* permutation=*/ NULL);
    IGRAPH_ASSERT(igraph_vcount(&g1) == igraph_vcount(&g2));
    igraph_destroy(&g2);

    igraph_correlated_game(&g2, &g1, 0.0, 0.3, /* permutation=*/ NULL);
    IGRAPH_ASSERT(igraph_vcount(&g1) == igraph_vcount(&g2));
    igraph_destroy(&g2);

    igraph_correlated_game(&g2, &g1, 1.0, 0.3, /* permutation=*/ NULL);
    IGRAPH_ASSERT(igraph_vcount(&g1) == igraph_vcount(&g2));
    igraph_destroy(&g2);

    igraph_destroy(&g1);

    VERIFY_FINALLY_STACK();

    return 0;
}
