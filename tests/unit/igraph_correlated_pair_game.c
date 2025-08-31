/*
   igraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include <unistd.h>
#include "test_utilities.h"

int main(void) {
    igraph_t g1, g2;
    igraph_bool_t same;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_correlated_pair_game(&g1, &g2, 0, /*corr*/ 1, /*p*/ 0.99,
           /*directed*/ 1, /*permutation*/ NULL);
    print_graph_canon(&g1);
    print_graph_canon(&g2);
    igraph_destroy(&g1);
    igraph_destroy(&g2);

    igraph_correlated_pair_game(&g1, &g2, 10, /*corr*/ 1, /*p*/ 0.5,
           /*directed*/ 1, /*permutation*/ NULL);
    igraph_is_same_graph(&g1, &g2, &same);
    IGRAPH_ASSERT(same);

    igraph_destroy(&g1);
    igraph_destroy(&g2);

    VERIFY_FINALLY_STACK();
    return 0;
}
