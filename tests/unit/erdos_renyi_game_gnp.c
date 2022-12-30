/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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
    igraph_rng_seed(igraph_rng_default(), 137);

    for (igraph_integer_t size=2; size < 5; size++) {
        for (igraph_integer_t i=0; i < 100; i++) {
            igraph_t g;
            igraph_bool_t simple;

            igraph_erdos_renyi_game_gnp(&g, size, 0.5, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);

            igraph_is_simple(&g, &simple);
            if (! simple) {
                printf("Erdos-Renyi GNP graph is not simple! size=%" IGRAPH_PRId ", i=%" IGRAPH_PRId ".\n",
                       size, i);
                print_graph(&g);
            }
            IGRAPH_ASSERT(simple);

            igraph_destroy(&g);
        }
    }

    for (igraph_integer_t size=2; size < 5; size++) {
        for (igraph_integer_t i=0; i < 100; i++) {
            igraph_t g;
            igraph_bool_t simple;

            igraph_erdos_renyi_game_gnp(&g, size, 0.5, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

            igraph_is_simple(&g, &simple);
            if (! simple) {
                printf("Erdos-Renyi GNP graph is not simple! size=%" IGRAPH_PRId ", i=%" IGRAPH_PRId ".\n",
                       size, i);
                print_graph(&g);
            }
            IGRAPH_ASSERT(simple);

            igraph_destroy(&g);
        }
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
