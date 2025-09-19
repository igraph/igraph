/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "bench.h"

int main(void) {
    igraph_t g;
    igraph_graph_list_t res;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_graph_list_init(&res, 0);

    igraph_empty(&g, 1000, IGRAPH_UNDIRECTED);
    BENCH(" 1 Decompose graph with 1000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_empty(&g, 10000, IGRAPH_UNDIRECTED);
    BENCH(" 2 Decompose graph with 10000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_empty(&g, 100000, IGRAPH_UNDIRECTED);
    BENCH(" 3 Decompose graph with 100000 isolated vertices",
          igraph_decompose(&g, &res, IGRAPH_WEAK, -1, -1);
         );
    igraph_destroy(&g);
    igraph_graph_list_clear(&res);

    igraph_graph_list_destroy(&res);

    return 0;
}
