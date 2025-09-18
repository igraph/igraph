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
    igraph_vector_int_list_t res;
    igraph_int_t res_int;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_vector_int_list_init(&res, 0);

    igraph_erdos_renyi_game_gnm(&g, 100, 3000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    BENCH(" 1 Cliques in random graph with 100 vertices and 3000 edges",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 0, -1);
         );
    igraph_destroy(&g);
    igraph_vector_int_list_clear(&res);

    igraph_erdos_renyi_game_gnm(&g, 200, 10000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    BENCH(" 2 Cliques in random graph with 200 vertices and 10000 edges, up to size 5",
          igraph_cliques(&g, &res, /* min_size= */ 0, /* max_size= */ 5, -1);
         );
    igraph_vector_int_list_clear(&res);
    BENCH(" 3 Clique number of the same graph with 200 vertices and 10000 edges",
          igraph_clique_number(&g, &res_int);
         );
    igraph_vector_int_list_clear(&res);
    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&res);

    return 0;
}
