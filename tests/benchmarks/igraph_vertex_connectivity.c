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
    igraph_int_t vconn;

    igraph_rng_seed(igraph_rng_default(), 54);
    BENCH_INIT();

    igraph_grg_game(&g, 100, 0.2, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);

    BENCH(" 1 Vertex connectivity of geometric random graph, n=100, r=0.2",
          igraph_vertex_connectivity(&g, &vconn, /* checks = */ false);
         );

    igraph_destroy(&g);

    igraph_grg_game(&g, 200, 0.141, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);

    BENCH(" 2 Vertex connectivity of geometric random graph, n=200, r=0.141",
          igraph_vertex_connectivity(&g, &vconn, /* checks = */ false);
         );

    igraph_destroy(&g);

    igraph_grg_game(&g, 400, 0.1, /* torus = */ false, /* x = */ NULL, /* y = */ NULL);

    BENCH(" 3 Vertex connectivity of geometric random graph, n=400, r=0.1",
          igraph_vertex_connectivity(&g, &vconn, /* checks = */ false);
         );

    igraph_destroy(&g);

    return 0;
}
