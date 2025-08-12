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
#include "bench.h"

void bench_erdos_renyi(int nodes, int edges, int rewirings, int rep) {
   igraph_t graph;

   igraph_erdos_renyi_game_gnm(&graph, nodes, edges, false, false, false);

   char name[200];
   sprintf(name, "Erdős-Rényi graph rewire %d nodes, %d edges, %d swaps, %dx",
           nodes, edges, rewirings, rep);

   BENCH(name, REPEAT(igraph_rewire(&graph, rewirings, IGRAPH_REWIRING_SIMPLE), rep));

   igraph_destroy(&graph);
}

int main(void) {
   igraph_rng_seed(igraph_rng_default(), 137);
   BENCH_INIT();

   bench_erdos_renyi(/* nodes */100, /* edges */3000, /* rewirings */30000, /* rep */500);
   return 0;
}
