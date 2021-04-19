/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "bench.h"

#define N 6000
#define M 2000000

int main() {

    igraph_t g;
    igraph_vector_t trans;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, N, M,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init(&trans, igraph_vcount(&g));

    BENCH(" 1 Transitivity GNM",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN);
         );

    igraph_destroy(&g);
    igraph_barabasi_game(&g, N, /*power=*/ 1, M / N, /*outseq=*/ 0,
                         /*outpref=*/ 0, /*A=*/ 1, IGRAPH_UNDIRECTED,
                         IGRAPH_BARABASI_PSUMTREE, /*start_from=*/ 0);

    BENCH(" 2 Transitivity preferential attachment",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN);
         );

    igraph_destroy(&g);
    igraph_vector_destroy(&trans);

    return 0;
}
