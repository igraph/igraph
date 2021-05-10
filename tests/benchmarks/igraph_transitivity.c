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
    igraph_vs_t all_vertices;
    igraph_real_t avg_trans, global_trans;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, N, M,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_init(&trans, igraph_vcount(&g));
    igraph_vs_seq(&all_vertices, 0, igraph_vcount(&g) - 1);

    BENCH(" 1 Local transitivity, all vertices method, GNM",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 2 Local transitivity, subset method, GNM",
          igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 3 Average local transitivity GNM",
          igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 4 Global transitivity GNM",
          igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    igraph_barabasi_game(&g, N, /*power=*/ 1, M / N, /*outseq=*/ 0,
                         /*outpref=*/ 0, /*A=*/ 1, IGRAPH_UNDIRECTED,
                         IGRAPH_BARABASI_PSUMTREE, /*start_from=*/ 0);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_seq(&all_vertices, 0, igraph_vcount(&g) - 1);

    BENCH(" 5 Local transitivity, all vertices method, Barabasi",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 6 Local transitivity, subset method, Barabasi",
          igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 7 Average local transitivity, Barabasi",
          igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 8 Global transitivity, Barabasi",
          igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 500, 2000,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_seq(&all_vertices, 0, igraph_vcount(&g) - 1);

#define REPS 1000

    BENCH(" 9 Local transitivity, all vertices method, small GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH("10 Local transitivity, subset method, small GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH("11 Average local transitivity, small GNM",
          REPEAT(igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH("12 Global transitivity, small GNM",
          REPEAT(igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    igraph_erdos_renyi_game(&g, IGRAPH_ERDOS_RENYI_GNM, 50, 300,
                            IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_seq(&all_vertices, 0, igraph_vcount(&g) - 1);

#undef REPS

#define REPS 10000

    BENCH("13 Local transitivity, all vertices method, tiny GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH("14 Local transitivity, subset method, tiny GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH("15 Average local transitivity, tiny GNM",
          REPEAT(igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH("16 Global transitivity, tiny GNM",
          REPEAT(igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

#undef REPS

    igraph_vector_destroy(&trans);

    return 0;
}
