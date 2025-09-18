/*
   igraph library.
   Copyright (C) 2013-2024  The igraph development team <igraph@igraph.org>

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
    igraph_vector_t trans;
    igraph_vs_t all_vertices;
    igraph_real_t avg_trans, global_trans, tri_count;

    igraph_rng_seed(igraph_rng_default(), 42);
    BENCH_INIT();

#define N 6000
#define M 2000000

    igraph_erdos_renyi_game_gnm(&g, N, M, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_init(&trans, igraph_vcount(&g));
    igraph_vs_range(&all_vertices, 0, igraph_vcount(&g));

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

    BENCH(" 5 Count triangles per vertex, all vertices method, GNM",
          igraph_count_adjacent_triangles(&g, &trans, igraph_vss_all());
         );

    BENCH(" 6 Count triangles per vertex, subset method, GNM",
          igraph_count_adjacent_triangles(&g, &trans, all_vertices);
         );

    BENCH(" 7 Count all triangles, GNM",
          igraph_count_triangles(&g, &tri_count);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    printf("\n");

    igraph_barabasi_game(&g, N, /*power=*/ 1, M / N, /*outseq=*/ 0,
                         /*outpref=*/ 0, /*A=*/ 1, IGRAPH_UNDIRECTED,
                         IGRAPH_BARABASI_PSUMTREE, /*start_from=*/ 0);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_range(&all_vertices, 0, igraph_vcount(&g));

    BENCH(" 1 Local transitivity, all vertices method, Barabasi",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 2 Local transitivity, subset method, Barabasi",
          igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 3 Average local transitivity, Barabasi",
          igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 4 Global transitivity, Barabasi",
          igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN);
         );

    BENCH(" 5 Count triangles per vertex, all vertices method, Barabasi",
          igraph_count_adjacent_triangles(&g, &trans, igraph_vss_all());
          );

    BENCH(" 6 Count triangles per vertex, subset method, Barabasi",
          igraph_count_adjacent_triangles(&g, &trans, all_vertices);
          );

    BENCH(" 7 Count all triangles, Barabasi",
          igraph_count_triangles(&g, &tri_count);
          );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    printf("\n");

#undef N
#undef M

#define N 1000000
#define M 10000000

    igraph_erdos_renyi_game_gnm(&g, N, M, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_init(&trans, igraph_vcount(&g));
    igraph_vs_range(&all_vertices, 0, igraph_vcount(&g));

    BENCH(" 1 Local transitivity, all vertices method, large GNM",
          igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                                               IGRAPH_TRANSITIVITY_NAN);
          );

    BENCH(" 2 Local transitivity, subset method, large GNM",
          igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                                               IGRAPH_TRANSITIVITY_NAN);
          );

    BENCH(" 3 Average local transitivity, large GNM",
          igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN);
          );

    BENCH(" 4 Global transitivity, large GNM",
          igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN);
          );

    BENCH(" 5 Count triangles per vertex, all vertices method, large GNM",
          igraph_count_adjacent_triangles(&g, &trans, igraph_vss_all());
          );

    BENCH(" 6 Count triangles per vertex, subset method, large GNM",
          igraph_count_adjacent_triangles(&g, &trans, all_vertices);
          );

    BENCH(" 7 Count all triangles, large GNM",
          igraph_count_triangles(&g, &tri_count);
          );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

    printf("\n");

    igraph_erdos_renyi_game_gnm(&g, 500, 2000, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_range(&all_vertices, 0, igraph_vcount(&g));

#define REPS 1000

    BENCH(" 1 Local transitivity, all vertices method, small GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH(" 2 Local transitivity, subset method, small GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH(" 3 Average local transitivity, small GNM",
          REPEAT(igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH(" 4 Global transitivity, small GNM",
          REPEAT(igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH(" 5 Count triangles per vertex, all vertices method, small GNM",
          REPEAT(igraph_count_adjacent_triangles(&g, &trans, igraph_vss_all()), REPS);
         );

    BENCH(" 6 Count triangles per vertex, subset method, small GNM",
          REPEAT(igraph_count_adjacent_triangles(&g, &trans, all_vertices), REPS);
         );

    BENCH(" 7 Count all triangles, small GNM",
          REPEAT(igraph_count_triangles(&g, &tri_count), REPS);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

#undef REPS

    printf("\n");

    igraph_erdos_renyi_game_gnm(&g, 50, 300, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
    igraph_vector_resize(&trans, igraph_vcount(&g));
    igraph_vs_range(&all_vertices, 0, igraph_vcount(&g));

#define REPS 10000

    BENCH(" 1 Local transitivity, all vertices method, tiny GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, igraph_vss_all(),
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH(" 2 Local transitivity, subset method, tiny GNM",
          REPEAT(igraph_transitivity_local_undirected(&g, &trans, all_vertices,
                  IGRAPH_TRANSITIVITY_NAN), REPS);
         );

    BENCH(" 3 Average local transitivity, tiny GNM",
          REPEAT(igraph_transitivity_avglocal_undirected(&g, &avg_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH(" 4 Global transitivity, tiny GNM",
          REPEAT(igraph_transitivity_undirected(&g, &global_trans, IGRAPH_TRANSITIVITY_NAN),
                 REPS);
         );

    BENCH(" 5 Count triangles per vertex, all vertices method, tiny GNM",
          REPEAT(igraph_count_adjacent_triangles(&g, &trans, igraph_vss_all()), REPS);
         );

    BENCH(" 6 Count triangles per vertex, subset method, tiny GNM",
          REPEAT(igraph_count_adjacent_triangles(&g, &trans, all_vertices), REPS);
         );

    BENCH(" 7 Count all triangles, tiny GNM",
          REPEAT(igraph_count_triangles(&g, &tri_count), REPS);
         );

    igraph_vs_destroy(&all_vertices);
    igraph_destroy(&g);

#undef REPS

    igraph_vector_destroy(&trans);

    return 0;
}
