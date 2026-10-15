/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/* Self-contained random DAG generator, mirroring the one used in
 * tests/unit/igraph_minimum_path_cover.c (kept separate on purpose, per
 * the benchmark harness convention of standalone benchmark mains). */
static igraph_error_t random_dag(igraph_t *g, igraph_int_t n, igraph_int_t m) {
    igraph_vector_int_t perm, edges;
    igraph_vector_bool_t used;
    igraph_int_t i, added, attempts;
    igraph_int_t max_edges = n > 1 ? n * (n - 1) / 2 : 0;
    igraph_int_t max_attempts;

    if (m > max_edges) {
        m = max_edges;
    }
    max_attempts = (max_edges + 10) * 50;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&perm, n);
    for (i = 0; i < n; i++) {
        VECTOR(perm)[i] = i;
    }
    for (i = n - 1; i > 0; i--) {
        igraph_int_t j = RNG_INTEGER(0, i);
        igraph_int_t tmp = VECTOR(perm)[i];
        VECTOR(perm)[i] = VECTOR(perm)[j];
        VECTOR(perm)[j] = tmp;
    }

    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&used, n > 0 ? n * n : 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    added = 0;
    for (attempts = 0; added < m && attempts < max_attempts; attempts++) {
        igraph_int_t a = RNG_INTEGER(0, n - 1);
        igraph_int_t b = RNG_INTEGER(0, n - 1);
        if (VECTOR(perm)[a] >= VECTOR(perm)[b] || VECTOR(used)[a * n + b]) {
            continue;
        }
        VECTOR(used)[a * n + b] = true;
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, a));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, b));
        added++;
    }

    IGRAPH_CHECK(igraph_create(g, &edges, n, IGRAPH_DIRECTED));

    igraph_vector_int_destroy(&edges);
    igraph_vector_bool_destroy(&used);
    igraph_vector_int_destroy(&perm);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_t g;
    igraph_vector_int_list_t cover;
    igraph_int_t width;
    igraph_int_t sizes[] = {100, 500, 1000};
    int i;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_int_list_init(&cover, 0);

    random_dag(&g, 300, 1200);
    BENCH(" 1 Naive reduction + naive DFS solve, random DAG n=300 m=1200",
          igraph_minimum_path_cover(&g, &cover, &width, NULL,
              IGRAPH_MPC_REDUCTION_NAIVE, IGRAPH_MPC_SOLVER_NAIVE_DFS,
              NULL, NULL, NULL);
         );
    BENCH(" 2 Greedy reduction + naive DFS solve, random DAG n=300 m=1200",
          igraph_minimum_path_cover(&g, &cover, &width, NULL,
              IGRAPH_MPC_REDUCTION_GREEDY, IGRAPH_MPC_SOLVER_NAIVE_DFS,
              NULL, NULL, NULL);
         );
    BENCH(" 3 Greedy reduction + maxflow-reduction solve, random DAG n=300 m=1200",
          igraph_minimum_path_cover(&g, &cover, &width, NULL,
              IGRAPH_MPC_REDUCTION_GREEDY, IGRAPH_MPC_SOLVER_MAXFLOW_REDUCTION,
              NULL, NULL, NULL);
         );
    igraph_destroy(&g);

    for (i = 0; i < 3; i++) {
        igraph_int_t n = sizes[i];
        char label[128];

        random_dag(&g, n, n * 4);
        snprintf(label, sizeof(label),
                 " %d Greedy reduction + naive DFS solve, random DAG n=%" IGRAPH_PRId " m=%" IGRAPH_PRId,
                 4 + i, n, n * 4);
        BENCH(label,
              igraph_minimum_path_cover(&g, &cover, &width, NULL,
                  IGRAPH_MPC_REDUCTION_GREEDY, IGRAPH_MPC_SOLVER_NAIVE_DFS,
                  NULL, NULL, NULL);
             );
        igraph_destroy(&g);
    }

    igraph_vector_int_list_destroy(&cover);

    return 0;
}
