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

#include "test_utilities.h"

/* -------------------------------------------------------------------- */
/* Small DAG generators, ported from PerformanceMPC's graph.cpp          */
/* -------------------------------------------------------------------- */

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

static igraph_error_t binary_tree(igraph_t *g, igraph_int_t depth, igraph_bool_t reverse) {
    igraph_int_t n = (((igraph_int_t) 1) << depth) - 1;
    igraph_vector_int_t edges;
    igraph_int_t i;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    for (i = 0; i < n; i++) {
        igraph_int_t child;
        for (child = 2 * i + 1; child <= 2 * i + 2; child++) {
            if (child >= n) {
                continue;
            }
            if (reverse) {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, child));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            } else {
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
                IGRAPH_CHECK(igraph_vector_int_push_back(&edges, child));
            }
        }
    }
    IGRAPH_CHECK(igraph_create(g, &edges, n, IGRAPH_DIRECTED));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t random_x_partite(igraph_t *g, igraph_int_t x, igraph_int_t n, igraph_int_t m) {
    igraph_int_t total = x * n;
    igraph_int_t max_edges = x > 1 ? (x - 1) * n * n : 0;
    igraph_vector_int_t edges;
    igraph_vector_bool_t used;
    igraph_int_t added, attempts, max_attempts;

    if (m > max_edges) {
        m = max_edges;
    }
    max_attempts = (max_edges + 10) * 50;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&used, max_edges > 0 ? max_edges : 1);

    added = 0;
    for (attempts = 0; added < m && attempts < max_attempts; attempts++) {
        igraph_int_t layer = RNG_INTEGER(0, x - 2);
        igraph_int_t a_local = RNG_INTEGER(0, n - 1);
        igraph_int_t b_local = RNG_INTEGER(0, n - 1);
        igraph_int_t idx = (layer * n + a_local) * n + b_local;
        if (VECTOR(used)[idx]) {
            continue;
        }
        VECTOR(used)[idx] = true;
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, layer * n + a_local));
        IGRAPH_CHECK(igraph_vector_int_push_back(&edges, (layer + 1) * n + b_local));
        added++;
    }

    IGRAPH_CHECK(igraph_create(g, &edges, total, IGRAPH_DIRECTED));

    igraph_vector_bool_destroy(&used);
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Validity checkers, ported from PerformanceMPC's naive.cpp/antichain.cpp */
/* -------------------------------------------------------------------- */

static igraph_error_t check_valid_cover(
        const igraph_t *g, const igraph_vector_int_list_t *cover, igraph_bool_t *valid) {

    igraph_int_t n = igraph_vcount(g);
    igraph_int_t no_of_paths = igraph_vector_int_list_size(cover);
    igraph_vector_bool_t seen;
    igraph_int_t i, j;

    *valid = true;
    IGRAPH_VECTOR_BOOL_INIT_FINALLY(&seen, n);

    for (i = 0; i < no_of_paths && *valid; i++) {
        igraph_vector_int_t *path = igraph_vector_int_list_get_ptr(cover, i);
        igraph_int_t plen = igraph_vector_int_size(path);

        if (plen == 0) {
            *valid = false;
            break;
        }
        for (j = 0; j < plen; j++) {
            igraph_int_t v = VECTOR(*path)[j];
            if (v < 0 || v >= n) {
                *valid = false;
                break;
            }
            /* Paths in this covering-style formulation may legitimately
             * share vertices at merge points (this is a path COVER, not a
             * vertex-disjoint partition), so a repeated visit is fine. */
            VECTOR(seen)[v] = true;
            if (j > 0) {
                igraph_bool_t adjacent;
                IGRAPH_CHECK(igraph_are_adjacent(g, VECTOR(*path)[j - 1], v, &adjacent));
                if (!adjacent) {
                    *valid = false;
                    break;
                }
            }
        }
    }

    if (*valid) {
        for (i = 0; i < n; i++) {
            if (!VECTOR(seen)[i]) {
                *valid = false;
                break;
            }
        }
    }

    igraph_vector_bool_destroy(&seen);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static igraph_error_t check_is_antichain(
        const igraph_t *g, const igraph_vector_int_t *antichain, igraph_bool_t *valid) {

    igraph_int_t k = igraph_vector_int_size(antichain);
    igraph_int_t i, j;

    *valid = true;

    for (i = 0; i < k && *valid; i++) {
        igraph_vector_int_t reach;
        igraph_int_t a = VECTOR(*antichain)[i];

        IGRAPH_VECTOR_INT_INIT_FINALLY(&reach, 0);
        IGRAPH_CHECK(igraph_subcomponent(g, &reach, a, IGRAPH_OUT));
        for (j = 0; j < k; j++) {
            if (i == j) {
                continue;
            }
            if (igraph_vector_int_contains(&reach, VECTOR(*antichain)[j])) {
                *valid = false;
                break;
            }
        }
        igraph_vector_int_destroy(&reach);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Central cross-validation routine, mirrors PerformanceMPC's test_all() */
/* -------------------------------------------------------------------- */

static const igraph_mpc_reduction_t reductions[] = {
    IGRAPH_MPC_REDUCTION_NAIVE, IGRAPH_MPC_REDUCTION_GREEDY
};
static const igraph_mpc_solver_t solvers[] = {
    IGRAPH_MPC_SOLVER_NAIVE_DFS, IGRAPH_MPC_SOLVER_MAXFLOW_REDUCTION
};

static igraph_error_t test_all(const igraph_t *g, igraph_int_t expected_width) {
    igraph_int_t widths[4];
    igraph_int_t combo = 0;
    igraph_int_t ri, si;

    for (ri = 0; ri < 2; ri++) {
        for (si = 0; si < 2; si++, combo++) {
            igraph_vector_int_list_t cover, chain_cover;
            igraph_t flow_network;
            igraph_vector_int_t flow_network_flow, flow_network_demand;
            igraph_vector_int_t antichain;
            igraph_bool_t valid;

            IGRAPH_CHECK(igraph_vector_int_list_init(&cover, 0));
            IGRAPH_FINALLY(igraph_vector_int_list_destroy, &cover);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&flow_network_flow, 0);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&flow_network_demand, 0);

            IGRAPH_CHECK(igraph_minimum_path_cover(
                    g, &cover, &widths[combo], NULL, reductions[ri], solvers[si],
                    &flow_network, &flow_network_flow, &flow_network_demand));
            IGRAPH_FINALLY(igraph_destroy, &flow_network);

            if (expected_width >= 0) {
                IGRAPH_ASSERT(widths[combo] == expected_width);
            }
            if (combo > 0) {
                IGRAPH_ASSERT(widths[combo] == widths[0]);
            }

            IGRAPH_CHECK(check_valid_cover(g, &cover, &valid));
            IGRAPH_ASSERT(valid);

            IGRAPH_VECTOR_INT_INIT_FINALLY(&antichain, 0);
            IGRAPH_CHECK(igraph_maximum_antichain(
                    &flow_network, &flow_network_flow, &flow_network_demand, &antichain));
            IGRAPH_ASSERT(igraph_vector_int_size(&antichain) == widths[combo]);
            IGRAPH_CHECK(check_is_antichain(g, &antichain, &valid));
            IGRAPH_ASSERT(valid);
            igraph_vector_int_destroy(&antichain);
            IGRAPH_FINALLY_CLEAN(1);

            IGRAPH_CHECK(igraph_vector_int_list_init(&chain_cover, 0));
            IGRAPH_FINALLY(igraph_vector_int_list_destroy, &chain_cover);
            IGRAPH_CHECK(igraph_minimum_chain_cover(g, &cover, &chain_cover));
            IGRAPH_ASSERT(igraph_vector_int_list_size(&chain_cover) == widths[combo]);
            igraph_vector_int_list_destroy(&chain_cover);
            IGRAPH_FINALLY_CLEAN(1);

            igraph_destroy(&flow_network);
            igraph_vector_int_destroy(&flow_network_demand);
            igraph_vector_int_destroy(&flow_network_flow);
            igraph_vector_int_list_destroy(&cover);
            IGRAPH_FINALLY_CLEAN(4);
        }
    }

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Fixtures, transcribed 0-indexed from PerformanceMPC's test/data/test{1..4} */
/* -------------------------------------------------------------------- */

static igraph_error_t test_fixtures(void) {
    igraph_t g;

    igraph_small(&g, 2, IGRAPH_DIRECTED, 0, 1, -1);
    IGRAPH_CHECK(test_all(&g, 1));
    igraph_destroy(&g);

    igraph_small(&g, 1, IGRAPH_DIRECTED, -1);
    IGRAPH_CHECK(test_all(&g, 1));
    igraph_destroy(&g);

    igraph_small(&g, 2, IGRAPH_DIRECTED, -1);
    IGRAPH_CHECK(test_all(&g, 2));
    igraph_destroy(&g);

    igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,3, 2,3, 3,4, 3,5, -1);
    IGRAPH_CHECK(test_all(&g, 2));
    igraph_destroy(&g);

    return IGRAPH_SUCCESS;
}

static igraph_error_t test_random_dags(void) {
    igraph_int_t sizes[] = {5, 10, 25, 50};
    igraph_int_t si, seed;

    for (si = 0; si < 4; si++) {
        igraph_int_t n = sizes[si];
        igraph_int_t factors[] = {1, 2, 4};
        igraph_int_t fi;
        for (fi = 0; fi < 3; fi++) {
            for (seed = 0; seed < 2; seed++) {
                igraph_t g;
                igraph_rng_seed(igraph_rng_default(), 1000 + seed);
                IGRAPH_CHECK(random_dag(&g, n, n * factors[fi]));
                IGRAPH_CHECK(test_all(&g, -1));
                igraph_destroy(&g);
            }
        }
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t test_binary_trees(void) {
    igraph_int_t depth;

    for (depth = 2; depth <= 5; depth++) {
        igraph_t g;

        IGRAPH_CHECK(binary_tree(&g, depth, false));
        IGRAPH_CHECK(test_all(&g, -1));
        igraph_destroy(&g);

        IGRAPH_CHECK(binary_tree(&g, depth, true));
        IGRAPH_CHECK(test_all(&g, -1));
        igraph_destroy(&g);
    }

    return IGRAPH_SUCCESS;
}

static igraph_error_t test_x_partite(void) {
    igraph_int_t x;

    igraph_rng_seed(igraph_rng_default(), 2024);
    for (x = 2; x <= 3; x++) {
        igraph_t g;

        IGRAPH_CHECK(random_x_partite(&g, x, 5, 10));
        IGRAPH_CHECK(test_all(&g, -1));
        igraph_destroy(&g);
    }

    return IGRAPH_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* Negative/validation tests                                             */
/* -------------------------------------------------------------------- */

static igraph_error_t test_errors(void) {
    igraph_t g;
    igraph_vector_int_list_t cover;
    igraph_vector_int_t bad_weights, good_weights;

    IGRAPH_CHECK(igraph_vector_int_list_init(&cover, 0));

    /* undirected graph -> IGRAPH_EINVAL */
    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0,1, 1,2, -1);
    CHECK_ERROR(igraph_minimum_path_cover(
            &g, &cover, NULL, NULL, IGRAPH_MPC_REDUCTION_NAIVE,
            IGRAPH_MPC_SOLVER_NAIVE_DFS, NULL, NULL, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);

    /* cyclic graph -> IGRAPH_EINVAL */
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, 2,0, -1);
    CHECK_ERROR(igraph_minimum_path_cover(
            &g, &cover, NULL, NULL, IGRAPH_MPC_REDUCTION_NAIVE,
            IGRAPH_MPC_SOLVER_NAIVE_DFS, NULL, NULL, NULL), IGRAPH_EINVAL);
    igraph_destroy(&g);

    /* mismatched vertex_weights size -> IGRAPH_EINVAL */
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0,1, 1,2, -1);
    igraph_vector_int_init(&bad_weights, 1);
    CHECK_ERROR(igraph_minimum_path_cover(
            &g, &cover, NULL, &bad_weights, IGRAPH_MPC_REDUCTION_NAIVE,
            IGRAPH_MPC_SOLVER_NAIVE_DFS, NULL, NULL, NULL), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&bad_weights);

    /* non-uniform vertex_weights with the greedy reduction -> IGRAPH_UNIMPLEMENTED */
    igraph_vector_int_init_int(&good_weights, 3, 1, 2, 1);
    CHECK_ERROR(igraph_minimum_path_cover(
            &g, &cover, NULL, &good_weights, IGRAPH_MPC_REDUCTION_GREEDY,
            IGRAPH_MPC_SOLVER_NAIVE_DFS, NULL, NULL, NULL), IGRAPH_UNIMPLEMENTED);
    igraph_vector_int_destroy(&good_weights);

    igraph_destroy(&g);

    igraph_vector_int_list_destroy(&cover);

    return IGRAPH_SUCCESS;
}

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 42);

    RUN_TEST(test_errors);
    RUN_TEST(test_fixtures);
    RUN_TEST(test_random_dags);
    RUN_TEST(test_binary_trees);
    RUN_TEST(test_x_partite);

    VERIFY_FINALLY_STACK();

    return 0;
}
