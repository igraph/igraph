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

/* Exact-set correctness tests for igraph_simple_cycles(), added for #2911.
 *
 * Every cycle is described canonically by the sorted list of its edge IDs,
 * and the list of cycles is sorted. Two results are equal only if they contain
 * exactly the same cycles, each exactly once.
 *
 * - The results are compared with an independent brute-force enumerator that
 *   tests every subset of the edges of a small graph. It shares no code with
 *   Johnson's algorithm.
 * - On larger graphs, the result must not change when the vertices are
 *   permuted: igraph_permute_vertices() keeps edge IDs, so the canonical
 *   descriptions must be identical. This idea is adapted from
 *   check_cycles_permuted() in PR #2916 by @GenieTim.
 * - Explicit cases: the #2911 graph under all vertex permutations, multigraphs,
 *   and early termination through the callback and through max_results.
 */

#include <igraph.h>

#include "test_utilities.h"

#define ORACLE_MAX_EDGES 12

static igraph_int_t failures = 0;

static const char *mode_name(const igraph_t *graph, igraph_neimode_t mode) {
    if (!igraph_is_directed(graph)) {
        return "undirected";
    }
    switch (mode) {
    case IGRAPH_OUT: return "directed OUT";
    case IGRAPH_IN:  return "directed IN";
    default:         return "directed ALL";
    }
}

/* Sorts the edge IDs of each cycle, then sorts the list of cycles. */
static void canonicalize(igraph_vector_int_list_t *cycles) {
    const igraph_int_t n = igraph_vector_int_list_size(cycles);
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_vector_int_sort(igraph_vector_int_list_get_ptr(cycles, i));
    }
    igraph_vector_int_list_sort(cycles, &igraph_vector_int_lex_cmp);
}

static igraph_bool_t cycle_lists_equal(const igraph_vector_int_list_t *a,
                                       const igraph_vector_int_list_t *b) {
    const igraph_int_t n = igraph_vector_int_list_size(a);
    if (igraph_vector_int_list_size(b) != n) {
        return false;
    }
    for (igraph_int_t i = 0; i < n; i++) {
        if (!igraph_vector_int_all_e(igraph_vector_int_list_get_ptr(a, i),
                                     igraph_vector_int_list_get_ptr(b, i))) {
            return false;
        }
    }
    return true;
}

static void report_failure(const char *what, const igraph_t *graph, igraph_neimode_t mode,
                           igraph_int_t min_len, igraph_int_t max_len,
                           const igraph_vector_int_list_t *expected,
                           const igraph_vector_int_list_t *found) {
    failures++;
    if (failures > 5) {
        return;
    }
    printf("FAILURE (%s): %s, min_cycle_length = %" IGRAPH_PRId
           ", max_cycle_length = %" IGRAPH_PRId "\n",
           what, mode_name(graph, mode), min_len, max_len);
    print_graph(graph);
    if (expected) {
        printf("Expected cycles (edge IDs):\n");
        print_vector_int_list(expected);
    }
    if (found) {
        printf("Found cycles (edge IDs):\n");
        print_vector_int_list(found);
    }
    printf("\n");
}

/* Brute-force enumeration of all simple cycles as sorted edge-ID sets.
 *
 * A non-empty edge subset is a simple cycle if it is connected and
 *  - directed graphs, OUT or IN mode: each of its vertices has in-degree 1
 *    and out-degree 1 within the subset;
 *  - undirected graphs or ALL mode: each of its vertices has degree 2 within
 *    the subset, counting a self-loop twice.
 * OUT and IN mode find the same edge sets, traversed in opposite directions.
 * A single self-loop is a cycle of length 1, and two parallel edges (or, in
 * ALL mode, a mutual pair) form a cycle of length 2. The length of a cycle
 * is its number of edges. Negative limits mean no limit. */
static void brute_force_cycles(const igraph_t *graph, igraph_neimode_t mode,
                               igraph_int_t min_len, igraph_int_t max_len,
                               igraph_vector_int_list_t *result) {
    const igraph_int_t n = igraph_vcount(graph);
    const igraph_int_t m = igraph_ecount(graph);
    const igraph_bool_t directed = igraph_is_directed(graph) && mode != IGRAPH_ALL;
    igraph_vector_int_t indeg, outdeg, parent, cycle;

    IGRAPH_ASSERT(m <= ORACLE_MAX_EDGES);

    igraph_vector_int_list_clear(result);
    igraph_vector_int_init(&indeg, n);
    igraph_vector_int_init(&outdeg, n);
    igraph_vector_int_init(&parent, n);
    igraph_vector_int_init(&cycle, 0);

    for (igraph_uint_t mask = 1; mask < ((igraph_uint_t) 1 << m); mask++) {
        igraph_bool_t ok = true;
        igraph_int_t components = 0;

        igraph_vector_int_clear(&cycle);
        for (igraph_int_t e = 0; e < m; e++) {
            if (mask & ((igraph_uint_t) 1 << e)) {
                igraph_vector_int_push_back(&cycle, e);
            }
        }
        const igraph_int_t len = igraph_vector_int_size(&cycle);
        if (max_len >= 0 && len > max_len) {
            continue;
        }
        if (len < min_len) {
            continue;
        }

        igraph_vector_int_null(&indeg);
        igraph_vector_int_null(&outdeg);
        for (igraph_int_t v = 0; v < n; v++) {
            VECTOR(parent)[v] = v;
        }

        for (igraph_int_t i = 0; i < len; i++) {
            const igraph_int_t e = VECTOR(cycle)[i];
            igraph_int_t u = IGRAPH_FROM(graph, e), v = IGRAPH_TO(graph, e);
            VECTOR(outdeg)[u]++;
            VECTOR(indeg)[v]++;
            /* union-find, without balancing; the graphs are tiny */
            while (VECTOR(parent)[u] != u) {
                u = VECTOR(parent)[u];
            }
            while (VECTOR(parent)[v] != v) {
                v = VECTOR(parent)[v];
            }
            VECTOR(parent)[u] = v;
        }

        for (igraph_int_t v = 0; v < n && ok; v++) {
            const igraph_int_t in = VECTOR(indeg)[v], out = VECTOR(outdeg)[v];
            if (in == 0 && out == 0) {
                continue; /* not in the subset */
            }
            if (directed) {
                ok = (in == 1 && out == 1);
            } else {
                ok = (in + out == 2);
            }
            if (VECTOR(parent)[v] == v) {
                components++;
            }
        }
        if (ok && components == 1) {
            igraph_vector_int_list_push_back_copy(result, &cycle);
        }
    }

    igraph_vector_int_destroy(&cycle);
    igraph_vector_int_destroy(&parent);
    igraph_vector_int_destroy(&outdeg);
    igraph_vector_int_destroy(&indeg);

    canonicalize(result);
}

/* Checks that 'vertices' and 'edges' describe the same closed walk in the
 * graph, that it follows the edge directions required by 'mode', and that no
 * vertex is repeated. */
static igraph_bool_t is_valid_cycle(const igraph_t *graph, igraph_neimode_t mode,
                                    const igraph_vector_int_t *vertices,
                                    const igraph_vector_int_t *edges) {
    const igraph_int_t len = igraph_vector_int_size(edges);
    const igraph_bool_t directed = igraph_is_directed(graph) && mode != IGRAPH_ALL;
    igraph_bitset_t seen;
    igraph_bool_t ok = true;

    if (len == 0 || igraph_vector_int_size(vertices) != len) {
        return false;
    }

    igraph_bitset_init(&seen, igraph_vcount(graph));
    for (igraph_int_t i = 0; i < len && ok; i++) {
        const igraph_int_t v = VECTOR(*vertices)[i];
        const igraph_int_t w = VECTOR(*vertices)[(i + 1) % len];
        const igraph_int_t e = VECTOR(*edges)[i];
        const igraph_int_t from = IGRAPH_FROM(graph, e), to = IGRAPH_TO(graph, e);

        if (IGRAPH_BIT_TEST(seen, v)) {
            ok = false;
        }
        IGRAPH_BIT_SET(seen, v);

        if (!directed) {
            ok = ok && ((from == v && to == w) || (from == w && to == v));
        } else if (mode == IGRAPH_OUT) {
            ok = ok && (from == v && to == w);
        } else {
            ok = ok && (from == w && to == v);
        }
    }
    igraph_bitset_destroy(&seen);

    return ok;
}

/* Runs igraph_simple_cycles(), checks that each cycle is valid and within the
 * length limits, and stores the canonical edge-set description in 'result'. */
static void find_cycles_canonical(const igraph_t *graph, igraph_neimode_t mode,
                                  igraph_int_t min_len, igraph_int_t max_len,
                                  igraph_vector_int_list_t *result) {
    igraph_vector_int_list_t vertices;

    igraph_vector_int_list_init(&vertices, 0);
    igraph_simple_cycles(graph, &vertices, result, mode, min_len, max_len, IGRAPH_UNLIMITED);

    const igraph_int_t n = igraph_vector_int_list_size(result);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&vertices) == n);
    for (igraph_int_t i = 0; i < n; i++) {
        const igraph_vector_int_t *edges = igraph_vector_int_list_get_ptr(result, i);
        const igraph_int_t len = igraph_vector_int_size(edges);
        if (!is_valid_cycle(graph, mode, igraph_vector_int_list_get_ptr(&vertices, i), edges) ||
            (max_len >= 0 && len > max_len) || len < min_len) {
            report_failure("invalid cycle", graph, mode, min_len, max_len, NULL, result);
            break;
        }
    }
    igraph_vector_int_list_destroy(&vertices);

    canonicalize(result);
}

/* Compares igraph_simple_cycles() with the brute-force enumerator. */
static void check_against_oracle(const igraph_t *graph, igraph_neimode_t mode,
                                 igraph_int_t min_len, igraph_int_t max_len) {
    igraph_vector_int_list_t expected, found;

    igraph_vector_int_list_init(&expected, 0);
    igraph_vector_int_list_init(&found, 0);

    brute_force_cycles(graph, mode, min_len, max_len, &expected);
    find_cycles_canonical(graph, mode, min_len, max_len, &found);

    if (!cycle_lists_equal(&expected, &found)) {
        report_failure("differs from brute force", graph, mode, min_len, max_len, &expected, &found);
    }

    igraph_vector_int_list_destroy(&found);
    igraph_vector_int_list_destroy(&expected);
}

/* Checks that the canonical result is the same for 'graph' and for
 * 'perm_count' random permutations of its vertices, plus the reversal. */
static void check_permutation_invariance(const igraph_t *graph, igraph_neimode_t mode,
                                         igraph_int_t min_len, igraph_int_t max_len,
                                         igraph_int_t perm_count) {
    const igraph_int_t n = igraph_vcount(graph);
    igraph_vector_int_list_t base, permuted_cycles;
    igraph_vector_int_t perm;
    igraph_t permuted;

    igraph_vector_int_list_init(&base, 0);
    igraph_vector_int_list_init(&permuted_cycles, 0);
    igraph_vector_int_init(&perm, n);

    find_cycles_canonical(graph, mode, min_len, max_len, &base);

    for (igraph_int_t p = 0; p <= perm_count; p++) {
        if (p == 0) {
            for (igraph_int_t i = 0; i < n; i++) {
                VECTOR(perm)[i] = n - 1 - i;
            }
        } else {
            for (igraph_int_t i = 0; i < n; i++) {
                VECTOR(perm)[i] = i;
            }
            igraph_vector_int_shuffle(&perm);
        }
        igraph_permute_vertices(graph, &permuted, &perm);
        find_cycles_canonical(&permuted, mode, min_len, max_len, &permuted_cycles);
        if (!cycle_lists_equal(&base, &permuted_cycles)) {
            report_failure("not invariant under vertex permutation", &permuted, mode,
                           min_len, max_len, &base, &permuted_cycles);
        }
        igraph_destroy(&permuted);
    }

    igraph_vector_int_destroy(&perm);
    igraph_vector_int_list_destroy(&permuted_cycles);
    igraph_vector_int_list_destroy(&base);
}

/* Random multigraph with self-loops; 'kind' selects the directedness and mode
 * so that all four combinations are covered in turn. */
static void random_graph(igraph_t *graph, igraph_int_t n, igraph_int_t m, int kind,
                         igraph_neimode_t *mode) {
    static const igraph_neimode_t modes[4] = { IGRAPH_ALL, IGRAPH_OUT, IGRAPH_IN, IGRAPH_ALL };
    igraph_erdos_renyi_game_gnm(graph, n, m, /* directed = */ kind != 0,
                                IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, /* edge_labeled = */ false);
    *mode = modes[kind];
}

/* Cycle length limits to try: none, a minimum, a maximum, both. */
static void length_limits(igraph_int_t t, igraph_int_t *min_len, igraph_int_t *max_len) {
    *min_len = (t % 7 == 3) ? 3 : -1;
    *max_len = (t % 5 == 2) ? 4 : (t % 11 == 6 ? 2 : -1);
}

static void test_oracle_random(void) {
    const igraph_int_t trials = 800;

    for (igraph_int_t t = 0; t < trials; t++) {
        igraph_t g;
        igraph_neimode_t mode;
        igraph_int_t min_len, max_len;
        const igraph_int_t n = RNG_INTEGER(1, 8);
        igraph_int_t m = RNG_INTEGER(n / 2, 2 * n);
        if (m > ORACLE_MAX_EDGES) {
            m = ORACLE_MAX_EDGES;
        }

        random_graph(&g, n, m, (int) (t % 4), &mode);
        length_limits(t, &min_len, &max_len);
        check_against_oracle(&g, mode, min_len, max_len);
        igraph_destroy(&g);
    }
}

static void test_permutations_random(void) {
    const igraph_int_t trials = 200;

    for (igraph_int_t t = 0; t < trials; t++) {
        igraph_t g;
        igraph_neimode_t mode;
        igraph_int_t min_len, max_len;
        const int kind = (int) (t % 4);
        const igraph_int_t n = RNG_INTEGER(5, 30);
        /* Keep undirected / ALL graphs sparser: they have many more cycles. */
        const igraph_int_t m = (kind == 1 || kind == 2) ?
                               RNG_INTEGER(n / 2, 2 * n) : RNG_INTEGER(n / 2, n + 6);

        random_graph(&g, n, m, kind, &mode);
        length_limits(t, &min_len, &max_len);
        check_permutation_invariance(&g, mode, min_len, max_len, 4);
        igraph_destroy(&g);
    }
}

/* The graph from #2911 under every vertex permutation. The only cycle is
 * 2 <-> 5, i.e. edges 3 and 6. */
static void test_issue_2911(void) {
    igraph_t g, permuted;
    igraph_vector_int_t perm;
    igraph_vector_int_list_t found;
    const igraph_neimode_t modes[3] = { IGRAPH_OUT, IGRAPH_IN, IGRAPH_ALL };
    const igraph_int_t n = 6;

    igraph_small(&g, n, IGRAPH_DIRECTED,
                 0, 1, 0, 5, 1, 5, 2, 5, 4, 3, 4, 5, 5, 2, 5, 3,
                 -1);

    igraph_vector_int_list_init(&found, 0);
    find_cycles_canonical(&g, IGRAPH_OUT, -1, -1, &found);
    if (igraph_vector_int_list_size(&found) != 1 ||
        igraph_vector_int_size(igraph_vector_int_list_get_ptr(&found, 0)) != 2 ||
        VECTOR(*igraph_vector_int_list_get_ptr(&found, 0))[0] != 3 ||
        VECTOR(*igraph_vector_int_list_get_ptr(&found, 0))[1] != 6) {
        report_failure("issue #2911 graph", &g, IGRAPH_OUT, -1, -1, NULL, &found);
    }
    igraph_vector_int_list_destroy(&found);

    /* Enumerate all 6! permutations in lexicographic order. */
    igraph_vector_int_init_range(&perm, 0, n);
    for (;;) {
        igraph_permute_vertices(&g, &permuted, &perm);
        for (int k = 0; k < 3; k++) {
            check_against_oracle(&permuted, modes[k], -1, -1);
        }
        igraph_destroy(&permuted);

        /* next permutation */
        igraph_int_t i = n - 2;
        while (i >= 0 && VECTOR(perm)[i] > VECTOR(perm)[i + 1]) {
            i--;
        }
        if (i < 0) {
            break;
        }
        igraph_int_t j = n - 1;
        while (VECTOR(perm)[j] < VECTOR(perm)[i]) {
            j--;
        }
        igraph_vector_int_swap_elements(&perm, i, j);
        for (igraph_int_t a = i + 1, b = n - 1; a < b; a++, b--) {
            igraph_vector_int_swap_elements(&perm, a, b);
        }
    }

    igraph_vector_int_destroy(&perm);
    igraph_destroy(&g);
}

/* Checks the cycle count against a hand-computed value, and the cycles
 * themselves against the brute-force enumerator, for all four
 * directedness / mode combinations given. */
static void check_count_and_oracle(const igraph_t *graph, igraph_neimode_t mode,
                                   igraph_int_t expected_count) {
    igraph_vector_int_list_t found;

    igraph_vector_int_list_init(&found, 0);
    find_cycles_canonical(graph, mode, -1, -1, &found);
    if (igraph_vector_int_list_size(&found) != expected_count) {
        printf("Expected %" IGRAPH_PRId " cycles, found %" IGRAPH_PRId ".\n",
               expected_count, igraph_vector_int_list_size(&found));
        report_failure("wrong cycle count", graph, mode, -1, -1, NULL, &found);
    }
    igraph_vector_int_list_destroy(&found);

    check_against_oracle(graph, mode, -1, -1);
}

static void test_multigraphs(void) {
    igraph_t g;

    /* Triple edge 0 -> 1 and a single edge 1 -> 0: three 2-cycles when
     * directions are followed; in ALL mode, any two of the four edges. */
    igraph_small(&g, 2, IGRAPH_DIRECTED, 0, 1, 0, 1, 0, 1, 1, 0, -1);
    check_count_and_oracle(&g, IGRAPH_OUT, 3);
    check_count_and_oracle(&g, IGRAPH_IN, 3);
    check_count_and_oracle(&g, IGRAPH_ALL, 6);
    igraph_destroy(&g);

    /* Undirected triple edge: three 2-cycles. */
    igraph_small(&g, 2, IGRAPH_UNDIRECTED, 0, 1, 0, 1, 0, 1, -1);
    check_count_and_oracle(&g, IGRAPH_ALL, 3);
    igraph_destroy(&g);

    /* Mutual pairs on a directed triangle: in ALL mode each mutual pair is a
     * 2-cycle, and each of the 2^3 choices of one edge per pair is a 3-cycle. */
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0, 1, 1, 0, 1, 2, 2, 1, 2, 0, 0, 2, -1);
    check_count_and_oracle(&g, IGRAPH_OUT, 5);
    check_count_and_oracle(&g, IGRAPH_IN, 5);
    check_count_and_oracle(&g, IGRAPH_ALL, 3 + 8);
    igraph_destroy(&g);

    /* Multiple self-loops on one vertex, plus a self-loop elsewhere,
     * attached to a 2-cycle. Each self-loop is a cycle of its own. */
    igraph_small(&g, 3, IGRAPH_DIRECTED, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 2, -1);
    check_count_and_oracle(&g, IGRAPH_OUT, 5);
    check_count_and_oracle(&g, IGRAPH_IN, 5);
    check_count_and_oracle(&g, IGRAPH_ALL, 5);
    igraph_destroy(&g);

    igraph_small(&g, 3, IGRAPH_UNDIRECTED, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 2, 2, 2, -1);
    check_count_and_oracle(&g, IGRAPH_ALL, 5);
    igraph_destroy(&g);

    /* A directed ring of degree-2 vertices whose edges are all doubled, and a
     * second ring hanging off it: exercises the start-vertex skip with
     * parallel edges. */
    igraph_small(&g, 6, IGRAPH_DIRECTED,
                 0, 1, 0, 1, 1, 2, 2, 0, 2, 0,
                 3, 4, 4, 5, 5, 3, 3, 1,
                 -1);
    check_count_and_oracle(&g, IGRAPH_OUT, 4 + 1);
    check_count_and_oracle(&g, IGRAPH_IN, 4 + 1);
    igraph_destroy(&g);
}

typedef struct {
    igraph_int_t stop_after;
    igraph_int_t calls;
    igraph_vector_int_list_t edges;
} stop_state_t;

static igraph_error_t stop_after_k(const igraph_vector_int_t *vertices,
                                   const igraph_vector_int_t *edges,
                                   void *arg) {
    stop_state_t *st = (stop_state_t *) arg;
    IGRAPH_UNUSED(vertices);
    st->calls++;
    IGRAPH_ASSERT(st->calls <= st->stop_after); /* never called after IGRAPH_STOP */
    IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(&st->edges, edges));
    return st->calls == st->stop_after ? IGRAPH_STOP : IGRAPH_SUCCESS;
}

/* Every cycle in 'part' must be distinct and occur in 'all' (both canonical). */
static igraph_bool_t is_distinct_subset(const igraph_vector_int_list_t *part,
                                        const igraph_vector_int_list_t *all) {
    const igraph_int_t n = igraph_vector_int_list_size(part);
    for (igraph_int_t i = 0; i < n; i++) {
        const igraph_vector_int_t *c = igraph_vector_int_list_get_ptr(part, i);
        igraph_bool_t present = false;
        if (i > 0 && igraph_vector_int_all_e(c, igraph_vector_int_list_get_ptr(part, i - 1))) {
            return false;
        }
        for (igraph_int_t j = 0; j < igraph_vector_int_list_size(all) && !present; j++) {
            present = igraph_vector_int_all_e(c, igraph_vector_int_list_get_ptr(all, j));
        }
        if (!present) {
            return false;
        }
    }
    return true;
}

static void check_early_stop(const igraph_t *graph, igraph_neimode_t mode) {
    igraph_vector_int_list_t all, vertices, edges;

    igraph_vector_int_list_init(&all, 0);
    igraph_vector_int_list_init(&vertices, 0);
    igraph_vector_int_list_init(&edges, 0);

    find_cycles_canonical(graph, mode, -1, -1, &all);
    const igraph_int_t total = igraph_vector_int_list_size(&all);
    IGRAPH_ASSERT(total >= 3);

    const igraph_int_t ks[] = { 1, 2, total / 2, total - 1, total, total + 5 };
    for (size_t i = 0; i < sizeof(ks) / sizeof(ks[0]); i++) {
        const igraph_int_t k = ks[i];
        const igraph_int_t expected = k < total ? k : total;

        /* Callback returning IGRAPH_STOP after k cycles. */
        stop_state_t st;
        st.stop_after = k;
        st.calls = 0;
        igraph_vector_int_list_init(&st.edges, 0);
        /* IGRAPH_STOP from the callback must be converted to success. A leaked
         * IGRAPH_STOP would not trigger the error handler, so check it here. */
        IGRAPH_ASSERT(igraph_simple_cycles_callback(graph, mode, -1, -1, &stop_after_k, &st)
                      == IGRAPH_SUCCESS);
        IGRAPH_ASSERT(st.calls == expected);
        canonicalize(&st.edges);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&st.edges) == expected);
        IGRAPH_ASSERT(is_distinct_subset(&st.edges, &all));
        igraph_vector_int_list_destroy(&st.edges);

        /* max_results = k */
        igraph_simple_cycles(graph, &vertices, &edges, mode, -1, -1, k);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&vertices) == expected);
        IGRAPH_ASSERT(igraph_vector_int_list_size(&edges) == expected);
        for (igraph_int_t j = 0; j < expected; j++) {
            IGRAPH_ASSERT(is_valid_cycle(graph, mode,
                                         igraph_vector_int_list_get_ptr(&vertices, j),
                                         igraph_vector_int_list_get_ptr(&edges, j)));
        }
        canonicalize(&edges);
        IGRAPH_ASSERT(is_distinct_subset(&edges, &all));
    }

    igraph_vector_int_list_destroy(&edges);
    igraph_vector_int_list_destroy(&vertices);
    igraph_vector_int_list_destroy(&all);
}

static void test_early_stop(void) {
    igraph_t g;

    igraph_full(&g, 5, IGRAPH_DIRECTED, /* loops = */ false);
    check_early_stop(&g, IGRAPH_OUT);
    check_early_stop(&g, IGRAPH_IN);
    igraph_destroy(&g);

    igraph_wheel(&g, 8, IGRAPH_WHEEL_UNDIRECTED, 0);
    check_early_stop(&g, IGRAPH_ALL);
    igraph_destroy(&g);

    /* A multigraph with self-loops, in ALL mode. */
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0, 1, 1, 0, 1, 2, 2, 3, 3, 1, 3, 3, 2, 0, -1);
    check_early_stop(&g, IGRAPH_ALL);
    check_early_stop(&g, IGRAPH_OUT);
    igraph_destroy(&g);
}

/* Runs one group of tests and reports how many checks failed in it. */
#define RUN_SECTION(name) \
    do { \
        const igraph_int_t before = failures; \
        name(); \
        VERIFY_FINALLY_STACK(); \
        if (failures > before) { \
            printf(#name ": %" IGRAPH_PRId " failed checks.\n", failures - before); \
        } \
    } while (0)

int main(void) {
    igraph_rng_seed(igraph_rng_default(), 2911);

    RUN_SECTION(test_issue_2911);
    RUN_SECTION(test_multigraphs);
    RUN_SECTION(test_early_stop);
    RUN_SECTION(test_oracle_random);
    RUN_SECTION(test_permutations_random);

    if (failures > 0) {
        printf("%" IGRAPH_PRId " failed checks in total.\n", failures);
        return 1;
    }

    return 0;
}
