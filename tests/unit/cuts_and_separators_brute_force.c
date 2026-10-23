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

/* Compares igraph_all_st_cuts(), igraph_all_st_mincuts(),
 * igraph_minimum_size_separators() and igraph_all_minimal_st_separators()
 * with brute-force enumeration on small graphs. The checks are done on the
 * returned sets themselves, not on printed output, and use their own
 * connectivity code instead of other igraph cut or separator functions.
 * Vertex and edge sets are stored as bit masks. */

#include <igraph.h>

#include "test_utilities.h"

#include <stdint.h>
#include <stdlib.h>

#define MAX_M 40

typedef uint32_t vset;
typedef uint64_t eset;

typedef struct {
    int n, m;
    igraph_bool_t directed;
    int from[MAX_M], to[MAX_M];
} testgraph;

static void fail(const testgraph *g, const char *what) {
    fprintf(stderr, "FAILED: %s\n%s graph with %d vertices, edges:",
            what, g->directed ? "directed" : "undirected", g->n);
    for (int e = 0; e < g->m; e++) {
        fprintf(stderr, " %d-%d", g->from[e], g->to[e]);
    }
    fprintf(stderr, "\n");
    abort();
}

/* Vertices reachable from src through the vertices in 'allowed', ignoring
 * the edges in 'removed'. */
static vset reach(const testgraph *g, int src, vset allowed, eset removed) {
    vset seen = (vset) 1 << src, prev = 0;
    while (seen != prev) {
        prev = seen;
        for (int e = 0; e < g->m; e++) {
            vset u = (vset) 1 << g->from[e], v = (vset) 1 << g->to[e];
            if ((removed >> e) & 1) {
                continue;
            }
            if ((seen & u) && (allowed & v)) {
                seen |= v;
            }
            if (!g->directed && (seen & v) && (allowed & u)) {
                seen |= u;
            }
        }
    }
    return seen;
}

static vset all_vertices(const testgraph *g) {
    return ((vset) 1 << g->n) - 1;
}

static int cmp_eset(const void *a, const void *b) {
    eset x = *(const eset *) a, y = *(const eset *) b;
    return (x > y) - (x < y);
}

/* Bit masks of the sets in 'list', which must contain no repeated or
 * out-of-range elements. */
static size_t to_masks(const testgraph *g, const igraph_vector_int_list_t *list,
                       int limit, eset *masks) {
    size_t len = igraph_vector_int_list_size(list);
    for (size_t i = 0; i < len; i++) {
        const igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(list, i);
        masks[i] = 0;
        for (igraph_int_t j = 0; j < igraph_vector_int_size(v); j++) {
            igraph_int_t x = VECTOR(*v)[j];
            if (x < 0 || x >= limit || ((masks[i] >> x) & 1)) {
                fail(g, "invalid or repeated element in a returned set");
            }
            masks[i] |= (eset) 1 << x;
        }
    }
    return len;
}

/* The expected sets must be exactly the returned ones, each listed once. */
static void expect_sets(const testgraph *g, const char *what,
                        eset *expected, size_t n_expected,
                        eset *actual, size_t n_actual) {
    qsort(expected, n_expected, sizeof(eset), cmp_eset);
    qsort(actual, n_actual, sizeof(eset), cmp_eset);
    /* Brute force may find the same edge set from several vertex sets. */
    size_t k = 0;
    for (size_t i = 0; i < n_expected; i++) {
        if (k == 0 || expected[k - 1] != expected[i]) {
            expected[k++] = expected[i];
        }
    }
    if (k != n_actual) {
        fprintf(stderr, "expected %zu sets, got %zu\n", k, n_actual);
        fail(g, what);
    }
    for (size_t i = 0; i < k; i++) {
        if (expected[i] != actual[i]) {
            fprintf(stderr, "expected 0x%llx, got 0x%llx\n",
                    (unsigned long long) expected[i], (unsigned long long) actual[i]);
            fail(g, what);
        }
    }
}

/* Edges leaving the vertex set X, and their total capacity. */
static eset leaving(const testgraph *g, vset X) {
    eset result = 0;
    for (int e = 0; e < g->m; e++) {
        if (((X >> g->from[e]) & 1) && !((X >> g->to[e]) & 1)) {
            result |= (eset) 1 << e;
        }
    }
    return result;
}

static long weight(const testgraph *g, eset cut, const int *capacity) {
    long w = 0;
    for (int e = 0; e < g->m; e++) {
        w += ((cut >> e) & 1) ? (capacity ? capacity[e] : 1) : 0;
    }
    return w;
}

static igraph_bool_t is_minimal_cut(const testgraph *g, int s, int t, eset cut) {
    if ((reach(g, s, all_vertices(g), cut) >> t) & 1) {
        return false;
    }
    for (int e = 0; e < g->m; e++) {
        if (((cut >> e) & 1) && !((reach(g, s, all_vertices(g), cut & ~((eset) 1 << e)) >> t) & 1)) {
            return false;
        }
    }
    return true;
}

static void make_graph(const testgraph *g, igraph_t *graph) {
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, 0);
    for (int e = 0; e < g->m; e++) {
        igraph_vector_int_push_back(&edges, g->from[e]);
        igraph_vector_int_push_back(&edges, g->to[e]);
    }
    IGRAPH_ASSERT(igraph_create(graph, &edges, g->n, g->directed) == IGRAPH_SUCCESS);
    igraph_vector_int_destroy(&edges);
}

/* Checks igraph_all_st_cuts() (all minimal cuts) or igraph_all_st_mincuts()
 * (all minimum capacity cuts). The target must be reachable from the source. */
static void check_cuts(const testgraph *g, int s, int t, const int *capacity, igraph_bool_t min_only) {
    static eset expected[1 << 18], cuts_masks[1 << 16], parts_masks[1 << 16];
    igraph_t graph;
    igraph_vector_int_list_t cuts, parts;
    igraph_vector_t cap;
    igraph_real_t value;
    long best = -1;
    size_t n_expected = 0;

    for (vset X = 0; X <= all_vertices(g); X++) {
        if (!((X >> s) & 1) || ((X >> t) & 1)) {
            continue;
        }
        eset c = leaving(g, X);
        if (min_only) {
            long w = weight(g, c, capacity);
            if (best < 0 || w < best) {
                best = w;
                n_expected = 0;
            }
            if (w != best) {
                continue;
            }
        } else if (!is_minimal_cut(g, s, t, c)) {
            continue;
        }
        expected[n_expected++] = c;
    }

    make_graph(g, &graph);
    igraph_vector_int_list_init(&cuts, 0);
    igraph_vector_int_list_init(&parts, 0);
    igraph_vector_init(&cap, g->m);
    for (int e = 0; capacity && e < g->m; e++) {
        VECTOR(cap)[e] = capacity[e];
    }
    if (min_only) {
        IGRAPH_ASSERT(igraph_all_st_mincuts(&graph, &value, &cuts, &parts, s, t, capacity ? &cap : NULL) == IGRAPH_SUCCESS);
        if (value != (igraph_real_t) best) {
            fail(g, "all_st_mincuts returned the wrong cut value");
        }
    } else {
        IGRAPH_ASSERT(igraph_all_st_cuts(&graph, &cuts, &parts, s, t) == IGRAPH_SUCCESS);
    }

    size_t n_actual = to_masks(g, &cuts, g->m, cuts_masks);
    IGRAPH_ASSERT(to_masks(g, &parts, g->n, parts_masks) == n_actual);
    for (size_t i = 0; i < n_actual; i++) {
        vset X = (vset) parts_masks[i];
        /* Each cut must be the set of edges leaving its partition. */
        if (!((X >> s) & 1) || ((X >> t) & 1) || leaving(g, X) != cuts_masks[i]) {
            fail(g, "a cut does not match its partition");
        }
    }
    expect_sets(g, min_only ? "all_st_mincuts" : "all_st_cuts", expected, n_expected, cuts_masks, n_actual);

    igraph_vector_destroy(&cap);
    igraph_vector_int_list_destroy(&cuts);
    igraph_vector_int_list_destroy(&parts);
    igraph_destroy(&graph);
}

static int size_of(eset set) {
    int size = 0;
    for (; set; set &= set - 1) {
        size++;
    }
    return size;
}

static int components_without(const testgraph *g, vset S) {
    vset left = all_vertices(g) & ~S;
    int count = 0;
    for (int v = 0; v < g->n; v++) {
        if ((left >> v) & 1) {
            left &= ~reach(g, v, left, 0);
            count++;
        }
    }
    return count;
}

/* Is S, which does not contain s and t, a minimal (s,t) separator? */
static igraph_bool_t is_minimal_separator(const testgraph *g, vset S, int s, int t) {
    if (((S >> s) & 1) || ((S >> t) & 1) || ((reach(g, s, all_vertices(g) & ~S, 0) >> t) & 1)) {
        return false;
    }
    for (int v = 0; v < g->n; v++) {
        vset smaller = S & ~((vset) 1 << v);
        if (((S >> v) & 1) && !((reach(g, s, all_vertices(g) & ~smaller, 0) >> t) & 1)) {
            return false;
        }
    }
    return true;
}

/* Checks igraph_minimum_size_separators() (min_size) or
 * igraph_all_minimal_st_separators(). Like igraph, the latter is only compared
 * on non-empty separators, see https://github.com/igraph/igraph/issues/2517 */
static void check_separators(const testgraph *g, igraph_bool_t min_size) {
    static eset expected[1 << 10], actual[1 << 10];
    igraph_t graph;
    igraph_vector_int_list_t separators;
    size_t n_expected = 0;

    for (vset S = 1; S <= all_vertices(g); S++) {
        igraph_bool_t found = false;
        for (int s = 0; s < g->n && !found; s++) {
            for (int t = s + 1; t < g->n && !found; t++) {
                found = is_minimal_separator(g, S, s, t);
            }
        }
        if (found) {
            expected[n_expected++] = S;
        }
    }
    if (min_size) {
        /* A minimum size separator is a minimal separator of the smallest size.
         * There are none if the graph is disconnected. */
        int smallest = g->n;
        size_t k = 0;
        for (size_t i = 0; i < n_expected; i++) {
            smallest = size_of(expected[i]) < smallest ? size_of(expected[i]) : smallest;
        }
        for (size_t i = 0; i < n_expected; i++) {
            if (size_of(expected[i]) == smallest && components_without(g, 0) == 1) {
                expected[k++] = expected[i];
            }
        }
        n_expected = k;
    }

    make_graph(g, &graph);
    igraph_vector_int_list_init(&separators, 0);
    if (min_size) {
        IGRAPH_ASSERT(igraph_minimum_size_separators(&graph, &separators) == IGRAPH_SUCCESS);
    } else {
        IGRAPH_ASSERT(igraph_all_minimal_st_separators(&graph, &separators) == IGRAPH_SUCCESS);
    }
    size_t n_actual = to_masks(g, &separators, g->n, actual);
    expect_sets(g, min_size ? "minimum_size_separators" : "all_minimal_st_separators",
                expected, n_expected, actual, n_actual);

    igraph_vector_int_list_destroy(&separators);
    igraph_destroy(&graph);
}

static void random_graph(testgraph *g, int n, igraph_bool_t directed) {
    int percent = (int) RNG_INTEGER(15, 70);
    g->n = n;
    g->m = 0;
    g->directed = directed;
    for (int u = 0; u < n; u++) {
        for (int v = directed ? 0 : u + 1; v < n; v++) {
            if (u != v && (int) RNG_INTEGER(0, 99) < percent) {
                g->from[g->m] = u;
                g->to[g->m] = v;
                g->m++;
            }
        }
    }
}

static void from_edges(testgraph *g, int n, igraph_bool_t directed, const int (*edges)[2], int m) {
    g->n = n;
    g->m = m;
    g->directed = directed;
    for (int e = 0; e < m; e++) {
        g->from[e] = edges[e][0];
        g->to[e] = edges[e][1];
    }
}

int main(void) {
    testgraph g;
    int capacity[MAX_M];

    igraph_rng_seed(igraph_rng_default(), 20260521);

    /* Graph from https://github.com/igraph/igraph/issues/2532 */
    static const int edges_2532[][2] = {
        {1, 11}, {3, 13}, {4, 14}, {5, 15}, {6, 16}, {7, 17}, {8, 18}, {9, 19},
        {0, 8}, {18, 0}, {0, 9}, {19, 0}, {2, 1}, {11, 2}, {14, 1}, {11, 4},
        {2, 4}, {14, 2}, {2, 5}, {15, 2}, {14, 3}, {13, 4}, {18, 3}, {13, 8},
        {16, 5}, {15, 6}, {17, 5}, {15, 7}, {19, 6}, {16, 9}, {18, 7}, {17, 8}
    };
    from_edges(&g, 20, true, edges_2532, 32);
    check_cuts(&g, 2, 0, NULL, true);

    /* Graph from https://github.com/igraph/igraph/issues/2521 */
    static const int edges_2521[][2] = {
        {0, 8}, {0, 9}, {1, 2}, {1, 4}, {2, 4}, {2, 5},
        {3, 4}, {3, 8}, {5, 6}, {5, 7}, {6, 9}, {7, 8}
    };
    from_edges(&g, 10, false, edges_2521, 12);
    check_separators(&g, true);
    check_separators(&g, false);

    for (int round = 0; round < 12; round++) {
        int n = (int) RNG_INTEGER(2, 6);

        random_graph(&g, n, true);
        for (int e = 0; e < g.m; e++) {
            capacity[e] = (int) RNG_INTEGER(1, 4);
        }
        for (int s = 0; s < n; s++) {
            for (int t = 0; t < n; t++) {
                /* The result for an unreachable target is not documented, skip it. */
                if (s != t && ((reach(&g, s, all_vertices(&g), 0) >> t) & 1)) {
                    check_cuts(&g, s, t, NULL, false);
                    check_cuts(&g, s, t, NULL, true);
                    check_cuts(&g, s, t, capacity, true);
                }
            }
        }

        random_graph(&g, n, false);
        check_separators(&g, true);
        check_separators(&g, false);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
