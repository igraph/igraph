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

/* Checks the cut and separator enumeration functions against brute force on
 * small random graphs. The existing tests compare printed output with stored
 * results, which cannot show that the stored results are correct.
 *
 * Related issues: igraph/igraph#2532, igraph/igraph#2521, igraph/igraph#2536.
 */

#include <igraph.h>

#include "test_utilities.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define BIT(i) (((uint64_t) 1) << (i))
#define MAX_SETS (1 << 12)

static int compare_masks(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *) a, y = *(const uint64_t *) b;
    return (x > y) - (x < y);
}

static int popcount(uint64_t x) {
    int c = 0;
    for (; x; x &= x - 1) {
        c++;
    }
    return c;
}

static uint64_t to_mask(const igraph_vector_int_t *ids) {
    uint64_t mask = 0;
    for (igraph_int_t i = 0; i < igraph_vector_int_size(ids); i++) {
        mask |= BIT(VECTOR(*ids)[i]);
    }
    return mask;
}

/* Sorted masks of a list of sets. Duplicates are kept, so that comparing with
 * the duplicate-free brute-force result detects them. */
static size_t sorted_masks(const igraph_vector_int_list_t *list, uint64_t *masks) {
    size_t count = (size_t) igraph_vector_int_list_size(list);
    IGRAPH_ASSERT(count <= MAX_SETS);
    for (size_t i = 0; i < count; i++) {
        masks[i] = to_mask(igraph_vector_int_list_get_ptr(list, (igraph_int_t) i));
    }
    qsort(masks, count, sizeof(uint64_t), compare_masks);
    return count;
}

static void assert_same(const uint64_t *a, size_t na, const uint64_t *b, size_t nb) {
    IGRAPH_ASSERT(na == nb);
    for (size_t i = 0; i < na; i++) {
        IGRAPH_ASSERT(a[i] == b[i]);
    }
}

/* Edge cuts in directed graphs */

/* Is 'target' reachable from 'source' when the edges in 'removed' are deleted? */
static igraph_bool_t reachable(const igraph_t *graph, uint64_t removed,
                               igraph_int_t source, igraph_int_t target) {
    uint64_t seen = BIT(source), old = 0;
    while (seen != old) {
        old = seen;
        for (igraph_int_t e = 0; e < igraph_ecount(graph); e++) {
            if (!(removed & BIT(e)) && (seen & BIT(IGRAPH_FROM(graph, e)))) {
                seen |= BIT(IGRAPH_TO(graph, e));
            }
        }
    }
    return seen & BIT(target);
}

static igraph_real_t capacity_of(uint64_t cut, const igraph_vector_t *capacity) {
    igraph_real_t sum = 0;
    for (int e = 0; e < 64; e++) {
        if (cut & BIT(e)) {
            sum += capacity ? VECTOR(*capacity)[e] : 1;
        }
    }
    return sum;
}

/* Each cut must be the set of edges leaving its partition, which contains the
 * source but not the target. If 'value' is not negative, each cut must also
 * have exactly that capacity. */
static void check_partitions(const igraph_t *graph, igraph_int_t s, igraph_int_t t,
                             const igraph_vector_t *capacity, igraph_real_t value,
                             const igraph_vector_int_list_t *cuts,
                             const igraph_vector_int_list_t *parts) {
    IGRAPH_ASSERT(igraph_vector_int_list_size(cuts) == igraph_vector_int_list_size(parts));
    for (igraph_int_t i = 0; i < igraph_vector_int_list_size(cuts); i++) {
        uint64_t cut = to_mask(igraph_vector_int_list_get_ptr(cuts, i));
        uint64_t part = to_mask(igraph_vector_int_list_get_ptr(parts, i));
        uint64_t leaving = 0;
        for (igraph_int_t e = 0; e < igraph_ecount(graph); e++) {
            if ((part & BIT(IGRAPH_FROM(graph, e))) && !(part & BIT(IGRAPH_TO(graph, e)))) {
                leaving |= BIT(e);
            }
        }
        IGRAPH_ASSERT((part & BIT(s)) && !(part & BIT(t)));
        IGRAPH_ASSERT(leaving == cut);
        IGRAPH_ASSERT(value < 0 || capacity_of(cut, capacity) == value);
    }
}

/* Checks igraph_all_st_cuts() and igraph_all_st_mincuts() for one vertex pair.
 * With 'brute_force', the results are compared with all minimal cuts found by
 * trying every edge subset (at most 12 edges). */
static void check_pair(const igraph_t *graph, igraph_int_t s, igraph_int_t t,
                       const igraph_vector_t *capacity, igraph_bool_t brute_force) {
    igraph_int_t m = igraph_ecount(graph);
    igraph_vector_int_list_t cuts, parts;
    igraph_real_t value, flow, best = IGRAPH_INFINITY;
    static uint64_t all[MAX_SETS], min[MAX_SETS], found[MAX_SETS];
    size_t num_all = 0, num_min = 0, num_found;

    igraph_vector_int_list_init(&cuts, 0);
    igraph_vector_int_list_init(&parts, 0);
    IGRAPH_ASSERT(igraph_all_st_cuts(graph, &cuts, &parts, s, t) == IGRAPH_SUCCESS);
    check_partitions(graph, s, t, capacity, -1, &cuts, &parts);
    num_found = sorted_masks(&cuts, found);
    igraph_vector_int_list_destroy(&cuts);
    igraph_vector_int_list_destroy(&parts);

    if (brute_force) {
        IGRAPH_ASSERT(m <= 12);
        for (uint64_t mask = 0; mask < BIT(m); mask++) {
            igraph_bool_t minimal = !reachable(graph, mask, s, t);
            for (igraph_int_t e = 0; minimal && e < m; e++) {
                minimal = !(mask & BIT(e)) || reachable(graph, mask & ~BIT(e), s, t);
            }
            if (minimal) {
                all[num_all++] = mask;
                best = fmin(best, capacity_of(mask, capacity));
            }
        }
        assert_same(found, num_found, all, num_all);
        for (size_t i = 0; i < num_all; i++) {
            if (capacity_of(all[i], capacity) == best) {
                min[num_min++] = all[i];
            }
        }
    }

    igraph_vector_int_list_init(&cuts, 0);
    igraph_vector_int_list_init(&parts, 0);
    IGRAPH_ASSERT(igraph_all_st_mincuts(graph, &value, &cuts, &parts, s, t, capacity) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_st_mincut_value(graph, &flow, s, t, capacity) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(value == flow);
    check_partitions(graph, s, t, capacity, value, &cuts, &parts);
    if (brute_force) {
        IGRAPH_ASSERT(value == best);
        num_found = sorted_masks(&cuts, found);
        assert_same(found, num_found, min, num_min);
    }
    igraph_vector_int_list_destroy(&cuts);
    igraph_vector_int_list_destroy(&parts);
}

static void test_cuts(void) {
    igraph_t g;
    igraph_vector_t capacity;

    /* Graph from igraph/igraph#2532: igraph_all_st_mincuts() once returned
     * cuts that were not of minimum size. Too large for brute force. */
    igraph_small(&g, 20, IGRAPH_DIRECTED,
                 1, 11, 3, 13, 4, 14, 5, 15, 6, 16, 7, 17, 8, 18, 9, 19,
                 0, 8, 18, 0, 0, 9, 19, 0, 2, 1, 11, 2, 14, 1, 11, 4,
                 2, 4, 14, 2, 2, 5, 15, 2, 14, 3, 13, 4, 18, 3, 13, 8,
                 16, 5, 15, 6, 17, 5, 15, 7, 19, 6, 16, 9, 18, 7, 17, 8,
                 -1);
    check_pair(&g, 2, 0, NULL, false);
    igraph_destroy(&g);

    /* Random graphs with unit and with integer capacities. The last third may
     * have parallel edges and self-loops. */
    igraph_rng_seed(igraph_rng_default(), 20260);
    for (int trial = 0; trial < 36; trial++) {
        igraph_int_t n = RNG_INTEGER(3, 6);
        igraph_int_t m = RNG_INTEGER(n - 1, 12);
        igraph_vector_int_t edges;

        igraph_vector_int_init(&edges, 0);
        igraph_vector_init(&capacity, m);
        for (igraph_int_t e = 0; e < m; e++) {
            igraph_vector_int_push_back(&edges, RNG_INTEGER(0, n - 1));
            igraph_vector_int_push_back(&edges, RNG_INTEGER(0, n - 1));
            VECTOR(capacity)[e] = RNG_INTEGER(1, 3);
        }
        igraph_create(&g, &edges, n, IGRAPH_DIRECTED);
        if (trial < 24) {
            igraph_simplify(&g, true, true, NULL);
            igraph_vector_resize(&capacity, igraph_ecount(&g));
        }
        for (igraph_int_t s = 0; s < n; s++) {
            for (igraph_int_t t = 0; t < n; t++) {
                if (s != t && reachable(&g, 0, s, t)) {
                    check_pair(&g, s, t, NULL, true);
                    check_pair(&g, s, t, &capacity, true);
                }
            }
        }
        igraph_vector_destroy(&capacity);
        igraph_vector_int_destroy(&edges);
        igraph_destroy(&g);
    }
}

/* Vertex separators in undirected graphs */

static uint64_t adj[16];

/* Vertices reachable from 'v' without using the vertices in 'removed'. */
static uint64_t component(igraph_int_t n, uint64_t removed, igraph_int_t v) {
    uint64_t seen = BIT(v), old = 0;
    while (seen != old) {
        old = seen;
        for (igraph_int_t u = 0; u < n; u++) {
            if (seen & BIT(u)) {
                seen |= adj[u] & ~removed;
            }
        }
    }
    return seen;
}

/* Does removing 'mask' leave at least two components? If 'full' is set, at
 * least two of them must also be adjacent to every vertex of 'mask'; such a
 * set is exactly a minimal (a, b)-separator for some a and b. */
static igraph_bool_t separates(igraph_int_t n, uint64_t mask, igraph_bool_t full) {
    uint64_t rest = (BIT(n) - 1) & ~mask;
    int count = 0;
    while (rest) {
        igraph_int_t v = 0;
        uint64_t comp, nbrs = 0;
        while (!(rest & BIT(v))) {
            v++;
        }
        comp = component(n, mask, v);
        for (igraph_int_t u = 0; u < n; u++) {
            if (comp & BIT(u)) {
                nbrs |= adj[u];
            }
        }
        count += !full || (nbrs & mask) == mask;
        rest &= ~comp;
    }
    return count >= 2;
}

static void check_separators(const igraph_t *graph) {
    igraph_int_t n = igraph_vcount(graph);
    igraph_vector_int_list_t list;
    static uint64_t expected[MAX_SETS], found[MAX_SETS];
    size_t num_expected = 0, num_found;
    int min_size = 99;

    for (igraph_int_t e = 0; e < igraph_ecount(graph); e++) {
        adj[IGRAPH_FROM(graph, e)] |= BIT(IGRAPH_TO(graph, e));
        adj[IGRAPH_TO(graph, e)] |= BIT(IGRAPH_FROM(graph, e));
    }
    igraph_bool_t connected = component(n, 0, 0) == BIT(n) - 1;

    /* Minimum size separators; none are expected in disconnected or complete graphs. */
    for (uint64_t mask = 0; connected && mask < BIT(n); mask++) {
        if (separates(n, mask, false) && popcount(mask) < min_size) {
            min_size = popcount(mask);
        }
    }
    for (uint64_t mask = 0; connected && mask < BIT(n); mask++) {
        if (separates(n, mask, false) && popcount(mask) == min_size) {
            expected[num_expected++] = mask;
        }
    }
    igraph_vector_int_list_init(&list, 0);
    IGRAPH_ASSERT(igraph_minimum_size_separators(graph, &list) == IGRAPH_SUCCESS);
    num_found = sorted_masks(&list, found);
    assert_same(found, num_found, expected, num_expected);
    igraph_vector_int_list_destroy(&list);

    /* Minimal separators. Only connected graphs are checked, because the
     * handling of the empty separator differs. */
    num_expected = 0;
    for (uint64_t mask = 0; mask < BIT(n); mask++) {
        if (separates(n, mask, true)) {
            expected[num_expected++] = mask;
        }
    }
    igraph_vector_int_list_init(&list, 0);
    IGRAPH_ASSERT(igraph_all_minimal_st_separators(graph, &list) == IGRAPH_SUCCESS);
    num_found = sorted_masks(&list, found);
    if (connected) {
        assert_same(found, num_found, expected, num_expected);
    }
    igraph_vector_int_list_destroy(&list);

    for (igraph_int_t u = 0; u < n; u++) {
        adj[u] = 0;
    }
}

static void test_separators(void) {
    igraph_t g;

    /* Graph from igraph/igraph#2521: igraph_minimum_size_separators() once
     * returned sets that were not separators. */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 8, 0, 9, 1, 2, 1, 4, 2, 4, 2, 5, 3, 4, 3, 8, 5, 6, 5, 7, 6, 9, 7, 8,
                 -1);
    check_separators(&g);
    igraph_destroy(&g);

    /* Random graphs, connected and not. */
    igraph_rng_seed(igraph_rng_default(), 20261);
    for (int trial = 0; trial < 60; trial++) {
        igraph_int_t n = RNG_INTEGER(4, 9);
        igraph_vector_int_t edges;

        igraph_vector_int_init(&edges, 0);
        for (igraph_int_t u = 0; u < n; u++) {
            for (igraph_int_t v = u + 1; v < n; v++) {
                if (RNG_INTEGER(0, 9) < 4) {
                    igraph_vector_int_push_back(&edges, u);
                    igraph_vector_int_push_back(&edges, v);
                }
            }
        }
        igraph_create(&g, &edges, n, IGRAPH_UNDIRECTED);
        check_separators(&g);
        igraph_destroy(&g);
        igraph_vector_int_destroy(&edges);
    }
}

int main(void) {
    test_cuts();
    VERIFY_FINALLY_STACK();

    test_separators();
    VERIFY_FINALLY_STACK();

    return 0;
}
