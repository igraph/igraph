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

/* Tests that igraph_betweenness_subset() and igraph_edge_betweenness_subset()
 * count each qualifying pair of vertices exactly once, also when the source
 * and target sets overlap only partially. See
 * https://github.com/igraph/igraph/issues/2359 */

#include <igraph.h>
#include "test_utilities.h"

/* Brute-force reference implementation. It lists all shortest paths from each
 * source with igraph_get_all_shortest_paths(), and for each qualifying pair
 * (s, t), adds 1 / (number of shortest s-t paths) to every inner vertex and
 * every edge of each shortest s-t path.
 *
 * When paths are undirected, the pair {s, t} qualifies if one of its vertices
 * is a source and the other is a target. Each such pair is counted once. */
void brute_force_subset_betweenness(
        const igraph_t *graph, const igraph_vector_t *weights,
        const igraph_vector_bool_t *is_source, const igraph_vector_bool_t *is_target,
        igraph_bool_t directed,
        igraph_vector_t *vres, igraph_vector_t *eres) {

    igraph_int_t vcount = igraph_vcount(graph);
    igraph_neimode_t mode = (directed && igraph_is_directed(graph)) ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_int_list_t vpaths, epaths;
    igraph_vector_int_t path_count;

    igraph_vector_int_list_init(&vpaths, 0);
    igraph_vector_int_list_init(&epaths, 0);
    igraph_vector_int_init(&path_count, vcount);

    igraph_vector_resize(vres, vcount);
    igraph_vector_null(vres);
    igraph_vector_resize(eres, igraph_ecount(graph));
    igraph_vector_null(eres);

    for (igraph_int_t s = 0; s < vcount; s++) {
        igraph_int_t npaths;

        if (!VECTOR(*is_source)[s]) {
            continue;
        }

        igraph_get_all_shortest_paths(graph, weights, &vpaths, &epaths, NULL,
                                      s, igraph_vss_all(), mode);
        npaths = igraph_vector_int_list_size(&vpaths);

        /* Count the shortest paths to each vertex. */
        igraph_vector_int_null(&path_count);
        for (igraph_int_t i = 0; i < npaths; i++) {
            const igraph_vector_int_t *vpath = igraph_vector_int_list_get_ptr(&vpaths, i);
            VECTOR(path_count)[igraph_vector_int_tail(vpath)]++;
        }

        for (igraph_int_t i = 0; i < npaths; i++) {
            const igraph_vector_int_t *vpath = igraph_vector_int_list_get_ptr(&vpaths, i);
            const igraph_vector_int_t *epath = igraph_vector_int_list_get_ptr(&epaths, i);
            igraph_int_t len = igraph_vector_int_size(vpath);
            igraph_int_t t = igraph_vector_int_tail(vpath);
            igraph_real_t frac;

            if (t == s || !VECTOR(*is_target)[t]) {
                continue;
            }
            /* For undirected paths, a pair whose vertices are both sources
             * and both targets is reached from both of its ends. Count it
             * only from its smaller end. */
            if (mode == IGRAPH_ALL && VECTOR(*is_source)[t] && VECTOR(*is_target)[s] && t < s) {
                continue;
            }

            frac = 1.0 / VECTOR(path_count)[t];
            for (igraph_int_t j = 1; j < len - 1; j++) {
                VECTOR(*vres)[VECTOR(*vpath)[j]] += frac;
            }
            for (igraph_int_t j = 0; j < len - 1; j++) {
                VECTOR(*eres)[VECTOR(*epath)[j]] += frac;
            }
        }
    }

    igraph_vector_int_destroy(&path_count);
    igraph_vector_int_list_destroy(&epaths);
    igraph_vector_int_list_destroy(&vpaths);
}

/* Converts a boolean membership vector to a vertex selector. */
void bool_to_vs(const igraph_vector_bool_t *mask, igraph_vector_int_t *ids, igraph_vs_t *vs) {
    igraph_int_t n = igraph_vector_bool_size(mask);
    igraph_vector_int_clear(ids);
    for (igraph_int_t i = 0; i < n; i++) {
        if (VECTOR(*mask)[i]) {
            igraph_vector_int_push_back(ids, i);
        }
    }
    igraph_vs_vector(vs, ids);
}

/* Compares both subset betweenness functions with the brute-force reference
 * on random graphs and random source and target sets. */
void random_test(igraph_bool_t directed_graph, igraph_bool_t directed, igraph_bool_t weighted) {
    igraph_t g;
    igraph_vector_t weights, vres, eres, vres_bf, eres_bf;
    igraph_vector_bool_t is_source, is_target;
    igraph_vector_int_t source_ids, target_ids;
    igraph_vs_t sources, targets;

    igraph_vector_init(&weights, 0);
    igraph_vector_init(&vres, 0);
    igraph_vector_init(&eres, 0);
    igraph_vector_init(&vres_bf, 0);
    igraph_vector_init(&eres_bf, 0);
    igraph_vector_bool_init(&is_source, 0);
    igraph_vector_bool_init(&is_target, 0);
    igraph_vector_int_init(&source_ids, 0);
    igraph_vector_int_init(&target_ids, 0);

    for (int iter = 0; iter < 100; iter++) {
        igraph_int_t n = RNG_INTEGER(2, 10);
        igraph_int_t m = RNG_INTEGER(0, 2 * n);

        igraph_erdos_renyi_game_gnm(&g, n, m, directed_graph, IGRAPH_MULTI_SW, IGRAPH_EDGE_UNLABELED);

        /* Small integer weights make ties between path lengths exact. */
        igraph_vector_resize(&weights, m);
        for (igraph_int_t i = 0; i < m; i++) {
            VECTOR(weights)[i] = RNG_INTEGER(1, 3);
        }

        igraph_vector_bool_resize(&is_source, n);
        igraph_vector_bool_resize(&is_target, n);
        for (igraph_int_t i = 0; i < n; i++) {
            VECTOR(is_source)[i] = RNG_BOOL();
            VECTOR(is_target)[i] = RNG_BOOL();
        }
        bool_to_vs(&is_source, &source_ids, &sources);
        bool_to_vs(&is_target, &target_ids, &targets);

        igraph_betweenness_subset(&g, weighted ? &weights : NULL, &vres,
                                  sources, targets, igraph_vss_all(), directed, false);
        igraph_edge_betweenness_subset(&g, weighted ? &weights : NULL, &eres,
                                       sources, targets, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                       directed, false);
        brute_force_subset_betweenness(&g, weighted ? &weights : NULL,
                                       &is_source, &is_target, directed,
                                       &vres_bf, &eres_bf);

        if (!igraph_vector_all_almost_e(&vres, &vres_bf, 1e-10) ||
            !igraph_vector_all_almost_e(&eres, &eres_bf, 1e-10)) {
            printf("Mismatch with brute-force result.\n");
            print_graph(&g);
            if (weighted) {
                printf("Weights: ");
                print_vector(&weights);
            }
            printf("Sources: ");
            print_vector_int(&source_ids);
            printf("Targets: ");
            print_vector_int(&target_ids);
            printf("Directed paths: %s\n", directed ? "true" : "false");
            printf("Vertex betweenness:   ");
            print_vector(&vres);
            printf("Brute-force result:   ");
            print_vector(&vres_bf);
            printf("Edge betweenness:     ");
            print_vector(&eres);
            printf("Brute-force result:   ");
            print_vector(&eres_bf);
            IGRAPH_FATAL("Subset betweenness does not match brute-force result.");
        }

        igraph_vs_destroy(&sources);
        igraph_vs_destroy(&targets);
        igraph_destroy(&g);
    }

    igraph_vector_int_destroy(&target_ids);
    igraph_vector_int_destroy(&source_ids);
    igraph_vector_bool_destroy(&is_target);
    igraph_vector_bool_destroy(&is_source);
    igraph_vector_destroy(&eres_bf);
    igraph_vector_destroy(&vres_bf);
    igraph_vector_destroy(&eres);
    igraph_vector_destroy(&vres);
    igraph_vector_destroy(&weights);
}

void print_both(const igraph_t *g, igraph_vs_t sources, igraph_vs_t targets, igraph_bool_t directed) {
    igraph_vector_t res;
    igraph_vector_init(&res, 0);

    igraph_betweenness_subset(g, NULL, &res, sources, targets, igraph_vss_all(), directed, false);
    printf("Vertices: ");
    print_vector(&res);

    igraph_edge_betweenness_subset(g, NULL, &res, sources, targets,
                                   igraph_ess_all(IGRAPH_EDGEORDER_ID), directed, false);
    printf("Edges:    ");
    print_vector(&res);

    igraph_vector_destroy(&res);
}

int main(void) {
    igraph_t g;
    igraph_vector_int_t v13;

    igraph_rng_seed(igraph_rng_default(), 42);

    igraph_vector_int_init_int(&v13, 2, 1, 3);

    /* Path 0-1-2 and isolated vertex 3, sources {0, 1}, targets {2, 3}.
     * The qualifying pairs are {0,2}, {0,3}, {1,2} and {1,3}. Only {0,2} has
     * an inner vertex, 1. Edge 0-1 is used by {0,2}; edge 1-2 is used by
     * {0,2} and {1,2}. Expected: vertices ( 0 1 0 0 ), edges ( 1 2 ). */
    printf("Path of length 2 with an isolated vertex, disjoint sources and targets:\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, -1);
    print_both(&g, igraph_vss_range(0, 2), igraph_vss_range(2, 4), IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    /* Path 0-1-2-3, sources {0, 1}, targets {1, 3}. The qualifying pairs are
     * {0,1}, {0,3} and {1,3}; {1,1} is not a pair. Vertex 1 is inside the path
     * of {0,3}; vertex 2 is inside the paths of {0,3} and {1,3}. Edge 0-1 is
     * used by {0,1} and {0,3}; edges 1-2 and 2-3 by {0,3} and {1,3}.
     * Expected: vertices ( 0 1 2 0 ), edges ( 2 2 2 ). */
    printf("\nPath of length 3, partially overlapping sources and targets:\n");
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, -1);
    print_both(&g, igraph_vss_range(0, 2), igraph_vss_vector(&v13), IGRAPH_UNDIRECTED);
    igraph_destroy(&g);

    /* The same directed path, with directed paths ignored, must give the same
     * result as the undirected path. With directed paths considered, the
     * qualifying ordered pairs are (0,1), (0,3) and (1,3), which happen to
     * give the same result. After swapping sources and targets, the
     * qualifying ordered pairs are (1,0), (3,0) and (3,1), none of which is
     * connected by a directed path. */
    printf("\nDirected path of length 3, directed paths ignored:\n");
    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, -1);
    print_both(&g, igraph_vss_range(0, 2), igraph_vss_vector(&v13), IGRAPH_UNDIRECTED);
    printf("\nDirected path of length 3, directed paths considered:\n");
    print_both(&g, igraph_vss_range(0, 2), igraph_vss_vector(&v13), IGRAPH_DIRECTED);
    printf("\nDirected path of length 3, sources and targets swapped:\n");
    print_both(&g, igraph_vss_vector(&v13), igraph_vss_range(0, 2), IGRAPH_DIRECTED);
    igraph_destroy(&g);

    igraph_vector_int_destroy(&v13);

    random_test(IGRAPH_UNDIRECTED, IGRAPH_UNDIRECTED, false);
    random_test(IGRAPH_UNDIRECTED, IGRAPH_UNDIRECTED, true);
    random_test(IGRAPH_DIRECTED, IGRAPH_UNDIRECTED, false);
    random_test(IGRAPH_DIRECTED, IGRAPH_UNDIRECTED, true);
    random_test(IGRAPH_DIRECTED, IGRAPH_DIRECTED, false);
    random_test(IGRAPH_DIRECTED, IGRAPH_DIRECTED, true);

    VERIFY_FINALLY_STACK();

    return 0;
}
