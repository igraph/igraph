//igraph/src/centrality/betweenness.c 500
//
#include "igraph_error.h"
#include "igraph_vector_list.h"
#include "test_utilities.h"


 // TODO also do \function igraph_betweenness_subset
 // TODO do the cutoff (lot of work)

igraph_error_t igraph_direct_betweenness_cutoff(
        const igraph_t *graph, igraph_vector_t *res,
        const igraph_vs_t vids, igraph_bool_t directed,
        const igraph_vs_t sources, const igraph_vs_t targets,
        const igraph_vector_t *weights, igraph_real_t cutoff) {
    igraph_vector_int_list_t vertices;
    igraph_vit_t vit, vit_sources, vit_targets, vit_sources2;
    igraph_neimode_t mode;
    igraph_integer_t no_of_selected_nodes;
    igraph_integer_t i_res;
    igraph_vs_t targets_rest;
    igraph_vector_int_t targets_rest_vec;
    igraph_integer_t no_of_sources;
    igraph_integer_t no_of_targets;

    IGRAPH_CHECK(igraph_vs_size(graph, &vids, &no_of_selected_nodes));
    IGRAPH_CHECK(igraph_vs_size(graph, &sources, &no_of_sources));

    if (!igraph_is_directed(graph)) {
        directed = false;
    }

    if (directed) {
        mode = IGRAPH_OUT;
    } else {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit_sources));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_sources);

    IGRAPH_CHECK(igraph_vit_create(graph, sources, &vit_sources2));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_sources2);

    IGRAPH_CHECK(igraph_vit_create(graph, targets, &vit_targets));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit_targets);

    igraph_vector_resize(res, no_of_selected_nodes);
    igraph_vector_null(res);

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vertices, 0);

    if (directed) {
        //FINALLY
        igraph_vs_copy(&targets_rest, &targets);
    }

    for (IGRAPH_VIT_RESET(vit_sources); !IGRAPH_VIT_END(vit_sources); IGRAPH_VIT_NEXT(vit_sources)) {
        if (!directed) {
            igraph_bool_t source_in_targets = false;
            for (IGRAPH_VIT_RESET(vit_targets); !IGRAPH_VIT_END(vit_targets); IGRAPH_VIT_NEXT(vit_targets)) {
                if (IGRAPH_VIT_GET(vit_targets) == IGRAPH_VIT_GET(vit_sources)) {
                    source_in_targets = true;
                    break;
                }
            }
            if (source_in_targets) {
                IGRAPH_VECTOR_INT_INIT_FINALLY(&targets_rest_vec, 0);
                for (IGRAPH_VIT_RESET(vit_targets); !IGRAPH_VIT_END(vit_targets); IGRAPH_VIT_NEXT(vit_targets)) {
                    igraph_bool_t target_in_sources = false;
                    if (IGRAPH_VIT_GET(vit_targets) < IGRAPH_VIT_GET(vit_sources)) {
                        for (IGRAPH_VIT_RESET(vit_sources2); !IGRAPH_VIT_END(vit_sources2); IGRAPH_VIT_NEXT(vit_sources2)) {
                            if (IGRAPH_VIT_GET(vit_targets) == IGRAPH_VIT_GET(vit_sources2)) {
                                target_in_sources = true;
                                break;
                            }
                        }
                    }
                    if (!target_in_sources) {
                        igraph_vector_int_push_back(&targets_rest_vec, IGRAPH_VIT_GET(vit_targets));
                    }
                }
                igraph_vs_vector_copy(&targets_rest, &targets_rest_vec);
                igraph_vector_int_destroy(&targets_rest_vec);
                IGRAPH_FINALLY_CLEAN(1);
            } else {
                igraph_vs_copy(&targets_rest, &targets);
            }
        }
        /*
        printf("targets: ");
        igraph_vector_int_t tmp;
        igraph_vector_int_init(&tmp, 0);
        igraph_vs_as_vector(graph, targets_rest, &tmp);
        igraph_vector_int_print(&tmp);
        igraph_vector_int_destroy(&tmp);
        */
        IGRAPH_CHECK(igraph_vs_size(graph, &targets_rest, &no_of_targets));
        if (cutoff == -1) {
            if (!weights) {
                /* when there is more than one shortest path between two vertices,
                 * all of them will be returned. */
                /*   The vectors are
                 *   ordered according to their target vertex: first the shortest paths to
                 *   vertex 0, then to vertex 1, etc.
                 *   */
                igraph_get_all_shortest_paths(graph, &vertices, NULL, NULL, IGRAPH_VIT_GET(vit_sources), targets_rest, mode);
            } else {
                igraph_get_all_shortest_paths_dijkstra(graph, &vertices, NULL, NULL, IGRAPH_VIT_GET(vit_sources), targets_rest, weights, mode);
            }
        } else {
            IGRAPH_ERROR("Cutoffs are not implemented, set it to -1", IGRAPH_EINVAL);
        }
        igraph_integer_t no_of_paths = igraph_vector_int_list_size(&vertices);
        for (igraph_integer_t i = 0; i < no_of_paths; i++) {
            igraph_vector_int_t *path = igraph_vector_int_list_get_ptr(&vertices, i);
            igraph_integer_t no_of_path_nodes = igraph_vector_int_size(path);
            igraph_integer_t no_of_geodesics_with_these_endpoints = 0;
            for (igraph_integer_t k = 0; k < no_of_paths; k++) {
                igraph_vector_int_t *path2 = igraph_vector_int_list_get_ptr(&vertices, k);
                igraph_integer_t no_of_path2_nodes = igraph_vector_int_size(path2);
                if (VECTOR(*path)[0] == VECTOR(*path2)[0] && VECTOR(*path)[no_of_path_nodes - 1] == VECTOR(*path2)[no_of_path2_nodes -1]) {
                    no_of_geodesics_with_these_endpoints++;
                }
            }
            for (igraph_integer_t j = 1; j < no_of_path_nodes - 1; j++) { //don't count the endpoints
                for (IGRAPH_VIT_RESET(vit), i_res = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i_res ++) {
                    if (VECTOR(*path)[j] == IGRAPH_VIT_GET(vit)) {
                        VECTOR(*res)[i_res] += 1.0 / no_of_geodesics_with_these_endpoints;
                    }
                }
            }
        }
        if (!directed) {
            igraph_vs_destroy(&targets_rest);
        }
    }

    igraph_vector_int_list_destroy(&vertices);
    igraph_vit_destroy(&vit);
    igraph_vit_destroy(&vit_sources);
    igraph_vit_destroy(&vit_sources2);
    igraph_vit_destroy(&vit_targets);

    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

void compare_directed_subset(igraph_t *graph, igraph_bool_t directed, igraph_vs_t sources, igraph_vs_t targets) {
    igraph_vector_t test_res, res;
    igraph_vector_init(&res, 0);
    igraph_vector_init(&test_res, 0);
    igraph_betweenness_subset(graph, &res,
                       igraph_vss_all(), directed,
                       sources, targets,
                       NULL);
    //igraph_vector_print(&res);

    igraph_direct_betweenness_cutoff(graph, &test_res,
                       igraph_vss_all(), directed,
                       sources, targets,
                       NULL, -1);
    //igraph_vector_print(&test_res);

    for (igraph_integer_t i = 0; i < igraph_vector_size(&res); i++) {
        if (igraph_cmp_epsilon(VECTOR(res)[i], VECTOR(test_res)[i], 0.01)) {
            printf("index: %ld, value igraph: %f, value direct: %f\n", i, VECTOR(res)[i], VECTOR(test_res)[i]);
        }
    }

    igraph_vector_destroy(&res);
    igraph_vector_destroy(&test_res);
}

void compare_directed(igraph_t *graph, igraph_bool_t directed) {
    printf("all sources and targets:\n");
    compare_directed_subset(graph, directed, igraph_vss_all(), igraph_vss_all());
    printf("two sources, two different targets:\n");
    compare_directed_subset(graph, directed, igraph_vss_range(0, 2), igraph_vss_range(2,4));
    printf("two sources, same two targets:\n");
    compare_directed_subset(graph, directed, igraph_vss_range(0,2), igraph_vss_range(0,2));
    printf("one source, one different target:\n");
    compare_directed_subset(graph, directed, igraph_vss_range(0, 1), igraph_vss_range(2,3));
}

void compare(igraph_t *graph) {
    printf("directed test:\n");
    compare_directed(graph, true);
    printf("undirected test:\n");
    compare_directed(graph, false);
}

int main() {
    igraph_t graph;
    printf("\nundirected line\n");
    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0,1, 1,2, -1);
    compare(&graph);
    igraph_destroy(&graph);

    printf("\nundirected star\n");
    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0,4, 1,4, 2,4, 3,4, -1);
    compare(&graph);
    igraph_destroy(&graph);

    printf("five vertices, 0 connected to 1 and 2, undirected\n");
    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0,1, 0,2, -1);
    compare(&graph);
    igraph_destroy(&graph);

    printf("\nbarabasi game, undirected:\n");
    igraph_barabasi_game(/* graph= */    &graph,
                                         /* n= */        100,
                                         /* power= */    1,
                                         /* m= */        3,
                                         /* outseq= */   0,
                                         /* outpref= */  0,
                                         /* A= */        1,
                                         /* directed= */ 0,
                                         /* algo= */     IGRAPH_BARABASI_BAG,
                                         /* start_from= */ 0);

    compare(&graph);
    igraph_destroy(&graph);

    printf("\ndirected graph:\n");
    igraph_erdos_renyi_game_gnp(&graph, 100, 0.3, IGRAPH_DIRECTED, /*loops*/true);
    compare(&graph);
    igraph_destroy(&graph);

    printf("\nundirected graph:\n");
    igraph_erdos_renyi_game_gnp(&graph, 100, 0.3, IGRAPH_UNDIRECTED, /*loops*/true);
    compare(&graph);
    igraph_destroy(&graph);
}
