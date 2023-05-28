//igraph/src/centrality/betweenness.c 500
//
#include "igraph_vector_list.h"
#include "test_utilities.h"


 // TODO also do \function igraph_betweenness_subset
 // TODO do the cutoff

//TODO  * The betweenness centrality of a vertex is the number of geodesics
// * going through it. If there are more than one geodesic between two
// * vertices, the value of these geodesics are weighted by one over the
// * number of geodesics.
igraph_error_t igraph_direct_betweenness_cutoff(
        const igraph_t *graph, igraph_vector_t *res,
        const igraph_vs_t vids, igraph_bool_t directed,
        const igraph_vector_t *weights, igraph_real_t cutoff) {
    igraph_vector_int_list_t vertices;
    igraph_vit_t vit;
    igraph_vs_t vs_rest;
    igraph_neimode_t mode;
    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_selected_nodes;
    igraph_integer_t i_res;

    IGRAPH_CHECK(igraph_vs_size(graph, &vids, &no_of_selected_nodes));

    if (directed) {
        mode = IGRAPH_OUT;
    } else {
        mode = IGRAPH_ALL;
    }

    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);

    igraph_vector_resize(res, no_of_selected_nodes);
    igraph_vector_null(res);

    IGRAPH_VECTOR_INT_LIST_INIT_FINALLY(&vertices, 0);

    for (igraph_integer_t i_all = 0; i_all < no_of_nodes; i_all ++) {
        igraph_vs_range(&vs_rest, i_all, no_of_nodes);
        if (cutoff == -1) {
            if (!weights) {
                /* when there is more than one shortest path between two vertices,
                 * all of them will be returned. */
                /*   The vectors are
                 *   ordered according to their target vertex: first the shortest paths to
                 *   vertex 0, then to vertex 1, etc.
                 *   */
                igraph_get_all_shortest_paths(graph, &vertices, NULL, NULL, i_all, vs_rest, mode);
            } else {
                igraph_get_all_shortest_paths_dijkstra(graph, &vertices, NULL, NULL, i_all, vs_rest, weights, mode);
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
    }

    igraph_vector_int_list_destroy(&vertices);
    igraph_vit_destroy(&vit);

    return IGRAPH_SUCCESS;
}

void compare(igraph_t *graph) {
    igraph_vector_t test_res, res;
    igraph_vector_init(&res, 0);
    igraph_vector_init(&test_res, 0);
    igraph_betweenness(graph, &res,
                       igraph_vss_all(), true,
                       NULL);
    //igraph_vector_print(&res);

    igraph_direct_betweenness_cutoff(graph, &test_res,
                       igraph_vss_all(), true,
                       NULL, -1);
    //igraph_vector_print(&test_res);

    for (igraph_integer_t i = 0; i < igraph_vector_size(&res); i++) {
        if (igraph_cmp_epsilon(VECTOR(res)[i], VECTOR(test_res)[i], 0.01)) {
            printf("index: %ld, values: %f, %f", i, VECTOR(res)[i], VECTOR(test_res)[i]);
        }
    }

    igraph_vector_destroy(&res);
    igraph_vector_destroy(&test_res);
}


int main() {
    igraph_t graph;

    igraph_small(&graph, 5, IGRAPH_UNDIRECTED, 0,1, 0,2, -1);
    compare(&graph);
    igraph_destroy(&graph);

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
