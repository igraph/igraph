/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

typedef struct igraph_lazy_adjlist2_t {
    const igraph_t *graph;
    igraph_integer_t length;
    igraph_integer_t data_length;
    igraph_vector_int_t *adjs;
    igraph_integer_t *data;
    igraph_integer_t next_data;
    igraph_neimode_t mode;
    igraph_loops_t loops;
    igraph_multiple_t multiple;
} igraph_lazy_adjlist2_t;
igraph_error_t igraph_lazy_adjlist2_init(const igraph_t *graph,
                             igraph_lazy_adjlist2_t *al,
                             igraph_neimode_t mode,
                             igraph_loops_t loops,
                             igraph_multiple_t multiple) {
    if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
        IGRAPH_ERROR("Cannor create lazy adjacency list view", IGRAPH_EINVMODE);
    }

    if (!igraph_is_directed(graph)) {
        mode = IGRAPH_ALL;
    }

    al->mode = mode;
    al->loops = loops;
    al->multiple = multiple;
    al->graph = graph;

    al->length = igraph_vcount(graph);
    al->data_length = igraph_ecount(graph) * 2;
    al->adjs = IGRAPH_CALLOC(al->length, igraph_vector_int_t);
    al->data = IGRAPH_CALLOC(al->data_length, igraph_integer_t);
    al->next_data = 0;

    if (al->adjs == NULL || al->data == NULL) {
        IGRAPH_ERROR("Cannot create lazy adjacency list view", IGRAPH_ENOMEM);
    }

    return IGRAPH_SUCCESS;
}

void igraph_lazy_adjlist2_destroy(igraph_lazy_adjlist2_t *al) {
    IGRAPH_FREE(al->adjs);
    IGRAPH_FREE(al->data);
}

igraph_vector_int_t *igraph_i_lazy_adjlist2_get_real(igraph_lazy_adjlist2_t *al,
        igraph_integer_t pno) {
    igraph_integer_t no = pno;
    igraph_error_t ret;

    if (al->adjs[no].stor_begin == NULL) {
        igraph_vector_int_view(&al->adjs[no], al->data + al->next_data, al->data_length - al->next_data);
        /* igraph_neighbors calls a resize. Since we're handling our own
         * vectors, we can't have our stor_begin be realloced. A realloc should only
         * happen when the vector is too small to handle the neighbors, which
         * should never happen, because the stor_end is set to one past the
         * end of the latest integer.
         */
        ret = igraph_neighbors(al->graph, &al->adjs[no], no, al->mode);
        if (ret != IGRAPH_SUCCESS) {
            igraph_error("", IGRAPH_FILE_BASENAME, __LINE__, ret);
            return NULL;
        }
        al->next_data += igraph_vector_int_size(&(al->adjs[no]));
    }

    return &al->adjs[no];
}

igraph_error_t igraph_neighborhood_adjl2(const igraph_t *graph, igraph_vector_int_list_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t *neis;
    igraph_vector_int_t tmp;
    igraph_lazy_adjlist2_t adjlist2;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_int_list_clear(res);

    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        IGRAPH_CHECK(igraph_lazy_adjlist2_init(graph, &adjlist2, mode, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    } else {
        IGRAPH_CHECK(igraph_lazy_adjlist2_init(graph, &adjlist2, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    }
    IGRAPH_FINALLY(igraph_lazy_adjlist2_destroy, &adjlist2);
    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_int_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, node));
        }
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            neis = igraph_i_lazy_adjlist2_get_real(&adjlist2, actnode);
            //IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, &tmp));
    }

    igraph_lazy_adjlist2_destroy(&adjlist2);
    igraph_vector_int_destroy(&tmp);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_neighborhood_adjl(const igraph_t *graph, igraph_vector_int_list_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t *neis;
    igraph_vector_int_t tmp;
    igraph_lazy_adjlist_t adjlist;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_int_list_clear(res);

    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    } else {
        IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    }
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_int_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, node));
        }
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            neis = igraph_lazy_adjlist_get(&adjlist, actnode);
            //IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, &tmp));
    }

    igraph_lazy_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&tmp);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

igraph_error_t igraph_neighborhood_adj(const igraph_t *graph, igraph_vector_int_list_t *res,
                        igraph_vs_t vids, igraph_integer_t order,
                        igraph_neimode_t mode, igraph_integer_t mindist) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_dqueue_int_t q;
    igraph_vit_t vit;
    igraph_integer_t i, j;
    igraph_integer_t *added;
    igraph_vector_int_t *neis;
    igraph_vector_int_t tmp;
    igraph_adjlist_t adjlist;

    if (order < 0) {
        IGRAPH_ERROR("Negative order in neighborhood size", IGRAPH_EINVAL);
    }

    if (mindist < 0 || mindist > order) {
        IGRAPH_ERROR("Minimum distance should be between zero and order",
                     IGRAPH_EINVAL);
    }

    added = IGRAPH_CALLOC(no_of_nodes, igraph_integer_t);
    if (added == 0) {
        IGRAPH_ERROR("Cannot calculate neighborhood size", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);
    IGRAPH_DQUEUE_INT_INIT_FINALLY(&q, 100);
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_vector_int_list_reserve(res, IGRAPH_VIT_SIZE(vit)));
    igraph_vector_int_list_clear(res);

    if (!igraph_is_directed(graph) || mode == IGRAPH_ALL) {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    }
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    for (i = 0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
        igraph_integer_t node = IGRAPH_VIT_GET(vit);
        added[node] = i + 1;
        igraph_vector_int_clear(&tmp);
        if (mindist == 0) {
            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, node));
        }
        if (order > 0) {
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, node));
            IGRAPH_CHECK(igraph_dqueue_int_push(&q, 0));
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t actnode = igraph_dqueue_int_pop(&q);
            igraph_integer_t actdist = igraph_dqueue_int_pop(&q);
            igraph_integer_t n;
            neis = igraph_adjlist_get(&adjlist, actnode);
            //IGRAPH_CHECK(igraph_neighbors(graph, &neis, actnode, mode));
            n = igraph_vector_int_size(neis);

            if (actdist < order - 1) {
                /* we add them to the q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, nei));
                        IGRAPH_CHECK(igraph_dqueue_int_push(&q, actdist + 1));
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            } else {
                /* we just count them but don't add them to q */
                for (j = 0; j < n; j++) {
                    igraph_integer_t nei = VECTOR(*neis)[j];
                    if (added[nei] != i + 1) {
                        added[nei] = i + 1;
                        if (actdist + 1 >= mindist) {
                            IGRAPH_CHECK(igraph_vector_int_push_back(&tmp, nei));
                        }
                    }
                }
            }

        } /* while q not empty */

        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, &tmp));
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&tmp);
    igraph_vit_destroy(&vit);
    igraph_dqueue_int_destroy(&q);
    IGRAPH_FREE(added);
    IGRAPH_FINALLY_CLEAN(5);

    return IGRAPH_SUCCESS;
}

void    do_benchmark(igraph_t *g, igraph_vs_t vids, igraph_integer_t repeat)
{
    igraph_vector_int_list_t result_orig;
    igraph_vector_int_list_t result_adj;
    igraph_vector_int_list_t result_adjl;
    igraph_vector_int_list_t result_adjl2;

    igraph_vector_int_list_init(&result_orig, 0);
    igraph_vector_int_list_init(&result_adj, 0);
    igraph_vector_int_list_init(&result_adjl, 0);
    igraph_vector_int_list_init(&result_adjl2, 0);

    for (igraph_integer_t order = 1; order <= 3; order++) {
        printf("order %" IGRAPH_PRId ":\n", order);
        BENCH("Original function:",
                REPEAT(igraph_neighborhood(g, &result_orig, vids, order,
                        /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));

        BENCH("Using a lazy adjlist:",
                REPEAT(igraph_neighborhood_adjl(g, &result_adjl, vids, order,
                        /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));

        BENCH("Using lazy adjlist 2:",
                REPEAT(igraph_neighborhood_adjl2(g, &result_adjl2, vids, order,
                        /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));

        BENCH("Using adjlist:",
                REPEAT(igraph_neighborhood_adj(g, &result_adj, vids, order,
                        /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));
        for (igraph_integer_t i = 0; i <= order; i++) {
            IGRAPH_ASSERT(!igraph_vector_int_lex_cmp(&VECTOR(result_orig)[i], &VECTOR(result_adj)[i]));
            IGRAPH_ASSERT(!igraph_vector_int_lex_cmp(&VECTOR(result_adj)[i], &VECTOR(result_adjl)[i]));
            IGRAPH_ASSERT(!igraph_vector_int_lex_cmp(&VECTOR(result_adjl)[i], &VECTOR(result_adjl2)[i]));
        }
    }
    igraph_vector_int_list_destroy(&result_orig);
    igraph_vector_int_list_destroy(&result_adj);
    igraph_vector_int_list_destroy(&result_adjl);
    igraph_vector_int_list_destroy(&result_adjl2);
}

int main(void) {
    igraph_t g_full, g_ring, g_er;
    igraph_vs_t vids_all, vids_50, vids_5000, vids_200;

    igraph_vs_all(&vids_all);
    igraph_vs_range(&vids_50, 0, 50);
    igraph_vs_range(&vids_5000, 0, 5000);
    igraph_vs_range(&vids_200, 0, 200);

    igraph_full(&g_full, 500, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_ring(&g_ring, 50000, IGRAPH_UNDIRECTED, /* mutual */ 0, /*circular*/ 0);
    igraph_erdos_renyi_game(&g_er, IGRAPH_ERDOS_RENYI_GNM, 2000, 20000, IGRAPH_UNDIRECTED, /*loops*/ 0);

    printf("Select all vertices:\n\n");
    printf("Full graph:\n");
    do_benchmark(&g_full, vids_all, 10);
    printf("\nRing graph:\n");
    do_benchmark(&g_ring, vids_all, 40);
    printf("\nRandom graph:\n");
    do_benchmark(&g_er, vids_all, 15);

    printf("\n\nSelect 10%% of vertices:\n\n");
    printf("Full graph:\n");
    do_benchmark(&g_full, vids_50, 100);
    printf("\nRing graph:\n");
    do_benchmark(&g_ring, vids_5000, 400);
    printf("\nRandom graph:\n");
    do_benchmark(&g_er, vids_200, 150);

    igraph_destroy(&g_full);
    igraph_destroy(&g_ring);
    igraph_destroy(&g_er);
    igraph_vs_destroy(&vids_all);
    igraph_vs_destroy(&vids_50);
    igraph_vs_destroy(&vids_5000);
    igraph_vs_destroy(&vids_200);
}
