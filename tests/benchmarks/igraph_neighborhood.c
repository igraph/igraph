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

void    do_benchmark(igraph_t *g, igraph_vs_t vids, igraph_integer_t order, igraph_integer_t repeat)
{
    igraph_vector_int_list_t result_orig;
    igraph_vector_int_list_t result_adj;

    igraph_vector_int_list_init(&result_orig, 0);
    igraph_vector_int_list_init(&result_adj, 0);
    BENCH("Original function:",
          REPEAT(igraph_neighborhood(g, &result_orig, vids, order,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));

    BENCH("Using adjlist:",
          REPEAT(igraph_neighborhood_adj(g, &result_adj, vids, order,
                  /*mode*/ IGRAPH_ALL, /*mindist*/ 0), repeat));
    for (igraph_integer_t i = 0; i <= order; i++) {
        IGRAPH_ASSERT(!igraph_vector_int_lex_cmp(&VECTOR(result_orig)[i], &VECTOR(result_adj)[i]));
    }
    igraph_vector_int_list_destroy(&result_orig);
    igraph_vector_int_list_destroy(&result_adj);
}

int main() {
    igraph_t g_full, g_ring, g_er;
    igraph_vs_t vids_all;

    igraph_vs_all(&vids_all);

    igraph_full(&g_full, 500, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_ring(&g_ring, 50000, IGRAPH_UNDIRECTED, /* mutual */ 0, /*circular*/ 0);
    igraph_erdos_renyi_game(&g_er, IGRAPH_ERDOS_RENYI_GNM, 2000, 20000, IGRAPH_UNDIRECTED, /*loops*/ 0);

    printf("Full graph:\n");
    do_benchmark(&g_full, vids_all, 2, 1);
    printf("Ring graph:\n");
    do_benchmark(&g_ring, vids_all, 2, 1);
    printf("Random graph:\n");
    do_benchmark(&g_er, vids_all, 2, 1);

    igraph_destroy(&g_full);
    igraph_destroy(&g_ring);
}
