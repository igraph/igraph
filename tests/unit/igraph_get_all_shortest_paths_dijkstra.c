/* igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

void check_nrgeo(const igraph_t *graph, igraph_vs_t vs,
                 const igraph_vector_int_list_t *paths,
                 const igraph_vector_int_t *nrgeo) {
    igraph_int_t i, n;
    igraph_vector_int_t nrgeo2, *path;
    igraph_vit_t vit;

    n = igraph_vcount(graph);
    igraph_vector_int_init(&nrgeo2, n);
    if (igraph_vector_int_size(nrgeo) != n) {
        printf("nrgeo vector length must be %" IGRAPH_PRId ", was %" IGRAPH_PRId, n, igraph_vector_int_size(nrgeo));
        return;
    }

    n = igraph_vector_int_list_size(paths);
    for (i = 0; i < n; i++) {
        path = igraph_vector_int_list_get_ptr(paths, i);
        if (path == 0) {
            printf("Null path found in result vector at index %" IGRAPH_PRId "\n", i);
            return;
        }
        if (igraph_vector_int_size(path) == 0) {
            printf("Empty path found in result vector at index %" IGRAPH_PRId "\n", i);
            return;
        }
        VECTOR(nrgeo2)[igraph_vector_int_tail(path)] += 1;
    }

    igraph_vit_create(graph, vs, &vit);
    for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
        igraph_int_t node = IGRAPH_VIT_GET(vit);
        if (VECTOR(*nrgeo)[node] - VECTOR(nrgeo2)[node]) {
            printf("nrgeo[%" IGRAPH_PRId "] invalid, observed = %" IGRAPH_PRId ", expected = %" IGRAPH_PRId "\n",
                   node, VECTOR(*nrgeo)[node], VECTOR(nrgeo2)[node]);
        }
    }
    igraph_vit_destroy(&vit);

    igraph_vector_int_destroy(&nrgeo2);
}

void print_and_destroy_items(igraph_vector_int_list_t* vec) {
    igraph_int_t i;

    for (i = 0; i < igraph_vector_int_list_size(vec); i++) {
        igraph_vector_int_print(igraph_vector_int_list_get_ptr(vec, i));
    }

    igraph_vector_int_list_clear(vec);
}

int main(void) {

    igraph_t g;
    igraph_vector_int_list_t vertices, edges;

    igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1, 1 };
    igraph_real_t weights2[] = { 0, 2, 1, 0, 5, 2, 1, 1, 0, 2, 2, 8, 1, 1, 3, 1, 1, 4, 2, 1 };
    igraph_int_t dim[] = { 4, 4 };

    igraph_vector_t weights_vec;
    igraph_vector_int_t nrgeo;
    igraph_vector_int_t dim_vec;
    igraph_vs_t vs;

    igraph_vector_int_init(&nrgeo, 0);

    /* Simple ring graph without weights */

    igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);

    igraph_vector_int_list_init(&vertices, 0);
    igraph_vector_int_list_init(&edges, 0);
    igraph_vs_vector_small(&vs, 1, 3, 4, 5, 2, 1,  -1);

    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ &vertices, /*edges=*/ &edges,  /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ NULL, /*mode=*/ IGRAPH_OUT);
    check_nrgeo(&g, vs, &vertices, &nrgeo);
    print_and_destroy_items(&vertices);
    print_and_destroy_items(&edges);

    /* Same ring, but with weights */

    weights_vec = igraph_vector_view(weights, sizeof(weights) / sizeof(weights[0]));
    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ &vertices, /*edges=*/ NULL, /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ &weights_vec, /*mode=*/ IGRAPH_OUT);
    check_nrgeo(&g, vs, &vertices, &nrgeo);
    print_and_destroy_items(&vertices);

    /* we are now testing the combination of vertices == NULL and edges != NUL */

    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ NULL, /*edges=*/ &edges, /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ &weights_vec, /*mode=*/ IGRAPH_OUT);
    print_and_destroy_items(&edges);

    igraph_destroy(&g);

    /* More complicated example */

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 1, 0, 2, 0, 3,   1, 2, 1, 4, 1, 5,
                 2, 3, 2, 6,         3, 2, 3, 6,
                 4, 5, 4, 7,         5, 6, 5, 8, 5, 9,
                 7, 5, 7, 8,         8, 9,
                 5, 2,
                 2, 1,
                 -1);

    weights_vec = igraph_vector_view(weights2, sizeof(weights2) / sizeof(weights2[0]));
    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ &vertices, /*edges=*/ &edges, /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ &weights_vec, /*mode=*/ IGRAPH_OUT);

    check_nrgeo(&g, vs, &vertices, &nrgeo);

    /* Sort the paths in a deterministic manner to avoid problems with
     * different qsort() implementations on different platforms */
    igraph_vector_int_list_sort(&vertices, igraph_vector_int_colex_cmp);
    igraph_vector_int_list_sort(&edges, igraph_vector_int_colex_cmp);
    print_and_destroy_items(&vertices);
    print_and_destroy_items(&edges);

    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    /* Regular lattice with some heavyweight edges */
    dim_vec = igraph_vector_int_view(dim, sizeof(dim) / sizeof(dim[0]));
    igraph_square_lattice(&g, &dim_vec, 1, 0, 0, 0);
    igraph_vs_vector_small(&vs, 3, 12, 15, -1);
    igraph_vector_init(&weights_vec, 24);
    igraph_vector_fill(&weights_vec, 1);
    VECTOR(weights_vec)[2] = 100;
    VECTOR(weights_vec)[8] = 100; /* 1-->2, 4-->8 */
    igraph_get_all_shortest_paths_dijkstra(
                &g,
                /*vertices=*/ 0, /*edges=*/ 0, /*nrgeo=*/ &nrgeo,
                /*from=*/ 0, /*to=*/ vs,
                /*weights=*/ &weights_vec, /*mode=*/ IGRAPH_OUT);
    igraph_vector_destroy(&weights_vec);
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);

    printf("%" IGRAPH_PRId " ", VECTOR(nrgeo)[3]);
    printf("%" IGRAPH_PRId " ", VECTOR(nrgeo)[12]);
    printf("%" IGRAPH_PRId "\n", VECTOR(nrgeo)[15]);

    igraph_vector_int_list_destroy(&vertices);
    igraph_vector_int_list_destroy(&edges);
    igraph_vector_int_destroy(&nrgeo);

    VERIFY_FINALLY_STACK();
    return 0;
}
