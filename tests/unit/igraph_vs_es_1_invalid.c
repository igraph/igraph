/* igraph library.
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

/* Single-vertex and single-edge selectors must reject out-of-range IDs,
 * including negative ones and IGRAPH_INTEGER_MAX. */

int main(void) {
    igraph_t g;
    igraph_vit_t vit;
    igraph_eit_t eit;
    igraph_vector_int_t res, vertices, edges;
    igraph_vector_t weights;
    const igraph_int_t bad_ids[] = { -1, -2, 4, IGRAPH_INTEGER_MAX };
    const size_t bad_count = sizeof(bad_ids) / sizeof(bad_ids[0]);

    igraph_vector_int_init(&res, 0);
    igraph_vector_int_init(&vertices, 0);
    igraph_vector_int_init(&edges, 0);
    igraph_vector_init_int(&weights, 3, 1, 2, 3);

    /* Path graph with 4 vertices and 3 edges. */
    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, -1);

    printf("Valid IDs:\n");
    IGRAPH_ASSERT(igraph_vit_create(&g, igraph_vss_1(3), &vit) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(IGRAPH_VIT_SIZE(vit) == 1 && IGRAPH_VIT_GET(vit) == 3);
    igraph_vit_destroy(&vit);
    IGRAPH_ASSERT(igraph_eit_create(&g, igraph_ess_1(2), &eit) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(IGRAPH_EIT_SIZE(eit) == 1 && IGRAPH_EIT_GET(eit) == 2);
    igraph_eit_destroy(&eit);
    igraph_degree(&g, &res, igraph_vss_1(1), IGRAPH_ALL, IGRAPH_LOOPS);
    print_vector_int(&res);
    igraph_get_shortest_path(&g, NULL, &vertices, &edges, 0, 3, IGRAPH_ALL);
    print_vector_int(&vertices);
    igraph_get_widest_path(&g, &vertices, &edges, 0, 3, &weights, IGRAPH_ALL);
    print_vector_int(&vertices);

    printf("Invalid IDs:\n");
    for (size_t i = 0; i < bad_count; i++) {
        igraph_int_t id = bad_ids[i];

        CHECK_ERROR(igraph_vit_create(&g, igraph_vss_1(id), &vit), IGRAPH_EINVVID);
        CHECK_ERROR(igraph_degree(&g, &res, igraph_vss_1(id), IGRAPH_ALL, IGRAPH_LOOPS), IGRAPH_EINVVID);

        CHECK_ERROR(igraph_get_shortest_path(&g, NULL, &vertices, &edges, 0, id, IGRAPH_ALL), IGRAPH_EINVVID);
        CHECK_ERROR(igraph_get_shortest_path(&g, &weights, &vertices, &edges, 0, id, IGRAPH_ALL), IGRAPH_EINVVID);
        CHECK_ERROR(igraph_get_shortest_path_dijkstra(&g, &vertices, &edges, 0, id, &weights, IGRAPH_ALL), IGRAPH_EINVVID);
        CHECK_ERROR(igraph_get_shortest_path_bellman_ford(&g, &vertices, &edges, 0, id, &weights, IGRAPH_ALL), IGRAPH_EINVVID);
        CHECK_ERROR(igraph_get_widest_path(&g, &vertices, &edges, 0, id, &weights, IGRAPH_ALL), IGRAPH_EINVVID);

        CHECK_ERROR(igraph_eit_create(&g, igraph_ess_1(id), &eit), IGRAPH_EINVEID);
        CHECK_ERROR(igraph_delete_edges(&g, igraph_ess_1(id)), IGRAPH_EINVEID);
        IGRAPH_ASSERT(igraph_ecount(&g) == 3);

        printf("%" IGRAPH_PRId " rejected\n", id);
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&edges);
    igraph_vector_int_destroy(&vertices);
    igraph_vector_int_destroy(&res);

    VERIFY_FINALLY_STACK();

    return 0;
}
