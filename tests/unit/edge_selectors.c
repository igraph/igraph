/*
   igraph library.
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
#include "test_utilities.h"

void check(igraph_t *graph, igraph_es_t *es) {
    igraph_eit_t eit;
    igraph_int_t edge;
    IGRAPH_ASSERT(igraph_eit_create(graph, *es, &eit) == IGRAPH_SUCCESS);
    for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
        edge = IGRAPH_EIT_GET(eit);
        printf("%" IGRAPH_PRId " %" IGRAPH_PRId "\n", IGRAPH_FROM(graph, edge), IGRAPH_TO(graph, edge));
    }
    igraph_eit_destroy(&eit);
}

int main(void) {
    igraph_t g, g_no_vertices, g_no_edges;
    igraph_es_t es, es_copy;
    igraph_vector_int_t v, check_as_vector;
    igraph_eit_t eit;

    igraph_small(&g, 5, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, -1);
    igraph_small(&g_no_vertices, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_no_edges, 5, IGRAPH_UNDIRECTED, -1);

    printf("Checking es_vector:\n");
    igraph_vector_int_init_int(&v, 3, 2, 3, 4);
    igraph_es_vector(&es, &v);
    check(&g, &es);
    CHECK_ERROR(igraph_eit_create(&g_no_edges, es, &eit), IGRAPH_EINVEID);
    CHECK_ERROR(igraph_eit_create(&g_no_vertices, es, &eit), IGRAPH_EINVEID);

    printf("es_vector_copy\n");
    igraph_es_vector_copy(&es, &v);
    check(&g, &es);

    printf("es_copy\n");
    igraph_es_copy(&es_copy, &es);
    check(&g, &es);
    igraph_es_destroy(&es_copy);

    printf("es_as_vector\n");
    igraph_vector_int_clear(&v);
    igraph_es_as_vector(&g, es, &v);
    igraph_vector_int_print(&v);
    igraph_es_destroy(&es);
    igraph_vector_int_destroy(&v);

    printf("es_vector with negative entry should fail.\n");
    igraph_vector_int_init_int(&v, 3, -2, 3, 4);
    igraph_es_vector(&es, &v);
    CHECK_ERROR(igraph_eit_create(&g, es, &eit), IGRAPH_EINVEID);
    igraph_vector_int_destroy(&v);

    printf("Checking es_range:\n");
    igraph_es_range(&es, 2, 5);
    check(&g, &es);
    CHECK_ERROR(igraph_eit_create(&g_no_edges, es, &eit), IGRAPH_EINVEID);
    CHECK_ERROR(igraph_eit_create(&g_no_vertices, es, &eit), IGRAPH_EINVEID);

    printf("Checking eit_as_vector using seq:\n");
    igraph_eit_create(&g, es, &eit);
    igraph_vector_int_init_int(&check_as_vector, 0);
    igraph_eit_as_vector(&eit, &check_as_vector);
    igraph_vector_int_print(&check_as_vector);
    igraph_vector_int_destroy(&check_as_vector);

    printf("Checking ess_range using es_range parameters:\n");
    es = igraph_ess_range(2, 5);
    check(&g, &es);
    CHECK_ERROR(igraph_eit_create(&g_no_edges, es, &eit), IGRAPH_EINVEID);
    CHECK_ERROR(igraph_eit_create(&g_no_vertices, es, &eit), IGRAPH_EINVEID);

    printf("Checking whether ess_range accepts an empty range.\n");
    es = igraph_ess_range(2, 2);
    check(&g, &es);
    CHECK_ERROR(igraph_eit_create(&g_no_edges, es, &eit), IGRAPH_EINVEID);
    CHECK_ERROR(igraph_eit_create(&g_no_vertices, es, &eit), IGRAPH_EINVEID);

    printf("Checking es_path:\n");
    igraph_vector_int_init_int(&v, 3, 4, 3, 2);
    igraph_es_path(&es, &v, /*directed*/0);
    check(&g, &es);
    CHECK_ERROR(igraph_eit_create(&g_no_vertices, es, &eit), IGRAPH_EINVVID);
    CHECK_ERROR(igraph_eit_create(&g_no_edges, es, &eit), IGRAPH_EINVAL);
    igraph_es_destroy(&es);

    igraph_es_path(&es, &v, /*directed*/1);
    CHECK_ERROR(igraph_eit_create(&g, es, &eit), IGRAPH_EINVAL);
    igraph_vector_int_destroy(&v);
    igraph_es_destroy(&es);

    printf("es_path with negative entry should fail.\n");
    igraph_vector_int_init_int(&v, 3, -4, 3, 2);
    igraph_es_path(&es, &v, /*directed*/0);
    CHECK_ERROR(igraph_eit_create(&g, es, &eit), IGRAPH_EINVVID);

    printf("Checking es_type.\n");
    IGRAPH_ASSERT(igraph_es_type(&es) == IGRAPH_ES_PATH);
    igraph_es_destroy(&es);

    printf("Checking es_none.\n");
    igraph_es_none(&es);
    check(&g, &es);

    printf("es_pairs:\n");
    igraph_vector_int_destroy(&v);
    igraph_vector_int_init_int(&v, 4, 1,1, 2,0);
    igraph_es_pairs(&es, &v, IGRAPH_DIRECTED);
    check(&g, &es);
    igraph_es_destroy(&es);

    printf("es_pairs for nonexistent pair should fail.\n");
    igraph_vector_int_destroy(&v);
    igraph_vector_int_init_int(&v, 4, 6,1, 2,0);
    igraph_es_pairs(&es, &v, IGRAPH_DIRECTED);
    CHECK_ERROR(igraph_eit_create(&g, es, &eit), IGRAPH_EINVVID);

    igraph_vector_int_destroy(&v);
    igraph_destroy(&g);
    igraph_destroy(&g_no_vertices);
    igraph_destroy(&g_no_edges);

    VERIFY_FINALLY_STACK();
    return 0;
}
