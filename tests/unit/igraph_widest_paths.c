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

#include "test_utilities.inc"

int main() {

    igraph_integer_t n, m;
    igraph_t g;
    igraph_matrix_t res1, res2;
    igraph_vs_t from, to;
    igraph_vector_t w;
    igraph_integer_t i;

    igraph_vector_ptr_t vertices;
    igraph_vector_ptr_t edges;
    igraph_vector_int_t predecessors;
    igraph_vector_int_t inbound_edges;

    igraph_vector_int_t vertices2;
    igraph_vector_int_t edges2;


    /* ==================================================================== */
    /* 1. null weight vector                                                */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, -1);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_seq(&from, 0, n-1);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    CHECK_ERROR(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ NULL,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ NULL,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ NULL, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 2, /* weights */ NULL, /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res1);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 2. Number of weights don't match number of edges                     */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, -1);
    igraph_vector_init_real(&w, m+3, -1.0, IGRAPH_NAN, 3.0, 2.0, 5.0, 1.0);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_seq(&from, 0, n-1);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    CHECK_ERROR(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ NULL,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ NULL,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ NULL, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 2, /* weights */ NULL, /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 3. NaN values in weights                                             */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 3, -1);
    igraph_vector_init_real(&w, m, -1.0, IGRAPH_NAN, 3.0);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_seq(&from, 0, n-1);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    CHECK_ERROR(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 2, /* weights */ &w, /* mode */ IGRAPH_OUT), IGRAPH_EINVAL);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);

    /* ==================================================================== */
    /* 4. Empty graph                                                       */

    n = 0;
    m = 0;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, -1);
    igraph_vector_init(&w, m);

    igraph_matrix_init(&res1, 0, 0);
    igraph_matrix_init(&res2, 0, 0);
    igraph_vs_none(&from);
    igraph_vs_none(&to);


    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 5. 1 node graph                                                      */
    printf("=== 5. Testing 1 Node Graph ===\n");

    n = 1;
    m = 0;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, -1);
    igraph_vector_init(&w, m);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_seq(&from, 0, n-1);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 6. Unrechable Nodes                                                  */
    printf("\n=== 6. Testing Unreachable Nodes ===\n");

    n = 4;
    m = 4;
    igraph_small(&g, n, IGRAPH_DIRECTED, 0, 2, 0, 1, 1, 2, 3, 2, -1);
    igraph_vector_init_real(&w, m, 1.0, 2.0, 3.0, 5.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 7. Self Loops                                                        */
    printf("\n=== 7. Testing Self Loops ===\n");

    n = 3;
    m = 6;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 2, 0, 1, 1, 2, 0, 0, 1, 1, 2, 2, -1);
    igraph_vector_init_real(&w, m, 1.0, 2.0, 3.0, 5.0, 5.0, 1.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);
    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 8. Multiple Edges                                                    */
    printf("\n=== 8. Testing Multiple Edges ===\n");

    n = 2;
    m = 8;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 0, 1, 0, 1, 0, 1,
                                        0, 1, 0, 1, 0, 1, 0, 1, -1);
    igraph_vector_init_real(&w, m, 2.0, 2.0, 2.0, 10.0, 2.0, 2.0, 2.0, 2.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);
    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 9. Directed Graphs                                                   */
    printf("\n=== 9. Testing Directed Graphs ===\n");

    n = 2;
    m = 6;
    igraph_small(&g, n, IGRAPH_DIRECTED, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, -1);
    igraph_vector_init_real(&w, m, 100.0, 200.0, 100.0, 200.0, 1.0, 200.0 );

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);
    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 10. Mode                                                             */
    printf("\n=== 10. Testing Mode ===\n");

    n = 2;
    m = 6;
    igraph_small(&g, n, IGRAPH_DIRECTED, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, -1);
    igraph_vector_init_real(&w, m, 100.0, 200.0, 100.0, 200.0, 1.0, 200.0 );

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);
    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_IN) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_IN) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_IN,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_IN)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 11. Multiple Widest Paths                                            */
    printf("\n=== 11. Testing Multiple Widest Paths ===\n");

    n = 8;
    m = 9;
    igraph_small(&g, n, IGRAPH_DIRECTED, 0, 1, 1, 2, 2, 7, 0, 3, 3, 4, 4, 7,
                                            0, 5, 5, 6, 6, 7, -1);
    igraph_vector_init_real(&w, m, 100.0, 10.0, 200.0, 10.0, 15.0, 35.0, 3300.0, 10.0, 10.0 );

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);
    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    printf("node: 7\n");
    igraph_vector_int_t *vertex_path = VECTOR(vertices)[7];
    igraph_vector_int_t *edge_path = VECTOR(edges)[7];

    print_vector_int(vertex_path);
    print_vector_int(edge_path);
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 12. 5 Node Simple Graph                                              */
    printf("\n=== 12. Testing 5 Node Simple Graph ===\n");

    n = 5;
    m = 5;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 3,
                 -1);
    igraph_vector_init_real(&w, m, 8.0, 6.0, 10.0, 7.0, 5.0);

    igraph_matrix_init(&res1, 5, 5);
    igraph_matrix_init(&res2, 5, 5);
    igraph_vs_seq(&from, 0, 4);
    igraph_vs_seq(&to, 0, 4);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(&g, &res1, from, to, &w, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(&g, &res2, from, to, &w, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_matrix_format(&res2, stdout, "%f");

    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));

    igraph_vector_destroy(&w);
    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);
    igraph_destroy(&g);

    /* ==================================================================== */
    /* 13. 7 Node Wikipedia Graph                                           */
    printf("\n=== 13. Testing 7 Node Wikipedia Graph ===\n");


    n = 7;
    m = 11;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 1, 2, 1, 3, 2, 4, 2, 5,
                 3, 4, 3, 6, 4, 5, 4, 6, 5, 6,
                 -1);
    igraph_vector_init_real(&w, m, 15.0, 53.0, 40.0, 46.0, 31.0,
                            17.0, 3.0, 11.0, 29.0, 8.0, 40.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 3);
    igraph_vs_seq(&to, 0, n-1);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(&g, &res1, from, to, &w, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(&g, &res2, from, to, &w, IGRAPH_OUT) == IGRAPH_SUCCESS);
    print_matrix_format(&res2, stdout, "%f");

    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);


    /* === igraph_get_widest_paths === */

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vs_seq(&to, 0, n-1);
    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                                /* from */ 3, /* to */ to, /* weights */ &w,
                                /* mode */ IGRAPH_OUT,
                                /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges
                                ) == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vs_destroy(&to);
    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    /* === igraph_get_widest_path === */

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                                /* from */ 3, /* to */ 6, /* weights */ &w,
                                /* mode */ IGRAPH_OUT ) == IGRAPH_SUCCESS);
    printf("get_widest_path()\n");
    print_vector_int(&vertices2);
    print_vector_int(&edges2);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 14. Negative Widths                                                  */
    printf("\n=== 14. Testing Negative Weights ===\n");

    n = 4;
    m = 4;
    igraph_small(&g, n, IGRAPH_DIRECTED, 0, 1, 1, 2, 0, 3, 3, 2, -1);
    igraph_vector_init_real(&w, m, 2.0, -3.0, -5.0, 6.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);


    /* ==================================================================== */
    /* 15. Disconnected Graph                                               */
    printf("\n=== 15. Testing Disconnected Graphs ===\n");

    n = 5;
    m = 4;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 0, 2, 3, 4, -1);
    igraph_vector_init_real(&w, m, 2.0, 5.0, 5.0, 7.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_seq(&from, 0, n-1);
    igraph_vs_seq(&to, 0, n-1);

    igraph_vector_ptr_init(&vertices, n);
    igraph_vector_ptr_init(&edges, n);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        VECTOR(vertices)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(vertices)[i], 0);
        VECTOR(edges)[i] = calloc(1, sizeof(igraph_vector_int_t));
        igraph_vector_init(VECTOR(edges)[i], 0);
    }
    igraph_vector_int_init(&predecessors, 0);
    igraph_vector_int_init(&inbound_edges, 0);

    igraph_vector_int_init(&vertices2, 0);
    igraph_vector_int_init(&edges2, 0);

    IGRAPH_ASSERT(igraph_widest_paths_dijkstra(/* graph */ &g, /* res */ &res1,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_paths_floyd_warshall(/* graph */ &g, /* res */ &res2,
                            /* from */ from, /* to */ to,
                            /* weights */ &w,
                            /* mode */ IGRAPH_OUT) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_matrix_all_e(&res1, &res2));
    print_matrix_format(&res1, stdout, "%f");

    IGRAPH_ASSERT(igraph_get_widest_paths(/* graph */ &g, /* vertices */ &vertices, /* edges */ &edges,
                    /* from */ 0, /* to */ to, /* weights */ &w, /* mode */ IGRAPH_OUT,
                    /* predecessors */ &predecessors, /* inbound_edges */ &inbound_edges)
                    == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(/* graph */ &g, /* vertices */ &vertices2, /* edges */ &edges2,
                    /* from */ 0, /* to */ 0, /* weights */ &w, /* mode */ IGRAPH_OUT)
                    == IGRAPH_SUCCESS);
    for (i = 0; i < n; i++) {
        printf("node: %" IGRAPH_PRId "\n", i);
        igraph_vector_int_t *vertex_path = VECTOR(vertices)[i];
        igraph_vector_int_t *edge_path = VECTOR(edges)[i];
        print_vector_int(vertex_path);
        print_vector_int(edge_path);
    }
    printf("predecessors:\n");
    print_vector_int(&predecessors);
    printf("inbound_edges:\n");
    print_vector_int(&inbound_edges);

    igraph_vector_int_destroy(&edges2);
    igraph_vector_int_destroy(&vertices2);

    igraph_vector_int_destroy(&inbound_edges);
    igraph_vector_int_destroy(&predecessors);
    for (i = 0; i < igraph_vector_ptr_size(&vertices); i++) {
        igraph_vector_int_destroy(VECTOR(vertices)[i]);
        free(VECTOR(vertices)[i]);
        igraph_vector_int_destroy(VECTOR(edges)[i]);
        free(VECTOR(edges)[i]);
    }
    igraph_vector_ptr_destroy(&edges);
    igraph_vector_ptr_destroy(&vertices);

    igraph_vs_destroy(&to);
    igraph_vs_destroy(&from);
    igraph_matrix_destroy(&res2);
    igraph_matrix_destroy(&res1);

    igraph_vector_destroy(&w);
    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
