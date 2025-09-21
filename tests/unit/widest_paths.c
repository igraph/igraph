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

/* initialise structures used in "get_widest_path(s)" algorithms */
void init_vertices_and_edges(igraph_int_t n,
                            igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                            igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                            igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {
    igraph_vector_int_list_init(vertices, n);
    igraph_vector_int_list_init(edges, n);

    igraph_vector_int_init(parents, 0);
    igraph_vector_int_init(inbound_edges, 0);

    igraph_vector_int_init(vertices2, 0);
    igraph_vector_int_init(edges2, 0);
}

/* destroy structures used in "get_widest_path(s)" algorithms, ignores if NULL*/
void destroy_vertices_and_edges(igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                            igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                            igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {
    if (edges2)         igraph_vector_int_destroy(edges2);
    if (vertices2)      igraph_vector_int_destroy(vertices2);

    if (inbound_edges)  igraph_vector_int_destroy(inbound_edges);
    if (parents)   igraph_vector_int_destroy(parents);

    if (edges) igraph_vector_int_list_destroy(edges);
    if (vertices) igraph_vector_int_list_destroy(vertices);
}

/* destroy all structures, ignores if NULL */
void destroy_all(igraph_t *g, igraph_vector_t *w, igraph_matrix_t *res1, igraph_matrix_t *res2,
                    igraph_vs_t *from, igraph_vs_t *to,
                    igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                    igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                    igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {

    destroy_vertices_and_edges(vertices, edges, parents, inbound_edges, vertices2, edges2);

    if (to)     igraph_vs_destroy(to);
    if (from)   igraph_vs_destroy(from);
    if (res2)   igraph_matrix_destroy(res2);
    if (res1)   igraph_matrix_destroy(res1);

    if (w)      igraph_vector_destroy(w);
    if (g)      igraph_destroy(g);
}

/* confirm that all widest paths algorithms fail */
void check_invalid_input(igraph_t *g, igraph_vector_t *w, igraph_matrix_t *res,
                        igraph_int_t source, igraph_int_t destination,
                        igraph_vs_t *from, igraph_vs_t *to,
                        igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                        igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                        igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {
    CHECK_ERROR(igraph_widest_path_widths_dijkstra(g, res, *from, *to, w, IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_widest_path_widths_floyd_warshall(g, res, *from, *to, w, IGRAPH_OUT), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_paths(g, vertices, edges, source, *to, w, IGRAPH_OUT,
                parents, inbound_edges), IGRAPH_EINVAL);
    CHECK_ERROR(igraph_get_widest_path(g, vertices2, edges2, source, destination, w, IGRAPH_OUT), IGRAPH_EINVAL);
}

/* runs the width finding algorithms and assert they succeed */
void run_widest_paths(igraph_t *g, igraph_vector_t *w, igraph_matrix_t *res1, igraph_matrix_t *res2,
                    igraph_vs_t *from, igraph_vs_t *to, igraph_neimode_t mode) {
    IGRAPH_ASSERT(igraph_widest_path_widths_dijkstra(g, res1, *from, *to, w, mode) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_widest_path_widths_floyd_warshall(g, res2, *from, *to, w, mode) == IGRAPH_SUCCESS);
}

/* runs the path finding algorithms and asserts they succeed */
void run_get_widest_paths(igraph_t *g, igraph_vector_t *w,
                        igraph_int_t source, igraph_int_t destination,
                        igraph_vs_t *to, igraph_neimode_t mode,
                        igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                        igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                        igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {
    IGRAPH_ASSERT(igraph_get_widest_paths(g, vertices, edges, source, *to, w, mode,
                    parents, inbound_edges) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_get_widest_path(g, vertices2, edges2, source, destination,
                    w, mode) == IGRAPH_SUCCESS);
}

/* print results of just the matrices */
void check_and_print_matrices(igraph_matrix_t *res1, igraph_matrix_t *res2) {
    IGRAPH_ASSERT(igraph_matrix_all_e(res1, res2));
    print_matrix_format(res1, stdout, "%f");
}

/* prints results of all widest path algorithms */
void print_results(igraph_int_t n, igraph_matrix_t *res1, igraph_matrix_t *res2,
                    igraph_vector_int_list_t *vertices, igraph_vector_int_list_t *edges,
                    igraph_vector_int_t *parents, igraph_vector_int_t *inbound_edges,
                    igraph_vector_int_t *vertices2, igraph_vector_int_t *edges2) {
    igraph_int_t i;

    check_and_print_matrices(res1, res2);
    printf("\n");

    for (i = 0; i < n; i++) {
        printf("path to node %" IGRAPH_PRId ":\n", i);
        igraph_vector_int_t *vertex_path = igraph_vector_int_list_get_ptr(vertices, i);
        igraph_vector_int_t *edge_path = igraph_vector_int_list_get_ptr(edges, i);
        printf("  vertices: ");
        print_vector_int(vertex_path);
        printf("  edges:    ");
        print_vector_int(edge_path);
    }
    printf("\n");

    printf("parents:  ");
    print_vector_int(parents);
    printf("inbound_edges: ");
    print_vector_int(inbound_edges);
    printf("\n");

    printf("vertex path:   ");
    print_vector_int(vertices2);
    printf("edge path:     ");
    print_vector_int(edges2);
    printf("\n");
}


int main(void) {

    igraph_int_t n, m;
    igraph_t g;
    igraph_matrix_t res1, res2;
    igraph_vs_t from, to;
    igraph_vector_t w;

    igraph_vector_int_list_t vertices;
    igraph_vector_int_list_t edges;
    igraph_vector_int_t parents;
    igraph_vector_int_t inbound_edges;

    igraph_vector_int_t vertices2;
    igraph_vector_int_t edges2;


    /* ==================================================================== */
    /* 1. null weight vector                                                */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_range(&from, 0, n);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    check_invalid_input(&g, NULL, &res1, 0, 2, &from, &to, &vertices, &edges,
                        &parents, &inbound_edges, &vertices2, &edges2);

    destroy_vertices_and_edges(&vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, NULL, &res1, NULL, &from, &to, NULL, NULL, NULL, NULL, NULL, NULL);

    /* ==================================================================== */
    /* 2. Number of weights don't match number of edges                     */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_vector_init_real(&w, m+3, -1.0, 2.0, 3.0, 2.0, 5.0, 1.0);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_range(&from, 0, n);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    check_invalid_input(&g, &w, &res1, 0, 2, &from, &to, &vertices, &edges,
                        &parents, &inbound_edges, &vertices2, &edges2);

    destroy_vertices_and_edges(&vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, NULL, &from, &to, NULL, NULL, NULL, NULL, NULL, NULL);

    /* ==================================================================== */
    /* 3. NaN values in weights                                             */

    n = 3;
    m = 3;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 2, 0, -1);
    igraph_vector_init_real(&w, m, -1.0, IGRAPH_NAN, 3.0);

    igraph_matrix_init(&res1, n, n);
    igraph_vs_range(&from, 0, n);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    check_invalid_input(&g, &w, &res1, 0, 2, &from, &to, &vertices, &edges,
                        &parents, &inbound_edges, &vertices2, &edges2);

    destroy_vertices_and_edges(&vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, NULL, &from, &to, NULL, NULL, NULL, NULL, NULL, NULL);

    /* ==================================================================== */
    /* 4. Empty graph                                                       */

    n = 0;
    m = 0;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, -1);
    igraph_vector_init(&w, m);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_none(&from);
    igraph_vs_none(&to);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    /* Should successfully run with nothing occuring */

    destroy_all(&g, &w, &res1, &res2, &from, &to, NULL, NULL, NULL, NULL, NULL, NULL);

    /* ==================================================================== */
    /* 5. 1 node graph                                                      */
    printf("=== 5. Testing 1 Node Graph ===\n");

    n = 1;
    m = 0;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, -1);
    igraph_vector_init(&w, m);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_all(&from);
    igraph_vs_all(&to);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 0, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 6. Unreachable Nodes                                                 */
    printf("\n=== 6. Testing Unreachable Nodes ===\n");

    n = 4;
    m = 4;
    igraph_small(&g, n, IGRAPH_DIRECTED, 0, 2, 0, 1, 1, 2, 3, 2, -1);
    igraph_vector_init_real(&w, m, 1.0, 2.0, 3.0, 5.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 3, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 7. Self Loops                                                        */
    printf("\n=== 7. Testing Self Loops ===\n");

    n = 3;
    m = 6;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 2, 0, 1, 1, 2, 0, 0, 1, 1, 2, 2, -1);
    igraph_vector_init_real(&w, m, 1.0, 2.0, 3.0, 5.0, 5.0, 1.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_all(&from);
    igraph_vs_all(&to);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 2, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    /* Self loops should be effectively ignored, since the widest path to
       yourself is always infinity. */
    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 8. Multiple Edges                                                    */
    printf("\n=== 8. Testing Multiple Edges ===\n");

    n = 2;
    m = 8;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 0, 1, 0, 1, 0, 1,
                                        0, 1, 0, 1, 0, 1, 0, 1, -1);
    igraph_vector_init_real(&w, m, 2.0, 2.0, 2.0, 10.0, 2.0, 2.0, 2.0, 2.0);

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 1, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 9. Directed Graphs                                                   */
    printf("\n=== 9. Testing Directed Graphs ===\n");

    n = 2;
    m = 6;
    igraph_small(&g, n, IGRAPH_DIRECTED, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, -1);
    igraph_vector_init_real(&w, m, 100.0, 200.0, 100.0, 200.0, 1.0, 200.0 );

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 1, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 10. Mode                                                             */
    printf("\n=== 10. Testing Mode ===\n");
    /* This is same as test 9, just with the mode reversed. */

    n = 2;
    m = 6;
    igraph_small(&g, n, IGRAPH_DIRECTED, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, -1);
    igraph_vector_init_real(&w, m, 100.0, 200.0, 100.0, 200.0, 1.0, 200.0 );

    igraph_matrix_init(&res1, 1, n);
    igraph_matrix_init(&res2, 1, n);
    igraph_vs_1(&from, 0);
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_IN);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 1, &to, IGRAPH_IN,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


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
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 7, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    check_and_print_matrices(&res1, &res2);

    printf("\npath to node 7:\n");
    igraph_vector_int_t *vertex_path = igraph_vector_int_list_get_ptr(&vertices, 7);
    igraph_vector_int_t *edge_path = igraph_vector_int_list_get_ptr(&edges, 7);

    printf("  vertices: ");
    print_vector_int(vertex_path);
    printf("  edges:    ");
    print_vector_int(edge_path);
    printf("parents:  ");
    print_vector_int(&parents);
    printf("inbound_edges: ");
    print_vector_int(&inbound_edges);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);

    /* ==================================================================== */
    /* 12. 5 Node Simple Graph                                              */
    printf("\n=== 12. Testing 5 Node Simple Graph ===\n");

    n = 5;
    m = 5;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0, 1, 1, 2, 2, 3, 3, 4, 0, 3,
                 -1);
    igraph_vector_init_real(&w, m, 8.0, 6.0, 10.0, 7.0, 5.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_all(&from);
    igraph_vs_all(&to);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 4, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);

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
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 3, /* destination */ 6, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


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
    igraph_vs_range(&to, 0, n);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 2, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);


    /* ==================================================================== */
    /* 15. Disconnected Graph                                               */
    printf("\n=== 15. Testing Disconnected Graphs ===\n");

    n = 5;
    m = 4;
    igraph_small(&g, n, IGRAPH_UNDIRECTED, 0, 1, 1, 2, 0, 2, 3, 4, -1);
    igraph_vector_init_real(&w, m, 2.0, 5.0, 5.0, 7.0);

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_all(&from);
    igraph_vs_all(&to);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_OUT);
    run_get_widest_paths(&g, &w, /* source */ 0, /* destination */ 4, &to, IGRAPH_OUT,
                            &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);

    /* ==================================================================== */
    /* 16. Infinite weights                                                 */

    n = 10;
    m = 6;
    igraph_small(&g, n, IGRAPH_UNDIRECTED,
                 0,1, 1,2, 3,4, 4,5, 6,7, 8,9,
                 -1);
    igraph_vector_init_real(&w, m,
                            0.5, IGRAPH_INFINITY, /* 0,1, 1,2 */
                            1.5, -IGRAPH_INFINITY, /* 3,4, 4,5 */
                            IGRAPH_INFINITY, /* 6,7 */
                            -IGRAPH_INFINITY); /* 8,9 */

    igraph_matrix_init(&res1, n, n);
    igraph_matrix_init(&res2, n, n);
    igraph_vs_all(&from);
    igraph_vs_all(&to);

    init_vertices_and_edges(n, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    /* Note: Since edges with negative infinite weight are ignored, there is
     * no path from 3 to 5. */
    run_widest_paths(&g, &w, &res1, &res2, &from, &to, IGRAPH_ALL);
    run_get_widest_paths(&g, &w, /* source */ 3, /* destination */ 5, &to, IGRAPH_ALL,
                         &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    print_results(n, &res1, &res2, &vertices, &edges, &parents, &inbound_edges, &vertices2, &edges2);

    destroy_all(&g, &w, &res1, &res2, &from, &to, &vertices, &edges, &parents, &inbound_edges,
                &vertices2, &edges2);

    VERIFY_FINALLY_STACK();

    return 0;
}
