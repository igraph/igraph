/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <stdio.h>

int main() {

    igraph_t g;
    FILE *karate, *neural;
    igraph_real_t res;
    igraph_vector_t types;
    igraph_vector_t degree, outdegree, indegree;

    igraph_real_t football_types[] = {
        7, 0, 2, 3, 7, 3, 2, 8, 8, 7, 3, 10, 6, 2, 6, 2, 7, 9, 6, 1, 9, 8, 8, 7, 10, 0, 6, 9,
        11, 1, 1, 6, 2, 0, 6, 1, 5, 0, 6, 2, 3, 7, 5, 6, 4, 0, 11, 2, 4, 11, 10, 8, 3, 11, 6,
        1, 9, 4, 11, 10, 2, 6, 9, 10, 2, 9, 4, 11, 8, 10, 9, 6, 3, 11, 3, 4, 9, 8, 8, 1, 5, 3,
        5, 11, 3, 6, 4, 9, 11, 0, 5, 4, 4, 7, 1, 9, 9, 10, 3, 6, 2, 1, 3, 0, 7, 0, 2, 3, 8, 0,
        4, 8, 4, 9, 11
    };

    karate = fopen("karate.gml", "r");
    igraph_read_graph_gml(&g, karate);
    fclose(karate);

    igraph_vector_init(&types, 0);
    igraph_degree(&g, &types, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);

    igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 0);
    printf("%.5f\n", res);

    igraph_destroy(&g);

    /*---------------------*/

    neural = fopen("celegansneural.gml", "r");
    igraph_read_graph_gml(&g, neural);
    fclose(neural);

    igraph_degree(&g, &types, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);

    igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 1);
    printf("%.5f\n", res);
    igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 0);
    printf("%.5f\n", res);

    igraph_destroy(&g);
    igraph_vector_destroy(&types);

    /*---------------------*/

    karate = fopen("karate.gml", "r");
    igraph_read_graph_gml(&g, karate);
    fclose(karate);

    igraph_vector_init(&degree, 0);
    igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
    igraph_vector_add_constant(&degree, -1);

    igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 0);
    printf("%.5f\n", res);

    igraph_destroy(&g);

    /*---------------------*/

    neural = fopen("celegansneural.gml", "r");
    igraph_read_graph_gml(&g, neural);
    fclose(neural);

    igraph_degree(&g, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
    igraph_vector_add_constant(&degree, -1);

    igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 1);
    printf("%.5f\n", res);
    igraph_assortativity(&g, &degree, 0, &res, /*directed=*/ 0);
    printf("%.5f\n", res);

    igraph_vector_destroy(&degree);

    /*---------------------*/

    igraph_vector_init(&indegree, 0);
    igraph_vector_init(&outdegree, 0);
    igraph_degree(&g, &indegree, igraph_vss_all(), IGRAPH_IN, /*loops=*/ 1);
    igraph_degree(&g, &outdegree, igraph_vss_all(), IGRAPH_OUT, /*loops=*/ 1);
    igraph_vector_add_constant(&indegree, -1);
    igraph_vector_add_constant(&outdegree, -1);

    igraph_assortativity(&g, &outdegree, &indegree, &res, /*directed=*/ 1);
    printf("%.5f\n", res);

    igraph_vector_destroy(&indegree);
    igraph_vector_destroy(&outdegree);

    /*---------------------*/

    igraph_assortativity_degree(&g, &res, /*directed=*/ 1);
    printf("%.5f\n", res);

    igraph_destroy(&g);

    /*---------------------*/

    karate = fopen("karate.gml", "r");
    igraph_read_graph_gml(&g, karate);
    fclose(karate);

    igraph_assortativity_degree(&g, &res, /*directed=*/ 1);
    printf("%.5f\n", res);

    igraph_destroy(&g);

    /*---------------------*/

    igraph_small(&g, sizeof(football_types) / sizeof(igraph_real_t),
                 IGRAPH_UNDIRECTED,
                 0, 1, 2, 3, 0, 4, 4, 5, 3, 5, 2, 6, 6, 7, 7, 8, 8, 9, 0, 9, 4, 9, 5, 10, 10, 11, 5, 11,
                 3, 11, 12, 13, 2, 13, 2, 14, 12, 14, 14, 15, 13, 15, 2, 15, 4, 16, 9, 16, 0, 16,
                 16, 17, 12, 17, 12, 18, 18, 19, 17, 20, 20, 21, 8, 21, 7, 21, 9, 22, 7, 22, 21,
                 22, 8, 22, 22, 23, 9, 23, 4, 23, 16, 23, 0, 23, 11, 24, 24, 25, 1, 25, 3, 26, 12,
                 26, 14, 26, 26, 27, 17, 27, 1, 27, 17, 27, 4, 28, 11, 28, 24, 28, 19, 29, 29,
                 30, 19, 30, 18, 31, 31, 32, 21, 32, 15, 32, 13, 32, 6, 32, 0, 33, 1, 33, 25, 33,
                 19, 33, 31, 34, 26, 34, 12, 34, 18, 34, 34, 35, 0, 35, 29, 35, 19, 35, 30, 35,
                 18, 36, 12, 36, 20, 36, 19, 36, 36, 37, 1, 37, 25, 37, 33, 37, 18, 38, 16, 38,
                 28, 38, 26, 38, 14, 38, 12, 38, 38, 39, 6, 39, 32, 39, 13, 39, 15, 39, 7, 40, 3,
                 40, 40, 41, 8, 41, 4, 41, 23, 41, 9, 41, 0, 41, 16, 41, 34, 42, 29, 42, 18, 42,
                 26, 42, 42, 43, 36, 43, 26, 43, 31, 43, 38, 43, 12, 43, 14, 43, 19, 44, 35, 44,
                 30, 44, 44, 45, 13, 45, 33, 45, 1, 45, 37, 45, 25, 45, 21, 46, 46, 47, 22, 47,
                 6, 47, 15, 47, 2, 47, 39, 47, 32, 47, 44, 48, 48, 49, 32, 49, 46, 49, 30, 50,
                 24, 50, 11, 50, 28, 50, 50, 51, 40, 51, 8, 51, 22, 51, 21, 51, 3, 52, 40, 52, 5,
                 52, 52, 53, 25, 53, 48, 53, 49, 53, 46, 53, 39, 54, 31, 54, 38, 54, 14, 54, 34,
                 54, 18, 54, 54, 55, 31, 55, 6, 55, 35, 55, 29, 55, 19, 55, 30, 55, 27, 56, 56,
                 57, 1, 57, 42, 57, 44, 57, 48, 57, 3, 58, 6, 58, 17, 58, 36, 58, 36, 59, 58, 59,
                 59, 60, 10, 60, 39, 60, 6, 60, 47, 60, 13, 60, 15, 60, 2, 60, 43, 61, 47, 61,
                 54, 61, 18, 61, 26, 61, 31, 61, 34, 61, 61, 62, 20, 62, 45, 62, 17, 62, 27, 62,
                 56, 62, 27, 63, 58, 63, 59, 63, 42, 63, 63, 64, 9, 64, 32, 64, 60, 64, 2, 64, 6,
                 64, 47, 64, 13, 64, 0, 65, 27, 65, 17, 65, 63, 65, 56, 65, 20, 65, 65, 66, 59,
                 66, 24, 66, 44, 66, 48, 66, 16, 67, 41, 67, 46, 67, 53, 67, 49, 67, 67, 68, 15,
                 68, 50, 68, 21, 68, 51, 68, 7, 68, 22, 68, 8, 68, 4, 69, 24, 69, 28, 69, 50, 69,
                 11, 69, 69, 70, 43, 70, 65, 70, 20, 70, 56, 70, 62, 70, 27, 70, 60, 71, 18, 71,
                 14, 71, 34, 71, 54, 71, 38, 71, 61, 71, 31, 71, 71, 72, 2, 72, 10, 72, 3, 72,
                 40, 72, 52, 72, 7, 73, 49, 73, 53, 73, 67, 73, 46, 73, 73, 74, 2, 74, 72, 74, 5,
                 74, 10, 74, 52, 74, 3, 74, 40, 74, 20, 75, 66, 75, 48, 75, 57, 75, 44, 75, 75,
                 76, 27, 76, 59, 76, 20, 76, 70, 76, 66, 76, 56, 76, 62, 76, 73, 77, 22, 77, 7,
                 77, 51, 77, 21, 77, 8, 77, 77, 78, 23, 78, 50, 78, 28, 78, 22, 78, 8, 78, 68,
                 78, 7, 78, 51, 78, 31, 79, 43, 79, 30, 79, 19, 79, 29, 79, 35, 79, 55, 79, 79,
                 80, 37, 80, 29, 80, 16, 81, 5, 81, 40, 81, 10, 81, 72, 81, 3, 81, 81, 82, 74,
                 82, 39, 82, 77, 82, 80, 82, 30, 82, 29, 82, 7, 82, 53, 83, 81, 83, 69, 83, 73,
                 83, 46, 83, 67, 83, 49, 83, 83, 84, 24, 84, 49, 84, 52, 84, 3, 84, 74, 84, 10,
                 84, 81, 84, 5, 84, 3, 84, 6, 85, 14, 85, 38, 85, 43, 85, 80, 85, 12, 85, 26, 85,
                 31, 85, 44, 86, 53, 86, 75, 86, 57, 86, 48, 86, 80, 86, 66, 86, 86, 87, 17, 87,
                 62, 87, 56, 87, 24, 87, 20, 87, 65, 87, 49, 88, 58, 88, 83, 88, 69, 88, 46, 88,
                 53, 88, 73, 88, 67, 88, 88, 89, 1, 89, 37, 89, 25, 89, 33, 89, 55, 89, 45, 89,
                 5, 90, 8, 90, 23, 90, 0, 90, 11, 90, 50, 90, 24, 90, 69, 90, 28, 90, 29, 91, 48,
                 91, 66, 91, 69, 91, 44, 91, 86, 91, 57, 91, 80, 91, 91, 92, 35, 92, 15, 92, 86,
                 92, 48, 92, 57, 92, 61, 92, 66, 92, 75, 92, 0, 93, 23, 93, 80, 93, 16, 93, 4,
                 93, 82, 93, 91, 93, 41, 93, 9, 93, 34, 94, 19, 94, 55, 94, 79, 94, 80, 94, 29,
                 94, 30, 94, 82, 94, 35, 94, 70, 95, 69, 95, 76, 95, 62, 95, 56, 95, 27, 95, 17,
                 95, 87, 95, 37, 95, 48, 96, 17, 96, 76, 96, 27, 96, 56, 96, 65, 96, 20, 96, 87,
                 96, 5, 97, 86, 97, 58, 97, 11, 97, 59, 97, 63, 97, 97, 98, 77, 98, 48, 98, 84,
                 98, 40, 98, 10, 98, 5, 98, 52, 98, 81, 98, 89, 99, 34, 99, 14, 99, 85, 99, 54,
                 99, 18, 99, 31, 99, 61, 99, 71, 99, 14, 99, 99, 100, 82, 100, 13, 100, 2, 100,
                 15, 100, 32, 100, 64, 100, 47, 100, 39, 100, 6, 100, 51, 101, 30, 101, 94,
                 101, 1, 101, 79, 101, 58, 101, 19, 101, 55, 101, 35, 101, 29, 101, 100, 102,
                 74, 102, 52, 102, 98, 102, 72, 102, 40, 102, 10, 102, 3, 102, 102, 103, 33,
                 103, 45, 103, 25, 103, 89, 103, 37, 103, 1, 103, 70, 103, 72, 104, 11, 104,
                 0, 104, 93, 104, 67, 104, 41, 104, 16, 104, 87, 104, 23, 104, 4, 104, 9, 104,
                 89, 105, 103, 105, 33, 105, 62, 105, 37, 105, 45, 105, 1, 105, 80, 105, 25,
                 105, 25, 106, 56, 106, 92, 106, 2, 106, 13, 106, 32, 106, 60, 106, 6, 106,
                 64, 106, 15, 106, 39, 106, 88, 107, 75, 107, 98, 107, 102, 107, 72, 107, 40,
                 107, 81, 107, 5, 107, 10, 107, 84, 107, 4, 108, 9, 108, 7, 108, 51, 108, 77,
                 108, 21, 108, 78, 108, 22, 108, 68, 108, 79, 109, 30, 109, 63, 109, 1, 109,
                 33, 109, 103, 109, 105, 109, 45, 109, 25, 109, 89, 109, 37, 109, 67, 110,
                 13, 110, 24, 110, 80, 110, 88, 110, 49, 110, 73, 110, 46, 110, 83, 110, 53,
                 110, 23, 111, 64, 111, 46, 111, 78, 111, 8, 111, 21, 111, 51, 111, 7, 111,
                 108, 111, 68, 111, 77, 111, 52, 112, 96, 112, 97, 112, 57, 112, 66, 112, 63,
                 112, 44, 112, 92, 112, 75, 112, 91, 112, 28, 113, 20, 113, 95, 113, 59, 113,
                 70, 113, 17, 113, 87, 113, 76, 113, 65, 113, 96, 113, 83, 114, 88, 114, 110,
                 114, 53, 114, 49, 114, 73, 114, 46, 114, 67, 114, 58, 114, 15, 114, 104, 114,
                 -1);
    igraph_simplify(&g, /*multiple=*/ 1, /*loops=*/ 1, /*edge_comb=*/ 0);
    igraph_vector_view(&types, football_types,
                       sizeof(football_types) / sizeof(igraph_real_t));
    igraph_assortativity_nominal(&g, &types, &res, /*directed=*/ 0);
    printf("%.5f\n", res);

    igraph_destroy(&g);

    return 0;
}

