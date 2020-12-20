/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2013  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.inc"

void test_density(const igraph_t *graph, igraph_bool_t loops) {
    igraph_real_t density;

    if (igraph_density(graph, &density, loops)) {
        printf("FAILED!\n");
        return;
    }

    if (igraph_is_nan(density)) {
        printf("nan\n");
    } else {
        printf("%.4f\n", density);
    }
}

int main() {

    igraph_t g;
    igraph_vector_t v;

    igraph_vector_init(&v, 0);

    /* Test graphs with no vertices and no edges */
    igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with one vertex and no edges */
    igraph_create(&g, &v, 1, IGRAPH_UNDIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with one vertex and a loop edge */
    igraph_vector_resize(&v, 2);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 0;
    igraph_create(&g, &v, 1, IGRAPH_UNDIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with one vertex and two loop edges */
    igraph_vector_resize(&v, 4);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 0;
    VECTOR(v)[2] = 0;
    VECTOR(v)[3] = 0;
    igraph_create(&g, &v, 1, IGRAPH_UNDIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with two vertices and one edge between them */
    igraph_vector_resize(&v, 2);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    igraph_create(&g, &v, 2, IGRAPH_UNDIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with two vertices, one edge between them and a loop on one
     * of them */
    igraph_vector_resize(&v, 4);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 1;
    igraph_create(&g, &v, 2, IGRAPH_UNDIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Test graphs with two vertices, one edge between them and a loop on both
     * of them */
    igraph_vector_resize(&v, 6);
    VECTOR(v)[0] = 0;
    VECTOR(v)[1] = 1;
    VECTOR(v)[2] = 1;
    VECTOR(v)[3] = 1;
    VECTOR(v)[4] = 0;
    VECTOR(v)[5] = 0;
    igraph_create(&g, &v, 2, IGRAPH_UNDIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    igraph_create(&g, &v, 1, IGRAPH_DIRECTED);
    test_density(&g, 1);
    igraph_destroy(&g);
    printf("======\n");

    /* Zachary karate club graph */
    igraph_famous(&g, "zachary");
    test_density(&g, 0);
    test_density(&g, 1);
    igraph_destroy(&g);

    igraph_vector_destroy(&v);

    VERIFY_FINALLY_STACK();

    return 0;
}
