/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int check_convex_hull(igraph_matrix_t* coords) {
    igraph_vector_t result;
    igraph_matrix_t resmat;
    long int i;

    /* Testing with index output mode */
    igraph_vector_init(&result, 1);
    if (igraph_convex_hull(coords, &result, 0)) {
        return 1;
    }

    for (i = 0; i < igraph_vector_size(&result); i++) {
        printf("%ld ", (long)VECTOR(result)[i]);
    }
    printf("\n");
    igraph_vector_destroy(&result);

    /* Testing with coordinate output mode */
    igraph_matrix_init(&resmat, 0, 0);
    if (igraph_convex_hull(coords, 0, &resmat)) {
        return 1;
    }

    for (i = 0; i < igraph_matrix_nrow(&resmat); i++) {
        printf("%.3f %.3f ", MATRIX(resmat, i, 0), MATRIX(resmat, i, 1));
    }
    printf("\n");

    igraph_matrix_destroy(&resmat);

    return 0;
}

int test_simple() {
    igraph_real_t coords_array[][2] = {
        {3, 2}, {5, 1}, {4, 4}, {6, 4}, {4, 3},
        {2, 5}, {1, 3}, {2, 4}, {6, 3}, {9, 2}
    };
    igraph_matrix_t coords;
    int i, result;

    printf("test_simple\n");

    igraph_matrix_init(&coords, 10, 2);
    for (i = 0; i < 20; i++) {
        MATRIX(coords, i / 2, i % 2) = coords_array[i / 2][i % 2];
    }
    result = check_convex_hull(&coords);
    igraph_matrix_destroy(&coords);
    return result;
}

int test_collinear() {
    igraph_real_t coords_array[][2] =
    {{3, 2}, {5, 1}, {7, 0}, {9, -1}, {11, -2}};
    igraph_matrix_t coords;
    int i, result;

    printf("test_collinear\n");

    igraph_matrix_init(&coords, 5, 2);
    for (i = 0; i < 10; i++) {
        MATRIX(coords, i / 2, i % 2) = coords_array[i / 2][i % 2];
    }
    result = check_convex_hull(&coords);
    igraph_matrix_destroy(&coords);
    return result;
}

int test_degenerate() {
    igraph_matrix_t coords;
    int result;

    printf("test_degenerate\n");

    igraph_matrix_init(&coords, 2, 2);
    MATRIX(coords, 0, 0) = 3;
    MATRIX(coords, 0, 1) = 2;
    MATRIX(coords, 1, 0) = 5;
    MATRIX(coords, 1, 1) = 1;
    result = check_convex_hull(&coords);

    igraph_matrix_resize(&coords, 1, 2);
    MATRIX(coords, 0, 0) = 3;
    MATRIX(coords, 0, 1) = 2;
    result = check_convex_hull(&coords);

    igraph_matrix_resize(&coords, 0, 2);
    result = check_convex_hull(&coords);

    igraph_matrix_destroy(&coords);
    return result;
}

int test_bug_805() {
    igraph_real_t coords_array[][2] = {
        {0, 0}, {1, 0}, {0.707, 0.707}, {0, 1}, {-0.707, 0.707}, {-1, 0},
        {-0.707, -0.707}, {0, -1}, {0.707, -0.707}, {2, 0}, {1.414, 1.414}, {0, 2},
        {-1.414, 1.414}, {-2, 0}, {-1.414, -1.414}, {0, -2}, {1.414, -1.414}, {3, 0},
        {2.121, 2.121}, {0, 3}, {-2.121, 2.121}, {-3, 0}, {-2.121, -2.121}, {0, -3},
        {2.121, -2.121}, {4, 0}, {2.828, 2.828}, {0, 4}, {-2.828, 2.828}, {-4, 0},
        {-2.828, -2.828}, {0, -4}, {2.828, -2.828}
    };
    igraph_matrix_t coords;
    int i, result;

    printf("test_bug_805\n");

    igraph_matrix_init(&coords, 33, 2);
    for (i = 0; i < 66; i++) {
        MATRIX(coords, i / 2, i % 2) = coords_array[i / 2][i % 2];
    }
    result = check_convex_hull(&coords);
    igraph_matrix_destroy(&coords);
    return result;
}

int test_bug_1115() {
    igraph_real_t coords_array[][2] = {
        {37, 52}, {49, 49}, {52, 64}, {20, 26}, {40, 30}, {21, 47}, {17, 63}, {31, 62},
        {52, 33}, {51, 21}, {42, 41}, {31, 32}, {5, 25}, {12, 42}, {36, 16}, {52, 41},
        {27, 23}, {17, 33}, {13, 13}, {57, 58}, {62, 42}, {42, 57}, {16, 57}, {8, 52},
        {7, 38}, {27, 68}, {30, 48}, {43, 67}, {58, 48}, {58, 27}, {37, 69}, {38, 46},
        {46, 10}, {61, 33}, {62, 63}, {63, 69}, {32, 22}, {45, 35}, {59, 15}, {5, 6},
        {10, 17}, {21, 10}, {5, 64}, {30, 15}, {39, 10}, {32, 39}, {25, 32}, {25, 55},
        {48, 28}, {56, 37}, {30, 40}
    };
    igraph_matrix_t coords;
    int i, result;

    printf("test_bug_1115\n");

    igraph_matrix_init(&coords, 51, 2);
    for (i = 0; i < 102; i++) {
        MATRIX(coords, i / 2, i % 2) = coords_array[i / 2][i % 2];
    }
    result = check_convex_hull(&coords);
    igraph_matrix_destroy(&coords);
    return result;
}

int main() {
    int result;

    result = test_simple();
    if (result != 0) {
        return result;
    }

    result = test_collinear();
    if (result != 0) {
        return result;
    }

    result = test_degenerate();
    if (result != 0) {
        return result;
    }

    result = test_bug_805();
    if (result != 0) {
        return result;
    }

    result = test_bug_1115();
    if (result != 0) {
        return result;
    }

    return 0;
}
