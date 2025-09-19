/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#define TOSTR1(x) #x
#define TOSTR(x) TOSTR1(x)

void bench(int m, int n, int rep) {
    igraph_matrix_t mat;

    igraph_matrix_init(&mat, m, n);
    for (igraph_int_t j=0; j < n; j++) {
        for (igraph_int_t i=0; i < m; i++) {
            MATRIX(mat, i, j) = RNG_UNIF(-1, 1);
        }
    }

    char name[200];
    sprintf(name, "Transpose %5d x %5d, %dx", m, n, rep);
    BENCH( name, REPEAT(igraph_matrix_transpose(&mat), rep) );

    igraph_matrix_destroy(&mat);
}

int main(void) {
    BENCH_INIT();

    bench(30, 30, 100000);
    bench(100, 100, 10000);
    bench(1000, 1000, 100);
    bench(1024, 1024, 100); /* naive implementation has bad cache behaviour with power of 2 sizes */
    bench(1023, 1025, 100); /* non-symmetric */
    bench(3000, 3000, 10);
    /* skinny non-symmetric: */
    bench(100, 10000, 100);
    bench(10000, 100, 100);

    return 0;
}
