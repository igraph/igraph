/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

int main() {
    igraph_matrix_t m;
    igraph_vector_t x, y;

    igraph_vector_init_real(&x, 3, 1.0, 2.0, 3.0);
    igraph_vector_init_real(&y, 4, 4.0, 5.0, 6.0, 7.0);

    igraph_matrix_init(&m, 4, 3);
    MATRIX(m, 0, 0) = 1;
    MATRIX(m, 0, 1) = 2;
    MATRIX(m, 0, 2) = 3;
    MATRIX(m, 1, 0) = 2;
    MATRIX(m, 1, 1) = 3;
    MATRIX(m, 1, 2) = 4;
    MATRIX(m, 2, 0) = 3;
    MATRIX(m, 2, 1) = 4;
    MATRIX(m, 2, 2) = 5;
    MATRIX(m, 3, 0) = 4;
    MATRIX(m, 3, 1) = 5;
    MATRIX(m, 3, 2) = 6;

    igraph_blas_dgemv(0, 2, &m, &x, 3, &y);
    igraph_vector_print(&y);

    igraph_vector_destroy(&x);
    igraph_vector_destroy(&y);
    igraph_matrix_destroy(&m);

    return 0;
}

