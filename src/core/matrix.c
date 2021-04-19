/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_types.h"
#include "igraph_matrix.h"

#define BASE_IGRAPH_REAL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_IGRAPH_REAL

#define BASE_INT
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_INT

#define BASE_LONG
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_LONG

#define BASE_CHAR
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_CHAR

#define BASE_BOOL
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_BOOL

#define BASE_COMPLEX
#include "igraph_pmt.h"
#include "matrix.pmt"
#include "igraph_pmt_off.h"
#undef BASE_COMPLEX

#ifndef USING_R
int igraph_matrix_complex_print(const igraph_matrix_complex_t *m) {

    long int nr = igraph_matrix_complex_nrow(m);
    long int nc = igraph_matrix_complex_ncol(m);
    long int i, j;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            igraph_complex_t z = MATRIX(*m, i, j);
            if (j != 0) {
                putchar(' ');
            }
            printf("%g%+gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
        }
        printf("\n");
    }

    return 0;
}
#endif

int igraph_matrix_complex_fprint(const igraph_matrix_complex_t *m,
                                 FILE *file) {

    long int nr = igraph_matrix_complex_nrow(m);
    long int nc = igraph_matrix_complex_ncol(m);
    long int i, j;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            igraph_complex_t z = MATRIX(*m, i, j);
            if (j != 0) {
                fputc(' ', file);
            }
            fprintf(file, "%g%+gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
        }
        fprintf(file, "\n");
    }

    return 0;
}

int igraph_matrix_complex_real(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *real) {
    long int nrow = igraph_matrix_complex_nrow(v);
    long int ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_real(&v->data, &real->data));
    return 0;
}

int igraph_matrix_complex_imag(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *imag) {
    long int nrow = igraph_matrix_complex_nrow(v);
    long int ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_imag(&v->data, &imag->data));
    return 0;
}

int igraph_matrix_complex_realimag(const igraph_matrix_complex_t *v,
                                   igraph_matrix_t *real,
                                   igraph_matrix_t *imag) {
    long int nrow = igraph_matrix_complex_nrow(v);
    long int ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_realimag(&v->data, &real->data,
                 &imag->data));
    return 0;
}

int igraph_matrix_complex_create(igraph_matrix_complex_t *v,
                                 const igraph_matrix_t *real,
                                 const igraph_matrix_t *imag) {
    IGRAPH_CHECK(igraph_vector_complex_create(&v->data, &real->data,
                 &imag->data));
    return 0;
}

int igraph_matrix_complex_create_polar(igraph_matrix_complex_t *v,
                                       const igraph_matrix_t *r,
                                       const igraph_matrix_t *theta) {
    IGRAPH_CHECK(igraph_vector_complex_create_polar(&v->data, &r->data,
                 &theta->data));
    return 0;
}

igraph_bool_t igraph_matrix_all_e_tol(const igraph_matrix_t *lhs,
                                      const igraph_matrix_t *rhs,
                                      igraph_real_t tol) {
    return igraph_vector_e_tol(&lhs->data, &rhs->data, tol);
}

int igraph_matrix_zapsmall(igraph_matrix_t *m, igraph_real_t tol) {
    return igraph_vector_zapsmall(&m->data, tol);
}
