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
igraph_error_t igraph_matrix_complex_print(const igraph_matrix_complex_t *m) {

    igraph_integer_t nr = igraph_matrix_complex_nrow(m);
    igraph_integer_t nc = igraph_matrix_complex_ncol(m);
    igraph_integer_t i, j;
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

    return IGRAPH_SUCCESS;
}
#endif

igraph_error_t igraph_matrix_complex_fprint(const igraph_matrix_complex_t *m,
                                 FILE *file) {

    igraph_integer_t nr = igraph_matrix_complex_nrow(m);
    igraph_integer_t nc = igraph_matrix_complex_ncol(m);
    igraph_integer_t i, j;
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

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_real
 * \brief Gives the real part of a complex matrix.
 *
 * \param v Pointer to a complex matrix.
 * \param real Pointer to an initialized matrix. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_real(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *real) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(v);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_real(&v->data, &real->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_imag
 * \brief Gives the imaginary part of a complex matrix.
 *
 * \param v Pointer to a complex matrix.
 * \param imag Pointer to an initialized matrix. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_imag(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *imag) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(v);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_imag(&v->data, &imag->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_realimag
 * \brief Gives the real and imaginary parts of a complex matrix.
 *
 * \param v Pointer to a complex matrix.
 * \param real Pointer to an initialized matrix. The real part will be stored here.
 * \param imag Pointer to an initialized matrix. The imaginary part will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_realimag(const igraph_matrix_complex_t *v,
                                   igraph_matrix_t *real,
                                   igraph_matrix_t *imag) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(v);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(v);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_realimag(&v->data, &real->data,
                 &imag->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_create
 * \brief Creates a complex matrix from a real and imaginary part.
 *
 * \param v Pointer to an uninitialized complex matrix.
 * \param real Pointer to the real part of the complex matrix.
 * \param imag Pointer to the imaginary part of the complex matrix.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_create(igraph_matrix_complex_t *v,
                                 const igraph_matrix_t *real,
                                 const igraph_matrix_t *imag) {
    igraph_integer_t nrowr = igraph_matrix_nrow(real);
    igraph_integer_t ncolr = igraph_matrix_ncol(real);
    igraph_integer_t nrowi = igraph_matrix_nrow(imag);
    igraph_integer_t ncoli = igraph_matrix_ncol(imag);

    if (nrowr != nrowi) {
        IGRAPH_ERRORF("Number of rows in real matrix (%" IGRAPH_PRId
                ") not equal to number of rows in imaginary matrix (%"
                IGRAPH_PRId ").", IGRAPH_EINVAL, nrowr, nrowi);
    }
    if (ncolr != ncoli) {
        IGRAPH_ERRORF("Number of columns in real matrix (%" IGRAPH_PRId
                ") not equal to number of columns in imaginary matrix (%"
                IGRAPH_PRId ").", IGRAPH_EINVAL, ncolr, ncoli);
    }

    igraph_matrix_complex_init(v, nrowr, ncolr);

    for (igraph_integer_t i = 0; i < nrowr * ncolr; i++) {
        VECTOR(v->data)[i] = igraph_complex(VECTOR(real->data)[i], VECTOR(imag->data)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_create_polar
 * \brief Creates a complex matrix from a magnitude and an angle.
 *
 * \param v Pointer to an uninitialized complex matrix.
 * \param r Pointer to the magnitude of the complex matrix.
 * \param theta Pointer to the angle of the complex matrix.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_create_polar(igraph_matrix_complex_t *v,
                                       const igraph_matrix_t *r,
                                       const igraph_matrix_t *theta) {
    igraph_integer_t nrowr = igraph_matrix_nrow(r);
    igraph_integer_t ncolr = igraph_matrix_ncol(r);
    igraph_integer_t nrowi = igraph_matrix_nrow(theta);
    igraph_integer_t ncoli = igraph_matrix_ncol(theta);

    if (nrowr != nrowi) {
        IGRAPH_ERRORF("Number of rows in magnitude matrix (%" IGRAPH_PRId
                ") not equal to number of rows in angle matrix (%"
                IGRAPH_PRId ").", IGRAPH_EINVAL, nrowr, nrowi);
    }
    if (ncolr != ncoli) {
        IGRAPH_ERRORF("Number of columns in magnitude matrix (%" IGRAPH_PRId
                ") not equal to number of columns in angle matrix (%"
                IGRAPH_PRId ").", IGRAPH_EINVAL, ncolr, ncoli);
    }

    igraph_matrix_complex_init(v, nrowr, ncolr);

    for (igraph_integer_t i = 0; i < nrowr * ncolr; i++) {
        VECTOR(v->data)[i] = igraph_complex_polar(VECTOR(r->data)[i], VECTOR(theta->data)[i]);
    }
    return IGRAPH_SUCCESS;
}

igraph_bool_t igraph_matrix_all_e_tol(const igraph_matrix_t *lhs,
                                      const igraph_matrix_t *rhs,
                                      igraph_real_t tol) {
    return igraph_vector_e_tol(&lhs->data, &rhs->data, tol);
}

igraph_error_t igraph_matrix_zapsmall(igraph_matrix_t *m, igraph_real_t tol) {
    return igraph_vector_zapsmall(&m->data, tol);
}
