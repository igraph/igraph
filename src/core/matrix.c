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

#include "igraph_matrix.h"
#include "igraph_types.h"

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

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_real
 * \brief Gives the real part of a complex matrix.
 *
 * \param m Pointer to a complex matrix.
 * \param real Pointer to an initialized matrix. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_real(const igraph_matrix_complex_t *m,
                               igraph_matrix_t *real) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(m);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(m);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_real(&m->data, &real->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_imag
 * \brief Gives the imaginary part of a complex matrix.
 *
 * \param m Pointer to a complex matrix.
 * \param imag Pointer to an initialized matrix. The result will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_imag(const igraph_matrix_complex_t *m,
                               igraph_matrix_t *imag) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(m);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(m);
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_imag(&m->data, &imag->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_realimag
 * \brief Gives the real and imaginary parts of a complex matrix.
 *
 * \param m Pointer to a complex matrix.
 * \param real Pointer to an initialized matrix. The real part will be stored here.
 * \param imag Pointer to an initialized matrix. The imaginary part will be stored here.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_realimag(const igraph_matrix_complex_t *m,
                                   igraph_matrix_t *real,
                                   igraph_matrix_t *imag) {
    igraph_integer_t nrow = igraph_matrix_complex_nrow(m);
    igraph_integer_t ncol = igraph_matrix_complex_ncol(m);
    IGRAPH_CHECK(igraph_matrix_resize(real, nrow, ncol));
    IGRAPH_CHECK(igraph_matrix_resize(imag, nrow, ncol));
    IGRAPH_CHECK(igraph_vector_complex_realimag(&m->data, &real->data,
                 &imag->data));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_create
 * \brief Creates a complex matrix from a real and imaginary part.
 *
 * \param m Pointer to an uninitialized complex matrix.
 * \param real Pointer to the real part of the complex matrix.
 * \param imag Pointer to the imaginary part of the complex matrix.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_create(igraph_matrix_complex_t *m,
                                 const igraph_matrix_t *real,
                                 const igraph_matrix_t *imag) {
    igraph_integer_t nrowr = igraph_matrix_nrow(real);
    igraph_integer_t ncolr = igraph_matrix_ncol(real);
    igraph_integer_t nrowi = igraph_matrix_nrow(imag);
    igraph_integer_t ncoli = igraph_matrix_ncol(imag);

    if (nrowr != nrowi || ncolr != ncoli) {
        IGRAPH_ERRORF("Dimensions of real (%" IGRAPH_PRId " by %" IGRAPH_PRId ") and "
                "imaginary (%" IGRAPH_PRId " by %" IGRAPH_PRId ") matrices must match.",
                IGRAPH_EINVAL, nrowr, ncolr, nrowi, ncoli);
    }

    IGRAPH_CHECK(igraph_matrix_complex_init(m, nrowr, ncolr));

    for (igraph_integer_t i = 0; i < nrowr * ncolr; i++) {
        VECTOR(m->data)[i] = igraph_complex(VECTOR(real->data)[i], VECTOR(imag->data)[i]);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup matrix
 * \function igraph_matrix_complex_create_polar
 * \brief Creates a complex matrix from a magnitude and an angle.
 *
 * \param m Pointer to an uninitialized complex matrix.
 * \param r Pointer to a real matrix containing magnitudes.
 * \param theta Pointer to a real matrix containing arguments (phase angles).
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the
 * number of elements in the matrix.
 */

igraph_error_t igraph_matrix_complex_create_polar(igraph_matrix_complex_t *m,
                                       const igraph_matrix_t *r,
                                       const igraph_matrix_t *theta) {
    igraph_integer_t nrowr = igraph_matrix_nrow(r);
    igraph_integer_t ncolr = igraph_matrix_ncol(r);
    igraph_integer_t nrowt = igraph_matrix_nrow(theta);
    igraph_integer_t ncolt = igraph_matrix_ncol(theta);

    if (nrowr != nrowt || ncolr != ncolt) {
        IGRAPH_ERRORF("Dimensions of magnitude (%" IGRAPH_PRId " by %" IGRAPH_PRId ") and "
                "angle (%" IGRAPH_PRId " by %" IGRAPH_PRId ") matrices must match.",
                IGRAPH_EINVAL, nrowr, ncolr, nrowt, ncolt);
    }

    IGRAPH_CHECK(igraph_matrix_complex_init(m, nrowr, ncolr));

    for (igraph_integer_t i = 0; i < nrowr * ncolr; i++) {
        VECTOR(m->data)[i] = igraph_complex_polar(VECTOR(r->data)[i], VECTOR(theta->data)[i]);
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_matrix_complex_all_almost_e
 * \brief Are all elements almost equal?
 *
 * Checks if the elements of two complex matrices are equal within a relative tolerance.
 *
 * \param lhs The first matrix.
 * \param rhs The second matrix.
 * \param eps Relative tolerance, see \ref igraph_complex_almost_equals() for details.
 * \return True if the two matrices are almost equal, false if there is at least
 *     one differing element or if the matrices are not of the same dimensions.
 */
igraph_bool_t igraph_matrix_complex_all_almost_e(igraph_matrix_complex_t *lhs,
                                                  igraph_matrix_complex_t *rhs,
                                                  igraph_real_t eps) {
    return lhs->ncol == rhs->ncol && lhs->nrow == rhs->nrow &&
            igraph_vector_complex_all_almost_e(&lhs->data, &rhs->data, eps);
}

/**
 * Deprecated in favour of \ref igraph_matrix_all_almost_e() which uses
 * relative tolerances. Will be removed in 0.11.
 *
 * Checks if two matrices are equal within an absolute tolerance.
 */
igraph_bool_t igraph_matrix_all_e_tol(const igraph_matrix_t *lhs,
                                      const igraph_matrix_t *rhs,
                                      igraph_real_t tol) {
    return lhs->ncol == rhs->ncol && lhs->nrow == rhs->nrow &&
            igraph_vector_e_tol(&lhs->data, &rhs->data, tol);
}


/**
 * \function igraph_matrix_all_almost_e
 * \brief Are all elements almost equal?
 *
 * Checks if the elements of two matrices are equal within a relative tolerance.
 *
 * \param lhs The first matrix.
 * \param rhs The second matrix.
 * \param eps Relative tolerance, see \ref igraph_almost_equals() for details.
 * \return True if the two matrices are almost equal, false if there is at least
 *     one differing element or if the matrices are not of the same dimensions.
 */
igraph_bool_t igraph_matrix_all_almost_e(const igraph_matrix_t *lhs,
                                         const igraph_matrix_t *rhs,
                                         igraph_real_t eps) {
    return lhs->ncol == rhs->ncol && lhs->nrow == rhs->nrow &&
            igraph_vector_all_almost_e(&lhs->data, &rhs->data, eps);
}


/**
 * \function igraph_matrix_zapsmall
 * \brief Replaces small elements of a matrix by exact zeros.
 *
 * Matrix elements which are smaller in magnitude than the given absolute
 * tolerance will be replaced by exact zeros. The default tolerance
 * corresponds to two-thirds of the representable digits of \type igraph_real_t,
 * i.e. <code>DBL_EPSILON^(2/3)</code> which is approximately <code>10^-10</code>.
 *
 * \param m   The matrix to process, it will be changed in-place.
 * \param tol Tolerance value. Numbers smaller than this in magnitude will
 *            be replaced by zeros. Pass in zero to use the default tolerance.
 *            Must not be negative.
 * \return Error code.
 *
 * \sa \ref igraph_matrix_all_almost_e() and \ref igraph_almost_equals() to
 * perform comparisons with relative tolerances.
 */
igraph_error_t igraph_matrix_zapsmall(igraph_matrix_t *m, igraph_real_t tol) {
    return igraph_vector_zapsmall(&m->data, tol);
}

/**
 * \function igraph_matrix_complex_zapsmall
 * \brief Replaces small elements of a complex matrix by exact zeros.
 *
 * Similarly to \ref igraph_matrix_zapsmall(), small elements will be replaced
 * by zeros. The operation is performed separately on the real and imaginary
 * parts of the numbers. This way, complex numbers with a large real part and
 * tiny imaginary part will effectively be transformed to real numbers.
 * The default tolerance
 * corresponds to two-thirds of the representable digits of \type igraph_real_t,
 * i.e. <code>DBL_EPSILON^(2/3)</code> which is approximately <code>10^-10</code>.
 *
 * \param m   The matrix to process, it will be changed in-place.
 * \param tol Tolerance value. Real and imaginary parts smaller than this in
 *            magnitude will be replaced by zeros. Pass in zero to use the default
 *            tolerance. Must not be negative.
 * \return Error code.
 *
 * \sa \ref igraph_matrix_complex_all_almost_e() and
 * \ref igraph_complex_almost_equals() to perform comparisons with relative
 * tolerances.
 */
igraph_error_t igraph_matrix_complex_zapsmall(igraph_matrix_complex_t *m, igraph_real_t tol) {
    return igraph_vector_complex_zapsmall(&m->data, tol);
}
