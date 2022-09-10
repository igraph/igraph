/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
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

#include "igraph_blas.h"
#include "linalg/blas_internal.h"

#include <limits.h>

/**
 * \function igraph_blas_dgemv
 * \brief Matrix-vector multiplication using BLAS, vector version.
 *
 * This function is a somewhat more user-friendly interface to
 * the \c dgemv function in BLAS. \c dgemv performs the operation
 * y = alpha*A*x + beta*y, where x and y are vectors and A is an
 * appropriately sized matrix (symmetric or non-symmetric).
 *
 * \param transpose whether to transpose the matrix \p A
 * \param alpha     the constant \p alpha
 * \param a         the matrix \p A
 * \param x         the vector \p x
 * \param beta      the constant \p beta
 * \param y         the vector \p y (which will be modified in-place)
 *
 * Time complexity: O(nk) if the matrix is of size n x k
 *
 * \return \c IGRAPH_EOVERFLOW if the matrix is too large for BLAS,
 *         \c IGRAPH_SUCCESS otherwise.
 * \sa \ref igraph_blas_dgemv_array if you have arrays instead of
 *     vectors.
 *
 * \example examples/simple/blas.c
 */
igraph_error_t igraph_blas_dgemv(igraph_bool_t transpose, igraph_real_t alpha,
                       const igraph_matrix_t *a, const igraph_vector_t *x,
                       igraph_real_t beta, igraph_vector_t *y) {
    char trans = transpose ? 'T' : 'N';
    int m, n;
    int inc = 1;

    if (igraph_matrix_nrow(a) > INT_MAX || igraph_matrix_ncol(a) > INT_MAX) {
        IGRAPH_ERROR("Matrix too large for BLAS", IGRAPH_EOVERFLOW);
    }

    m = (int) igraph_matrix_nrow(a);
    n = (int) igraph_matrix_ncol(a);

    IGRAPH_ASSERT(igraph_vector_size(x) == transpose ? m : n);
    IGRAPH_ASSERT(igraph_vector_size(y) == transpose ? n : m);

#ifdef HAVE_GFORTRAN
    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 VECTOR(*x), &inc, &beta, VECTOR(*y), &inc, /* trans_len = */ 1);
#else
    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 VECTOR(*x), &inc, &beta, VECTOR(*y), &inc);
#endif

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_blas_dgemm
 * \brief Matrix-matrix multiplication using BLAS.
 *
 * This function is a somewhat more user-friendly interface to
 * the \c dgemm function in BLAS. \c dgemm calculates
 * alpha*a*b + beta*c, where a, b and c are matrices, of which a and b
 * can be transposed.
 *
 * \param transpose_a whether to transpose the matrix \p a
 * \param transpose_b whether to transpose the matrix \p b
 * \param alpha       the constant \c alpha
 * \param a           the matrix \c a
 * \param b           the matrix \c b
 * \param beta        the constant \c beta
 * \param c           the matrix \c c. The result will also be stored here.
 *                    If beta is zero, c will be resized to fit the result.
 *
 * Time complexity: O(n m k) where matrix a is of size n × k, and matrix b is of
 * size k × m.
 *
 * \return \c IGRAPH_EOVERFLOW if the matrix is too large for BLAS,
 *         \c IGRAPH_EINVAL if the matrices have incompatible sizes,
 *         \c IGRAPH_SUCCESS otherwise.
 *
 * \example examples/simple/blas_dgemm.c
 */
igraph_error_t igraph_blas_dgemm(igraph_bool_t transpose_a, igraph_bool_t transpose_b,
        igraph_real_t alpha, const igraph_matrix_t *a, const igraph_matrix_t *b,
        igraph_real_t beta, igraph_matrix_t *c) {
    char trans_a = transpose_a ? 'T' : 'N';
    char trans_b = transpose_b ? 'T' : 'N';
    int m, n, k, lda, ldb, ldc;
    igraph_integer_t nrow_oa = transpose_a ? igraph_matrix_ncol(a) : igraph_matrix_nrow(a);
    igraph_integer_t ncol_oa = transpose_a ? igraph_matrix_nrow(a) : igraph_matrix_ncol(a);
    igraph_integer_t nrow_ob = transpose_b ? igraph_matrix_ncol(b) : igraph_matrix_nrow(b);
    igraph_integer_t ncol_ob = transpose_b ? igraph_matrix_nrow(b) : igraph_matrix_ncol(b);

    if (ncol_oa != nrow_ob) {
        IGRAPH_ERRORF("%" IGRAPH_PRId "-by-%" IGRAPH_PRId " and %" IGRAPH_PRId "-by-%" IGRAPH_PRId
               " matrices cannot be multiplied, incompatible dimensions.", IGRAPH_EINVAL,
               nrow_oa, ncol_oa, nrow_ob, ncol_ob);
    }
    if (beta != 0 && (ncol_oa != igraph_matrix_ncol(c) || nrow_oa != igraph_matrix_nrow(c))) {
        IGRAPH_ERRORF("%" IGRAPH_PRId "-by-%" IGRAPH_PRId " and %" IGRAPH_PRId "-by-%" IGRAPH_PRId
               " matrices cannot be added, incompatible dimensions.", IGRAPH_EINVAL,
               nrow_oa, ncol_ob, igraph_matrix_nrow(c), igraph_matrix_ncol(c));
    }
    if (nrow_oa > INT_MAX || ncol_oa > INT_MAX) {
        IGRAPH_ERROR("Matrix A too large for BLAS.", IGRAPH_EOVERFLOW);
    }
    if (ncol_ob > INT_MAX) {
        IGRAPH_ERROR("Matrix B too large for BLAS.", IGRAPH_EOVERFLOW);
    }
    if (beta == 0) {
        IGRAPH_CHECK(igraph_matrix_resize(c, nrow_oa, ncol_ob));
    }

    m = (int) nrow_oa;
    k = (int) ncol_oa;
    n = (int) ncol_ob;
    lda = (int) igraph_matrix_nrow(a);
    ldb = (int) igraph_matrix_nrow(b);
    ldc = (int) igraph_matrix_nrow(c);


#ifdef HAVE_GFORTRAN
    igraphdgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, VECTOR(a->data),
                 &lda, VECTOR(b->data), &ldb, &beta, VECTOR(c->data), &ldc,
                 /*trans_a_len*/ 1, /*trans_b_len*/ 1);
#else
    igraphdgemm_(&trans_a, &trans_b, &m, &n, &k, &alpha, VECTOR(a->data),
                 &lda, VECTOR(b->data), &ldb, &beta, VECTOR(c->data), &ldc);
#endif

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_blas_dgemv_array
 * \brief Matrix-vector multiplication using BLAS, array version.
 *
 * This function is a somewhat more user-friendly interface to
 * the \c dgemv function in BLAS. \c dgemv performs the operation
 * y = alpha*A*x + beta*y, where x and y are vectors and A is an
 * appropriately sized matrix (symmetric or non-symmetric).
 *
 * \param transpose whether to transpose the matrix \p A
 * \param alpha     the constant \p alpha
 * \param a         the matrix \p A
 * \param x         the vector \p x as a regular C array
 * \param beta      the constant \p beta
 * \param y         the vector \p y as a regular C array
 *                  (which will be modified in-place)
 *
 * Time complexity: O(nk) if the matrix is of size n x k
 *
 * \return \c IGRAPH_EOVERFLOW if the matrix is too large for BLAS,
 *         \c IGRAPH_SUCCESS otherwise.
 *
 * \sa \ref igraph_blas_dgemv if you have vectors instead of
 *     arrays.
 */
igraph_error_t igraph_blas_dgemv_array(igraph_bool_t transpose, igraph_real_t alpha,
                             const igraph_matrix_t* a, const igraph_real_t* x,
                             igraph_real_t beta, igraph_real_t* y) {
    char trans = transpose ? 'T' : 'N';
    int m, n;
    int inc = 1;

    if (igraph_matrix_nrow(a) > INT_MAX || igraph_matrix_ncol(a) > INT_MAX) {
        IGRAPH_ERROR("Matrix too large for BLAS", IGRAPH_EOVERFLOW);
    }

    m = (int) igraph_matrix_nrow(a);
    n = (int) igraph_matrix_ncol(a);

#ifdef HAVE_GFORTRAN
    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 (igraph_real_t*)x, &inc, &beta, y, &inc, /* trans_len = */ 1);
#else
    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 (igraph_real_t*)x, &inc, &beta, y, &inc);
#endif

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_blas_dnrm2
 * \brief Euclidean norm of a vector.
 *
 * \param v The vector.
 * \return Real value, the norm of \p v.
 *
 * Time complexity: O(n) where n is the length of the vector.
 */
igraph_real_t igraph_blas_dnrm2(const igraph_vector_t *v) {
    if (igraph_vector_size(v) > INT_MAX) {
        IGRAPH_ERROR("Vector too large for BLAS", IGRAPH_EOVERFLOW);
    }

    int n = (int) igraph_vector_size(v);
    int one = 1;
    return igraphdnrm2_(&n, VECTOR(*v), &one);
}

/**
 * \function igraph_blas_ddot
 * \brief Dot product of two vectors.
 *
 * \param v1 The first vector.
 * \param v2 The second vector.
 * \param res Pointer to a real, the result will be stored here.
 *
 * Time complexity: O(n) where n is the length of the vectors.
 *
 * \example examples/simple/blas.c
 */
igraph_error_t igraph_blas_ddot(const igraph_vector_t *v1, const igraph_vector_t *v2,
                       igraph_real_t *res) {

    if (igraph_vector_size(v1) > INT_MAX) {
        IGRAPH_ERROR("Vector too large for BLAS", IGRAPH_EOVERFLOW);
    }

    int n = (int) igraph_vector_size(v1);
    int one = 1;

    if (igraph_vector_size(v2) != n) {
        IGRAPH_ERROR("Dot product of vectors with different dimensions.",
                     IGRAPH_EINVAL);
    }

    *res = igraphddot_(&n, VECTOR(*v1), &one, VECTOR(*v2), &one);

    return IGRAPH_SUCCESS;
}
