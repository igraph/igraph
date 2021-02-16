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

#include "igraph_error.h"
#include "igraph_blas.h"

#include "linalg/blas_internal.h"

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
 * \sa \ref igraph_blas_dgemv_array if you have arrays instead of
 *     vectors.
 *
 * \example examples/simple/blas.c
 */
void igraph_blas_dgemv(igraph_bool_t transpose, igraph_real_t alpha,
                       const igraph_matrix_t* a, const igraph_vector_t* x,
                       igraph_real_t beta, igraph_vector_t* y) {
    char trans = transpose ? 'T' : 'N';
    int m, n;
    int inc = 1;

    m = (int) igraph_matrix_nrow(a);
    n = (int) igraph_matrix_ncol(a);

    IGRAPH_ASSERT(igraph_vector_size(x) == transpose ? m : n);
    IGRAPH_ASSERT(igraph_vector_size(y) == transpose ? n : m);

    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 VECTOR(*x), &inc, &beta, VECTOR(*y), &inc);
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
 * \sa \ref igraph_blas_dgemv if you have vectors instead of
 *     arrays.
 */
void igraph_blas_dgemv_array(igraph_bool_t transpose, igraph_real_t alpha,
                             const igraph_matrix_t* a, const igraph_real_t* x,
                             igraph_real_t beta, igraph_real_t* y) {
    char trans = transpose ? 'T' : 'N';
    int m, n;
    int inc = 1;

    m = (int) igraph_matrix_nrow(a);
    n = (int) igraph_matrix_ncol(a);

    igraphdgemv_(&trans, &m, &n, &alpha, VECTOR(a->data), &m,
                 (igraph_real_t*)x, &inc, &beta, y, &inc);
}

igraph_real_t igraph_blas_dnrm2(const igraph_vector_t *v) {
    int n = igraph_vector_size(v);
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
int igraph_blas_ddot(const igraph_vector_t *v1, const igraph_vector_t *v2,
                       igraph_real_t *res) {

    int n = igraph_vector_size(v1);
    int one = 1;

    if (igraph_vector_size(v2) != n) {
        IGRAPH_ERROR("Dot product of vectors with different dimensions.",
                     IGRAPH_EINVAL);
    }

    *res = igraphddot_(&n, VECTOR(*v1), &one, VECTOR(*v2), &one);

    return 0;
}
