/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_lapack.h"

#include "linalg/lapack_internal.h"

/**
 * \function igraph_lapack_dgetrf
 * \brief LU factorization of a general M-by-N matrix.
 *
 * The factorization has the form
 *      A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if m > n), and U is upper
 * triangular (upper trapezoidal if m &lt; n).
 * \param a The input/output matrix. On entry, the M-by-N matrix to be
 *      factored. On exit, the factors L and U from the factorization
 *      A = P * L * U; the unit diagonal elements of L are not
 *      stored.
 * \param ipiv An integer vector, the pivot indices are stored here,
 *      unless it is a null pointer. Row \c i of the matrix was
 *      interchanged with row <code>ipiv[i]</code>.
 * \param info LAPACK error code. Zero on successful exit. If its value is
 *      a positive number i, it indicates that U(i,i) is exactly zero.
 *      The factorization has been
 *      completed, but the factor U is exactly singular, and division
 *      by zero will occur if it is used to solve a system of
 *      equations. If LAPACK returns an error, i.e. a negative info
 *      value, then an igraph error is generated as well.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_lapack_dgetrf(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
                         int *info) {
    int m = (int) igraph_matrix_nrow(a);
    int n = (int) igraph_matrix_ncol(a);
    int lda = m > 0 ? m : 1;
    igraph_vector_int_t *myipiv = ipiv, vipiv;

    if (!ipiv) {
        IGRAPH_CHECK(igraph_vector_int_init(&vipiv, m < n ? m : n));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &vipiv);
        myipiv = &vipiv;
    } else {
        IGRAPH_CHECK(igraph_vector_int_resize(ipiv, m < n ? m : n));
    }

    igraphdgetrf_(&m, &n, VECTOR(a->data), &lda, VECTOR(*myipiv), info);

    if (*info > 0) {
        IGRAPH_WARNING("LU: factor is exactly singular.");
    } else if (*info < 0) {
        switch (*info) {
        case -1:
            IGRAPH_ERROR("Invalid number of rows.", IGRAPH_ELAPACK);
            break;
        case -2:
            IGRAPH_ERROR("Invalid number of columns.", IGRAPH_ELAPACK);
            break;
        case -3:
            IGRAPH_ERROR("Invalid input matrix.", IGRAPH_ELAPACK);
            break;
        case -4:
            IGRAPH_ERROR("Invalid LDA parameter.", IGRAPH_ELAPACK);
            break;
        case -5:
            IGRAPH_ERROR("Invalid pivot vector.", IGRAPH_ELAPACK);
            break;
        case -6:
            IGRAPH_ERROR("Invalid info argument.", IGRAPH_ELAPACK);
            break;
        default:
            IGRAPH_ERROR("Unknown LAPACK error.", IGRAPH_ELAPACK);
            break;
        }
    }

    if (!ipiv) {
        igraph_vector_int_destroy(&vipiv);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_lapack_dgetrs
 * \brief Solve general system of linear equations using LU factorization.
 *
 * This function calls LAPACK to solve a system of linear equations
 *      A * X = B  or  A' * X = B
 * with a general N-by-N matrix A using the LU factorization
 * computed by \ref igraph_lapack_dgetrf.
 * \param transpose Logical scalar, whether to transpose the input
 *      matrix.
 * \param a A matrix containing the L and U factors from the
 *      factorization A = P*L*U. L is expected to be unitriangular,
 *      diagonal entries are those of U. If A is singular, no warning or
 *      error wil be given and random output will be returned.
 * \param ipiv An integer vector, the pivot indices from \ref
 *      igraph_lapack_dgetrf() must be given here. Row \c i of A was
 *      interchanged with row <code>ipiv[i]</code>.
 * \param b The right hand side matrix must be given here. The solution
            will also be placed here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_lapack_dgetrs(igraph_bool_t transpose, const igraph_matrix_t *a,
                         const igraph_vector_int_t *ipiv, igraph_matrix_t *b) {
    char trans = transpose ? 'T' : 'N';
    int n = (int) igraph_matrix_nrow(a);
    int nrhs = (int) igraph_matrix_ncol(b);
    int lda = n > 0 ? n : 1;
    int ldb = n > 0 ? n : 1;
    int info;

    if (n != igraph_matrix_ncol(a)) {
        IGRAPH_ERROR("Cannot LU solve matrix.", IGRAPH_NONSQUARE);
    }
    if (n != igraph_matrix_nrow(b)) {
        IGRAPH_ERROR("Cannot LU solve matrix, RHS of wrong size.", IGRAPH_EINVAL);
    }
    if (igraph_vector_int_size(ipiv) > 0) {
        igraph_integer_t min, max;
        igraph_vector_int_minmax(ipiv, &min, &max);
        if (max > n || min < 1) {
            IGRAPH_ERROR("Pivot index out of range.", IGRAPH_EINVAL);
        }
    }
    if (igraph_vector_int_size(ipiv) != n) {
        IGRAPH_ERROR("Pivot vector length must match number of matrix rows.", IGRAPH_EINVAL);
    }
    igraphdgetrs_(&trans, &n, &nrhs, VECTOR(a->data), &lda, VECTOR(*ipiv),
                  VECTOR(b->data), &ldb, &info);

    if (info < 0) {
        switch (info) {
        case -1:
            IGRAPH_ERROR("Invalid transpose argument.", IGRAPH_ELAPACK);
            break;
        case -2:
            IGRAPH_ERROR("Invalid number of rows/columns.", IGRAPH_ELAPACK);
            break;
        case -3:
            IGRAPH_ERROR("Invalid number of RHS vectors.", IGRAPH_ELAPACK);
            break;
        case -4:
            IGRAPH_ERROR("Invalid LU matrix.", IGRAPH_ELAPACK);
            break;
        case -5:
            IGRAPH_ERROR("Invalid LDA parameter.", IGRAPH_ELAPACK);
            break;
        case -6:
            IGRAPH_ERROR("Invalid pivot vector.", IGRAPH_ELAPACK);
            break;
        case -7:
            IGRAPH_ERROR("Invalid RHS matrix.", IGRAPH_ELAPACK);
            break;
        case -8:
            IGRAPH_ERROR("Invalid LDB parameter.", IGRAPH_ELAPACK);
            break;
        case -9:
            IGRAPH_ERROR("Invalid info argument.", IGRAPH_ELAPACK);
            break;
        default:
            IGRAPH_ERROR("Unknown LAPACK error.", IGRAPH_ELAPACK);
            break;
        }
    }

    return 0;
}

/**
 * \function igraph_lapack_dgesv
 * Solve system of linear equations with LU factorization
 *
 * This function computes the solution to a real system of linear
 * equations A * X = B, where A is an N-by-N matrix and X and B are
 * N-by-NRHS matrices.
 *
 * </para><para>The LU decomposition with partial pivoting and row
 * interchanges is used to factor A as
 *    A = P * L * U,
 * where P is a permutation matrix, L is unit lower triangular, and U is
 * upper triangular.  The factored form of A is then used to solve the
 * system of equations A * X = B.
 * \param a Matrix. On entry the N-by-N coefficient matrix, on exit,
 *        the factors L and U from the factorization A=P*L*U; the unit
 *        diagonal elements of L are not stored.
 * \param ipiv An integer vector or a null pointer. If not a null
 *        pointer, then the pivot indices that define the permutation
 *        matrix P, are stored here. Row i of the matrix was
 *        interchanged with row IPIV(i).
 * \param b Matrix, on entry the right hand side matrix should be
 *        stored here. On exit, if there was no error, and the info
 *        argument is zero, then it contains the solution matrix X.
 * \param info The LAPACK info code. If it is positive, then
 *        U(info,info) is exactly zero. In this case the factorization
 *        has been completed, but the factor U is exactly
 *        singular, so the solution could not be computed.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_lapack_dgesv.c
 */

int igraph_lapack_dgesv(igraph_matrix_t *a, igraph_vector_int_t *ipiv,
                        igraph_matrix_t *b, int *info) {

    int n = (int) igraph_matrix_nrow(a);
    int nrhs = (int) igraph_matrix_ncol(b);
    int lda = n > 0 ? n : 1;
    int ldb = n > 0 ? n : 1;
    igraph_vector_int_t *myipiv = ipiv, vipiv;

    if (n != igraph_matrix_ncol(a)) {
        IGRAPH_ERROR("Cannot LU solve matrix.", IGRAPH_NONSQUARE);
    }
    if (n != igraph_matrix_nrow(b)) {
        IGRAPH_ERROR("Cannot LU solve matrix, RHS of wrong size.", IGRAPH_EINVAL);
    }

    if (!ipiv) {
        IGRAPH_CHECK(igraph_vector_int_init(&vipiv, n));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &vipiv);
        myipiv = &vipiv;
    }

    igraphdgesv_(&n, &nrhs, VECTOR(a->data), &lda, VECTOR(*myipiv),
                 VECTOR(b->data), &ldb, info);

    if (*info > 0) {
        IGRAPH_WARNING("LU: factor is exactly singular.");
    } else if (*info < 0) {
        switch (*info) {
        case -1:
            IGRAPH_ERROR("Invalid number of rows/column.", IGRAPH_ELAPACK);
            break;
        case -2:
            IGRAPH_ERROR("Invalid number of RHS vectors.", IGRAPH_ELAPACK);
            break;
        case -3:
            IGRAPH_ERROR("Invalid input matrix.", IGRAPH_ELAPACK);
            break;
        case -4:
            IGRAPH_ERROR("Invalid LDA parameter.", IGRAPH_ELAPACK);
            break;
        case -5:
            IGRAPH_ERROR("Invalid pivot vector.", IGRAPH_ELAPACK);
            break;
        case -6:
            IGRAPH_ERROR("Invalid RHS matrix.", IGRAPH_ELAPACK);
            break;
        case -7:
            IGRAPH_ERROR("Invalid LDB parameter.", IGRAPH_ELAPACK);
            break;
        case -8:
            IGRAPH_ERROR("Invalid info argument.", IGRAPH_ELAPACK);
            break;
        default:
            IGRAPH_ERROR("Unknown LAPACK error.", IGRAPH_ELAPACK);
            break;
        }
    }

    if (!ipiv) {
        igraph_vector_int_destroy(&vipiv);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_lapack_dsyevr
 * Selected eigenvalues and optionally eigenvectors of a symmetric matrix
 *
 * Calls the DSYEVR LAPACK function to compute selected eigenvalues
 * and, optionally, eigenvectors of a real symmetric matrix A.
 * Eigenvalues and eigenvectors can be selected by specifying either
 * a range of values or a range of indices for the desired eigenvalues.
 *
 * </para><para>See more in the LAPACK documentation.
 * \param A Matrix, on entry it contains the symmetric input
 *        matrix. Only the leading N-by-N upper triangular part is
 *        used for the computation.
 * \param which Constant that gives which eigenvalues (and possibly
 *        the corresponding eigenvectors) to calculate. Possible
 *        values are \c IGRAPH_LAPACK_DSYEV_ALL, all eigenvalues;
 *        \c IGRAPH_LAPACK_DSYEV_INTERVAL, all eigenvalues in the
 *        half-open interval (vl,vu];
 *        \c IGRAPH_LAPACK_DSYEV_SELECT, the il-th through iu-th
 *        eigenvalues.
 * \param vl If \p which is \c IGRAPH_LAPACK_DSYEV_INTERVAL, then
 *        this is the lower bound of the interval to be searched for
 *        eigenvalues. See also the \p vestimate argument.
 * \param vu If \p which is \c IGRAPH_LAPACK_DSYEV_INTERVAL, then
 *        this is the upper bound of the interval to be searched for
 *        eigenvalues. See also the \p vestimate argument.
 * \param vestimate An upper bound for the number of eigenvalues in
 *        the (vl,vu] interval, if \p which is \c
 *        IGRAPH_LAPACK_DSYEV_INTERVAL. Memory is allocated only for
 *        the given number of eigenvalues (and eigenvectors), so this
 *        upper bound must be correct.
 * \param il The index of the smallest eigenvalue to return, if \p
 *        which is \c IGRAPH_LAPACK_DSYEV_SELECT.
 * \param iu The index of the largets eigenvalue to return, if \p
 *        which is \c IGRAPH_LAPACK_DSYEV_SELECT.
 * \param abstol The absolute error tolerance for the eigevalues. An
 *        approximate eigenvalue is accepted as converged when it is
 *        determined to lie in an interval [a,b] of width less than or
 *        equal to abstol + EPS * max(|a|,|b|), where EPS is the
 *        machine precision.
 * \param values An initialized vector, the eigenvalues are stored
 *        here, unless it is a null pointer. It will be resized as
 *        needed.
 * \param vectors An initialized matrix, the eigenvectors are stored
 *        in its columns, unless it is a null pointer. It will be
 *        resized as needed.
 * \param support An integer vector. If not a null pointer, then it
 *        will be resized to (2*max(1,M)) (M is a the total number of
 *        eigenvalues found). Then the support of the eigenvectors in
 *        \p vectors is stored here, i.e., the indices
 *        indicating the nonzero elements in \p vectors.
 *        The i-th eigenvector is nonzero only in elements
 *        support(2*i-1) through support(2*i).
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_lapack_dsyevr.c
 */

int igraph_lapack_dsyevr(const igraph_matrix_t *A,
                         igraph_lapack_dsyev_which_t which,
                         igraph_real_t vl, igraph_real_t vu, int vestimate,
                         int il, int iu, igraph_real_t abstol,
                         igraph_vector_t *values, igraph_matrix_t *vectors,
                         igraph_vector_int_t *support) {

    igraph_matrix_t Acopy;
    char jobz = vectors ? 'V' : 'N', range, uplo = 'U';
    int n = (int) igraph_matrix_nrow(A), lda = n, ldz = n;
    int m, info;
    igraph_vector_t *myvalues = values, vvalues;
    igraph_vector_int_t *mysupport = support, vsupport;
    igraph_vector_t work;
    igraph_vector_int_t iwork;
    int lwork = -1, liwork = -1;

    if (n != igraph_matrix_ncol(A)) {
        IGRAPH_ERROR("Cannot find eigenvalues/vectors.", IGRAPH_NONSQUARE);
    }
    if (which == IGRAPH_LAPACK_DSYEV_INTERVAL &&
        (vestimate < 1 || vestimate > n)) {
        IGRAPH_ERROR("Estimated (upper bound) number of eigenvalues must be "
                     "between 1 and n.", IGRAPH_EINVAL);
    }
    if (which == IGRAPH_LAPACK_DSYEV_SELECT && iu - il < 0) {
        IGRAPH_ERROR("Invalid 'il' and/or 'iu' values.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
    IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);

    IGRAPH_VECTOR_INIT_FINALLY(&work, 1);
    IGRAPH_CHECK(igraph_vector_int_init(&iwork, 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &iwork);

    if (!values) {
        IGRAPH_VECTOR_INIT_FINALLY(&vvalues, 0);
        myvalues = &vvalues;
    }
    if (!support) {
        IGRAPH_CHECK(igraph_vector_int_init(&vsupport, 0));
        IGRAPH_FINALLY(igraph_vector_int_destroy, &vsupport);
        mysupport = &vsupport;
    }

    IGRAPH_CHECK(igraph_vector_resize(myvalues, n));

    switch (which) {
    case IGRAPH_LAPACK_DSYEV_ALL:
        range = 'A';
        IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2 * n));
        if (vectors) {
            IGRAPH_CHECK(igraph_matrix_resize(vectors, n, n));
        }
        break;
    case IGRAPH_LAPACK_DSYEV_INTERVAL:
        range = 'V';
        IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2 * vestimate));
        if (vectors) {
            IGRAPH_CHECK(igraph_matrix_resize(vectors, n, vestimate));
        }
        break;
    case IGRAPH_LAPACK_DSYEV_SELECT:
        range = 'I';
        IGRAPH_CHECK(igraph_vector_int_resize(mysupport, 2 * (iu - il + 1)));
        if (vectors) {
            IGRAPH_CHECK(igraph_matrix_resize(vectors, n, iu - il + 1));
        }
        break;
    }

    igraphdsyevr_(&jobz, &range, &uplo, &n, &MATRIX(Acopy, 0, 0), &lda,
                  &vl, &vu, &il, &iu, &abstol, &m, VECTOR(*myvalues),
                  vectors ? &MATRIX(*vectors, 0, 0) : 0, &ldz, VECTOR(*mysupport),
                  VECTOR(work), &lwork, VECTOR(iwork), &liwork, &info);

    if (info != 0) {
        IGRAPH_ERROR("Invalid argument to dsyevr in workspace query.", IGRAPH_EINVAL);
    }

    lwork = (int) VECTOR(work)[0];
    liwork = VECTOR(iwork)[0];
    IGRAPH_CHECK(igraph_vector_resize(&work, lwork));
    IGRAPH_CHECK(igraph_vector_int_resize(&iwork, liwork));

    igraphdsyevr_(&jobz, &range, &uplo, &n, &MATRIX(Acopy, 0, 0), &lda,
                  &vl, &vu, &il, &iu, &abstol, &m, VECTOR(*myvalues),
                  vectors ? &MATRIX(*vectors, 0, 0) : 0, &ldz, VECTOR(*mysupport),
                  VECTOR(work), &lwork, VECTOR(iwork), &liwork, &info);

    if (info != 0) {
        IGRAPH_ERROR("Invalid argument to dsyevr in calculation.", IGRAPH_EINVAL);
    }

    if (values) {
        IGRAPH_CHECK(igraph_vector_resize(values, m));
    }
    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, m));
    }
    if (support) {
        IGRAPH_CHECK(igraph_vector_int_resize(support, m));
    }

    if (!support) {
        igraph_vector_int_destroy(&vsupport);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!values) {
        igraph_vector_destroy(&vvalues);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&iwork);
    igraph_vector_destroy(&work);
    igraph_matrix_destroy(&Acopy);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_lapack_dgeev
 * Eigenvalues and optionally eigenvectors of a non-symmetric matrix
 *
 * This function calls LAPACK to compute, for an N-by-N real
 * nonsymmetric matrix A, the eigenvalues and, optionally, the left
 * and/or right eigenvectors.
 *
 * </para><para>
 * The right eigenvector v(j) of A satisfies
 *                    A * v(j) = lambda(j) * v(j)
 * where lambda(j) is its eigenvalue.
 * The left eigenvector u(j) of A satisfies
 *                u(j)**H * A = lambda(j) * u(j)**H
 * where u(j)**H denotes the conjugate transpose of u(j).
 *
 * </para><para>
 * The computed eigenvectors are normalized to have Euclidean norm
 * equal to 1 and largest component real.
 *
 * \param A matrix. On entry it contains the N-by-N input matrix.
 * \param valuesreal Pointer to an initialized vector, or a null
 *        pointer. If not a null pointer, then the real parts of the
 *        eigenvalues are stored here. The vector will be resized as
 *        needed.
 * \param valuesimag Pointer to an initialized vector, or a null
 *        pointer. If not a null pointer, then the imaginary parts of
 *        the eigenvalues are stored here. The vector will be resized
 *        as needed.
 * \param vectorsleft Pointer to an initialized matrix, or a null
 *        pointer. If not a null pointer, then the left eigenvectors
 *        are stored in the columns of the matrix. The matrix will be
 *        resized as needed.
 * \param vectorsright Pointer to an initialized matrix, or a null
 *        pointer. If not a null pointer, then the right eigenvectors
 *        are stored in the columns of the matrix. The matrix will be
 *        resized as needed.
 * \param info This argument is used for two purposes. As an input
 *        argument it gives whether an igraph error should be
 *        generated if the QR algorithm fails to compute all
 *        eigenvalues. If \p info is non-zero, then an error is
 *        generated, otherwise only a warning is given.
 *        On exit it contains the LAPACK error code.
 *        Zero means successful exit.
 *        A negative values means that some of the arguments had an
 *        illegal value, this always triggers an igraph error. An i
 *        positive  value means that the QR algorithm failed to
 *        compute all the eigenvalues, and no eigenvectors have been
 *        computed; element i+1:N of \p valuesreal and \p valuesimag
 *        contain eigenvalues which have converged. This case only
 *        generates an igraph error, if \p info was non-zero on entry.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \example examples/simple/igraph_lapack_dgeev.c
 */

int igraph_lapack_dgeev(const igraph_matrix_t *A,
                        igraph_vector_t *valuesreal,
                        igraph_vector_t *valuesimag,
                        igraph_matrix_t *vectorsleft,
                        igraph_matrix_t *vectorsright,
                        int *info) {

    char jobvl = vectorsleft  ? 'V' : 'N';
    char jobvr = vectorsright ? 'V' : 'N';
    int n = (int) igraph_matrix_nrow(A);
    int lda = n, ldvl = n, ldvr = n, lwork = -1;
    igraph_vector_t work;
    igraph_vector_t *myreal = valuesreal, *myimag = valuesimag, vreal, vimag;
    igraph_matrix_t Acopy;
    int error = *info;

    if (igraph_matrix_ncol(A) != n) {
        IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev).", IGRAPH_NONSQUARE);
    }

    IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
    IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);

    IGRAPH_VECTOR_INIT_FINALLY(&work, 1);

    if (!valuesreal) {
        IGRAPH_VECTOR_INIT_FINALLY(&vreal, n);
        myreal = &vreal;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(myreal, n));
    }
    if (!valuesimag) {
        IGRAPH_VECTOR_INIT_FINALLY(&vimag, n);
        myimag = &vimag;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(myimag, n));
    }
    if (vectorsleft) {
        IGRAPH_CHECK(igraph_matrix_resize(vectorsleft, n, n));
    }
    if (vectorsright) {
        IGRAPH_CHECK(igraph_matrix_resize(vectorsright, n, n));
    }

    igraphdgeev_(&jobvl, &jobvr, &n, &MATRIX(Acopy, 0, 0), &lda,
                 VECTOR(*myreal), VECTOR(*myimag),
                 vectorsleft  ? &MATRIX(*vectorsleft, 0, 0) : 0, &ldvl,
                 vectorsright ? &MATRIX(*vectorsright, 0, 0) : 0, &ldvr,
                 VECTOR(work), &lwork, info);

    lwork = (int) VECTOR(work)[0];
    IGRAPH_CHECK(igraph_vector_resize(&work, lwork));

    igraphdgeev_(&jobvl, &jobvr, &n, &MATRIX(Acopy, 0, 0), &lda,
                 VECTOR(*myreal), VECTOR(*myimag),
                 vectorsleft  ? &MATRIX(*vectorsleft, 0, 0) : 0, &ldvl,
                 vectorsright ? &MATRIX(*vectorsright, 0, 0) : 0, &ldvr,
                 VECTOR(work), &lwork, info);

    if (*info < 0) {
        IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev).", IGRAPH_ELAPACK);
    } else if (*info > 0) {
        if (error) {
            IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev).", IGRAPH_ELAPACK);
        } else {
            IGRAPH_WARNING("Cannot calculate eigenvalues (dgeev).");
        }
    }

    if (!valuesimag) {
        igraph_vector_destroy(&vimag);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!valuesreal) {
        igraph_vector_destroy(&vreal);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&work);
    igraph_matrix_destroy(&Acopy);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_lapack_dgeevx
 * Eigenvalues/vectors of nonsymmetric matrices, expert mode
 *
 * This function calculates the eigenvalues and optionally the left
 * and/or right eigenvectors of a nonsymmetric N-by-N real matrix.
 *
 * </para><para>
 * Optionally also, it computes a balancing transformation to improve
 * the conditioning of the eigenvalues and eigenvectors (\p ilo, \p ihi,
 * \p scale, and \p abnrm), reciprocal condition numbers for the
 * eigenvalues (\p rconde), and reciprocal condition numbers for the
 * right eigenvectors (\p rcondv).
 *
 * </para><para>
 * The right eigenvector v(j) of A satisfies
 *                   A * v(j) = lambda(j) * v(j)
 * where lambda(j) is its eigenvalue.
 * The left eigenvector u(j) of A satisfies
 *               u(j)^H * A = lambda(j) * u(j)^H
 * where u(j)^H denotes the conjugate transpose of u(j).
 *
 * </para><para>
 * The computed eigenvectors are normalized to have Euclidean norm
 * equal to 1 and largest component real.
 *
 * </para><para>
 * Balancing a matrix means permuting the rows and columns to make it
 * more nearly upper triangular, and applying a diagonal similarity
 * transformation D * A * D^(-1), where D is a diagonal matrix, to
 * make its rows and columns closer in norm and the condition numbers
 * of its eigenvalues and eigenvectors smaller.  The computed
 * reciprocal condition numbers correspond to the balanced matrix.
 * Permuting rows and columns will not change the condition numbers
 * (in exact arithmetic) but diagonal scaling will.  For further
 * explanation of balancing, see section 4.10.2 of the LAPACK
 * Users' Guide.
 *
 * \param balance Scalar that indicated, whether the input matrix
 *   should be balanced. Possible values:
 *   \clist
 *     \cli IGRAPH_LAPACK_DGEEVX_BALANCE_NONE
 *          no not diagonally scale or permute.
 *     \cli IGRAPH_LAPACK_DGEEVX_BALANCE_PERM
 *          perform permutations to make the matrix more nearly upper
 *          triangular. Do not diagonally scale.
 *     \cli IGRAPH_LAPACK_DGEEVX_BALANCE_SCALE
 *          diagonally scale the matrix, i.e. replace A by
 *          D*A*D^(-1), where D is a diagonal matrix, chosen to make
 *          the rows and columns of A more equal in norm. Do not
 *          permute.
 *     \cli IGRAPH_LAPACK_DGEEVX_BALANCE_BOTH
 *          both diagonally scale and permute A.
 *   \endclist
 * \param A The input matrix, must be square.
 * \param valuesreal An initialized vector, or a NULL pointer. If not
 *   a NULL pointer, then the real parts of the eigenvalues are stored
 *   here. The vector will be resized, as needed.
 * \param valuesimag An initialized vector, or a NULL pointer. If not
 *   a NULL pointer, then the imaginary parts of the eigenvalues are stored
 *   here. The vector will be resized, as needed.
 * \param vectorsleft An initialized matrix or a NULL pointer. If not
 *   a null pointer, then the left eigenvectors are stored here. The
 *   order corresponds to the eigenvalues and the eigenvectors are
 *   stored in a compressed form. If the j-th eigenvalue is real then
 *   column j contains the corresponding eigenvector. If the j-th and
 *   (j+1)-th eigenvalues form a complex conjugate pair, then the j-th
 *   and (j+1)-th columns contain their corresponding eigenvectors.
 * \param vectorsright An initialized matrix or a NULL pointer. If not
 *   a null pointer, then the right eigenvectors are stored here. The
 *   format is the same, as for the \p vectorsleft argument.
 * \param ilo
 * \param ihi \p ilo and \p ihi are integer values determined when A was
 *   balanced.  The balanced A(i,j) = 0 if I>J and
 *   J=1,...,ilo-1 or I=ihi+1,...,N.
 * \param scale Pointer to an initialized vector or a NULL pointer. If
 *   not a NULL pointer, then details of the permutations and scaling
 *   factors applied when balancing \p A, are stored here.
 *   If P(j) is the index of the row and column
 *   interchanged with row and column j, and D(j) is the scaling
 *   factor applied to row and column j, then
 *   \clist
 *      \cli scale(J) = P(J),    for J = 1,...,ilo-1
 *      \cli scale(J) = D(J),    for J = ilo,...,ihi
 *      \cli scale(J) = P(J)     for J = ihi+1,...,N.
 *   \endclist
 *   The order in which the interchanges are made is N to \p ihi+1,
 *   then 1 to \p ilo-1.
 * \param abnrm Pointer to a real variable, the one-norm of the
 *   balanced matrix is stored here. (The one-norm is the maximum of
 *   the sum of absolute values of elements in any column.)
 * \param rconde An initialized vector or a NULL pointer. If not a
 *   null pointer, then the reciprocal condition numbers of the
 *   eigenvalues are stored here.
 * \param rcondv An initialized vector or a NULL pointer. If not a
 *   null pointer, then the reciprocal condition numbers of the right
 *   eigenvectors are stored here.
 * \param info This argument is used for two purposes. As an input
 *        argument it gives whether an igraph error should be
 *        generated if the QR algorithm fails to compute all
 *        eigenvalues. If \p info is non-zero, then an error is
 *        generated, otherwise only a warning is given.
 *        On exit it contains the LAPACK error code.
 *        Zero means successful exit.
 *        A negative values means that some of the arguments had an
 *        illegal value, this always triggers an igraph error. An i
 *        positive  value means that the QR algorithm failed to
 *        compute all the eigenvalues, and no eigenvectors have been
 *        computed; element i+1:N of \p valuesreal and \p valuesimag
 *        contain eigenvalues which have converged. This case only
 *        generated an igraph error, if \p info was non-zero on entry.
 * \return Error code.
 *
 * Time complexity: TODO
 *
 * \example examples/simple/igraph_lapack_dgeevx.c
 */

int igraph_lapack_dgeevx(igraph_lapack_dgeevx_balance_t balance,
                         const igraph_matrix_t *A,
                         igraph_vector_t *valuesreal,
                         igraph_vector_t *valuesimag,
                         igraph_matrix_t *vectorsleft,
                         igraph_matrix_t *vectorsright,
                         int *ilo, int *ihi, igraph_vector_t *scale,
                         igraph_real_t *abnrm,
                         igraph_vector_t *rconde,
                         igraph_vector_t *rcondv,
                         int *info) {

    char balanc;
    char jobvl = vectorsleft  ? 'V' : 'N';
    char jobvr = vectorsright ? 'V' : 'N';
    char sense;
    int n = (int) igraph_matrix_nrow(A);
    int lda = n, ldvl = n, ldvr = n, lwork = -1;
    igraph_vector_t work;
    igraph_vector_int_t iwork;
    igraph_matrix_t Acopy;
    int error = *info;
    igraph_vector_t *myreal = valuesreal, *myimag = valuesimag, vreal, vimag;
    igraph_vector_t *myscale = scale, vscale;

    if (igraph_matrix_ncol(A) != n) {
        IGRAPH_ERROR("Cannot calculate eigenvalues (dgeevx).", IGRAPH_NONSQUARE);
    }

    switch (balance) {
    case IGRAPH_LAPACK_DGEEVX_BALANCE_NONE:
        balanc = 'N';
        break;
    case IGRAPH_LAPACK_DGEEVX_BALANCE_PERM:
        balanc = 'P';
        break;
    case IGRAPH_LAPACK_DGEEVX_BALANCE_SCALE:
        balanc = 'S';
        break;
    case IGRAPH_LAPACK_DGEEVX_BALANCE_BOTH:
        balanc = 'B';
        break;
    default:
        IGRAPH_ERROR("Invalid 'balance' argument.", IGRAPH_EINVAL);
        break;
    }

    if (!rconde && !rcondv) {
        sense = 'N';
    } else if (rconde && !rcondv) {
        sense = 'E';
    } else if (!rconde && rcondv) {
        sense = 'V';
    } else {
        sense = 'B';
    }

    IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
    IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);

    IGRAPH_VECTOR_INIT_FINALLY(&work, 1);
    IGRAPH_CHECK(igraph_vector_int_init(&iwork, n));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &iwork);

    if (!valuesreal) {
        IGRAPH_VECTOR_INIT_FINALLY(&vreal, n);
        myreal = &vreal;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(myreal, n));
    }
    if (!valuesimag) {
        IGRAPH_VECTOR_INIT_FINALLY(&vimag, n);
        myimag = &vimag;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(myimag, n));
    }
    if (!scale) {
        IGRAPH_VECTOR_INIT_FINALLY(&vscale, n);
        myscale = &vscale;
    } else {
        IGRAPH_CHECK(igraph_vector_resize(scale, n));
    }
    if (vectorsleft) {
        IGRAPH_CHECK(igraph_matrix_resize(vectorsleft, n, n));
    }
    if (vectorsright) {
        IGRAPH_CHECK(igraph_matrix_resize(vectorsright, n, n));
    }

    igraphdgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, &MATRIX(Acopy, 0, 0),
                  &lda, VECTOR(*myreal), VECTOR(*myimag),
                  vectorsleft  ? &MATRIX(*vectorsleft, 0, 0) : 0, &ldvl,
                  vectorsright ? &MATRIX(*vectorsright, 0, 0) : 0, &ldvr,
                  ilo, ihi, VECTOR(*myscale), abnrm,
                  rconde ? VECTOR(*rconde) : 0,
                  rcondv ? VECTOR(*rcondv) : 0,
                  VECTOR(work), &lwork, VECTOR(iwork), info);

    lwork = (int) VECTOR(work)[0];
    IGRAPH_CHECK(igraph_vector_resize(&work, lwork));

    igraphdgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, &MATRIX(Acopy, 0, 0),
                  &lda, VECTOR(*myreal), VECTOR(*myimag),
                  vectorsleft  ? &MATRIX(*vectorsleft, 0, 0) : 0, &ldvl,
                  vectorsright ? &MATRIX(*vectorsright, 0, 0) : 0, &ldvr,
                  ilo, ihi, VECTOR(*myscale), abnrm,
                  rconde ? VECTOR(*rconde) : 0,
                  rcondv ? VECTOR(*rcondv) : 0,
                  VECTOR(work), &lwork, VECTOR(iwork), info);

    if (*info < 0) {
        IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev).", IGRAPH_ELAPACK);
    } else if (*info > 0) {
        if (error) {
            IGRAPH_ERROR("Cannot calculate eigenvalues (dgeev).", IGRAPH_ELAPACK);
        } else {
            IGRAPH_WARNING("Cannot calculate eigenvalues (dgeev).");
        }
    }

    if (!scale) {
        igraph_vector_destroy(&vscale);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (!valuesimag) {
        igraph_vector_destroy(&vimag);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (!valuesreal) {
        igraph_vector_destroy(&vreal);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_int_destroy(&iwork);
    igraph_vector_destroy(&work);
    igraph_matrix_destroy(&Acopy);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

int igraph_lapack_dgehrd(const igraph_matrix_t *A,
                         int ilo, int ihi,
                         igraph_matrix_t *result) {

    int n = (int) igraph_matrix_nrow(A);
    int lda = n;
    int lwork = -1;
    igraph_vector_t work;
    igraph_real_t optwork;
    igraph_vector_t tau;
    igraph_matrix_t Acopy;
    int info = 0;
    int i;

    if (igraph_matrix_ncol(A) != n) {
        IGRAPH_ERROR("Hessenberg reduction failed.", IGRAPH_NONSQUARE);
    }

    if (ilo < 1 || ihi > n || ilo > ihi) {
        IGRAPH_ERROR("Invalid `ilo' and/or `ihi'.", IGRAPH_EINVAL);
    }

    if (n <= 1) {
        IGRAPH_CHECK(igraph_matrix_update(result, A));
        return 0;
    }

    IGRAPH_CHECK(igraph_matrix_copy(&Acopy, A));
    IGRAPH_FINALLY(igraph_matrix_destroy, &Acopy);
    IGRAPH_VECTOR_INIT_FINALLY(&tau, n - 1);

    igraphdgehrd_(&n, &ilo, &ihi, &MATRIX(Acopy, 0, 0), &lda, VECTOR(tau),
                  &optwork, &lwork, &info);

    if (info != 0) {
        IGRAPH_ERROR("Internal Hessenberg transformation error.",
                     IGRAPH_EINTERNAL);
    }

    lwork = (int) optwork;
    IGRAPH_VECTOR_INIT_FINALLY(&work, lwork);

    igraphdgehrd_(&n, &ilo, &ihi, &MATRIX(Acopy, 0, 0), &lda, VECTOR(tau),
                  VECTOR(work), &lwork, &info);

    if (info != 0) {
        IGRAPH_ERROR("Internal Hessenberg transformation error.",
                     IGRAPH_EINTERNAL);
    }

    igraph_vector_destroy(&work);
    igraph_vector_destroy(&tau);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_matrix_update(result, &Acopy));

    igraph_matrix_destroy(&Acopy);
    IGRAPH_FINALLY_CLEAN(1);

    for (i = 0; i < n - 2; i++) {
        int j;
        for (j = i + 2; j < n; j++) {
            MATRIX(*result, j, i) = 0.0;
        }
    }

    return 0;
}
