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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_matrix.h"

#ifndef IGRAPH_ARPACK_H
#define IGRAPH_ARPACK_H

#include "igraph_decls.h"

__BEGIN_DECLS

/**
 * \section about_arpack ARPACK interface in igraph
 *
 * <para>
 * ARPACK is a library for solving large scale eigenvalue problems.
 * The package is designed to compute a few eigenvalues and corresponding
 * eigenvectors of a general \c n by \c n matrix \c A. It is
 * most appropriate for large sparse or structured matrices \c A where
 * structured means that a matrix-vector product <code>w &lt;- Av</code> requires
 * order \c n rather than the usual order <code>n^2</code> floating point
 * operations. Please see
 * http://www.caam.rice.edu/software/ARPACK/ for details.
 * </para>
 *
 * <para>
 * The eigenvalue calculation in ARPACK (in the simplest
 * case) involves the calculation of the \c Av product where \c A
 * is the matrix we work with and \c v is an arbitrary vector. A
 * user-defined function of type \ref igraph_arpack_function_t
 * is expected to perform this product. If the product can be done
 * efficiently, e.g. if the matrix is sparse, then ARPACK is usually
 * able to calculate the eigenvalues very quickly.
 * </para>
 *
 * <para>In igraph, eigenvalue/eigenvector calculations usually
 * involve the following steps:
 * \olist
 *   \oli Initialization of an \ref igraph_arpack_options_t data
 *        structure using \ref igraph_arpack_options_init.
 *   \oli Setting some options in the initialized \ref
 *        igraph_arpack_options_t object.
 *   \oli Defining a function of type \ref igraph_arpack_function_t.
 *        The input of this function is a vector, and the output
 *        should be the output matrix multiplied by the input vector.
 *   \oli Calling \ref igraph_arpack_rssolve() (is the matrix is
 *        symmetric), or \ref igraph_arpack_rnsolve().
 * \endolist
 * The \ref igraph_arpack_options_t object can be used multiple
 * times.
 * </para>
 *
 * <para>
 * If we have many eigenvalue problems to solve, then it might worth
 * to create an \ref igraph_arpack_storage_t object, and initialize it
 * via \ref igraph_arpack_storage_init(). This structure contains all
 * memory needed for ARPACK (with the given upper limit regerding to
 * the size of the eigenvalue problem). Then many problems can be
 * solved using the same \ref igraph_arpack_storage_t object, without
 * always reallocating the required memory.
 * The \ref igraph_arpack_storage_t object needs to be destroyed by
 * calling \ref igraph_arpack_storage_destroy() on it, when it is not
 * needed any more.
 * </para>
 *
 * <para>
 * igraph does not contain all
 * ARPACK routines, only the ones dealing with symmetric and
 * non-symmetric eigenvalue problems using double precision real
 * numbers.
 * </para>
 *
 */

/**
 * \struct igraph_arpack_options_t
 * \brief Options for ARPACK
 *
 * This data structure contains the options of thee ARPACK eigenvalue
 * solver routines. It must be initialized by calling \ref
 * igraph_arpack_options_init() on it. Then it can be used for
 * multiple ARPACK calls, as the ARPACK solvers do not modify it.
 *
 * Input options:
 * \member bmat Character. Whether to solve a standard ('I') ot a
 *    generalized problem ('B').
 * \member n Dimension of the eigenproblem.
 * \member which Specifies which eigenvalues/vectors to
 *    compute. Possible values for symmetric matrices:
 *    \clist \cli LA
 *                Compute \c nev largest (algebraic) eigenvalues.
 *           \cli SA
 *                Compute \c nev smallest (algebraic) eigenvalues.
 *           \cli LM
 *                Compute \c nev largest (in magnitude) eigenvalues.
 *           \cli SM
 *                Compute \c nev smallest (in magnitude) eigenvalues.
 *           \cli BE
 *                Compute \c nev eigenvalues, half from each end of
 *                   the spectrum. When \c nev is odd, compute one
 *                   more from the high en than from the low
 *                   end. \endclist
 *    Possible values for non-symmetric matrices:
 *    \clist \cli LM
 *                Compute \c nev largest (in magnitude) eigenvalues.
 *           \cli SM
 *                Compute \c nev smallest (in magnitude) eigenvalues.
 *           \cli LR
 *                Compute \c nev eigenvalues of largest real part.
 *           \cli SR
 *                Compute \c nev eigenvalues of smallest real part.
 *           \cli LI
 *                Compute \c nev eigenvalues of largest imaginary part.
 *           \cli SI
 *                Compute \c nev eigenvalues of smallest imaginary
 *                    part. \endclist
 * \member nev The number of eigenvalues to be computed.
 * \member tol Stopping criterion: the relative accuracy
 *    of the Ritz value is considered acceptable if its error is less
 *    than \c tol times its estimated value. If this is set to zero
 *    then machine precision is used.
 * \member ncv Number of Lanczos vectors to be generated. Setting this
 *    to zero means that \ref igraph_arpack_rssolve and \ref igraph_arpack_rnsolve
 *    will determine a suitable value for \c ncv automatically.
 * \member ldv Numberic scalar. It should be set to
 *    zero in the current igraph implementation.
 * \member ishift Either zero or one. If zero then the shifts are
 *    provided by the user via reverse communication. If one then exact
 *    shifts with respect to the reduced tridiagonal matrix \c T.
 *    Please always set this to one.
 * \member mxiter Maximum number of Arnoldi update iterations allowed.
 * \member nb Blocksize to be used in the recurrence. Please always
 *    leave this on the default value, one.
 * \member mode The type of the eigenproblem to be solved.
 *    Possible values if the input matrix is symmetric:
 *    \olist
 *      \oli A*x=lambda*x, A is symmetric.
 *      \oli A*x=lambda*M*x, A is
 *       symmetric, M is symmetric positive definite.
 *      \oli K*x=lambda*M*x, K is
 *        symmetric, M is symmetric positive semi-definite.
 *      \oli K*x=lambda*KG*x, K is
 *       symmetric positive semi-definite, KG is symmetric
 *       indefinite.
 *     \oli A*x=lambda*M*x, A is
 *       symmetric, M is symmetric positive
 *       semi-definite. (Cayley transformed mode.) \endolist
 *    Please note that only \c mode ==1 was tested and other values
 *    might not work properly.
 *    Possible values if the input matrix is not symmetric:
 *    \olist
 *     \oli A*x=lambda*x.
 *     \oli A*x=lambda*M*x, M is
 *       symmetric positive definite.
 *     \oli A*x=lambda*M*x, M is
 *       symmetric semi-definite.
 *     \oli A*x=lambda*M*x, M is
 *           symmetric semi-definite. \endolist
 *     Please note that only \c mode == 1 was tested and other values
 *     might not work properly.
 * \member start Whether to use the supplied starting vector (1), or
 *    use a random starting vector (0). The starting vector must be
 *    supplied in the first column of the \c vectors argument of the
 *    \ref igraph_arpack_rssolve() of \ref igraph_arpack_rnsolve() call.
 *
 * Output options:
 * \member info Error flag of ARPACK. Possible values:
 *    \clist \cli 0
 *                Normal exit.
 *           \cli 1
 *                Maximum number of iterations taken.
 *           \cli 3
 *                No shifts could be applied during a cycle of the
 *         Implicitly restarted Arnoldi iteration. One possibility
 *         is to increase the size of \c ncv relative to \c
 *           nev. \endclist
 *    ARPACK can return other error flags as well, but these are
 *    converted to igraph errors, see \ref igraph_error_type_t.
 * \member ierr Error flag of the second ARPACK call (one eigenvalue
 *     computation usually involves two calls to ARPACK). This is
 *     always zero, as other error codes are converted to igraph errors.
 * \member noiter Number of Arnoldi iterations taken.
 * \member nconv Number of converged Ritz values. This
 *     represents the number of Ritz values that satisfy the
 *     convergence critetion.
 * \member numop Total number of matrix-vector multiplications.
 * \member numopb Not used currently.
 * \member numreo Total number of steps of re-orthogonalization.
 *
 * Internal options:
 * \member lworkl Do not modify this option.
 * \member sigma The shift for the shift-invert mode.
 * \member sigmai The imaginary part of the shift, for the
 *    non-symmetric or complex shift-invert mode.
 * \member iparam Do not modify this option.
 * \member ipntr Do not modify this option.
 *
 */

typedef struct igraph_arpack_options_t {
    /* INPUT */
    char bmat[1];         /* I-standard problem, G-generalized */
    int n;            /* Dimension of the eigenproblem */
    char which[2];        /* LA, SA, LM, SM, BE */
    int nev;                 /* Number of eigenvalues to be computed */
    igraph_real_t tol;        /* Stopping criterion */
    int ncv;          /* Number of columns in V */
    int ldv;          /* Leading dimension of V */
    int ishift;       /* 0-reverse comm., 1-exact with tridiagonal */
    int mxiter;              /* Maximum number of update iterations to take */
    int nb;           /* Block size on the recurrence, only 1 works */
    int mode;     /* The kind of problem to be solved (1-5)
                   1: A*x=l*x, A symmetric
                   2: A*x=l*M*x, A symm. M pos. def.
                   3: K*x = l*M*x, K symm., M pos. semidef.
                   4: K*x = l*KG*x, K s. pos. semidef. KG s. indef.
                   5: A*x = l*M*x, A symm., M symm. pos. semidef. */
    int start;        /* 0: random, 1: use the supplied vector */
    int lworkl;       /* Size of temporary storage, default is fine */
    igraph_real_t sigma;          /* The shift for modes 3,4,5 */
    igraph_real_t sigmai;     /* The imaginary part of shift for rnsolve */
    /* OUTPUT */
    int info;     /* What happened, see docs */
    int ierr;     /* What happened  in the dseupd call */
    int noiter;       /* The number of iterations taken */
    int nconv;
    int numop;        /* Number of OP*x operations */
    int numopb;       /* Number of B*x operations if BMAT='G' */
    int numreo;       /* Number of steps of re-orthogonalizations */
    /* INTERNAL */
    int iparam[11];
    int ipntr[14];
} igraph_arpack_options_t;

/**
 * \struct igraph_arpack_storage_t
 * \brief Storage for ARPACK
 *
 * Public members, do not modify them directly, these are considered
 * to be read-only.
 * \member maxn Maximum rank of matrix.
 * \member maxncv Maximum NCV.
 * \member maxldv Maximum LDV.
 *
 * These members are considered to be private:
 * \member workl Working memory.
 * \member workd Working memory.
 * \member d Memory for eigenvalues.
 * \member resid Memory for residuals.
 * \member ax Working memory.
 * \member select Working memory.
 * \member di Memory for eigenvalues, non-symmetric case only.
 * \member workev Working memory, non-symmetric case only.
 */

typedef struct igraph_arpack_storage_t {
    int maxn, maxncv, maxldv;
    igraph_real_t *v;
    igraph_real_t *workl;
    igraph_real_t *workd;
    igraph_real_t *d;
    igraph_real_t *resid;
    igraph_real_t *ax;
    int *select;
    igraph_real_t *di;        /* These two only for non-symmetric problems */
    igraph_real_t *workev;
} igraph_arpack_storage_t;

DECLDIR void igraph_arpack_options_init(igraph_arpack_options_t *o);

DECLDIR int igraph_arpack_storage_init(igraph_arpack_storage_t *s, long int maxn,
                                       long int maxncv, long int maxldv, igraph_bool_t symm);
DECLDIR void igraph_arpack_storage_destroy(igraph_arpack_storage_t *s);

/**
 * \typedef igraph_arpack_function_t
 * Type of the ARPACK callback function
 *
 * \param to Pointer to an \c igraph_real_t, the result of the
 *    matrix-vector product is expected to be stored here.
 * \param from Pointer to an \c igraph_real_t, the input matrix should
 *    be multiplied by the vector stored here.
 * \param n The length of the vector (which is the same as the order
 *    of the input matrix).
 * \param extra Extra argument to the matrix-vector calculation
 *    function. This is coming from the \ref igraph_arpack_rssolve()
 *    or \ref igraph_arpack_rnsolve() function.
 * \return Error code, if not zero, then the ARPACK solver considers
 *    this as an error, stops and calls the igraph error handler.
 */

typedef int igraph_arpack_function_t(igraph_real_t *to, const igraph_real_t *from,
                                     int n, void *extra);

DECLDIR int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
                                  igraph_arpack_options_t *options,
                                  igraph_arpack_storage_t *storage,
                                  igraph_vector_t *values, igraph_matrix_t *vectors);

DECLDIR int igraph_arpack_rnsolve(igraph_arpack_function_t *fun, void *extra,
                                  igraph_arpack_options_t *options,
                                  igraph_arpack_storage_t *storage,
                                  igraph_matrix_t *values, igraph_matrix_t *vectors);

DECLDIR int igraph_arpack_unpack_complex(igraph_matrix_t *vectors, igraph_matrix_t *values,
        long int nev);

__END_DECLS

#endif
