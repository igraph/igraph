/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 noet: */
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

#include "igraph_arpack.h"
#include "igraph_memory.h"

#include "linalg/arpack_internal.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* The ARPACK example file dssimp.f is used as a template */

static int igraph_i_arpack_err_dsaupd(int error) {
    switch (error) {
    case  1:      return IGRAPH_ARPACK_MAXIT;
    case  3:      return IGRAPH_ARPACK_NOSHIFT;
    case -1:      return IGRAPH_ARPACK_NPOS;
    case -2:      return IGRAPH_ARPACK_NEVNPOS;
    case -3:      return IGRAPH_ARPACK_NCVSMALL;
    case -4:      return IGRAPH_ARPACK_NONPOSI;
    case -5:      return IGRAPH_ARPACK_WHICHINV;
    case -6:      return IGRAPH_ARPACK_BMATINV;
    case -7:      return IGRAPH_ARPACK_WORKLSMALL;
    case -8:      return IGRAPH_ARPACK_TRIDERR;
    case -9:      return IGRAPH_ARPACK_ZEROSTART;
    case -10:     return IGRAPH_ARPACK_MODEINV;
    case -11:     return IGRAPH_ARPACK_MODEBMAT;
    case -12:     return IGRAPH_ARPACK_ISHIFT;
    case -13:     return IGRAPH_ARPACK_NEVBE;
    case -9999:   return IGRAPH_ARPACK_NOFACT;
    default:      return IGRAPH_ARPACK_UNKNOWN;
    }
}

static int igraph_i_arpack_err_dseupd(int error) {
    switch (error) {
    case -1:      return IGRAPH_ARPACK_NPOS;
    case -2:      return IGRAPH_ARPACK_NEVNPOS;
    case -3:      return IGRAPH_ARPACK_NCVSMALL;
    case -5:      return IGRAPH_ARPACK_WHICHINV;
    case -6:      return IGRAPH_ARPACK_BMATINV;
    case -7:      return IGRAPH_ARPACK_WORKLSMALL;
    case -8:      return IGRAPH_ARPACK_TRIDERR;
    case -9:      return IGRAPH_ARPACK_ZEROSTART;
    case -10:     return IGRAPH_ARPACK_MODEINV;
    case -11:     return IGRAPH_ARPACK_MODEBMAT;
    case -12:     return IGRAPH_ARPACK_NEVBE;
    case -14:     return IGRAPH_ARPACK_FAILED;
    case -15:     return IGRAPH_ARPACK_HOWMNY;
    case -16:     return IGRAPH_ARPACK_HOWMNYS;
    case -17:     return IGRAPH_ARPACK_EVDIFF;
    default:      return IGRAPH_ARPACK_UNKNOWN;
    }

}

static int igraph_i_arpack_err_dnaupd(int error) {
    switch (error) {
    case  1:      return IGRAPH_ARPACK_MAXIT;
    case  3:      return IGRAPH_ARPACK_NOSHIFT;
    case -1:      return IGRAPH_ARPACK_NPOS;
    case -2:      return IGRAPH_ARPACK_NEVNPOS;
    case -3:      return IGRAPH_ARPACK_NCVSMALL;
    case -4:      return IGRAPH_ARPACK_NONPOSI;
    case -5:      return IGRAPH_ARPACK_WHICHINV;
    case -6:      return IGRAPH_ARPACK_BMATINV;
    case -7:      return IGRAPH_ARPACK_WORKLSMALL;
    case -8:      return IGRAPH_ARPACK_TRIDERR;
    case -9:      return IGRAPH_ARPACK_ZEROSTART;
    case -10:     return IGRAPH_ARPACK_MODEINV;
    case -11:     return IGRAPH_ARPACK_MODEBMAT;
    case -12:     return IGRAPH_ARPACK_ISHIFT;
    case -9999:   return IGRAPH_ARPACK_NOFACT;
    default:      return IGRAPH_ARPACK_UNKNOWN;
    }
}

static int igraph_i_arpack_err_dneupd(int error) {
    switch (error) {
    case  1:      return IGRAPH_ARPACK_REORDER;
    case -1:      return IGRAPH_ARPACK_NPOS;
    case -2:      return IGRAPH_ARPACK_NEVNPOS;
    case -3:      return IGRAPH_ARPACK_NCVSMALL;
    case -5:      return IGRAPH_ARPACK_WHICHINV;
    case -6:      return IGRAPH_ARPACK_BMATINV;
    case -7:      return IGRAPH_ARPACK_WORKLSMALL;
    case -8:      return IGRAPH_ARPACK_SHUR;
    case -9:      return IGRAPH_ARPACK_LAPACK;
    case -10:     return IGRAPH_ARPACK_MODEINV;
    case -11:     return IGRAPH_ARPACK_MODEBMAT;
    case -12:     return IGRAPH_ARPACK_HOWMNYS;
    case -13:     return IGRAPH_ARPACK_HOWMNY;
    case -14:     return IGRAPH_ARPACK_FAILED;
    case -15:     return IGRAPH_ARPACK_EVDIFF;
    default:      return IGRAPH_ARPACK_UNKNOWN;
    }
}

/**
 * \function igraph_arpack_options_init
 * Initialize ARPACK options
 *
 * Initializes ARPACK options, set them to default values.
 * You can always pass the initialized \ref igraph_arpack_options_t
 * object to built-in igraph functions without any modification. The
 * built-in igraph functions modify the options to perform their
 * calculation, e.g. \ref igraph_pagerank() always searches for the
 * eigenvalue with the largest magnitude, regardless of the supplied
 * value.
 * </para><para>
 * If you want to implement your own function involving eigenvalue
 * calculation using ARPACK, however, you will likely need to set up
 * the fields for yourself.
 * \param o The \ref igraph_arpack_options_t object to initialize.
 *
 * Time complexity: O(1).
 */

void igraph_arpack_options_init(igraph_arpack_options_t *o) {
    o->bmat[0] = 'I';
    o->n = 0;         /* needs to be updated! */
    o->which[0] = 'X'; o->which[1] = 'X';
    o->nev = 1;
    o->tol = 0;
    o->ncv = 0;       /* 0 means "automatic" */
    o->ldv = o->n;        /* will be updated to (real) n */
    o->ishift = 1;
    o->mxiter = 3000;
    o->nb = 1;
    o->mode = 1;
    o->start = 0;
    o->lworkl = 0;
    o->sigma = 0;
    o->sigmai = 0;
    o->info = o->start;

    o->iparam[0] = o->ishift; o->iparam[1] = 0; o->iparam[2] = o->mxiter; o->iparam[3] = o->nb;
    o->iparam[4] = 0; o->iparam[5] = 0; o->iparam[6] = o->mode; o->iparam[7] = 0;
    o->iparam[8] = 0; o->iparam[9] = 0; o->iparam[10] = 0;
}

/**
 * \function igraph_arpack_storage_init
 * Initialize ARPACK storage
 *
 * You only need this function if you want to run multiple eigenvalue
 * calculations using ARPACK, and want to spare the memory
 * allocation/deallocation between each two runs. Otherwise it is safe
 * to supply a null pointer as the \c storage argument of both \ref
 * igraph_arpack_rssolve() and \ref igraph_arpack_rnsolve() to make
 * memory allocated and deallocated automatically.
 *
 * </para><para>Don't forget to call the \ref
 * igraph_arpack_storage_destroy() function on the storage object if
 * you don't need it any more.
 * \param s The \ref igraph_arpack_storage_t object to initialize.
 * \param maxn The maximum order of the matrices.
 * \param maxncv The maximum NCV parameter intended to use.
 * \param maxldv The maximum LDV parameter intended to use.
 * \param symm Whether symmetric or non-symmetric problems will be
 *    solved using this \ref igraph_arpack_storage_t. (You cannot use
 *    the same storage both with symmetric and non-symmetric solvers.)
 * \return Error code.
 *
 * Time complexity: O(maxncv*(maxldv+maxn)).
 */

int igraph_arpack_storage_init(igraph_arpack_storage_t *s, long int maxn,
                               long int maxncv, long int maxldv,
                               igraph_bool_t symm) {

    /* TODO: check arguments */
    s->maxn = (int) maxn;
    s->maxncv = (int) maxncv;
    s->maxldv = (int) maxldv;

#define CHECKMEM(x) \
    if (!x) { \
        IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

    s->v = IGRAPH_CALLOC(maxldv * maxncv, igraph_real_t); CHECKMEM(s->v);
    s->workd = IGRAPH_CALLOC(3 * maxn, igraph_real_t); CHECKMEM(s->workd);
    s->d = IGRAPH_CALLOC(2 * maxncv, igraph_real_t); CHECKMEM(s->d);
    s->resid = IGRAPH_CALLOC(maxn, igraph_real_t); CHECKMEM(s->resid);
    s->ax = IGRAPH_CALLOC(maxn, igraph_real_t); CHECKMEM(s->ax);
    s->select = IGRAPH_CALLOC(maxncv, int); CHECKMEM(s->select);

    if (symm) {
        s->workl = IGRAPH_CALLOC(maxncv * (maxncv + 8), igraph_real_t); CHECKMEM(s->workl);
        s->di = 0;
        s->workev = 0;
    } else {
        s->workl = IGRAPH_CALLOC(3 * maxncv * (maxncv + 2), igraph_real_t); CHECKMEM(s->workl);
        s->di = IGRAPH_CALLOC(2 * maxncv, igraph_real_t); CHECKMEM(s->di);
        s->workev = IGRAPH_CALLOC(3 * maxncv, igraph_real_t); CHECKMEM(s->workev);
        IGRAPH_FINALLY_CLEAN(2);
    }

#undef CHECKMEM

    IGRAPH_FINALLY_CLEAN(7);
    return 0;
}

/**
 * \function igraph_arpack_storage_destroy
 * Deallocate ARPACK storage
 *
 * \param s The \ref igraph_arpack_storage_t object for which the
 *    memory will be deallocated.
 *
 * Time complexity: operating system dependent.
 */

void igraph_arpack_storage_destroy(igraph_arpack_storage_t *s) {

    if (s->di) {
        IGRAPH_FREE(s->di);
    }
    if (s->workev) {
        IGRAPH_FREE(s->workev);
    }

    IGRAPH_FREE(s->workl);
    IGRAPH_FREE(s->select);
    IGRAPH_FREE(s->ax);
    IGRAPH_FREE(s->resid);
    IGRAPH_FREE(s->d);
    IGRAPH_FREE(s->workd);
    IGRAPH_FREE(s->v);
}

/**
 * "Solver" for 1x1 eigenvalue problems since ARPACK sometimes blows up with
 * these.
 */
static int igraph_i_arpack_rssolve_1x1(igraph_arpack_function_t *fun, void *extra,
                                       igraph_arpack_options_t* options,
                                       igraph_vector_t* values, igraph_matrix_t* vectors) {
    igraph_real_t a, b;
    int nev = options->nev;

    if (nev <= 0) {
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_NEVNPOS);
    }

    /* Probe the value in the matrix */
    a = 1;
    if (fun(&b, &a, 1, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }

    options->nconv = nev;

    if (values != 0) {
        IGRAPH_CHECK(igraph_vector_resize(values, 1));
        VECTOR(*values)[0] = b;
    }

    if (vectors != 0) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, 1, 1));
        MATRIX(*vectors, 0, 0) = 1;
    }

    return IGRAPH_SUCCESS;
}

/**
 * "Solver" for 1x1 eigenvalue problems since ARPACK sometimes blows up with
 * these.
 */
static int igraph_i_arpack_rnsolve_1x1(igraph_arpack_function_t *fun, void *extra,
                                       igraph_arpack_options_t* options,
                                       igraph_matrix_t* values, igraph_matrix_t* vectors) {
    igraph_real_t a, b;
    int nev = options->nev;

    if (nev <= 0) {
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_NEVNPOS);
    }

    /* Probe the value in the matrix */
    a = 1;
    if (fun(&b, &a, 1, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }

    options->nconv = nev;

    if (values != 0) {
        IGRAPH_CHECK(igraph_matrix_resize(values, 1, 2));
        MATRIX(*values, 0, 0) = b; MATRIX(*values, 0, 1) = 0;
    }

    if (vectors != 0) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, 1, 1));
        MATRIX(*vectors, 0, 0) = 1;
    }

    return IGRAPH_SUCCESS;
}

/**
 * "Solver" for 2x2 nonsymmetric eigenvalue problems since ARPACK sometimes
 * blows up with these.
 */
static int igraph_i_arpack_rnsolve_2x2(igraph_arpack_function_t *fun, void *extra,
                                       igraph_arpack_options_t* options, igraph_matrix_t* values,
                                       igraph_matrix_t* vectors) {
    igraph_real_t vec[2], mat[4];
    igraph_real_t a, b, c, d;
    igraph_real_t trace, det, tsq4_minus_d;
    igraph_complex_t eval1, eval2;
    igraph_complex_t evec1[2], evec2[2];
    igraph_bool_t swap_evals = 0;
    igraph_bool_t complex_evals = 0;
    int nev = options->nev;

    if (nev <= 0) {
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_NEVNPOS);
    }
    if (nev > 2) {
        nev = 2;
    }

    /* Probe the values in the matrix */
    vec[0] = 1; vec[1] = 0;
    if (fun(mat, vec, 2, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }
    vec[0] = 0; vec[1] = 1;
    if (fun(mat + 2, vec, 2, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }
    a = mat[0]; b = mat[2]; c = mat[1]; d = mat[3];

    /* Get the trace and the determinant */
    trace = a + d;
    det = a * d - b * c;
    tsq4_minus_d = trace * trace / 4 - det;

    /* Calculate the eigenvalues */
    complex_evals = tsq4_minus_d < 0;
    eval1 = igraph_complex_sqrt_real(tsq4_minus_d);
    if (complex_evals) {
        eval2 = igraph_complex_mul_real(eval1, -1);
    } else {
        /* to avoid having -0 in the imaginary part */
        eval2 = igraph_complex(-IGRAPH_REAL(eval1), 0);
    }
    eval1 = igraph_complex_add_real(eval1, trace / 2);
    eval2 = igraph_complex_add_real(eval2, trace / 2);

    if (c != 0) {
        evec1[0] = igraph_complex_sub_real(eval1, d);
        evec1[1] = igraph_complex(c, 0);
        evec2[0] = igraph_complex_sub_real(eval2, d);
        evec2[1] = igraph_complex(c, 0);
    } else if (b != 0) {
        evec1[0] = igraph_complex(b, 0);
        evec1[1] = igraph_complex_sub_real(eval1, a);
        evec2[0] = igraph_complex(b, 0);
        evec2[1] = igraph_complex_sub_real(eval2, a);
    } else {
        evec1[0] = igraph_complex(1, 0);
        evec1[1] = igraph_complex(0, 0);
        evec2[0] = igraph_complex(0, 0);
        evec2[1] = igraph_complex(1, 0);
    }

    /* Sometimes we have to swap eval1 with eval2 and evec1 with eval2;
     * determine whether we have to do it now */
    if (options->which[0] == 'S') {
        if (options->which[1] == 'M') {
            /* eval1 must be the one with the smallest magnitude */
            swap_evals = (igraph_complex_mod(eval1) > igraph_complex_mod(eval2));
        } else if (options->which[1] == 'R') {
            /* eval1 must be the one with the smallest real part */
            swap_evals = (IGRAPH_REAL(eval1) > IGRAPH_REAL(eval2));
        } else if (options->which[1] == 'I') {
            /* eval1 must be the one with the smallest imaginary part */
            swap_evals = (IGRAPH_IMAG(eval1) > IGRAPH_IMAG(eval2));
        } else {
            IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_WHICHINV);
        }
    } else if (options->which[0] == 'L') {
        if (options->which[1] == 'M') {
            /* eval1 must be the one with the largest magnitude */
            swap_evals = (igraph_complex_mod(eval1) < igraph_complex_mod(eval2));
        } else if (options->which[1] == 'R') {
            /* eval1 must be the one with the largest real part */
            swap_evals = (IGRAPH_REAL(eval1) < IGRAPH_REAL(eval2));
        } else if (options->which[1] == 'I') {
            /* eval1 must be the one with the largest imaginary part */
            swap_evals = (IGRAPH_IMAG(eval1) < IGRAPH_IMAG(eval2));
        } else {
            IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_WHICHINV);
        }
    } else if (options->which[0] == 'X' && options->which[1] == 'X') {
        /* No preference on the ordering of eigenvectors */
    } else {
        /* fprintf(stderr, "%c%c\n", options->which[0], options->which[1]); */
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_WHICHINV);
    }

    options->nconv = nev;

    if (swap_evals) {
        igraph_complex_t dummy;
        dummy = eval1; eval1 = eval2; eval2 = dummy;
        dummy = evec1[0]; evec1[0] = evec2[0]; evec2[0] = dummy;
        dummy = evec1[1]; evec1[1] = evec2[1]; evec2[1] = dummy;
    }

    if (complex_evals) {
        /* The eigenvalues are conjugate pairs, so we store only the
         * one with positive imaginary part */
        if (IGRAPH_IMAG(eval1) < 0) {
            eval1 = eval2;
            evec1[0] = evec2[0]; evec1[1] = evec2[1];
        }
    }

    if (values != 0) {
        IGRAPH_CHECK(igraph_matrix_resize(values, nev, 2));
        MATRIX(*values, 0, 0) = IGRAPH_REAL(eval1);
        MATRIX(*values, 0, 1) = IGRAPH_IMAG(eval1);
        if (nev > 1) {
            MATRIX(*values, 1, 0) = IGRAPH_REAL(eval2);
            MATRIX(*values, 1, 1) = IGRAPH_IMAG(eval2);
        }
    }

    if (vectors != 0) {
        if (complex_evals) {
            IGRAPH_CHECK(igraph_matrix_resize(vectors, 2, 2));
            MATRIX(*vectors, 0, 0) = IGRAPH_REAL(evec1[0]);
            MATRIX(*vectors, 1, 0) = IGRAPH_REAL(evec1[1]);
            MATRIX(*vectors, 0, 1) = IGRAPH_IMAG(evec1[0]);
            MATRIX(*vectors, 1, 1) = IGRAPH_IMAG(evec1[1]);
        } else {
            IGRAPH_CHECK(igraph_matrix_resize(vectors, 2, nev));
            MATRIX(*vectors, 0, 0) = IGRAPH_REAL(evec1[0]);
            MATRIX(*vectors, 1, 0) = IGRAPH_REAL(evec1[1]);
            if (nev > 1) {
                MATRIX(*vectors, 0, 1) = IGRAPH_REAL(evec2[0]);
                MATRIX(*vectors, 1, 1) = IGRAPH_REAL(evec2[1]);
            }
        }
    }

    return IGRAPH_SUCCESS;
}

/**
 * "Solver" for symmetric 2x2 eigenvalue problems since ARPACK sometimes blows
 * up with these.
 */
static int igraph_i_arpack_rssolve_2x2(igraph_arpack_function_t *fun, void *extra,
                                       igraph_arpack_options_t* options, igraph_vector_t* values,
                                       igraph_matrix_t* vectors) {
    igraph_real_t vec[2], mat[4];
    igraph_real_t a, b, c, d;
    igraph_real_t trace, det, tsq4_minus_d;
    igraph_real_t eval1, eval2;
    int nev = options->nev;

    if (nev <= 0) {
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_NEVNPOS);
    }
    if (nev > 2) {
        nev = 2;
    }

    /* Probe the values in the matrix */
    vec[0] = 1; vec[1] = 0;
    if (fun(mat, vec, 2, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }
    vec[0] = 0; vec[1] = 1;
    if (fun(mat + 2, vec, 2, extra)) {
        IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                     IGRAPH_ARPACK_PROD);
    }
    a = mat[0]; b = mat[2]; c = mat[1]; d = mat[3];

    /* Get the trace and the determinant */
    trace = a + d;
    det = a * d - b * c;
    tsq4_minus_d = trace * trace / 4 - det;

    if (tsq4_minus_d >= 0) {
        /* Both eigenvalues are real */
        eval1 = trace / 2 + sqrt(tsq4_minus_d);
        eval2 = trace / 2 - sqrt(tsq4_minus_d);
        if (c != 0) {
            mat[0] = eval1 - d; mat[2] = eval2 - d;
            mat[1] = c;       mat[3] = c;
        } else if (b != 0) {
            mat[0] = b;       mat[2] = b;
            mat[1] = eval1 - a; mat[3] = eval2 - a;
        } else {
            mat[0] = 1; mat[2] = 0;
            mat[1] = 0; mat[3] = 1;
        }
    } else {
        /* Both eigenvalues are complex. Should not happen with symmetric
         * matrices. */
        IGRAPH_ERROR("ARPACK error, 2x2 matrix is not symmetric", IGRAPH_EINVAL);
    }

    /* eval1 is always the larger eigenvalue. If we want the smaller
     * one, we have to swap eval1 with eval2 and also the columns of mat */
    if (options->which[0] == 'S') {
        trace = eval1; eval1 = eval2; eval2 = trace;
        trace = mat[0]; mat[0] = mat[2]; mat[2] = trace;
        trace = mat[1]; mat[1] = mat[3]; mat[3] = trace;
    } else if (options->which[0] == 'L' || options->which[0] == 'B') {
        /* Nothing to do here */
    } else if (options->which[0] == 'X' && options->which[1] == 'X') {
        /* No preference on the ordering of eigenvectors */
    } else {
        IGRAPH_ERROR("ARPACK error", IGRAPH_ARPACK_WHICHINV);
    }

    options->nconv = nev;

    if (values != 0) {
        IGRAPH_CHECK(igraph_vector_resize(values, nev));
        VECTOR(*values)[0] = eval1;
        if (nev > 1) {
            VECTOR(*values)[1] = eval2;
        }
    }

    if (vectors != 0) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, 2, nev));
        MATRIX(*vectors, 0, 0) = mat[0];
        MATRIX(*vectors, 1, 0) = mat[1];
        if (nev > 1) {
            MATRIX(*vectors, 0, 1) = mat[2];
            MATRIX(*vectors, 1, 1) = mat[3];
        }
    }

    return IGRAPH_SUCCESS;
}

int igraph_arpack_rssort(igraph_vector_t *values, igraph_matrix_t *vectors,
                         const igraph_arpack_options_t *options,
                         igraph_real_t *d, const igraph_real_t *v) {

    igraph_vector_t order;
    char sort[2];
    int apply = 1;
    unsigned int n = (unsigned int) options->n;
    int nconv = options->nconv;
    int nev = options->nev;
    unsigned int nans = (unsigned int) (nconv < nev ? nconv : nev);
    unsigned int i;

#define which(a,b) (options->which[0]==a && options->which[1]==b)

    if (which('L', 'A')) {
        sort[0] = 'S'; sort[1] = 'A';
    } else if (which('S', 'A')) {
        sort[0] = 'L'; sort[1] = 'A';
    } else if (which('L', 'M')) {
        sort[0] = 'S'; sort[1] = 'M';
    } else if (which('S', 'M')) {
        sort[0] = 'L'; sort[1] = 'M';
    } else if (which('B', 'E')) {
        sort[0] = 'L'; sort[1] = 'A';
    }

    IGRAPH_CHECK(igraph_vector_init_seq(&order, 0, nconv - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &order);
#ifdef HAVE_GFORTRAN
    igraphdsortr_(sort, &apply, &nconv, d, VECTOR(order), /*which_len=*/ 2);
#else
    igraphdsortr_(sort, &apply, &nconv, d, VECTOR(order));
#endif

    /* BE is special */
    if (which('B', 'E')) {
        int w = 0, l1 = 0, l2 = nev - 1;
        igraph_vector_t order2, d2;
        IGRAPH_VECTOR_INIT_FINALLY(&order2, nev);
        IGRAPH_VECTOR_INIT_FINALLY(&d2, nev);
        while (l1 <= l2) {
            VECTOR(order2)[w] = VECTOR(order)[l1];
            VECTOR(d2)[w] = d[l1];
            w++; l1++;
            if (l1 <= l2) {
                VECTOR(order2)[w] = VECTOR(order)[l2];
                VECTOR(d2)[w] = d[l2];
                w++; l2--;
            }
        }
        igraph_vector_update(&order, &order2);
        igraph_vector_copy_to(&d2, d);
        igraph_vector_destroy(&order2);
        igraph_vector_destroy(&d2);
        IGRAPH_FINALLY_CLEAN(2);
    }

#undef which

    /* Copy values */
    if (values) {
        IGRAPH_CHECK(igraph_vector_resize(values, nans));
        memcpy(VECTOR(*values), d, sizeof(igraph_real_t) * nans);
    }

    /* Reorder vectors */
    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, nans));
        for (i = 0; i < nans; i++) {
            unsigned int idx = (unsigned int) VECTOR(order)[i];
            const igraph_real_t *ptr = v + n * idx;
            memcpy(&MATRIX(*vectors, 0, i), ptr, sizeof(igraph_real_t) * n);
        }
    }

    igraph_vector_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int igraph_arpack_rnsort(igraph_matrix_t *values, igraph_matrix_t *vectors,
                         const igraph_arpack_options_t *options,
                         igraph_real_t *dr, igraph_real_t *di,
                         igraph_real_t *v) {

    igraph_vector_t order;
    char sort[2];
    int apply = 1;
    unsigned int n = (unsigned int) options->n;
    int nconv = options->nconv;
    int nev = options->nev;
    unsigned int nans = (unsigned int) (nconv < nev ? nconv : nev);
    unsigned int i;

#define which(a,b) (options->which[0]==a && options->which[1]==b)

    if (which('L', 'M')) {
        sort[0] = 'S'; sort[1] = 'M';
    } else if (which('S', 'M')) {
        sort[0] = 'L'; sort[1] = 'M';
    } else if (which('L', 'R')) {
        sort[0] = 'S'; sort[1] = 'R';
    } else if (which('S', 'R')) {
        sort[0] = 'L'; sort[1] = 'R';
    } else if (which('L', 'I')) {
        sort[0] = 'S'; sort[1] = 'I';
    } else if (which('S', 'I')) {
        sort[0] = 'L'; sort[1] = 'I';
    }

#undef which

    IGRAPH_CHECK(igraph_vector_init_seq(&order, 0, nconv - 1));
    IGRAPH_FINALLY(igraph_vector_destroy, &order);
#ifdef HAVE_GFORTRAN
    igraphdsortc_(sort, &apply, &nconv, dr, di, VECTOR(order), /*which_len=*/ 2);
#else
    igraphdsortc_(sort, &apply, &nconv, dr, di, VECTOR(order));
#endif

    if (values) {
        IGRAPH_CHECK(igraph_matrix_resize(values, nans, 2));
        memcpy(&MATRIX(*values, 0, 0), dr, sizeof(igraph_real_t) * nans);
        memcpy(&MATRIX(*values, 0, 1), di, sizeof(igraph_real_t) * nans);
    }

    if (vectors) {
        int nc = 0, nr = 0, ncol, vx = 0;
        for (i = 0; i < nans; i++) {
            if (di[i] == 0) {
                nr++;
            } else {
                nc++;
            }
        }
        ncol = (nc / 2) * 2 + (nc % 2) * 2 + nr;
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, ncol));

        for (i = 0; i < nans; i++) {
            unsigned int idx;

            idx = (unsigned int) VECTOR(order)[i];

            if (di[i] == 0) {
                /* real eigenvalue, single eigenvector */
                memcpy(&MATRIX(*vectors, 0, vx), v + n * idx, sizeof(igraph_real_t) * n);
                vx++;
            } else if (di[i] > 0) {
                /* complex eigenvalue, positive imaginary part encountered first.
                 * ARPACK stores its eigenvector directly in two consecutive columns.
                 * The complex conjugate pair of the eigenvalue (if any) will be in
                 * the next column and we will skip it because we advance 'i' below */
                memcpy(&MATRIX(*vectors, 0, vx), v + n * idx, sizeof(igraph_real_t) * 2 * n);
                vx += 2;
                i++;
            } else {
                /* complex eigenvalue, negative imaginary part encountered first.
                     * The positive one will be the next one, but we need to copy the
                     * eigenvector corresponding to the eigenvalue with the positive
                     * imaginary part. */
                idx = (unsigned int) VECTOR(order)[i + 1];
                memcpy(&MATRIX(*vectors, 0, vx), v + n * idx, sizeof(igraph_real_t) * 2 * n);
                vx += 2;
                i++;
            }
        }
    }

    igraph_vector_destroy(&order);
    IGRAPH_FINALLY_CLEAN(1);

    if (values) {
        /* Strive to include complex conjugate eigenvalue pairs in a way that the
         * positive imaginary part comes first */
        for (i = 0; i < nans; i++) {
            if (MATRIX(*values, i, 1) == 0) {
                /* Real eigenvalue, nothing to do */
            } else if (MATRIX(*values, i, 1) < 0) {
                /* Negative imaginary part came first; negate the imaginary part for
                 * this eigenvalue and the next one (which is the complex conjugate
                 * pair), and skip it */
                MATRIX(*values, i, 1) *= -1;
                i++;
                if (i < nans) {
                    MATRIX(*values, i, 1) *= -1;
                }
            } else {
                /* Positive imaginary part; skip the next eigenvalue, which is the
                 * complex conjugate pair */
                i++;
            }
        }
    }

    return 0;
}

/**
 * \function igraph_i_arpack_auto_ncv
 * \brief Tries to set up the value of \c ncv in an \c igraph_arpack_options_t
 *        automagically.
 */
static void igraph_i_arpack_auto_ncv(igraph_arpack_options_t* options) {
    /* This is similar to how Octave determines the value of ncv, with some
     * modifications. */
    int min_ncv = options->nev * 2 + 1;

    /* Use twice the number of desired eigenvectors plus one by default */
    options->ncv = min_ncv;
    /* ...but use at least 20 Lanczos vectors... */
    if (options->ncv < 20) {
        options->ncv = 20;
    }
    /* ...but having ncv close to n leads to some problems with small graphs
     * (example: PageRank of "A <--> C, D <--> E, B"), so we don't let it
     * to be larger than n / 2...
     */
    if (options->ncv > options->n / 2) {
        options->ncv = options->n / 2;
    }
    /* ...but we need at least min_ncv. */
    if (options->ncv < min_ncv) {
        options->ncv = min_ncv;
    }
    /* ...but at most n */
    if (options->ncv > options->n) {
        options->ncv = options->n;
    }
}

/**
 * \function igraph_i_arpack_report_no_convergence
 * \brief Prints a warning that informs the user that the ARPACK solver
 *        did not converge.
 */
static void igraph_i_arpack_report_no_convergence(const igraph_arpack_options_t* options) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "ARPACK solver failed to converge (%d iterations, "
             "%d/%d eigenvectors converged)", options->iparam[2],
             options->iparam[4], options->nev);
    IGRAPH_WARNING(buf);
}

/**
 * \function igraph_arpack_rssolve
 * \brief ARPACK solver for symmetric matrices
 *
 * This is the ARPACK solver for symmetric matrices. Please use
 * \ref igraph_arpack_rnsolve() for non-symmetric matrices.
 * \param fun Pointer to an \ref igraph_arpack_function_t object,
 *     the function that performs the matrix-vector multiplication.
 * \param extra An extra argument to be passed to \c fun.
 * \param options An \ref igraph_arpack_options_t object.
 * \param storage An \ref igraph_arpack_storage_t object, or a null
 *     pointer. In the latter case memory allocation and deallocation
 *     is performed automatically. Either this or the \p vectors argument
 *     must be non-null if the ARPACK iteration is started from a
 *     given starting vector. If both are given \p vectors take
 *     precedence.
 * \param values If not a null pointer, then it should be a pointer to an
 *     initialized vector. The eigenvalues will be stored here. The
 *     vector will be resized as needed.
 * \param vectors If not a null pointer, then it must be a pointer to
 *     an initialized matrix. The eigenvectors will be stored in the
 *     columns of the matrix. The matrix will be resized as needed.
 *     Either this or the \p vectors argument must be non-null if the
 *     ARPACK iteration is started from a given starting vector. If
 *     both are given \p vectors take precedence.
 * \return Error code.
 *
 * Time complexity: depends on the matrix-vector
 * multiplication. Usually a small number of iterations is enough, so
 * if the matrix is sparse and the matrix-vector multiplication can be
 * done in O(n) time (the number of vertices), then the eigenvalues
 * are found in O(n) time as well.
 */

int igraph_arpack_rssolve(igraph_arpack_function_t *fun, void *extra,
                          igraph_arpack_options_t *options,
                          igraph_arpack_storage_t *storage,
                          igraph_vector_t *values, igraph_matrix_t *vectors) {

    igraph_real_t *v, *workl, *workd, *d, *resid, *ax;
    igraph_bool_t free_them = 0;
    int *select, i;

    int ido = 0;
    int rvec = vectors || storage ? 1 : 0; /* calculate eigenvectors? */
    char *all = "All";

    int origldv = options->ldv, origlworkl = options->lworkl,
        orignev = options->nev, origncv = options->ncv;
    igraph_real_t origtol = options->tol;
    char origwhich[2];

    origwhich[0] = options->which[0];
    origwhich[1] = options->which[1];

    /* Special case for 1x1 and 2x2 matrices in mode 1 */
    if (options->mode == 1 && options->n == 1) {
        return igraph_i_arpack_rssolve_1x1(fun, extra, options, values, vectors);
    } else if (options->mode == 1 && options->n == 2) {
        return igraph_i_arpack_rssolve_2x2(fun, extra, options, values, vectors);
    }

    /* Brush up options if needed */
    if (options->ldv == 0) {
        options->ldv = options->n;
    }
    if (options->ncv == 0) {
        igraph_i_arpack_auto_ncv(options);
    }
    if (options->lworkl == 0) {
        options->lworkl = options->ncv * (options->ncv + 8);
    }
    if (options->which[0] == 'X') {
        options->which[0] = 'L';
        options->which[1] = 'M';
    }

    if (storage) {
        /* Storage provided */
        if (storage->maxn < options->n) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`n')", IGRAPH_EINVAL);
        }
        if (storage->maxncv < options->ncv) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`ncv')", IGRAPH_EINVAL);
        }
        if (storage->maxldv < options->ldv) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`ldv')", IGRAPH_EINVAL);
        }

        v      = storage->v;
        workl  = storage->workl;
        workd  = storage->workd;
        d      = storage->d;
        resid  = storage->resid;
        ax     = storage->ax;
        select = storage->select;

    } else {
        /* Storage not provided */
        free_them = 1;

#define CHECKMEM(x) \
    if (!x) { \
        IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

        v = IGRAPH_CALLOC(options->ldv * options->ncv, igraph_real_t); CHECKMEM(v);
        workl = IGRAPH_CALLOC(options->lworkl, igraph_real_t); CHECKMEM(workl);
        workd = IGRAPH_CALLOC(3 * options->n, igraph_real_t); CHECKMEM(workd);
        d = IGRAPH_CALLOC(2 * options->ncv, igraph_real_t); CHECKMEM(d);
        resid = IGRAPH_CALLOC(options->n, igraph_real_t); CHECKMEM(resid);
        ax = IGRAPH_CALLOC(options->n, igraph_real_t); CHECKMEM(ax);
        select = IGRAPH_CALLOC(options->ncv, int); CHECKMEM(select);

#undef CHECKMEM

    }

    /* Set final bits */
    options->bmat[0] = 'I';
    options->iparam[0] = options->ishift;
    options->iparam[1] = 0;   // not referenced
    options->iparam[2] = options->mxiter;
    options->iparam[3] = 1;   // currently dsaupd() works only for nb=1
    options->iparam[4] = 0;
    options->iparam[5] = 0;   // not referenced
    options->iparam[6] = options->mode;
    options->iparam[7] = 0;   // return value
    options->iparam[8] = 0;   // return value
    options->iparam[9] = 0;   // return value
    options->iparam[10] = 0;  // return value
    options->info = options->start;
    if (options->start) {
        if (!storage && !vectors) {
            IGRAPH_ERROR("Starting vector not given", IGRAPH_EINVAL);
        }
        if (vectors && (igraph_matrix_nrow(vectors) != options->n ||
                        igraph_matrix_ncol(vectors) != 1)) {
            IGRAPH_ERROR("Invalid starting vector size", IGRAPH_EINVAL);
        }
        if (vectors) {
            for (i = 0; i < options->n; i++) {
                resid[i] = MATRIX(*vectors, i, 0);
            }
        }
    }

    /* Ok, we have everything */
    while (1) {
#ifdef HAVE_GFORTRAN
        igraphdsaupd_(&ido, options->bmat, &options->n, options->which,
                      &options->nev, &options->tol,
                      resid, &options->ncv, v, &options->ldv,
                      options->iparam, options->ipntr,
                      workd, workl, &options->lworkl, &options->info,
                      /*bmat_len=*/ 1, /*which_len=*/ 2);
#else
        igraphdsaupd_(&ido, options->bmat, &options->n, options->which,
                      &options->nev, &options->tol,
                      resid, &options->ncv, v, &options->ldv,
                      options->iparam, options->ipntr,
                      workd, workl, &options->lworkl, &options->info);
#endif

        if (ido == -1 || ido == 1) {
            igraph_real_t *from = workd + options->ipntr[0] - 1;
            igraph_real_t *to = workd + options->ipntr[1] - 1;
            if (fun(to, from, options->n, extra) != 0) {
                IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                             IGRAPH_ARPACK_PROD);
            }

        } else {
            break;
        }
    }

    if (options->info == 1) {
        igraph_i_arpack_report_no_convergence(options);
    }
    if (options->info != 0) {
        IGRAPH_ERROR("ARPACK error", igraph_i_arpack_err_dsaupd(options->info));
    }

    options->ierr = 0;
#ifdef HAVE_GFORTRAN
    igraphdseupd_(&rvec, all, select, d, v, &options->ldv,
                  &options->sigma, options->bmat, &options->n,
                  options->which, &options->nev, &options->tol,
                  resid, &options->ncv, v, &options->ldv, options->iparam,
                  options->ipntr, workd, workl, &options->lworkl,
                  &options->ierr, /*howmny_len=*/ 1, /*bmat_len=*/ 1,
                  /*which_len=*/ 2);
#else
    igraphdseupd_(&rvec, all, select, d, v, &options->ldv,
                  &options->sigma, options->bmat, &options->n,
                  options->which, &options->nev, &options->tol,
                  resid, &options->ncv, v, &options->ldv, options->iparam,
                  options->ipntr, workd, workl, &options->lworkl,
                  &options->ierr);
#endif

    if (options->ierr != 0) {
        IGRAPH_ERROR("ARPACK error", igraph_i_arpack_err_dseupd(options->ierr));
    }

    /* Save the result */

    options->noiter = options->iparam[2];
    options->nconv = options->iparam[4];
    options->numop = options->iparam[8];
    options->numopb = options->iparam[9];
    options->numreo = options->iparam[10];

    if (options->nconv < options->nev) {
        IGRAPH_WARNING("Not enough eigenvalues/vectors in symmetric ARPACK "
                       "solver");
    }

    if (values || vectors) {
        IGRAPH_CHECK(igraph_arpack_rssort(values, vectors, options, d, v));
    }

    options->ldv = origldv;
    options->ncv = origncv;
    options->lworkl = origlworkl;
    options->which[0] = origwhich[0]; options->which[1] = origwhich[1];
    options->tol = origtol;
    options->nev = orignev;

    /* Clean up if needed */
    if (free_them) {
        IGRAPH_FREE(select);
        IGRAPH_FREE(ax);
        IGRAPH_FREE(resid);
        IGRAPH_FREE(d);
        IGRAPH_FREE(workd);
        IGRAPH_FREE(workl);
        IGRAPH_FREE(v);
        IGRAPH_FINALLY_CLEAN(7);
    }
    return 0;
}

/**
 * \function igraph_arpack_rnsolve
 * \brief ARPACK solver for non-symmetric matrices
 *
 * Please always consider calling \ref igraph_arpack_rssolve() if your
 * matrix is symmetric, it is much faster.
 * \ref igraph_arpack_rnsolve() for non-symmetric matrices.
 * </para><para>
 * Note that ARPACK is not called for 2x2 matrices as an exact algebraic
 * solution exists in these cases.
 *
 * \param fun Pointer to an \ref igraph_arpack_function_t object,
 *     the function that performs the matrix-vector multiplication.
 * \param extra An extra argument to be passed to \c fun.
 * \param options An \ref igraph_arpack_options_t object.
 * \param storage An \ref igraph_arpack_storage_t object, or a null
 *     pointer. In the latter case memory allocation and deallocation
 *     is performed automatically.
 * \param values If not a null pointer, then it should be a pointer to an
 *     initialized matrix. The (possibly complex) eigenvalues will be
 *     stored here. The matrix will have two columns, the first column
 *     contains the real, the second the imaginary parts of the
 *     eigenvalues.
 *     The matrix will be resized as needed.
 * \param vectors If not a null pointer, then it must be a pointer to
 *     an initialized matrix. The eigenvectors will be stored in the
 *     columns of the matrix. The matrix will be resized as needed.
 *     Note that real eigenvalues will have real eigenvectors in a single
 *     column in this matrix; however, complex eigenvalues come in conjugate
 *     pairs and the result matrix will store the eigenvector corresponding to
 *     the eigenvalue with \em positive imaginary part only. Since in this case
 *     the eigenvector is also complex, it will occupy \em two columns in the
 *     eigenvector matrix (the real and the imaginary parts, in this order).
 *     Caveat: if the eigenvalue vector returns only the eigenvalue with the
 *     \em negative imaginary part for a complex conjugate eigenvalue pair, the
 *     result vector will \em still store the eigenvector corresponding to the
 *     eigenvalue with the positive imaginary part (since this is how ARPACK
 *     works).
 * \return Error code.
 *
 * Time complexity: depends on the matrix-vector
 * multiplication. Usually a small number of iterations is enough, so
 * if the matrix is sparse and the matrix-vector multiplication can be
 * done in O(n) time (the number of vertices), then the eigenvalues
 * are found in O(n) time as well.
 */

int igraph_arpack_rnsolve(igraph_arpack_function_t *fun, void *extra,
                          igraph_arpack_options_t *options,
                          igraph_arpack_storage_t *storage,
                          igraph_matrix_t *values, igraph_matrix_t *vectors) {

    igraph_real_t *v, *workl, *workd, *dr, *di, *resid, *workev;
    igraph_bool_t free_them = 0;
    int *select, i;

    int ido = 0;
    int rvec = vectors || storage ? 1 : 0;
    char *all = "All";

    int origldv = options->ldv, origlworkl = options->lworkl,
        orignev = options->nev, origncv = options->ncv;
    igraph_real_t origtol = options->tol;
    int d_size;
    char origwhich[2];

    origwhich[0] = options->which[0];
    origwhich[1] = options->which[1];

    /* Special case for 1x1 and 2x2 matrices in mode 1 */
    if (options->mode == 1 && options->n == 1) {
        return igraph_i_arpack_rnsolve_1x1(fun, extra, options, values, vectors);
    } else if (options->mode == 1 && options->n == 2) {
        return igraph_i_arpack_rnsolve_2x2(fun, extra, options, values, vectors);
    }

    /* Brush up options if needed */
    if (options->ldv == 0) {
        options->ldv = options->n;
    }
    if (options->ncv == 0) {
        igraph_i_arpack_auto_ncv(options);
    }
    if (options->lworkl == 0) {
        options->lworkl = 3 * options->ncv * (options->ncv + 2);
    }
    if (options->which[0] == 'X') {
        options->which[0] = 'L';
        options->which[1] = 'M';
    }

    if (storage) {
        /* Storage provided */
        if (storage->maxn < options->n) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`n')", IGRAPH_EINVAL);
        }
        if (storage->maxncv < options->ncv) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`ncv')", IGRAPH_EINVAL);
        }
        if (storage->maxldv < options->ldv) {
            IGRAPH_ERROR("Not enough storage for ARPACK (`ldv')", IGRAPH_EINVAL);
        }

        v      = storage->v;
        workl  = storage->workl;
        workd  = storage->workd;
        workev = storage->workev;
        dr     = storage->d;
        di     = storage->di;
        d_size = options->n;
        resid  = storage->resid;
        select = storage->select;

    } else {
        /* Storage not provided */
        free_them = 1;

#define CHECKMEM(x) \
    if (!x) { \
        IGRAPH_ERROR("Cannot allocate memory for ARPACK", IGRAPH_ENOMEM); \
    } \
    IGRAPH_FINALLY(igraph_free, x);

        v = IGRAPH_CALLOC(options->n * options->ncv, igraph_real_t); CHECKMEM(v);
        workl = IGRAPH_CALLOC(options->lworkl, igraph_real_t); CHECKMEM(workl);
        workd = IGRAPH_CALLOC(3 * options->n, igraph_real_t); CHECKMEM(workd);
        d_size = 2 * options->nev + 1 > options->ncv ? 2 * options->nev + 1 : options->ncv;
        dr = IGRAPH_CALLOC(d_size, igraph_real_t); CHECKMEM(dr);
        di = IGRAPH_CALLOC(d_size, igraph_real_t); CHECKMEM(di);
        resid = IGRAPH_CALLOC(options->n, igraph_real_t); CHECKMEM(resid);
        select = IGRAPH_CALLOC(options->ncv, int); CHECKMEM(select);
        workev = IGRAPH_CALLOC(3 * options->ncv, igraph_real_t); CHECKMEM(workev);

#undef CHECKMEM

    }

    /* Set final bits */
    options->bmat[0] = 'I';
    options->iparam[0] = options->ishift;
    options->iparam[1] = 0;   // not referenced
    options->iparam[2] = options->mxiter;
    options->iparam[3] = 1;   // currently dnaupd() works only for nb=1
    options->iparam[4] = 0;
    options->iparam[5] = 0;   // not referenced
    options->iparam[6] = options->mode;
    options->iparam[7] = 0;   // return value
    options->iparam[8] = 0;   // return value
    options->iparam[9] = 0;   // return value
    options->iparam[10] = 0;  // return value
    options->info = options->start;
    if (options->start) {
        if (!storage && !vectors) {
            IGRAPH_ERROR("Starting vector not given", IGRAPH_EINVAL);
        }
        if (vectors && (igraph_matrix_nrow(vectors) != options->n ||
                        igraph_matrix_ncol(vectors) != 1)) {
            IGRAPH_ERROR("Invalid starting vector size", IGRAPH_EINVAL);
        }
        if (vectors) {
            for (i = 0; i < options->n; i++) {
                resid[i] = MATRIX(*vectors, i, 0);
            }
        }
    }

    /* Ok, we have everything */
    while (1) {
#ifdef HAVE_GFORTRAN
        igraphdnaupd_(&ido, options->bmat, &options->n, options->which,
                      &options->nev, &options->tol,
                      resid, &options->ncv, v, &options->ldv,
                      options->iparam, options->ipntr,
                      workd, workl, &options->lworkl, &options->info,
                      /*bmat_len=*/ 1, /*which_len=*/ 2);
#else
        igraphdnaupd_(&ido, options->bmat, &options->n, options->which,
                      &options->nev, &options->tol,
                      resid, &options->ncv, v, &options->ldv,
                      options->iparam, options->ipntr,
                      workd, workl, &options->lworkl, &options->info);
#endif

        if (ido == -1 || ido == 1) {
            igraph_real_t *from = workd + options->ipntr[0] - 1;
            igraph_real_t *to = workd + options->ipntr[1] - 1;
            if (fun(to, from, options->n, extra) != 0) {
                IGRAPH_ERROR("ARPACK error while evaluating matrix-vector product",
                             IGRAPH_ARPACK_PROD);
            }
        } else {
            break;
        }
    }

    if (options->info == 1) {
        igraph_i_arpack_report_no_convergence(options);
    }
    if (options->info != 0 && options->info != -9999) {
        IGRAPH_ERROR("ARPACK error", igraph_i_arpack_err_dnaupd(options->info));
    }

    options->ierr = 0;
#ifdef HAVE_GFORTRAN
    igraphdneupd_(&rvec, all, select, dr, di, v, &options->ldv,
                  &options->sigma, &options->sigmai, workev, options->bmat,
                  &options->n, options->which, &options->nev, &options->tol,
                  resid, &options->ncv, v, &options->ldv, options->iparam,
                  options->ipntr, workd, workl, &options->lworkl,
                  &options->ierr, /*howmny_len=*/ 1, /*bmat_len=*/ 1,
                  /*which_len=*/ 2);
#else
    igraphdneupd_(&rvec, all, select, dr, di, v, &options->ldv,
                  &options->sigma, &options->sigmai, workev, options->bmat,
                  &options->n, options->which, &options->nev, &options->tol,
                  resid, &options->ncv, v, &options->ldv, options->iparam,
                  options->ipntr, workd, workl, &options->lworkl,
                  &options->ierr);
#endif

    if (options->ierr != 0) {
        IGRAPH_ERROR("ARPACK error", igraph_i_arpack_err_dneupd(options->info));
    }

    /* Save the result */

    options->noiter = options->iparam[2];
    options->nconv = options->iparam[4];
    options->numop = options->iparam[8];
    options->numopb = options->iparam[9];
    options->numreo = options->iparam[10];

    if (options->nconv < options->nev) {
        IGRAPH_WARNING("Not enough eigenvalues/vectors in ARPACK "
                       "solver");
    }

    /* ARPACK might modify stuff in 'options' so reset everything that could
     * potentially get modified */
    options->ldv = origldv;
    options->ncv = origncv;
    options->lworkl = origlworkl;
    options->which[0] = origwhich[0]; options->which[1] = origwhich[1];
    options->tol = origtol;
    options->nev = orignev;

    if (values || vectors) {
        IGRAPH_CHECK(igraph_arpack_rnsort(values, vectors, options,
                                          dr, di, v));
    }

    /* Clean up if needed */
    if (free_them) {
        IGRAPH_FREE(workev);
        IGRAPH_FREE(select);
        IGRAPH_FREE(resid);
        IGRAPH_FREE(di);
        IGRAPH_FREE(dr);
        IGRAPH_FREE(workd);
        IGRAPH_FREE(workl);
        IGRAPH_FREE(v);
        IGRAPH_FINALLY_CLEAN(8);
    }
    return 0;
}

/**
 * \function igraph_arpack_unpack_complex
 * \brief Make the result of the non-symmetric ARPACK solver more readable
 *
 * This function works on the output of \ref igraph_arpack_rnsolve and
 * brushes it up a bit: it only keeps \p nev eigenvalues/vectors and
 * every eigenvector is stored in two columns of the \p vectors
 * matrix.
 *
 * </para><para>
 * The output of the non-symmetric ARPACK solver is somewhat hard to
 * parse, as real eigenvectors occupy only one column in the matrix,
 * and the complex conjugate eigenvectors are not stored at all
 * (usually). The other problem is that the solver might return more
 * eigenvalues than requested. The common use of this function is to
 * call it directly after \ref igraph_arpack_rnsolve with its \p
 * vectors and \p values argument and \c options->nev as \p nev.
 * This will add the vectors for eigenvalues with a negative imaginary
 * part and return all vectors as 2 columns, a real and imaginary part.
 * \param vectors The eigenvector matrix, as returned by \ref
 *   igraph_arpack_rnsolve. It will be resized, typically it will be
 *   larger.
 * \param values The eigenvalue matrix, as returned by \ref
 *   igraph_arpack_rnsolve. It will be resized, typically extra,
 *   unneeded rows (=eigenvalues) will be removed.
 * \param nev The number of eigenvalues/vectors to keep. Can be less
 *   or equal than the number originally requested from ARPACK.
 * \return Error code.
 *
 * Time complexity: linear in the number of elements in the \p vectors
 * matrix.
 */

int igraph_arpack_unpack_complex(igraph_matrix_t *vectors, igraph_matrix_t *values,
                                 long int nev) {

    long int nodes = igraph_matrix_nrow(vectors);
    long int no_evs = igraph_matrix_nrow(values);
    long int i, j;
    long int new_vector_pos;
    long int vector_pos;
    igraph_matrix_t new_vectors;

    /* Error checks */
    if (nev < 0) {
        IGRAPH_ERROR("`nev' cannot be negative", IGRAPH_EINVAL);
    }
    if (nev > no_evs) {
        IGRAPH_ERROR("`nev' too large, we don't have that many in `values'",
                     IGRAPH_EINVAL);
    }

    for (i = no_evs -1; i >= nev; i--) {
        IGRAPH_CHECK(igraph_matrix_remove_row(values, i));
    }

    IGRAPH_CHECK(igraph_matrix_init(&new_vectors, nodes, nev * 2));
    IGRAPH_FINALLY(igraph_matrix_destroy, &new_vectors);

    new_vector_pos = 0;
    vector_pos = 0;
    for (i = 0; i < nev && vector_pos < igraph_matrix_ncol(vectors); i++) {
        if (MATRIX(*values, i, 1) == 0) {
            /* Real eigenvalue */
            for (j = 0; j < nodes; j++) {
                MATRIX(new_vectors, j, new_vector_pos) = MATRIX(*vectors, j, vector_pos);
            }
            new_vector_pos += 2;
            vector_pos += 1;
        } else {
            /* complex eigenvalue */
            for (j = 0; j < nodes; j++) {
                MATRIX(new_vectors, j, new_vector_pos) = MATRIX(*vectors, j, vector_pos);
                MATRIX(new_vectors, j, new_vector_pos + 1) = MATRIX(*vectors, j, vector_pos + 1);
            }

            /* handle the conjugate */

            /* first check if the conjugate eigenvalue is there */
            i++;
            if (i >= nev) {
                break;
            }

            if (MATRIX(*values, i, 1) != -MATRIX(*values, i-1, 1)) {
                IGRAPH_ERROR("Complex eigenvalue not followed by its conjugate.", IGRAPH_EINVAL);
            }

            /* then copy and negate */
            for (j = 0; j < nodes; j++) {
                MATRIX(new_vectors, j, new_vector_pos + 2) = MATRIX(*vectors, j, vector_pos);
                MATRIX(new_vectors, j, new_vector_pos + 3) = -MATRIX(*vectors, j, vector_pos + 1);
            }
            new_vector_pos += 4;
            vector_pos += 2;
        }
    }
    igraph_matrix_destroy(vectors);
    IGRAPH_CHECK(igraph_matrix_copy(vectors, &new_vectors));
    igraph_matrix_destroy(&new_vectors);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
