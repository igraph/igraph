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

#include "igraph_eigen.h"

#include "igraph_qsort.h"
#include "igraph_blas.h"
#include "igraph_interface.h"
#include "igraph_adjlist.h"

#include <string.h>
#include <math.h>
#include <float.h>

static int igraph_i_eigen_arpackfun_to_mat(igraph_arpack_function_t *fun,
                                    int n, void *extra,
                                    igraph_matrix_t *res) {

    int i;
    igraph_vector_t v;

    IGRAPH_CHECK(igraph_matrix_init(res, n, n));
    IGRAPH_FINALLY(igraph_matrix_destroy, res);
    IGRAPH_VECTOR_INIT_FINALLY(&v, n);
    VECTOR(v)[0] = 1;
    IGRAPH_CHECK(fun(/*to=*/ &MATRIX(*res, 0, 0), /*from=*/ VECTOR(v), n,
                             extra));
    for (i = 1; i < n; i++) {
        VECTOR(v)[i - 1] = 0;
        VECTOR(v)[i  ] = 1;
        IGRAPH_CHECK(fun(/*to=*/ &MATRIX(*res, 0, i), /*from=*/ VECTOR(v), n,
                                 extra));
    }
    igraph_vector_destroy(&v);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_lm(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    igraph_matrix_t vec1, vec2;
    igraph_vector_t val1, val2;
    int n = (int) igraph_matrix_nrow(A);
    int p1 = 0, p2 = which->howmany - 1, pr = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&val1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&val2, 0);

    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_init(&vec1, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &vec1);
        IGRAPH_CHECK(igraph_matrix_init(&vec2, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &vec1);
    }

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ 1, /*iu=*/ which->howmany,
                                      /*abstol=*/ 1e-14, &val1,
                                      vectors ? &vec1 : 0,
                                      /*support=*/ 0));

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ n - which->howmany + 1, /*iu=*/ n,
                                      /*abstol=*/ 1e-14, &val2,
                                      vectors ? &vec2 : 0,
                                      /*support=*/ 0));

    if (values) {
        IGRAPH_CHECK(igraph_vector_resize(values, which->howmany));
    }
    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, which->howmany));
    }

    while (pr < which->howmany) {
        if (p2 < 0 || fabs(VECTOR(val1)[p1]) > fabs(VECTOR(val2)[p2])) {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val1)[p1];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec1, 0, p1),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p1++;
            pr++;
        } else {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val2)[p2];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec2, 0, p2),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p2--;
            pr++;
        }
    }


    if (vectors) {
        igraph_matrix_destroy(&vec2);
        igraph_matrix_destroy(&vec1);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_destroy(&val2);
    igraph_vector_destroy(&val1);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_sm(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    igraph_vector_t val;
    igraph_matrix_t vec;
    int i, w = 0, n = (int) igraph_matrix_nrow(A);
    igraph_real_t small;
    int p1, p2, pr = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&val, 0);

    if (vectors) {
        IGRAPH_MATRIX_INIT_FINALLY(&vec, 0, 0);
    }

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_ALL, /*vl=*/ 0,
                                      /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ 0, /*iu=*/ 0,
                                      /*abstol=*/ 1e-14, &val,
                                      vectors ? &vec : 0,
                                      /*support=*/ 0));

    /* Look for smallest value */
    small = fabs(VECTOR(val)[0]);
    for (i = 1; i < n; i++) {
        igraph_real_t v = fabs(VECTOR(val)[i]);
        if (v < small) {
            small = v;
            w = i;
        }
    }
    p1 = w - 1; p2 = w;

    if (values) {
        IGRAPH_CHECK(igraph_vector_resize(values, which->howmany));
    }
    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, which->howmany));
    }

    while (pr < which->howmany) {
        if (p2 == n - 1 || fabs(VECTOR(val)[p1]) < fabs(VECTOR(val)[p2])) {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val)[p1];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec, 0, p1),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p1--;
            pr++;
        } else {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val)[p2];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec, 0, p2),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p2++;
            pr++;
        }
    }

    if (vectors) {
        igraph_matrix_destroy(&vec);
        IGRAPH_FINALLY_CLEAN(1);
    }
    igraph_vector_destroy(&val);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_la(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    /* TODO: ordering? */

    int n = (int) igraph_matrix_nrow(A);
    int il = n - which->howmany + 1;
    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ il, /*iu=*/ n,
                                      /*abstol=*/ 1e-14, values, vectors,
                                      /*support=*/ 0));
    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_sa(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    /* TODO: ordering? */

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ 1, /*iu=*/ which->howmany,
                                      /*abstol=*/ 1e-14, values, vectors,
                                      /*support=*/ 0));

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_be(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    /* TODO: ordering? */

    igraph_matrix_t vec1, vec2;
    igraph_vector_t val1, val2;
    int n = (int) igraph_matrix_nrow(A);
    int p1 = 0, p2 = which->howmany / 2, pr = 0;

    IGRAPH_VECTOR_INIT_FINALLY(&val1, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&val2, 0);

    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_init(&vec1, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &vec1);
        IGRAPH_CHECK(igraph_matrix_init(&vec2, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &vec1);
    }

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ 1, /*iu=*/ (which->howmany) / 2,
                                      /*abstol=*/ 1e-14, &val1,
                                      vectors ? &vec1 : 0,
                                      /*support=*/ 0));

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ n - (which->howmany) / 2, /*iu=*/ n,
                                      /*abstol=*/ 1e-14, &val2,
                                      vectors ? &vec2 : 0,
                                      /*support=*/ 0));

    if (values) {
        IGRAPH_CHECK(igraph_vector_resize(values, which->howmany));
    }
    if (vectors) {
        IGRAPH_CHECK(igraph_matrix_resize(vectors, n, which->howmany));
    }

    while (pr < which->howmany) {
        if (pr % 2) {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val1)[p1];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec1, 0, p1),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p1++;
            pr++;
        } else {
            if (values) {
                VECTOR(*values)[pr] = VECTOR(val2)[p2];
            }
            if (vectors) {
                memcpy(&MATRIX(*vectors, 0, pr), &MATRIX(vec2, 0, p2),
                       sizeof(igraph_real_t) * (size_t) n);
            }
            p2--;
            pr++;
        }
    }

    if (vectors) {
        igraph_matrix_destroy(&vec2);
        igraph_matrix_destroy(&vec1);
        IGRAPH_FINALLY_CLEAN(2);
    }
    igraph_vector_destroy(&val2);
    igraph_vector_destroy(&val1);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_all(const igraph_matrix_t *A,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_ALL, /*vl=*/ 0,
                                      /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ 0, /*iu=*/ 0,
                                      /*abstol=*/ 1e-14, values, vectors,
                                      /*support=*/ 0));

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_iv(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_INTERVAL,
                                      /*vl=*/ which->vl, /*vu=*/ which->vu,
                                      /*vestimate=*/ which->vestimate,
                                      /*il=*/ 0, /*iu=*/ 0,
                                      /*abstol=*/ 1e-14, values, vectors,
                                      /*support=*/ 0));

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack_sel(const igraph_matrix_t *A,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    IGRAPH_CHECK(igraph_lapack_dsyevr(A, IGRAPH_LAPACK_DSYEV_SELECT,
                                      /*vl=*/ 0, /*vu=*/ 0, /*vestimate=*/ 0,
                                      /*il=*/ which->il, /*iu=*/ which->iu,
                                      /*abstol=*/ 1e-14, values, vectors,
                                      /*support=*/ 0));

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_lapack(const igraph_matrix_t *A,
        const igraph_sparsemat_t *sA,
        igraph_arpack_function_t *fun,
        int n, void *extra,
        const igraph_eigen_which_t *which,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    const igraph_matrix_t *myA = A;
    igraph_matrix_t mA;

    /* First we need to create a dense square matrix */

    if (A) {
        n = (int) igraph_matrix_nrow(A);
    } else if (sA) {
        n = (int) igraph_sparsemat_nrow(sA);
        IGRAPH_CHECK(igraph_matrix_init(&mA, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &mA);
        IGRAPH_CHECK(igraph_sparsemat_as_matrix(&mA, sA));
        myA = &mA;
    } else if (fun) {
        IGRAPH_CHECK(igraph_i_eigen_arpackfun_to_mat(fun, n, extra, &mA));
        IGRAPH_FINALLY(igraph_matrix_destroy, &mA);
        myA = &mA;
    }

    switch (which->pos) {
    case IGRAPH_EIGEN_LM:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_lm(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SM:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_sm(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_LA:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_la(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SA:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_sa(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_BE:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_be(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_ALL:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_all(myA,
                     values,
                     vectors));
        break;
    case IGRAPH_EIGEN_INTERVAL:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_iv(myA, which,
                     values,
                     vectors));
        break;
    case IGRAPH_EIGEN_SELECT:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack_sel(myA, which,
                     values,
                     vectors));
        break;
    default:
        /* This cannot happen */
        break;
    }

    if (!A) {
        igraph_matrix_destroy(&mA);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

typedef struct igraph_i_eigen_matrix_sym_arpack_data_t {
    const igraph_matrix_t *A;
    const igraph_sparsemat_t *sA;
} igraph_i_eigen_matrix_sym_arpack_data_t;

static int igraph_i_eigen_matrix_sym_arpack_cb(igraph_real_t *to,
                                        const igraph_real_t *from,
                                        int n, void *extra) {

    igraph_i_eigen_matrix_sym_arpack_data_t *data =
        (igraph_i_eigen_matrix_sym_arpack_data_t *) extra;

    if (data->A) {
        igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0,
                                               data->A, from, /*beta=*/ 0.0, to);
    } else { /* data->sA */
        igraph_vector_t vto, vfrom;
        igraph_vector_view(&vto, to, n);
        igraph_vector_view(&vfrom, from, n);
        igraph_vector_null(&vto);
        igraph_sparsemat_gaxpy(data->sA, &vfrom, &vto);
    }
    return 0;
}

static int igraph_i_eigen_matrix_symmetric_arpack_be(const igraph_matrix_t *A,
        const igraph_sparsemat_t *sA,
        igraph_arpack_function_t *fun,
        int n, void *extra,
        const igraph_eigen_which_t *which,
        igraph_arpack_options_t *options,
        igraph_arpack_storage_t *storage,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    igraph_vector_t tmpvalues, tmpvalues2;
    igraph_matrix_t tmpvectors, tmpvectors2;

    int low = (int) floor(which->howmany / 2.0), high = (int) ceil(which->howmany / 2.0);
    int l1, l2, w;

    igraph_i_eigen_matrix_sym_arpack_data_t myextra;
    myextra.A = A;
    myextra.sA = sA;

    if (low + high >= n) {
        IGRAPH_ERROR("Requested too many eigenvalues/vectors", IGRAPH_EINVAL);
    }

    if (!fun) {
        fun = igraph_i_eigen_matrix_sym_arpack_cb;
        extra = (void*) &myextra;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&tmpvalues, high);
    IGRAPH_MATRIX_INIT_FINALLY(&tmpvectors, n, high);
    IGRAPH_VECTOR_INIT_FINALLY(&tmpvalues2, low);
    IGRAPH_MATRIX_INIT_FINALLY(&tmpvectors2, n, low);

    options->n = n;
    options->nev = high;
    options->ncv = 2 * options->nev < n ? 2 * options->nev : n;
    options->which[0] = 'L'; options->which[1] = 'A';

    IGRAPH_CHECK(igraph_arpack_rssolve(fun, extra, options, storage,
                                       &tmpvalues, &tmpvectors));

    options->nev = low;
    options->ncv = 2 * options->nev < n ? 2 * options->nev : n;
    options->which[0] = 'S'; options->which[1] = 'A';

    IGRAPH_CHECK(igraph_arpack_rssolve(fun, extra, options, storage,
                                       &tmpvalues2, &tmpvectors2));

    IGRAPH_CHECK(igraph_vector_resize(values, low + high));
    IGRAPH_CHECK(igraph_matrix_resize(vectors, n, low + high));

    l1 = 0; l2 = 0; w = 0;
    while (w < which->howmany) {
        VECTOR(*values)[w] = VECTOR(tmpvalues)[l1];
        memcpy(&MATRIX(*vectors, 0, w), &MATRIX(tmpvectors, 0, l1),
               (size_t) n * sizeof(igraph_real_t));
        w++; l1++;
        if (w < which->howmany) {
            VECTOR(*values)[w] = VECTOR(tmpvalues2)[l2];
            memcpy(&MATRIX(*vectors, 0, w), &MATRIX(tmpvectors2, 0, l2),
                   (size_t) n * sizeof(igraph_real_t));
            w++; l2++;
        }
    }

    igraph_matrix_destroy(&tmpvectors2);
    igraph_vector_destroy(&tmpvalues2);
    igraph_matrix_destroy(&tmpvectors);
    igraph_vector_destroy(&tmpvalues);
    IGRAPH_FINALLY_CLEAN(4);

    return 0;
}

static int igraph_i_eigen_matrix_symmetric_arpack(const igraph_matrix_t *A,
        const igraph_sparsemat_t *sA,
        igraph_arpack_function_t *fun,
        int n, void *extra,
        const igraph_eigen_which_t *which,
        igraph_arpack_options_t *options,
        igraph_arpack_storage_t *storage,
        igraph_vector_t *values,
        igraph_matrix_t *vectors) {

    /* For ARPACK we need a matrix multiplication operation.
       This can be done in any format, so everything is fine,
       we don't have to convert. */

    igraph_i_eigen_matrix_sym_arpack_data_t myextra;

    myextra.A = A;
    myextra.sA = sA;

    if (!options) {
        IGRAPH_ERROR("`options' must be given for ARPACK algorithm",
                     IGRAPH_EINVAL);
    }

    if (which->pos == IGRAPH_EIGEN_BE) {
        return igraph_i_eigen_matrix_symmetric_arpack_be(A, sA, fun, n, extra,
                which, options, storage,
                values, vectors);
    } else {

        switch (which->pos) {
        case IGRAPH_EIGEN_LM:
            options->which[0] = 'L'; options->which[1] = 'M';
            options->nev = which->howmany;
            break;
        case IGRAPH_EIGEN_SM:
            options->which[0] = 'S'; options->which[1] = 'M';
            options->nev = which->howmany;
            break;
        case IGRAPH_EIGEN_LA:
            options->which[0] = 'L'; options->which[1] = 'A';
            options->nev = which->howmany;
            break;
        case IGRAPH_EIGEN_SA:
            options->which[0] = 'S'; options->which[1] = 'A';
            options->nev = which->howmany;
            break;
        case IGRAPH_EIGEN_ALL:
            options->which[0] = 'L'; options->which[1] = 'M';
            options->nev = n;
            break;
        case IGRAPH_EIGEN_INTERVAL:
            IGRAPH_ERROR("Interval of eigenvectors with ARPACK",
                         IGRAPH_UNIMPLEMENTED);
            /* TODO */
            break;
        case IGRAPH_EIGEN_SELECT:
            IGRAPH_ERROR("Selected eigenvalues with ARPACK",
                         IGRAPH_UNIMPLEMENTED);
            /* TODO */
            break;
        default:
            /* This cannot happen */
            break;
        }

        options->n = n;
        options->ncv = 2 * options->nev < n ? 2 * options->nev : n;

        if (!fun) {
            fun = igraph_i_eigen_matrix_sym_arpack_cb;
            extra = (void*) &myextra;
        }

        IGRAPH_CHECK(igraph_arpack_rssolve(fun, extra, options, storage,
                                           values, vectors));
        return 0;
    }
}

/* Get the eigenvalues and the eigenvectors from the compressed
   form. Order them according to the ordering criteria.
   Comparison functions for the reordering first */

typedef int (*igraph_i_eigen_matrix_lapack_cmp_t)(void*, const void*,
        const void *);

typedef struct igraph_i_eml_cmp_t {
    const igraph_vector_t *mag, *real, *imag;
} igraph_i_eml_cmp_t;

/* TODO: these should be defined in some header */

#define EPS        (DBL_EPSILON*100)
#define LESS(a,b)  ((a) < (b)-EPS)
#define MORE(a,b)  ((a) > (b)+EPS)
#define ZERO(a)    ((a) > -EPS && (a) < EPS)
#define NONZERO(a) ((a) < -EPS || (a) > EPS)

/* Largest magnitude. Ordering is according to
   1 Larger magnitude
   2 Real eigenvalues before complex ones
   3 Larger real part
   4 Larger imaginary part */

static int igraph_i_eigen_matrix_lapack_cmp_lm(void *extra, const void *a,
                                        const void *b) {
    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_m = VECTOR(*myextra->mag)[*aa];
    igraph_real_t b_m = VECTOR(*myextra->mag)[*bb];

    if (LESS(a_m, b_m)) {
        return 1;
    } else if (MORE(a_m, b_m)) {
        return -1;
    } else {
        igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
        igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
        igraph_real_t b_r = VECTOR(*myextra->real)[*bb];
        igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];
        if (ZERO(a_i)    && NONZERO(b_i))  {
            return -1;
        }
        if (NONZERO(a_i) && ZERO(b_i))     {
            return  1;
        }
        if (MORE(a_r, b_r)) {
            return -1;
        }
        if (LESS(a_r, b_r)) {
            return  1;
        }
        if (MORE(a_i, b_i)) {
            return -1;
        }
        if (LESS(a_i, b_i)) {
            return  1;
        }
    }
    return 0;
}

/* Smallest marginude. Ordering is according to
   1 Magnitude (smaller first)
   2 Complex eigenvalues before real ones
   3 Smaller real part
   4 Smaller imaginary part
   This ensures that lm has exactly the opposite order to sm */

static int igraph_i_eigen_matrix_lapack_cmp_sm(void *extra, const void *a,
                                        const void *b) {
    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_m = VECTOR(*myextra->mag)[*aa];
    igraph_real_t b_m = VECTOR(*myextra->mag)[*bb];

    if (MORE(a_m, b_m)) {
        return 1;
    } else if (LESS(a_m, b_m)) {
        return -1;
    } else {
        igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
        igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
        igraph_real_t b_r = VECTOR(*myextra->real)[*bb];
        igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];
        if (NONZERO(a_i) && ZERO(b_i))    {
            return -1;
        }
        if (ZERO(a_i)    && NONZERO(b_i)) {
            return  1;
        }
        if (LESS(a_r, b_r)) {
            return -1;
        }
        if (MORE(a_r, b_r)) {
            return  1;
        }
        if (LESS(a_i, b_i)) {
            return -1;
        }
        if (MORE(a_i, b_i)) {
            return  1;
        }
    }
    return 0;
}

/* Largest real part. Ordering is according to
   1 Larger real part
   2 Real eigenvalues come before complex ones
   3 Larger complex part */

static int igraph_i_eigen_matrix_lapack_cmp_lr(void *extra, const void *a,
                                        const void *b) {

    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
    igraph_real_t b_r = VECTOR(*myextra->real)[*bb];

    if (MORE(a_r, b_r)) {
        return -1;
    } else if (LESS(a_r, b_r)) {
        return 1;
    } else {
        igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
        igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];
        if (ZERO(a_i) && NONZERO(b_i)) {
            return -1;
        }
        if (NONZERO(a_i) && ZERO(b_i)) {
            return  1;
        }
        if (MORE(a_i, b_i)) {
            return -1;
        }
        if (LESS(a_i, b_i)) {
            return  1;
        }
    }

    return 0;
}

/* Largest real part. Ordering is according to
   1 Smaller real part
   2 Complex eigenvalues come before real ones
   3 Smaller complex part
   This is opposite to LR
*/

static int igraph_i_eigen_matrix_lapack_cmp_sr(void *extra, const void *a,
                                        const void *b) {

    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
    igraph_real_t b_r = VECTOR(*myextra->real)[*bb];

    if (LESS(a_r, b_r)) {
        return -1;
    } else if (MORE(a_r, b_r)) {
        return 1;
    } else {
        igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
        igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];
        if (NONZERO(a_i) && ZERO(b_i)) {
            return -1;
        }
        if (ZERO(a_i) && NONZERO(b_i)) {
            return  1;
        }
        if (LESS(a_i, b_i)) {
            return -1;
        }
        if (MORE(a_i, b_i)) {
            return  1;
        }
    }

    return 0;
}

/* Order:
   1 Larger imaginary part
   2 Real eigenvalues before complex ones
   3 Larger real part */

static int igraph_i_eigen_matrix_lapack_cmp_li(void *extra, const void *a,
                                        const void *b) {

    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
    igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];

    if (MORE(a_i, b_i)) {
        return -1;
    } else if (LESS(a_i, b_i)) {
        return 1;
    } else {
        igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
        igraph_real_t b_r = VECTOR(*myextra->real)[*bb];
        if (ZERO(a_i) && NONZERO(b_i)) {
            return -1;
        }
        if (NONZERO(a_i) && ZERO(b_i)) {
            return  1;
        }
        if (MORE(a_r, b_r)) {
            return -1;
        }
        if (LESS(a_r, b_r)) {
            return  1;
        }
    }

    return 0;
}

/* Order:
   1 Smaller imaginary part
   2 Complex eigenvalues before real ones
   3 Smaller real part
   Order is opposite to LI */

static int igraph_i_eigen_matrix_lapack_cmp_si(void *extra, const void *a,
                                        const void *b) {

    igraph_i_eml_cmp_t *myextra = (igraph_i_eml_cmp_t *) extra;
    int *aa = (int*) a, *bb = (int*) b;
    igraph_real_t a_i = VECTOR(*myextra->imag)[*aa];
    igraph_real_t b_i = VECTOR(*myextra->imag)[*bb];

    if (LESS(a_i, b_i)) {
        return -1;
    } else if (MORE(a_i, b_i)) {
        return 1;
    } else {
        igraph_real_t a_r = VECTOR(*myextra->real)[*aa];
        igraph_real_t b_r = VECTOR(*myextra->real)[*bb];
        if (NONZERO(a_i) && ZERO(b_i)) {
            return -1;
        }
        if (ZERO(a_i) && NONZERO(b_i)) {
            return  1;
        }
        if (LESS(a_r, b_r)) {
            return -1;
        }
        if (MORE(a_r, b_r)) {
            return  1;
        }
    }

    return 0;
}

#undef EPS
#undef LESS
#undef MORE
#undef ZERO
#undef NONZERO

#define INITMAG()                           \
    do {                                  \
        int i;                              \
        IGRAPH_VECTOR_INIT_FINALLY(&mag, nev);              \
        hasmag=1;                               \
        for (i=0; i<nev; i++) {                     \
            VECTOR(mag)[i] = VECTOR(*real)[i] * VECTOR(*real)[i] +        \
                             VECTOR(*imag)[i] * VECTOR(*imag)[i];                \
        }                                   \
    } while (0)

static int igraph_i_eigen_matrix_lapack_reorder(const igraph_vector_t *real,
        const igraph_vector_t *imag,
        const igraph_matrix_t *compressed,
        const igraph_eigen_which_t *which,
        igraph_vector_complex_t *values,
        igraph_matrix_complex_t *vectors) {
    igraph_vector_int_t idx;
    igraph_vector_t mag;
    igraph_bool_t hasmag = 0;
    int nev = (int) igraph_vector_size(real);
    int howmany = 0, start = 0;
    int i;
    igraph_i_eigen_matrix_lapack_cmp_t cmpfunc = 0;
    igraph_i_eml_cmp_t vextra;
    void *extra;

    vextra.mag = &mag;
    vextra.real = real;
    vextra.imag = imag;
    extra = &vextra;

    IGRAPH_CHECK(igraph_vector_int_init(&idx, nev));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &idx);

    switch (which->pos) {
    case IGRAPH_EIGEN_LM:
        INITMAG();
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_lm;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_ALL:
        INITMAG();
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_sm;
        howmany = nev;
        break;
    case IGRAPH_EIGEN_SM:
        INITMAG();
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_sm;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_LR:
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_lr;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_SR:
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_sr;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_SELECT:
        INITMAG();
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_sm;
        start = which->il - 1;
        howmany = which->iu - which->il + 1;
        break;
    case IGRAPH_EIGEN_LI:
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_li;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_SI:
        cmpfunc = igraph_i_eigen_matrix_lapack_cmp_si;
        howmany = which->howmany;
        break;
    case IGRAPH_EIGEN_INTERVAL:
    case IGRAPH_EIGEN_BE:
    default:
        IGRAPH_ERROR("Unimplemented eigenvalue ordering", IGRAPH_UNIMPLEMENTED);
        break;
    }

    for (i = 0; i < nev; i++) {
        VECTOR(idx)[i] = i;
    }

    igraph_qsort_r(VECTOR(idx), (size_t) nev, sizeof(VECTOR(idx)[0]), extra,
                   cmpfunc);

    if (hasmag) {
        igraph_vector_destroy(&mag);
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (values) {
        IGRAPH_CHECK(igraph_vector_complex_resize(values, howmany));
        for (i = 0; i < howmany; i++) {
            int x = VECTOR(idx)[start + i];
            VECTOR(*values)[i] = igraph_complex(VECTOR(*real)[x],
                                                VECTOR(*imag)[x]);
        }
    }

    if (vectors) {
        int n = (int) igraph_matrix_nrow(compressed);
        IGRAPH_CHECK(igraph_matrix_complex_resize(vectors, n, howmany));
        for (i = 0; i < howmany; i++) {
            int j, x = VECTOR(idx)[start + i];
            if (VECTOR(*imag)[x] == 0) {
                /* real eigenvalue */
                for (j = 0; j < n; j++) {
                    MATRIX(*vectors, j, i) = igraph_complex(MATRIX(*compressed, j, x),
                                                            0.0);
                }
            } else {
                /* complex eigenvalue */
                int neg = 1, co = 0;
                if (VECTOR(*imag)[x] < 0) {
                    neg = -1;
                    co = 1;
                }
                for (j = 0; j < n; j++) {
                    MATRIX(*vectors, j, i) =
                        igraph_complex(MATRIX(*compressed, j, x - co),
                                       neg * MATRIX(*compressed, j, x + 1 - co));
                }
            }
        }
    }

    igraph_vector_int_destroy(&idx);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_eigen_matrix_lapack_common(const igraph_matrix_t *A,
                                        const igraph_eigen_which_t *which,
                                        igraph_vector_complex_t *values,
                                        igraph_matrix_complex_t *vectors) {

    igraph_vector_t valuesreal, valuesimag;
    igraph_matrix_t vectorsright, *myvectors = vectors ? &vectorsright : 0;
    int n = (int) igraph_matrix_nrow(A);
    int info = 1;

    IGRAPH_VECTOR_INIT_FINALLY(&valuesreal, n);
    IGRAPH_VECTOR_INIT_FINALLY(&valuesimag, n);
    if (vectors) {
        IGRAPH_MATRIX_INIT_FINALLY(&vectorsright, n, n);
    }
    IGRAPH_CHECK(igraph_lapack_dgeev(A, &valuesreal, &valuesimag,
                                     /*vectorsleft=*/ 0, myvectors, &info));

    IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_reorder(&valuesreal,
                 &valuesimag,
                 myvectors, which, values,
                 vectors));

    if (vectors) {
        igraph_matrix_destroy(&vectorsright);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_vector_destroy(&valuesimag);
    igraph_vector_destroy(&valuesreal);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;

}

static int igraph_i_eigen_matrix_lapack_lm(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_sm(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_lr(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}


static int igraph_i_eigen_matrix_lapack_sr(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_li(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_si(const igraph_matrix_t *A,
                                    const igraph_eigen_which_t *which,
                                    igraph_vector_complex_t *values,
                                    igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_select(const igraph_matrix_t *A,
                                        const igraph_eigen_which_t *which,
                                        igraph_vector_complex_t *values,
                                        igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack_all(const igraph_matrix_t *A,
                                     const igraph_eigen_which_t *which,
                                     igraph_vector_complex_t *values,
                                     igraph_matrix_complex_t *vectors) {
    return igraph_i_eigen_matrix_lapack_common(A, which, values, vectors);
}

static int igraph_i_eigen_matrix_lapack(const igraph_matrix_t *A,
                                 const igraph_sparsemat_t *sA,
                                 igraph_arpack_function_t *fun,
                                 int n, void *extra,
                                 const igraph_eigen_which_t *which,
                                 igraph_vector_complex_t *values,
                                 igraph_matrix_complex_t *vectors) {

    const igraph_matrix_t *myA = A;
    igraph_matrix_t mA;

    /* We need to create a dense square matrix first */

    if (A) {
        n = (int) igraph_matrix_nrow(A);
    } else if (sA) {
        n = (int) igraph_sparsemat_nrow(sA);
        IGRAPH_CHECK(igraph_matrix_init(&mA, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &mA);
        IGRAPH_CHECK(igraph_sparsemat_as_matrix(&mA, sA));
        myA = &mA;
    } else if (fun) {
        IGRAPH_CHECK(igraph_i_eigen_arpackfun_to_mat(fun, n, extra, &mA));
        IGRAPH_FINALLY(igraph_matrix_destroy, &mA);
    }

    switch (which->pos) {
    case IGRAPH_EIGEN_LM:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_lm(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SM:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_sm(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_LR:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_lr(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SR:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_sr(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_LI:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_li(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SI:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_si(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_SELECT:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_select(myA, which,
                     values, vectors));
        break;
    case IGRAPH_EIGEN_ALL:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack_all(myA, which,
                     values,
                     vectors));
        break;
    default:
        /* This cannot happen */
        break;
    }

    if (!A) {
        igraph_matrix_destroy(&mA);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

static int igraph_i_eigen_checks(const igraph_matrix_t *A,
                          const igraph_sparsemat_t *sA,
                          igraph_arpack_function_t *fun, int n) {

    if ( (A ? 1 : 0) + (sA ? 1 : 0) + (fun ? 1 : 0) != 1) {
        IGRAPH_ERROR("Exactly one of 'A', 'sA' and 'fun' must be given",
                     IGRAPH_EINVAL);
    }

    if (A) {
        if (n != igraph_matrix_ncol(A) || n != igraph_matrix_nrow(A)) {
            IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
        }
    } else if (sA) {
        if (n != igraph_sparsemat_ncol(sA) || n != igraph_sparsemat_nrow(sA)) {
            IGRAPH_ERROR("Invalid matrix", IGRAPH_NONSQUARE);
        }
    }

    return 0;
}

/**
 * \function igraph_eigen_matrix_symmetric
 *
 * \example examples/simple/igraph_eigen_matrix_symmetric.c
 */

int igraph_eigen_matrix_symmetric(const igraph_matrix_t *A,
                                  const igraph_sparsemat_t *sA,
                                  igraph_arpack_function_t *fun, int n,
                                  void *extra,
                                  igraph_eigen_algorithm_t algorithm,
                                  const igraph_eigen_which_t *which,
                                  igraph_arpack_options_t *options,
                                  igraph_arpack_storage_t *storage,
                                  igraph_vector_t *values,
                                  igraph_matrix_t *vectors) {

    IGRAPH_CHECK(igraph_i_eigen_checks(A, sA, fun, n));

    if (which->pos != IGRAPH_EIGEN_LM &&
        which->pos != IGRAPH_EIGEN_SM &&
        which->pos != IGRAPH_EIGEN_LA &&
        which->pos != IGRAPH_EIGEN_SA &&
        which->pos != IGRAPH_EIGEN_BE &&
        which->pos != IGRAPH_EIGEN_ALL &&
        which->pos != IGRAPH_EIGEN_INTERVAL &&
        which->pos != IGRAPH_EIGEN_SELECT) {
        IGRAPH_ERROR("Invalid 'pos' position in 'which'", IGRAPH_EINVAL);
    }

    switch (algorithm) {
    case IGRAPH_EIGEN_AUTO:
        if (which->howmany == n || n < 100) {
            IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack(A, sA, fun, n,
                         extra, which,
                         values, vectors));
        } else {
            IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_arpack(A, sA, fun, n,
                         extra, which,
                         options, storage,
                         values, vectors));
        }
        break;
    case IGRAPH_EIGEN_LAPACK:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_lapack(A, sA, fun, n, extra,
                     which, values,
                     vectors));
        break;
    case IGRAPH_EIGEN_ARPACK:
        IGRAPH_CHECK(igraph_i_eigen_matrix_symmetric_arpack(A, sA, fun, n, extra,
                     which, options,
                     storage,
                     values, vectors));
        break;
    default:
        IGRAPH_ERROR("Unknown 'algorithm'", IGRAPH_EINVAL);
    }

    return 0;
}

/**
 * \function igraph_eigen_matrix
 *
 */

int igraph_eigen_matrix(const igraph_matrix_t *A,
                        const igraph_sparsemat_t *sA,
                        igraph_arpack_function_t *fun, int n,
                        void *extra,
                        igraph_eigen_algorithm_t algorithm,
                        const igraph_eigen_which_t *which,
                        igraph_arpack_options_t *options,
                        igraph_arpack_storage_t *storage,
                        igraph_vector_complex_t *values,
                        igraph_matrix_complex_t *vectors) {

    IGRAPH_UNUSED(options);
    IGRAPH_UNUSED(storage);

    IGRAPH_CHECK(igraph_i_eigen_checks(A, sA, fun, n));

    if (which->pos != IGRAPH_EIGEN_LM &&
        which->pos != IGRAPH_EIGEN_SM &&
        which->pos != IGRAPH_EIGEN_LR &&
        which->pos != IGRAPH_EIGEN_SR &&
        which->pos != IGRAPH_EIGEN_LI &&
        which->pos != IGRAPH_EIGEN_SI &&
        which->pos != IGRAPH_EIGEN_SELECT &&
        which->pos != IGRAPH_EIGEN_ALL) {
        IGRAPH_ERROR("Invalid 'pos' position in 'which'", IGRAPH_EINVAL);
    }

    switch (algorithm) {
    case IGRAPH_EIGEN_AUTO:
        IGRAPH_ERROR("'AUTO' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_LAPACK:
        IGRAPH_CHECK(igraph_i_eigen_matrix_lapack(A, sA, fun, n, extra, which,
                     values, vectors));
        /* TODO */
        break;
    case IGRAPH_EIGEN_ARPACK:
        IGRAPH_ERROR("'ARPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_COMP_AUTO:
        IGRAPH_ERROR("'COMP_AUTO' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_COMP_LAPACK:
        IGRAPH_ERROR("'COMP_LAPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_COMP_ARPACK:
        IGRAPH_ERROR("'COMP_ARPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    default:
        IGRAPH_ERROR("Unknown `algorithm'", IGRAPH_EINVAL);
    }

    return 0;
}

static int igraph_i_eigen_adjacency_arpack_sym_cb(igraph_real_t *to,
        const igraph_real_t *from,
        int n, void *extra) {
    igraph_adjlist_t *adjlist = (igraph_adjlist_t *) extra;
    igraph_vector_int_t *neis;
    int i, j, nlen;

    for (i = 0; i < n; i++) {
        neis = igraph_adjlist_get(adjlist, i);
        nlen = igraph_vector_int_size(neis);
        to[i] = 0.0;
        for (j = 0; j < nlen; j++) {
            int nei = VECTOR(*neis)[j];
            to[i] += from[nei];
        }
    }

    return 0;
}

static int igraph_i_eigen_adjacency_arpack(const igraph_t *graph,
                                    const igraph_eigen_which_t *which,
                                    igraph_arpack_options_t *options,
                                    igraph_arpack_storage_t* storage,
                                    igraph_vector_t *values,
                                    igraph_matrix_t *vectors,
                                    igraph_vector_complex_t *cmplxvalues,
                                    igraph_matrix_complex_t *cmplxvectors) {

    IGRAPH_UNUSED(cmplxvalues);
    IGRAPH_UNUSED(cmplxvectors);

    igraph_adjlist_t adjlist;
    void *extra = (void*) &adjlist;
    int n = igraph_vcount(graph);

    if (!options) {
        IGRAPH_ERROR("`options' must be given for ARPACK algorithm",
                     IGRAPH_EINVAL);
    }

    if (igraph_is_directed(graph)) {
        IGRAPH_ERROR("ARPACK adjacency eigensolver not implemented for "
                     "directed graphs", IGRAPH_UNIMPLEMENTED);
    }
    if (which->pos == IGRAPH_EIGEN_INTERVAL) {
        IGRAPH_ERROR("ARPACK adjacency eigensolver does not implement "
                     "`INTERNAL' eigenvalues", IGRAPH_UNIMPLEMENTED);
    }
    if (which->pos == IGRAPH_EIGEN_SELECT) {
        IGRAPH_ERROR("ARPACK adjacency eigensolver does not implement "
                     "`SELECT' eigenvalues", IGRAPH_UNIMPLEMENTED);
    }
    if (which->pos == IGRAPH_EIGEN_ALL) {
        IGRAPH_ERROR("ARPACK adjacency eigensolver does not implement "
                     "`ALL' eigenvalues", IGRAPH_UNIMPLEMENTED);
    }

    switch (which->pos) {
    case IGRAPH_EIGEN_LM:
        options->which[0] = 'L'; options->which[1] = 'M';
        options->nev = which->howmany;
        break;
    case IGRAPH_EIGEN_SM:
        options->which[0] = 'S'; options->which[1] = 'M';
        options->nev = which->howmany;
        break;
    case IGRAPH_EIGEN_LA:
        options->which[0] = 'L'; options->which[1] = 'A';
        options->nev = which->howmany;
        break;
    case IGRAPH_EIGEN_SA:
        options->which[0] = 'S'; options->which[1] = 'A';
        options->nev = which->howmany;
        break;
    case IGRAPH_EIGEN_ALL:
        options->which[0] = 'L'; options->which[1] = 'M';
        options->nev = n;
        break;
    case IGRAPH_EIGEN_BE:
        IGRAPH_ERROR("Eigenvectors from both ends with ARPACK",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_INTERVAL:
        IGRAPH_ERROR("Interval of eigenvectors with ARPACK",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_SELECT:
        IGRAPH_ERROR("Selected eigenvalues with ARPACK",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    default:
        /* This cannot happen */
        break;
    }

    options->n = n;
    options->ncv = 2 * options->nev < n ? 2 * options->nev : n;

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_IN, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_eigen_adjacency_arpack_sym_cb,
                                       extra, options, storage, values, vectors));

    igraph_adjlist_destroy(&adjlist);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_eigen_adjacency
 *
 */

int igraph_eigen_adjacency(const igraph_t *graph,
                           igraph_eigen_algorithm_t algorithm,
                           const igraph_eigen_which_t *which,
                           igraph_arpack_options_t *options,
                           igraph_arpack_storage_t *storage,
                           igraph_vector_t *values,
                           igraph_matrix_t *vectors,
                           igraph_vector_complex_t *cmplxvalues,
                           igraph_matrix_complex_t *cmplxvectors) {

    if (which->pos != IGRAPH_EIGEN_LM &&
        which->pos != IGRAPH_EIGEN_SM &&
        which->pos != IGRAPH_EIGEN_LA &&
        which->pos != IGRAPH_EIGEN_SA &&
        which->pos != IGRAPH_EIGEN_BE &&
        which->pos != IGRAPH_EIGEN_SELECT &&
        which->pos != IGRAPH_EIGEN_INTERVAL &&
        which->pos != IGRAPH_EIGEN_ALL) {
        IGRAPH_ERROR("Invalid 'pos' position in 'which'", IGRAPH_EINVAL);
    }

    switch (algorithm) {
    case IGRAPH_EIGEN_AUTO:
        IGRAPH_ERROR("'AUTO' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_LAPACK:
        IGRAPH_ERROR("'LAPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_ARPACK:
        IGRAPH_CHECK(igraph_i_eigen_adjacency_arpack(graph, which, options,
                     storage, values, vectors,
                     cmplxvalues,
                     cmplxvectors));
        break;
    case IGRAPH_EIGEN_COMP_AUTO:
        IGRAPH_ERROR("'COMP_AUTO' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_COMP_LAPACK:
        IGRAPH_ERROR("'COMP_LAPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    case IGRAPH_EIGEN_COMP_ARPACK:
        IGRAPH_ERROR("'COMP_ARPACK' algorithm not implemented yet",
                     IGRAPH_UNIMPLEMENTED);
        /* TODO */
        break;
    default:
        IGRAPH_ERROR("Unknown `algorithm'", IGRAPH_EINVAL);
    }


    return 0;
}

/**
 * \function igraph_eigen_laplacian
 *
 */

int igraph_eigen_laplacian(const igraph_t *graph,
                           igraph_eigen_algorithm_t algorithm,
                           const igraph_eigen_which_t *which,
                           igraph_arpack_options_t *options,
                           igraph_arpack_storage_t *storage,
                           igraph_vector_t *values,
                           igraph_matrix_t *vectors,
                           igraph_vector_complex_t *cmplxvalues,
                           igraph_matrix_complex_t *cmplxvectors) {

    IGRAPH_UNUSED(graph);
    IGRAPH_UNUSED(algorithm);
    IGRAPH_UNUSED(which);
    IGRAPH_UNUSED(options);
    IGRAPH_UNUSED(storage);
    IGRAPH_UNUSED(values);
    IGRAPH_UNUSED(vectors);
    IGRAPH_UNUSED(cmplxvalues);
    IGRAPH_UNUSED(cmplxvectors);

    /* TODO */

    IGRAPH_ERROR("'igraph_eigen_laplacian'", IGRAPH_UNIMPLEMENTED);
}
