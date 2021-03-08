/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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
#include <pthread.h>
#include <string.h>

#include "linalg/arpack_internal.h"

#include "test_utilities.inc"

/* Test whether ARPACK is thread-safe. We will create two threads,
   each calling a different ARPACK eigensolver. We will make sure that
   the ARPACK calls from the two threads overlap */

typedef struct thread_data_t {
    igraph_matrix_t *m;
    igraph_vector_t *result;
    pthread_cond_t *cond;
    pthread_mutex_t *mutex;
    int *steps, *othersteps;
} thread_data_t;

int arpack_mult(igraph_real_t *to, igraph_real_t *from, int n,
                igraph_matrix_t *matrix) {
    /* TODO */
    igraph_blas_dgemv_array(/*transpose=*/ 0, /*alpha=*/ 1.0, matrix,
                                           from, /*beta=*/ 0.0, to);

    return 0;
}

/* This is the function performed by each thread. It calls the
   low-level ARPACK symmetric eigensolver, step by step. After each
   step, it synchronizes with the other thread.

   The synchronization ensures that the two threads are using the
   thread-local variables at the same time. If they are really
   thread-local, then ARPACK still delivers the correct solution for
   the two matrices. Otherwise the result is undefined: maybe results
   will be incorrect, or the program will crash.

   This function is basically a simplified copy of igraph_arpack_rssolve.
*/

void *thread_function(void *arg) {
    thread_data_t *data = (thread_data_t*) arg;
    igraph_matrix_t *M = data->m;
    igraph_vector_t *result = data->result;
    pthread_cond_t *cond = data->cond;
    pthread_mutex_t *mutex = data->mutex;
    igraph_arpack_options_t options;
    igraph_real_t *v, *workl, *workd, *d, *resid, *ax;
    int *select;
    int ido = 0;
#if IGRAPH_THREAD_SAFE
    int rvec = 1;
    char *all = "All";
#endif
    int i;

    igraph_arpack_options_init(&options);
    options.n = igraph_matrix_nrow(M);
    options.ldv = options.n;
    options.nev = 1;
    options.ncv = 3;
    options.lworkl = options.ncv * (options.ncv + 8);
    options.which[0] = 'L';
    options.which[1] = 'M';

    options.iparam[0] = options.ishift;
    options.iparam[2] = options.mxiter;
    options.iparam[3] = options.nb;
    options.iparam[4] = 0;
    options.iparam[6] = options.mode;
    options.info = options.start;

    v = IGRAPH_CALLOC(options.ldv * options.ncv, igraph_real_t);
    workl = IGRAPH_CALLOC(options.lworkl, igraph_real_t);
    workd = IGRAPH_CALLOC(3 * options.n, igraph_real_t);
    d = IGRAPH_CALLOC(2 * options.ncv, igraph_real_t);
    resid = IGRAPH_CALLOC(options.n, igraph_real_t);
    ax = IGRAPH_CALLOC(options.n, igraph_real_t);
    select = IGRAPH_CALLOC(options.ncv, int);

    if (!v || !workl || !workd || !d || !resid || !ax || !select) {
        printf("Out of memory\n");
        return 0;
    }

    while (1) {
#if IGRAPH_THREAD_SAFE
        igraphdsaupd_(&ido, options.bmat, &options.n, options.which,
                      &options.nev, &options.tol, resid, &options.ncv, v,
                      &options.ldv, options.iparam, options.ipntr, workd,
                      workl, &options.lworkl, &options.info);
#endif

        if (ido == -1 || ido == 1) {

            igraph_real_t *from = workd + options.ipntr[0] - 1;
            igraph_real_t *to = workd + options.ipntr[1] - 1;
            arpack_mult(to, from, options.n, M);

        } else {
            break;
        }

        pthread_mutex_lock(mutex);
        *(data->steps) += 1;
        if ( *(data->othersteps) == *(data->steps) ) {
            pthread_cond_signal(cond);
        }

        while ( *(data->othersteps) < * (data->steps) && *(data->othersteps) != -1 ) {
            pthread_cond_wait(cond, mutex);
        }
        pthread_mutex_unlock(mutex);
    }

    pthread_mutex_lock(mutex);
    *data->steps = -1;
    pthread_cond_signal(cond);
    pthread_mutex_unlock(mutex);

    if (options.info != 0) {
        printf("ARPACK error\n");
        return 0;
    }

#if IGRAPH_THREAD_SAFE
    igraphdseupd_(&rvec, all, select, d, v, &options.ldv,
                  &options.sigma, options.bmat, &options.n,
                  options.which, &options.nev, &options.tol,
                  resid, &options.ncv, v, &options.ldv, options.iparam,
                  options.ipntr, workd, workl, &options.lworkl,
                  &options.ierr);
#endif

    if (options.ierr != 0) {
        printf("ARPACK error\n");
        return 0;
    }

    igraph_vector_resize(result, options.n);
    for (i = 0; i < options.n; i++) {
        VECTOR(*result)[i] = v[i];
    }

    free(v);
    free(workl);
    free(workd);
    free(d);
    free(resid);
    free(ax);
    free(select);

    return 0;
}

int main() {
    pthread_t thread_id1, thread_id2;
    void *exit_status1, *exit_status2;
    igraph_matrix_t m1, m2;
    igraph_vector_t result1, result2;
    pthread_cond_t steps_cond = PTHREAD_COND_INITIALIZER;
    pthread_mutex_t steps_mutex = PTHREAD_MUTEX_INITIALIZER;
    int steps1 = 0, steps2 = 0;
    thread_data_t
    data1 = { &m1, &result1, &steps_cond, &steps_mutex, &steps1, &steps2 },
    data2 = { &m2, &result2, &steps_cond, &steps_mutex, &steps2, &steps1 };
    int i, j;

    /* Skip if igraph is not thread safe */
    if (!IGRAPH_THREAD_SAFE) {
        return 77;
    }

    igraph_matrix_init(&m1, 10, 10);
    igraph_matrix_init(&m2, 10, 10);
    igraph_vector_init(&result1, igraph_matrix_nrow(&m1));
    igraph_vector_init(&result2, igraph_matrix_nrow(&m2));

    igraph_rng_seed(igraph_rng_default(), 42);

    for (i = 0; i < igraph_matrix_nrow(&m1); i++) {
        for (j = 0; j <= i; j++) {
            MATRIX(m1, i, j) = MATRIX(m1, j, i) =
                                   igraph_rng_get_integer(igraph_rng_default(), 0, 10);
        }
    }

    for (i = 0; i < igraph_matrix_nrow(&m2); i++) {
        for (j = 0; j <= i; j++) {
            MATRIX(m2, i, j) = MATRIX(m2, j, i) =
                                   igraph_rng_get_integer(igraph_rng_default(), 0, 10);
        }
    }

    pthread_create(&thread_id1, NULL, thread_function, (void *) &data1);
    pthread_create(&thread_id2, NULL, thread_function, (void *) &data2);

    pthread_join(thread_id1, &exit_status1);
    pthread_join(thread_id2, &exit_status2);

    igraph_matrix_print(&m1);
    igraph_vector_print(&result1);
    printf("---\n");
    igraph_matrix_print(&m2);
    igraph_vector_print(&result2);

    igraph_vector_destroy(&result1);
    igraph_vector_destroy(&result2);
    igraph_matrix_destroy(&m1);
    igraph_matrix_destroy(&m2);

    VERIFY_FINALLY_STACK();

    return 0;
}
