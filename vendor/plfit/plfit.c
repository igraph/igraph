/* vim:set ts=4 sw=4 sts=4 et: */
/* plfit.c
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "gss.h"
#include "lbfgs.h"
#include "platform.h"
#include "plfit.h"
#include "kolmogorov.h"
#include "sampling.h"
#include "hzeta.h"

/* #define PLFIT_DEBUG */

#define DATA_POINTS_CHECK \
    if (n <= 0) { \
        PLFIT_ERROR("no data points", PLFIT_EINVAL); \
    }

#define XMIN_CHECK_ZERO \
    if (xmin <= 0) { \
        PLFIT_ERROR("xmin must be greater than zero", PLFIT_EINVAL); \
    }
#define XMIN_CHECK_ONE \
    if (xmin < 1) { \
        PLFIT_ERROR("xmin must be at least 1", PLFIT_EINVAL); \
    }

static int plfit_i_resample_continuous(double* xs_head, size_t num_smaller,
        size_t n, double alpha, double xmin, size_t num_samples, mt_rng_t* rng,
        double* result);
static int plfit_i_resample_discrete(double* xs_head, size_t num_smaller,
        size_t n, double alpha, double xmin, size_t num_samples, mt_rng_t* rng,
        double* result);

static int double_comparator(const void *a, const void *b) {
    const double *da = (const double*)a;
    const double *db = (const double*)b;
    return (*da > *db) - (*da < *db);
}

static int plfit_i_copy_and_sort(double* xs, size_t n, double** result) {
    *result = (double*)malloc(sizeof(double) * n);
    if (*result == 0) {
        PLFIT_ERROR("cannot create sorted copy of input data", PLFIT_ENOMEM);
    }

    memcpy(*result, xs, sizeof(double) * n);
    qsort(*result, n, sizeof(double), double_comparator);

    return PLFIT_SUCCESS;
}

/**
 * Given an unsorted array of doubles, counts how many elements there are that
 * are smaller than a given value.
 *
 * \param  begin          pointer to the beginning of the array
 * \param  end            pointer to the first element after the end of the array
 * \param  xmin           the threshold value
 *
 * \return the nubmer of elements in the array that are smaller than the given
 *         value.
 */
static size_t count_smaller(double* begin, double* end, double xmin) {
    double* p;
    size_t counter = 0;

    for (p = begin; p < end; p++) {
        if (*p < xmin) {
            counter++;
        }
    }

    return counter;
}

/**
 * Given an unsorted array of doubles, return another array that contains the
 * elements that are smaller than a given value
 *
 * \param  begin          pointer to the beginning of the array
 * \param  end            pointer to the first element after the end of the array
 * \param  xmin           the threshold value
 * \param  result_length  if not \c NULL, the number of unique elements in the
 *                        given array is returned here
 *
 * \return pointer to the head of the new array or 0 if there is not enough
 * memory
 */
static double* extract_smaller(double* begin, double* end, double xmin,
        size_t* result_length) {
    size_t counter = count_smaller(begin, end, xmin);
    double *p, *result;

    result = calloc(counter, sizeof(double));
    if (result == 0)
        return 0;

    for (p = result; begin < end; begin++) {
        if (*begin < xmin) {
            *p = *begin;
            p++;
        }
    }

    if (result_length) {
        *result_length = counter;
    }

    return result;
}

/**
 * Given a sorted array of doubles, return another array that contains pointers
 * into the array for the start of each block of identical elements.
 *
 * \param  begin          pointer to the beginning of the array
 * \param  end            pointer to the first element after the end of the array
 * \param  result_length  if not \c NULL, the number of unique elements in the
 *                        given array is returned here
 *
 * \return pointer to the head of the new array or 0 if there is not enough
 * memory
 */
static double** unique_element_pointers(double* begin, double* end, size_t* result_length) {
    double* ptr = begin;
    double** result;
    double prev_x;
    size_t num_elts = 15;
    size_t used_elts = 0;

    /* Special case: empty array */
    if (begin == end) {
        result = calloc(1, sizeof(double*));
        if (result != 0) {
            result[0] = 0;
        }
        return result;
    }

    /* Allocate initial result array, including the guard element */
    result = calloc(num_elts+1, sizeof(double*));
    if (result == 0)
        return 0;

    prev_x = *begin;
    result[used_elts++] = begin;

    /* Process the input array */
    for (ptr = begin+1; ptr < end; ptr++) {
        if (*ptr == prev_x)
            continue;

        /* New block found */
        if (used_elts >= num_elts) {
            /* Array full; allocate a new chunk */
            num_elts = num_elts*2 + 1;
            result = realloc(result, sizeof(double*) * (num_elts+1));
            if (result == 0)
                return 0;
        }

        /* Store the new element */
        result[used_elts++] = ptr;
        prev_x = *ptr;
    }

    /* Calculate the result length */
    if (result_length != 0) {
        *result_length = used_elts;
    }

    /* Add the guard entry to the end of the result */
    result[used_elts++] = 0;

    return result;
}

static void plfit_i_perform_finite_size_correction(plfit_result_t* result, size_t n) {
    result->alpha = result->alpha * (n-1) / n + 1.0 / n;
}

/********** Continuous power law distribution fitting **********/

static void plfit_i_logsum_less_than_continuous(double* begin, double* end,
        double xmin, double* result, size_t* m) {
    double logsum = 0.0;
    size_t count = 0;

    for (; begin != end; begin++) {
        if (*begin >= xmin) {
            count++;
            logsum += log(*begin / xmin);
        }
    }

    *m = count;
    *result = logsum;
}

static double plfit_i_logsum_continuous(double* begin, double* end, double xmin) {
    double logsum = 0.0;
    for (; begin != end; begin++)
        logsum += log(*begin / xmin);
    return logsum;
}

static int plfit_i_estimate_alpha_continuous(double* xs, size_t n,
        double xmin, double* alpha) {
    double result;
    size_t m;

    XMIN_CHECK_ZERO;

    plfit_i_logsum_less_than_continuous(xs, xs+n, xmin, &result, &m);

    if (m == 0) {
        PLFIT_ERROR("no data point was larger than xmin", PLFIT_EINVAL);
    }

    *alpha = 1 + m / result;

    return PLFIT_SUCCESS;
}

static int plfit_i_estimate_alpha_continuous_sorted(double* xs, size_t n,
        double xmin, double* alpha) {
    double* end = xs+n;

    XMIN_CHECK_ZERO;

    for (; xs != end && *xs < xmin; xs++);
    if (xs == end) {
        PLFIT_ERROR("no data point was larger than xmin", PLFIT_EINVAL);
    }

    *alpha = 1 + (end-xs) / plfit_i_logsum_continuous(xs, end, xmin);

    return PLFIT_SUCCESS;
}

static int plfit_i_ks_test_continuous(double* xs, double* xs_end,
        const double alpha, const double xmin, double* D) {
    /* Assumption: xs is sorted and cut off at xmin so the first element is
     * always larger than or equal to xmin. */
    double result = 0, n;
    int m = 0;

    n = xs_end - xs;

    while (xs < xs_end) {
        double d = fabs(1-pow(xmin / *xs, alpha-1) - m / n);

        if (d > result)
            result = d;

        xs++; m++;
    }

    *D = result;

    return PLFIT_SUCCESS;
}

static int plfit_i_calculate_p_value_continuous(double* xs, size_t n,
        const plfit_continuous_options_t *options, plfit_bool_t xmin_fixed,
        plfit_result_t *result) {
    long int num_trials;
    long int successes = 0;
    double *xs_head;
    size_t num_smaller;
    plfit_continuous_options_t options_no_p_value = *options;
    int retval = PLFIT_SUCCESS;

    if (options->p_value_method == PLFIT_P_VALUE_SKIP) {
        result->p = NAN;
        return PLFIT_SUCCESS;
    }

    if (options->p_value_method == PLFIT_P_VALUE_APPROXIMATE) {
        num_smaller = count_smaller(xs, xs + n, result->xmin);
        result->p = plfit_ks_test_one_sample_p(result->D, n - num_smaller);
        return PLFIT_SUCCESS;
    }

    options_no_p_value.p_value_method = PLFIT_P_VALUE_SKIP;
    num_trials = (long int)(0.25 / options->p_value_precision / options->p_value_precision);
    if (num_trials <= 0) {
        PLFIT_ERROR("invalid p-value precision", PLFIT_EINVAL);
    }

    /* Extract the head of xs that contains elements smaller than xmin */
    xs_head = extract_smaller(xs, xs+n, result->xmin, &num_smaller);
    if (xs_head == 0)
        PLFIT_ERROR("cannot calculate exact p-value", PLFIT_ENOMEM);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        /* Parallel section starts here. If we are compiling using OpenMP, each
         * thread will use its own RNG that is seeded from the master RNG. If
         * we are compiling without OpenMP, there is only one thread and it uses
         * the master RNG. This section must be critical to ensure that only one
         * thread is using the master RNG at the same time. */
#ifdef _OPENMP
        mt_rng_t private_rng;
#endif
        mt_rng_t *p_rng;
        double *ys;
        long int i;
        plfit_result_t result_synthetic;

#ifdef _OPENMP
#pragma omp critical
        {
            p_rng = &private_rng;
            mt_init_from_rng(p_rng, options->rng);
        }
#else
        p_rng = options->rng;
#endif

        /* Allocate memory to sample into */
        ys = calloc(n, sizeof(double));
        if (ys == 0) {
            retval = PLFIT_ENOMEM;
        } else {
            /* The main for loop starts here. */
#ifdef _OPENMP
#pragma omp for reduction(+:successes)
#endif
            for (i = 0; i < num_trials; i++) {
                plfit_i_resample_continuous(xs_head, num_smaller, n, result->alpha,
                        result->xmin, n, p_rng, ys);
                if (xmin_fixed) {
                    plfit_estimate_alpha_continuous(ys, n, result->xmin,
                                &options_no_p_value, &result_synthetic);
                } else {
                    plfit_continuous(ys, n, &options_no_p_value, &result_synthetic);
                }
                if (result_synthetic.D > result->D)
                    successes++;
            }
            free(ys);
        }

        /* End of parallelized part */
    }

    free(xs_head);

    if (retval == PLFIT_SUCCESS) {
        result->p = successes / ((double)num_trials);
    } else {
        PLFIT_ERROR("cannot calculate exact p-value", retval);
    }

    return retval;
}

int plfit_log_likelihood_continuous(double* xs, size_t n, double alpha,
        double xmin, double* L) {
    double logsum, c;
    size_t m;

    if (alpha <= 1) {
        PLFIT_ERROR("alpha must be greater than one", PLFIT_EINVAL);
    }
    XMIN_CHECK_ZERO;

    c = (alpha - 1) / xmin;
    plfit_i_logsum_less_than_continuous(xs, xs+n, xmin, &logsum, &m);
    *L = -alpha * logsum + log(c) * m;

    return PLFIT_SUCCESS;
}

int plfit_estimate_alpha_continuous_sorted(double* xs, size_t n, double xmin,
        const plfit_continuous_options_t* options, plfit_result_t *result) {
    double *begin, *end;

    if (!options)
        options = &plfit_continuous_default_options;

    begin = xs;
    end = xs + n;
    while (begin < end && *begin < xmin)
        begin++;

    PLFIT_CHECK(plfit_i_estimate_alpha_continuous_sorted(begin, end-begin,
                xmin, &result->alpha));
    PLFIT_CHECK(plfit_i_ks_test_continuous(begin, end, result->alpha,
                xmin, &result->D));

    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, end-begin);
    result->xmin = xmin;

    PLFIT_CHECK(plfit_log_likelihood_continuous(begin, end-begin, result->alpha,
                result->xmin, &result->L));
    PLFIT_CHECK(plfit_i_calculate_p_value_continuous(xs, n, options, 1, result));

    return PLFIT_SUCCESS;
}

int plfit_estimate_alpha_continuous(double* xs, size_t n, double xmin,
        const plfit_continuous_options_t* options, plfit_result_t *result) {
    double *xs_copy;

    if (!options)
        options = &plfit_continuous_default_options;

    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &xs_copy));
    PLFIT_CHECK(plfit_estimate_alpha_continuous_sorted(xs_copy, n, xmin,
                options, result));
    free(xs_copy);

    return PLFIT_SUCCESS;
}

typedef struct {
    double *begin;        /**< Pointer to the beginning of the array holding the data */
    double *end;          /**< Pointer to after the end of the array holding the data */
    double **probes;      /**< Pointers to the elements of the array that will be probed */
    size_t num_probes;    /**< Number of probes */
    plfit_result_t last;  /**< Result of the last evaluation */
} plfit_continuous_xmin_opt_data_t;

static double plfit_i_continuous_xmin_opt_evaluate(void* instance, double x) {
    plfit_continuous_xmin_opt_data_t* data = (plfit_continuous_xmin_opt_data_t*)instance;
    double* begin = data->probes[(long int)x];

    data->last.xmin = *begin;

#ifdef PLFIT_DEBUG
    printf("Trying with probes[%ld] = %.4f\n", (long int)x, *begin);
#endif

    plfit_i_estimate_alpha_continuous_sorted(begin, data->end-begin, *begin,
            &data->last.alpha);
    plfit_i_ks_test_continuous(begin, data->end, data->last.alpha, *begin,
            &data->last.D);

    return data->last.D;
}

static int plfit_i_continuous_xmin_opt_progress(void* instance, double x, double fx,
        double min, double fmin, double left, double right, int k) {
#ifdef PLFIT_DEBUG
    printf("Iteration #%d: [%.4f; %.4f), x=%.4f, fx=%.4f, min=%.4f, fmin=%.4f\n",
            k, left, right, x, fx, min, fmin);
#endif

    /* Continue only if `left' and `right' point to different integers */
    return (int)left == (int)right;
}

static int plfit_i_continuous_xmin_opt_linear_scan(
        plfit_continuous_xmin_opt_data_t* opt_data, plfit_result_t* best_result,
        size_t* best_n) {
    /* i must be signed, otherwise OpenMP on Windows will complain as it
     * supports signed types only. ssize_t is a POSIX extension so it won't
     * work */
    ptrdiff_t i = 0; /* initialize to work around incorrect warning issued by clang 9.0 */
    plfit_result_t global_best_result;
    size_t global_best_n;

    /* Prepare some variables */
    global_best_n = 0;
    global_best_result.D = DBL_MAX;
    global_best_result.xmin = 0;
    global_best_result.alpha = 0;

    /* Due to the OpenMP parallelization, we do things as follows. Each
     * OpenMP thread will search for the best D-score on its own and store
     * the result in a private local_best_result variable. The end of the
     * parallel block contains a critical section that threads will enter
     * one by one and compare their private local_best_result with a
     * global_best that is shared among the threads.
     */
#ifdef _OPENMP
#pragma omp parallel shared(global_best_result, global_best_n) private(i) firstprivate(opt_data)
#endif
    {
        /* These variables are private since they are declared within the
         * parallel block */
        plfit_result_t local_best_result;
        plfit_continuous_xmin_opt_data_t local_opt_data = *opt_data;
        size_t local_best_n;

        /* Initialize the local_best_result and local_best_n variables */
        local_best_n = 0;
        local_best_result.D = DBL_MAX;
        local_best_result.xmin = 0;
        local_best_result.alpha = 0;
        local_best_result.p = 0;

        /* The range of the for loop below is divided among the threads.
         * nowait means that there will be no implicit barrier at the end
         * of the loop so threads that get there earlier can enter the
         * critical section without waiting for the others */
#ifdef _OPENMP
#pragma omp for nowait schedule(dynamic,10)
#endif
        for (i = 0; i < local_opt_data.num_probes-1; i++) {
            plfit_i_continuous_xmin_opt_evaluate(&local_opt_data, i);
            if (local_opt_data.last.D < local_best_result.D) {
#ifdef PLFIT_DEBUG
                printf("Found new local best at %g with D=%g\n",
                        local_opt_data.last.xmin, local_opt_data.last.D);
#endif
                local_best_result = local_opt_data.last;
                local_best_n = local_opt_data.end - local_opt_data.probes[i] + 1;
            }
        }

        /* Critical section that finds the global best result from the
         * local ones collected by each thread */
#ifdef _OPENMP
#pragma omp critical
#endif
        if (local_best_result.D < global_best_result.D) {
            global_best_result = local_best_result;
            global_best_n = local_best_n;
#ifdef PLFIT_DEBUG
            printf("Found new global best at %g with D=%g\n", global_best_result.xmin,
                    global_best_result.D);
#endif
        }
    }

    *best_result = global_best_result;
    *best_n = global_best_n;

#ifdef PLFIT_DEBUG
    printf("Returning global best: %g\n", best_result->xmin);
#endif

    return PLFIT_SUCCESS;
}

int plfit_continuous(double* xs, size_t n, const plfit_continuous_options_t* options,
        plfit_result_t* result) {
    gss_parameter_t gss_param;
    plfit_continuous_xmin_opt_data_t opt_data;
    plfit_result_t best_result = {
        /* alpha = */ NAN,
        /* xmin = */ NAN,
        /* L = */ NAN,
        /* D = */ NAN,
        /* p = */ NAN
    };

    int success;
    size_t i, best_n, num_uniques = 0;
    double x, *px, **uniques;

    DATA_POINTS_CHECK;

    /* Sane defaults */
    best_n = n;
    if (!options)
        options = &plfit_continuous_default_options;

    /* Make a copy of xs and sort it */
    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &opt_data.begin));
    opt_data.end = opt_data.begin + n;

    /* Create an array containing pointers to the unique elements of the input. From
     * each block of unique elements, we add the pointer to the first one. */
    uniques = unique_element_pointers(opt_data.begin, opt_data.end, &num_uniques);
    if (uniques == 0)
        PLFIT_ERROR("cannot fit continuous power-law", PLFIT_ENOMEM);

    /* We will now determine the best xmin that yields the lowest D-score. The
     * 'success' variable will denote whether the search procedure we tried was
     * successful. If it is false after having exhausted all options, we fall
     * back to a linear search. */
    success = 0;
    switch (options->xmin_method) {
        case PLFIT_GSS_OR_LINEAR:
            /* Try golden section search first. */
            if (num_uniques > 5) {
                opt_data.probes = uniques;
                opt_data.num_probes = num_uniques;
                gss_parameter_init(&gss_param);
                success = (gss(0, opt_data.num_probes-5, &x, 0,
                        plfit_i_continuous_xmin_opt_evaluate,
                        plfit_i_continuous_xmin_opt_progress, &opt_data, &gss_param) == 0);
                if (success) {
                    px = opt_data.probes[(int)x];
                    best_n = opt_data.end-px+1;
                    best_result = opt_data.last;
                }
            }
            break;

        case PLFIT_STRATIFIED_SAMPLING:
            if (num_uniques >= 50) {
                /* Try stratified sampling to narrow down the interval where the minimum
                 * is likely to reside. We check 10% of the unique items, distributed
                 * evenly, find the one with the lowest D-score, and then check the
                 * area around it more thoroughly. */
                const size_t subdivision_length = 10;
                size_t num_strata = num_uniques / subdivision_length;
                double **strata = calloc(num_strata, sizeof(double*));

                for (i = 0; i < num_strata; i++) {
                    strata[i] = uniques[i * subdivision_length];
                }

                opt_data.probes = strata;
                opt_data.num_probes = num_strata;
                plfit_i_continuous_xmin_opt_linear_scan(&opt_data, &best_result, &best_n);

                opt_data.num_probes = 0;
                for (i = 0; i < num_strata; i++) {
                    if (*strata[i] == best_result.xmin) {
                        /* Okay, scan more thoroughly from strata[i-1] to strata[i+1],
                         * which is from uniques[(i-1)*subdivision_length] to
                         * uniques[(i+1)*subdivision_length */
                        opt_data.probes = uniques + (i > 0 ? (i-1)*subdivision_length : 0);
                        opt_data.num_probes = 0;
                        if (i != 0)
                            opt_data.num_probes += subdivision_length;
                        if (i != num_strata-1)
                            opt_data.num_probes += subdivision_length;
                        break;
                    }
                }

                free(strata);
                if (opt_data.num_probes > 0) {
                    /* Do a strict linear scan in the subrange determined above */
                    plfit_i_continuous_xmin_opt_linear_scan(&opt_data,
                            &best_result, &best_n);
                    success = 1;
                } else {
                    /* This should not happen, but we handle it anyway */
                    success = 0;
                }
            }
            break;

        default:
            /* Just use the linear search */
            break;
    }

    if (!success) {
        /* More advanced search methods failed or were skipped; try linear search */
        opt_data.probes = uniques;
        opt_data.num_probes = num_uniques;
        plfit_i_continuous_xmin_opt_linear_scan(&opt_data, &best_result, &best_n);
        success = 1;
    }

    /* Get rid of the uniques array, we don't need it any more */
    free(uniques);

    /* Sort out the result */
    *result = best_result;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, best_n);

    PLFIT_CHECK(plfit_log_likelihood_continuous(opt_data.begin + n - best_n, best_n,
            result->alpha, result->xmin, &result->L));
    PLFIT_CHECK(plfit_i_calculate_p_value_continuous(opt_data.begin, n, options, 0, result));

    /* Get rid of the copied data as well */
    free(opt_data.begin);

    return PLFIT_SUCCESS;
}

/********** Discrete power law distribution fitting **********/

typedef struct {
    size_t m;
    double logsum;
    double xmin;
} plfit_i_estimate_alpha_discrete_data_t;

static double plfit_i_logsum_discrete(double* begin, double* end, double xmin) {
    double logsum = 0.0;
    for (; begin != end; begin++)
        logsum += log(*begin);
    return logsum;
}

static void plfit_i_logsum_less_than_discrete(double* begin, double* end, double xmin,
        double* logsum, size_t* m) {
    double result = 0.0;
    size_t count = 0;

    for (; begin != end; begin++) {
        if (*begin < xmin)
            continue;

        result += log(*begin);
        count++;
    }

    *logsum = result;
    *m = count;
}

static lbfgsfloatval_t plfit_i_estimate_alpha_discrete_lbfgs_evaluate(
        void* instance, const lbfgsfloatval_t* x,
        lbfgsfloatval_t* g, const int n,
        const lbfgsfloatval_t step) {
    plfit_i_estimate_alpha_discrete_data_t* data;
    lbfgsfloatval_t result;
    double dx = step;
    double huge = 1e10;     /* pseudo-infinity; apparently DBL_MAX does not work */
    double lnhzeta_x=NAN;
    double lnhzeta_deriv_x=NAN;

    data = (plfit_i_estimate_alpha_discrete_data_t*)instance;

#ifdef PLFIT_DEBUG
    printf("- Evaluating at %.4f (step = %.4f, xmin = %.4f)\n", *x, step, data->xmin);
#endif

    if (isnan(*x)) {
        g[0] = huge;
        return huge;
    }

    /* Find the delta X value to estimate the gradient */
    if (dx > 0.001 || dx == 0)
        dx = 0.001;
    else if (dx < -0.001)
        dx = -0.001;

    /* Is x[0] in its valid range? */
    if (x[0] <= 1.0) {
        /* The Hurwitz zeta function is infinite in this case */
        g[0] = (dx > 0) ? -huge : huge;
        return huge;
    }
    if (x[0] + dx <= 1.0) {
        g[0] = huge;
        result = x[0] * data->logsum + data->m * hsl_sf_lnhzeta(x[0], data->xmin);
    } else {
        hsl_sf_lnhzeta_deriv_tuple(x[0], data->xmin, &lnhzeta_x, &lnhzeta_deriv_x);
        g[0] = data->logsum + data->m * lnhzeta_deriv_x;
        result = x[0] * data->logsum + data->m * lnhzeta_x;
    }

#ifdef PLFIT_DEBUG
    printf("  - Gradient: %.4f\n", g[0]);
    printf("  - Result: %.4f\n", result);
#endif

    return result;
}

static int plfit_i_estimate_alpha_discrete_lbfgs_progress(void* instance,
        const lbfgsfloatval_t* x, const lbfgsfloatval_t* g,
        const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
        int n, int k, int ls) {
    return 0;
}

static int plfit_i_estimate_alpha_discrete_linear_scan(double* xs, size_t n,
        double xmin, double* alpha, const plfit_discrete_options_t* options,
        plfit_bool_t sorted) {
    double curr_alpha, best_alpha, L, L_max;
    double logsum;
    size_t m;

    XMIN_CHECK_ONE;
    if (options->alpha.min <= 1.0) {
        PLFIT_ERROR("alpha.min must be greater than 1.0", PLFIT_EINVAL);
    }
    if (options->alpha.max < options->alpha.min) {
        PLFIT_ERROR("alpha.max must be greater than alpha.min", PLFIT_EINVAL);
    }
    if (options->alpha.step <= 0) {
        PLFIT_ERROR("alpha.step must be positive", PLFIT_EINVAL);
    }

    if (sorted) {
        logsum = plfit_i_logsum_discrete(xs, xs+n, xmin);
        m = n;
    } else {
        plfit_i_logsum_less_than_discrete(xs, xs+n, xmin, &logsum, &m);
    }

    best_alpha = options->alpha.min; L_max = -DBL_MAX;
    for (curr_alpha = options->alpha.min; curr_alpha <= options->alpha.max;
            curr_alpha += options->alpha.step) {
        L = -curr_alpha * logsum - m * hsl_sf_lnhzeta(curr_alpha, xmin);
        if (L > L_max) {
            L_max = L;
            best_alpha = curr_alpha;
        }
    }

    *alpha = best_alpha;

    return PLFIT_SUCCESS;
}

static int plfit_i_estimate_alpha_discrete_lbfgs(double* xs, size_t n, double xmin,
        double* alpha, const plfit_discrete_options_t* options, plfit_bool_t sorted) {
    lbfgs_parameter_t param;
    lbfgsfloatval_t* variables;
    plfit_i_estimate_alpha_discrete_data_t data;
    int ret;

    XMIN_CHECK_ONE;

    /* Initialize algorithm parameters */
    lbfgs_parameter_init(&param);
    param.max_iterations = 0;   /* proceed until infinity */

    /* Set up context for optimization */
    data.xmin = xmin;
    if (sorted) {
        data.logsum = plfit_i_logsum_discrete(xs, xs+n, xmin);
        data.m = n;
    } else {
        plfit_i_logsum_less_than_discrete(xs, xs+n, xmin, &data.logsum, &data.m);
    }

    /* Allocate space for the single alpha variable */
    variables = lbfgs_malloc(1);
    variables[0] = 3.0;       /* initial guess */

    /* Optimization */
    ret = lbfgs(1, variables, /* ptr_fx = */ 0,
            plfit_i_estimate_alpha_discrete_lbfgs_evaluate,
            plfit_i_estimate_alpha_discrete_lbfgs_progress,
            &data, &param);

    if (ret < 0 &&
        ret != LBFGSERR_ROUNDING_ERROR &&
        ret != LBFGSERR_MAXIMUMLINESEARCH &&
        ret != LBFGSERR_MINIMUMSTEP &&
        ret != LBFGSERR_CANCELED) {
        char buf[4096];
        snprintf(buf, 4096, "L-BFGS optimization signaled an error (error code = %d)", ret);
        lbfgs_free(variables);
        PLFIT_ERROR(buf, PLFIT_FAILURE);
    }
    *alpha = variables[0];

    /* Deallocate the variable array */
    lbfgs_free(variables);

    return PLFIT_SUCCESS;
}

static int plfit_i_estimate_alpha_discrete_fast(double* xs, size_t n, double xmin,
        double* alpha, const plfit_discrete_options_t* options, plfit_bool_t sorted) {
    plfit_continuous_options_t cont_options;

    if (!options)
        options = &plfit_discrete_default_options;

    plfit_continuous_options_init(&cont_options);
    cont_options.finite_size_correction = options->finite_size_correction;

    XMIN_CHECK_ONE;

    if (sorted) {
        return plfit_i_estimate_alpha_continuous_sorted(xs, n, xmin-0.5, alpha);
    } else {
        return plfit_i_estimate_alpha_continuous(xs, n, xmin-0.5, alpha);
    }
}

static int plfit_i_estimate_alpha_discrete(double* xs, size_t n, double xmin,
        double* alpha, const plfit_discrete_options_t* options,
        plfit_bool_t sorted) {
    switch (options->alpha_method) {
        case PLFIT_LBFGS:
            PLFIT_CHECK(plfit_i_estimate_alpha_discrete_lbfgs(xs, n, xmin, alpha,
                        options, sorted));
            break;

        case PLFIT_LINEAR_SCAN:
            PLFIT_CHECK(plfit_i_estimate_alpha_discrete_linear_scan(xs, n, xmin,
                        alpha, options, sorted));
            break;

        case PLFIT_PRETEND_CONTINUOUS:
            PLFIT_CHECK(plfit_i_estimate_alpha_discrete_fast(xs, n, xmin,
                        alpha, options, sorted));
            break;

        default:
            PLFIT_ERROR("unknown optimization method specified", PLFIT_EINVAL);
    }

    return PLFIT_SUCCESS;
}

static int plfit_i_ks_test_discrete(double* xs, double* xs_end, const double alpha,
        const double xmin, double* D) {
    /* Assumption: xs is sorted and cut off at xmin so the first element is
     * always larger than or equal to xmin. */
    double result = 0, n, lnhzeta, x;
    int m = 0;

    n = xs_end - xs;
    lnhzeta = hsl_sf_lnhzeta(alpha, xmin);

    while (xs < xs_end) {
        double d;

        x = *xs;

        /* Re the next line: this used to be the following:
         *
         * fabs( 1 - hzeta(alpha, x) / hzeta(alpha, xmin) - m / n)
         *
         * However, using the Hurwitz zeta directly sometimes yields
         * underflows (see Github pull request #17 and related issues).
         * hzeta(alpha, x) / hzeta(alpha, xmin) can be replaced with
         * exp(lnhzeta(alpha, x) - lnhzeta(alpha, xmin)), but then
         * we have 1 - exp(something), which is better to calculate
         * with a dedicated expm1() function.
         */
        d = fabs( expm1( hsl_sf_lnhzeta(alpha, x) - lnhzeta ) + m / n);

        if (d > result)
            result = d;

        do {
            xs++; m++;
        } while (xs < xs_end && *xs == x);
    }

    *D = result;

    return PLFIT_SUCCESS;
}

static int plfit_i_calculate_p_value_discrete(double* xs, size_t n,
        const plfit_discrete_options_t* options, plfit_bool_t xmin_fixed,
        plfit_result_t *result) {
    long int num_trials;
    long int successes = 0;
    double *xs_head;
    size_t num_smaller;
    plfit_discrete_options_t options_no_p_value = *options;
    int retval = PLFIT_SUCCESS;

    if (options->p_value_method == PLFIT_P_VALUE_SKIP) {
        /* skipping p-value calculation */
        result->p = NAN;
        return PLFIT_SUCCESS;
    }

    if (options->p_value_method == PLFIT_P_VALUE_APPROXIMATE) {
        /* p-value approximation; most likely an upper bound */
        num_smaller = count_smaller(xs, xs + n, result->xmin);
        result->p = plfit_ks_test_one_sample_p(result->D, n - num_smaller);
        return PLFIT_SUCCESS;
    }

    options_no_p_value.p_value_method = PLFIT_P_VALUE_SKIP;
    num_trials = (long int)(0.25 / options->p_value_precision / options->p_value_precision);
    if (num_trials <= 0) {
        PLFIT_ERROR("invalid p-value precision", PLFIT_EINVAL);
    }

    /* Extract the head of xs that contains elements smaller than xmin */
    xs_head = extract_smaller(xs, xs+n, result->xmin, &num_smaller);
    if (xs_head == 0)
        PLFIT_ERROR("cannot calculate exact p-value", PLFIT_ENOMEM);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        /* Parallel section starts here. If we are compiling using OpenMP, each
         * thread will use its own RNG that is seeded from the master RNG. If
         * we are compiling without OpenMP, there is only one thread and it uses
         * the master RNG. This section must be critical to ensure that only one
         * thread is using the master RNG at the same time. */
#ifdef _OPENMP
        mt_rng_t private_rng;
#endif
        mt_rng_t *p_rng;
        double *ys;
        long int i;
        plfit_result_t result_synthetic;

#ifdef _OPENMP
#pragma omp critical
        {
            p_rng = &private_rng;
            mt_init_from_rng(p_rng, options->rng);
        }
#else
        p_rng = options->rng;
#endif

        /* Allocate memory to sample into */
        ys = calloc(n, sizeof(double));
        if (ys == 0) {
            retval = PLFIT_ENOMEM;
        } else {
            /* The main for loop starts here. */
#ifdef _OPENMP
#pragma omp for reduction(+:successes)
#endif
            for (i = 0; i < num_trials; i++) {
                plfit_i_resample_discrete(xs_head, num_smaller, n, result->alpha,
                        result->xmin, n, p_rng, ys);
                if (xmin_fixed) {
                    plfit_estimate_alpha_discrete(ys, n, result->xmin,
                                &options_no_p_value, &result_synthetic);
                } else {
                    plfit_discrete(ys, n, &options_no_p_value, &result_synthetic);
                }
                if (result_synthetic.D > result->D)
                    successes++;
            }

            free(ys);
        }

        /* End of parallelized part */
    }

    free(xs_head);

    if (retval == PLFIT_SUCCESS) {
        result->p = successes / ((double)num_trials);
    } else {
        PLFIT_ERROR("cannot calculate exact p-value", retval);
    }

    return retval;
}

int plfit_log_likelihood_discrete(double* xs, size_t n, double alpha, double xmin, double* L) {
    double result;
    size_t m;

    if (alpha <= 1) {
        PLFIT_ERROR("alpha must be greater than one", PLFIT_EINVAL);
    }
    XMIN_CHECK_ONE;

    plfit_i_logsum_less_than_discrete(xs, xs+n, xmin, &result, &m);
    result = - alpha * result - m * hsl_sf_lnhzeta(alpha, xmin);

    *L = result;

    return PLFIT_SUCCESS;
}

int plfit_estimate_alpha_discrete(double* xs, size_t n, double xmin,
        const plfit_discrete_options_t* options, plfit_result_t *result) {
    double *xs_copy, *begin, *end;

    if (!options)
        options = &plfit_discrete_default_options;

    /* Check the validity of the input parameters */
    DATA_POINTS_CHECK;
    if (options->alpha_method == PLFIT_LINEAR_SCAN) {
        if (options->alpha.min <= 1.0) {
            PLFIT_ERROR("alpha.min must be greater than 1.0", PLFIT_EINVAL);
        }
        if (options->alpha.max < options->alpha.min) {
            PLFIT_ERROR("alpha.max must be greater than alpha.min", PLFIT_EINVAL);
        }
        if (options->alpha.step <= 0) {
            PLFIT_ERROR("alpha.step must be positive", PLFIT_EINVAL);
        }
    }

    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &xs_copy));

    begin = xs_copy; end = xs_copy + n;
    while (begin < end && *begin < xmin)
        begin++;

    PLFIT_CHECK(plfit_i_estimate_alpha_discrete(begin, end-begin, xmin, &result->alpha,
                options, /* sorted = */ 1));
    PLFIT_CHECK(plfit_i_ks_test_discrete(begin, end, result->alpha, xmin, &result->D));

    result->xmin = xmin;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, end-begin);

    PLFIT_CHECK(plfit_log_likelihood_discrete(begin, end-begin, result->alpha,
                result->xmin, &result->L));
    PLFIT_CHECK(plfit_i_calculate_p_value_discrete(xs, n, options, 1, result));

    free(xs_copy);

    return PLFIT_SUCCESS;
}

int plfit_discrete(double* xs, size_t n, const plfit_discrete_options_t* options,
        plfit_result_t* result) {
    double curr_D, curr_alpha;
    plfit_result_t best_result;
    double *xs_copy, *px, *end, *end_xmin, prev_x;
    size_t best_n;
    int m;

    if (!options)
        options = &plfit_discrete_default_options;

    /* Check the validity of the input parameters */
    DATA_POINTS_CHECK;
    if (options->alpha_method == PLFIT_LINEAR_SCAN) {
        if (options->alpha.min <= 1.0) {
            PLFIT_ERROR("alpha.min must be greater than 1.0", PLFIT_EINVAL);
        }
        if (options->alpha.max < options->alpha.min) {
            PLFIT_ERROR("alpha.max must be greater than alpha.min", PLFIT_EINVAL);
        }
        if (options->alpha.step <= 0) {
            PLFIT_ERROR("alpha.step must be positive", PLFIT_EINVAL);
        }
    }

    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &xs_copy));

    best_result.D = DBL_MAX;
    best_result.xmin = 1;
    best_result.alpha = 1;
    best_n = 0;

    /* Make sure there are at least three distinct values if possible */
    px = xs_copy; end = px + n; end_xmin = end - 1; m = 0;
    prev_x = *end_xmin;
    while (*end_xmin == prev_x && end_xmin > px)
        end_xmin--;
    prev_x = *end_xmin;
    while (*end_xmin == prev_x && end_xmin > px)
        end_xmin--;

    prev_x = 0;
    while (px < end_xmin) {
        while (px < end_xmin && *px == prev_x) {
            px++; m++;
        }

        plfit_i_estimate_alpha_discrete(px, n-m, *px, &curr_alpha, options,
                /* sorted = */ 1);
        plfit_i_ks_test_discrete(px, end, curr_alpha, *px, &curr_D);

        if (curr_D < best_result.D) {
            best_result.alpha = curr_alpha;
            best_result.xmin = *px;
            best_result.D = curr_D;
            best_n = n-m;
        }

        prev_x = *px;
        px++; m++;
    }

    *result = best_result;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, best_n);

    PLFIT_CHECK(plfit_log_likelihood_discrete(xs_copy+(n-best_n), best_n,
                result->alpha, result->xmin, &result->L));
    PLFIT_CHECK(plfit_i_calculate_p_value_discrete(xs_copy, n, options, 0, result));

    free(xs_copy);

    return PLFIT_SUCCESS;
}

/***** resampling routines to generate synthetic replicates ****/

static int plfit_i_resample_continuous(double* xs_head, size_t num_smaller,
        size_t n, double alpha, double xmin, size_t num_samples, mt_rng_t* rng,
        double* result)
{
    size_t num_orig_samples, i;

    /* Calculate how many samples have to be drawn from xs_head */
    num_orig_samples = (size_t) plfit_rbinom(num_samples, num_smaller / (double)n, rng);

    /* Draw the samples from xs_head */
    for (i = 0; i < num_orig_samples; i++, result++) {
        *result = xs_head[(size_t)plfit_runif(0, num_smaller, rng)];
    }

    /* Draw the remaining samples from the fitted distribution */
    PLFIT_CHECK(plfit_rpareto_array(xmin, alpha-1, num_samples-num_orig_samples, rng,
            result));

    return PLFIT_SUCCESS;
}

int plfit_resample_continuous(double* xs, size_t n, double alpha, double xmin,
        size_t num_samples, mt_rng_t* rng, double* result) {
    double *xs_head;
    size_t num_smaller = 0;
    int retval;

    /* Extract the head of xs that contains elements smaller than xmin */
    xs_head = extract_smaller(xs, xs+n, xmin, &num_smaller);
    if (xs_head == 0)
        PLFIT_ERROR("cannot resample continuous dataset", PLFIT_ENOMEM);

    retval = plfit_i_resample_continuous(xs_head, num_smaller, n, alpha, xmin,
                num_samples, rng, result);

    /* Free xs_head; we don't need it any more */
    free(xs_head);

    return retval;
}

static int plfit_i_resample_discrete(double* xs_head, size_t num_smaller, size_t n,
        double alpha, double xmin, size_t num_samples, mt_rng_t* rng,
        double* result)
{
    size_t num_orig_samples, i;

    /* Calculate how many samples have to be drawn from xs_head */
    num_orig_samples = (size_t) plfit_rbinom(num_samples, num_smaller / (double)n, rng);

    /* Draw the samples from xs_head */
    for (i = 0; i < num_orig_samples; i++, result++) {
        *result = xs_head[(size_t)plfit_runif(0, num_smaller, rng)];
    }

    /* Draw the remaining samples from the fitted distribution */
    PLFIT_CHECK(plfit_rzeta_array((long int)xmin, alpha,
                num_samples-num_orig_samples, rng, result));

    return PLFIT_SUCCESS;
}

int plfit_resample_discrete(double* xs, size_t n, double alpha, double xmin,
        size_t num_samples, mt_rng_t* rng, double* result) {
    double *xs_head;
    size_t num_smaller = 0;
    int retval;

    /* Extract the head of xs that contains elements smaller than xmin */
    xs_head = extract_smaller(xs, xs+n, xmin, &num_smaller);
    if (xs_head == 0)
        PLFIT_ERROR("cannot resample discrete dataset", PLFIT_ENOMEM);

    retval = plfit_i_resample_discrete(xs_head, num_smaller, n, alpha, xmin,
                num_samples, rng, result);

    /* Free xs_head; we don't need it any more */
    free(xs_head);

    return retval;
}

/******** calculating the p-value of a fitted model only *******/

int plfit_calculate_p_value_continuous(double* xs, size_t n,
        const plfit_continuous_options_t* options, plfit_bool_t xmin_fixed,
        plfit_result_t *result) {
    double* xs_copy;

    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &xs_copy));
    PLFIT_CHECK(plfit_i_calculate_p_value_continuous(xs_copy, n, options,
                xmin_fixed, result));
    free(xs_copy);

    return PLFIT_SUCCESS;
}

int plfit_calculate_p_value_discrete(double* xs, size_t n,
        const plfit_discrete_options_t* options, plfit_bool_t xmin_fixed,
        plfit_result_t *result) {
    double* xs_copy;

    PLFIT_CHECK(plfit_i_copy_and_sort(xs, n, &xs_copy));
    PLFIT_CHECK(plfit_i_calculate_p_value_discrete(xs_copy, n, options,
                xmin_fixed, result));
    free(xs_copy);

    return PLFIT_SUCCESS;
}
