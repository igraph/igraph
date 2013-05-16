/* plfit.c
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "gss.h"
#include "lbfgs.h"
#include "platform.h"
#include "plfit.h"
#include "kolmogorov.h"
#include "zeta.h"

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

static int double_comparator(const void *a, const void *b) {
    const double *da = (const double*)a;
    const double *db = (const double*)b;
    return (*da > *db) - (*da < *db);
}

/**
 * Given a sorted array of doubles, return another array that contains pointers
 * into the array for the start of each block of identical elements.
 *
 * \param  begin          pointer to the beginning of the array
 * \param  end            pointer to the first element after the end of the array
 * \param  result_length  if not \c NULL, the number of unique elements in the
 *                        given array is returned here
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

void plfit_i_logsum_less_than_continuous(double* begin, double* end,
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

double plfit_i_logsum_continuous(double* begin, double* end, double xmin) {
    double logsum = 0.0;
    for (; begin != end; begin++)
        logsum += log(*begin / xmin);
    return logsum;
}

int plfit_i_estimate_alpha_continuous(double* xs, size_t n,
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

int plfit_i_estimate_alpha_continuous_sorted(double* xs, size_t n,
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

int plfit_estimate_alpha_continuous(double* xs, size_t n, double xmin,
        const plfit_continuous_options_t* options, plfit_result_t *result) {
    double *xs_copy;

	if (!options)
		options = &plfit_continuous_default_options;

    /* Make a copy of xs and sort it */
    xs_copy = (double*)malloc(sizeof(double) * n);
    memcpy(xs_copy, xs, sizeof(double) * n);
    qsort(xs_copy, n, sizeof(double), double_comparator);

    PLFIT_CHECK(plfit_estimate_alpha_continuous_sorted(xs_copy, n, xmin,
				options, result));

    free(xs_copy);

    return PLFIT_SUCCESS;
}

int plfit_estimate_alpha_continuous_sorted(double* xs, size_t n, double xmin,
        const plfit_continuous_options_t* options, plfit_result_t *result) {
    double* end;

	if (!options)
		options = &plfit_continuous_default_options;

	end = xs + n;
    while (xs < end && *xs < xmin)
        xs++;
    n = (size_t) (end - xs);

    PLFIT_CHECK(plfit_i_estimate_alpha_continuous_sorted(xs, n,
				xmin, &result->alpha));
    PLFIT_CHECK(plfit_i_ks_test_continuous(xs, end, result->alpha,
				xmin, &result->D));

    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, n);
    result->xmin = xmin;
    result->p = plfit_ks_test_one_sample_p(result->D, n);
    plfit_log_likelihood_continuous(xs, n, result->alpha, result->xmin, &result->L);

    return PLFIT_SUCCESS;
}

typedef struct {
	double *begin;        /**< Pointer to the beginning of the array holding the data */
	double *end;          /**< Pointer to after the end of the array holding the data */
	double **uniques;     /**< Pointers to unique elements of the input array */
	plfit_result_t last;  /**< Result of the last evaluation */
} plfit_continuous_xmin_opt_data_t;

double plfit_i_continuous_xmin_opt_evaluate(void* instance, double x) {
	plfit_continuous_xmin_opt_data_t* data = (plfit_continuous_xmin_opt_data_t*)instance;
	double* begin = data->uniques[(int)x];

	data->last.xmin = *begin;

#ifdef PLFIT_DEBUG
	printf("Trying with xmin = %.4f\n", *begin);
#endif

	plfit_i_estimate_alpha_continuous_sorted(begin, (size_t) (data->end-begin), *begin,
			&data->last.alpha);
	plfit_i_ks_test_continuous(begin, data->end, data->last.alpha, *begin,
			&data->last.D);

	return data->last.D;
}

int plfit_i_continuous_xmin_opt_progress(void* instance, double x, double fx,
		double min, double fmin, double left, double right, int k) {
#ifdef PLFIT_DEBUG
    printf("Iteration #%d: [%.4f; %.4f), x=%.4f, fx=%.4f, min=%.4f, fmin=%.4f\n",
            k, left, right, x, fx, min, fmin);
#endif

	/* Continue only if `left' and `right' point to different integers */
	return (int)left == (int)right;
}

int plfit_continuous(double* xs, size_t n, const plfit_continuous_options_t* options,
        plfit_result_t* result) {
	gss_parameter_t gss_param;
	plfit_continuous_xmin_opt_data_t opt_data;
	plfit_result_t best_result;
	int success;
	size_t i, best_n, num_uniques;
    double x, *px;

    DATA_POINTS_CHECK;

	if (!options)
		options = &plfit_continuous_default_options;

    /* Make a copy of xs and sort it */
    opt_data.begin = (double*)malloc(sizeof(double) * n);
    memcpy(opt_data.begin, xs, sizeof(double) * n);
    qsort(opt_data.begin, n, sizeof(double), double_comparator);
    opt_data.end = opt_data.begin + n;

    /* Create an array containing pointers to the unique elements of the input. From
     * each block of unique elements, we add the pointer to the first one. */
    opt_data.uniques = unique_element_pointers(opt_data.begin, opt_data.end,
			&num_uniques);
    if (opt_data.uniques == 0)
        return PLFIT_ENOMEM;

    /* We will now determine the best xmin that yields the lowest D-score.
	 * First we try a golden section search if needed. If that fails, we try
	 * a linear search.
     */
	if (options->xmin_method == PLFIT_GSS_OR_LINEAR && num_uniques > 5) {
		gss_parameter_init(&gss_param);
		success = (gss(0, num_uniques-5, &x, 0,
				plfit_i_continuous_xmin_opt_evaluate,
				plfit_i_continuous_xmin_opt_progress, &opt_data, &gss_param) == 0);
		best_result = opt_data.last;
		/* plfit_i_continuous_xmin_opt_evaluate will set opt_data.last to
		 * indicate the location of the optimum and the value of D */
	} else {
		success = 0;
	}

	if (success) {
		/* calculate best_n because we'll need it later. Luckily x indicates
		 * the index in opt_data.uniques that we have to look up in order to
		 * find the first element in the array that is included */
		px = opt_data.uniques[(int)x];
		best_n = (size_t) (opt_data.end-px+1);
	} else {
		/* GSS failed or skipped; try linear search */

		/* Prepare some variables */
		best_n = 0;
		best_result.D = DBL_MAX;
		best_result.xmin = 0;
		best_result.alpha = 0;
		
		for (i = 0; i < num_uniques-1; i++) {
			plfit_i_continuous_xmin_opt_evaluate(&opt_data, i);
			if (opt_data.last.D < best_result.D) {
				best_result = opt_data.last;
				best_n = (size_t) (opt_data.end - 
						   opt_data.uniques[i] + 1);
			}
		}
	}

    /* Get rid of the uniques array, we don't need it any more */
    free(opt_data.uniques);

    /* Sort out the result */
    *result = best_result;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, best_n);
    result->p = plfit_ks_test_one_sample_p(result->D, best_n);
    plfit_log_likelihood_continuous(opt_data.begin + n - best_n, best_n,
			result->alpha, result->xmin, &result->L);

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

double plfit_i_logsum_discrete(double* begin, double* end, double xmin) {
    double logsum = 0.0;
    for (; begin != end; begin++)
        logsum += log(*begin);
    return logsum;
}

void plfit_i_logsum_less_than_discrete(double* begin, double* end, double xmin,
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

lbfgsfloatval_t plfit_i_estimate_alpha_discrete_lbfgs_evaluate(
        void* instance, const lbfgsfloatval_t* x,
        lbfgsfloatval_t* g, const int n,
        const lbfgsfloatval_t step) {
    plfit_i_estimate_alpha_discrete_data_t* data;
    lbfgsfloatval_t result;
    double dx = step;
    double huge = 1e10;     /* pseudo-infinity; apparently DBL_MAX does not work */

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
	if (x[0] + dx <= 1.0)
		g[0] = huge;
	else
		g[0] = data->logsum + data->m *
			(log(gsl_sf_hzeta(x[0] + dx, data->xmin)) - log(gsl_sf_hzeta(x[0], data->xmin))) / dx;

    result = x[0] * data->logsum + data->m * log(gsl_sf_hzeta(x[0], data->xmin));

#ifdef PLFIT_DEBUG
    printf("  - Gradient: %.4f\n", g[0]);
    printf("  - Result: %.4f\n", result);
#endif

    return result;
}

int plfit_i_estimate_alpha_discrete_lbfgs_progress(void* instance,
        const lbfgsfloatval_t* x, const lbfgsfloatval_t* g,
        const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
        int n, int k, int ls) {
    return 0;
}

int plfit_i_estimate_alpha_discrete_linear_scan(double* xs, size_t n, double xmin,
        double* alpha, const plfit_discrete_options_t* options,
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
        L = -curr_alpha * logsum - m * log(gsl_sf_hzeta(curr_alpha, xmin));
        if (L > L_max) {
            L_max = L;
            best_alpha = curr_alpha;
        }
    }

    *alpha = best_alpha;

    return PLFIT_SUCCESS;
}

int plfit_i_estimate_alpha_discrete_lbfgs(double* xs, size_t n, double xmin,
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

int plfit_i_estimate_alpha_discrete_fast(double* xs, size_t n, double xmin,
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

int plfit_i_estimate_alpha_discrete(double* xs, size_t n, double xmin,
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
    double result = 0, n, hzeta, x;
    int m = 0;

    n = xs_end - xs;
    hzeta = gsl_sf_hzeta(alpha, xmin);

    while (xs < xs_end) {
        double d;

        x = *xs;
        d = fabs(1-(gsl_sf_hzeta(alpha, x) / hzeta) - m / n);

        if (d > result)
            result = d;

        do {
            xs++; m++;
        } while (xs < xs_end && *xs == x);
    }

    *D = result;

    return PLFIT_SUCCESS;
}

int plfit_log_likelihood_discrete(double* xs, size_t n, double alpha, double xmin, double* L) {
    double result;
    size_t m;

    if (alpha <= 1) {
        PLFIT_ERROR("alpha must be greater than one", PLFIT_EINVAL);
    }
    XMIN_CHECK_ONE;

    plfit_i_logsum_less_than_discrete(xs, xs+n, xmin, &result, &m);
    result = - alpha * result - m * log(gsl_sf_hzeta(alpha, xmin));

    *L = result;

    return PLFIT_SUCCESS;
}

int plfit_estimate_alpha_discrete(double* xs, size_t n, double xmin,
        const plfit_discrete_options_t* options, plfit_result_t *result) {
    double *xs_copy, *end;

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

    /* Make a copy of xs and sort it */
    xs_copy = (double*)malloc(sizeof(double) * n);
    memcpy(xs_copy, xs, sizeof(double) * n);
    qsort(xs_copy, n, sizeof(double), double_comparator);

    xs = xs_copy; end = xs_copy + n;
    while (xs < end && *xs < xmin)
        xs++;
    n = (size_t) (end - xs);

    PLFIT_CHECK(plfit_i_estimate_alpha_discrete(xs, n, xmin, &result->alpha,
				options, /* sorted = */ 1));
    PLFIT_CHECK(plfit_i_ks_test_discrete(xs, end, result->alpha, xmin, &result->D));

    result->xmin = xmin;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, n);
    result->p = plfit_ks_test_one_sample_p(result->D, n);
    plfit_log_likelihood_discrete(xs, n, result->alpha, result->xmin, &result->L);

    free(xs_copy);

    return PLFIT_SUCCESS;
}

int plfit_discrete(double* xs, size_t n, const plfit_discrete_options_t* options,
        plfit_result_t* result) {
    double curr_D, curr_alpha;
    plfit_result_t best_result;
    double *xs_copy, *px, *end, *end_xmin, prev_x;
	size_t best_n;
    size_t m;

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

    /* Make a copy of xs and sort it */
    xs_copy = (double*)malloc(sizeof(double) * n);
    memcpy(xs_copy, xs, sizeof(double) * n);
    qsort(xs_copy, n, sizeof(double), double_comparator);

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

	plfit_i_estimate_alpha_discrete(px, n - m, *px,
					&curr_alpha, options, /* sorted = */ 1);
        plfit_i_ks_test_discrete(px, end, curr_alpha, *px, &curr_D);

        if (curr_D < best_result.D) {
            best_result.alpha = curr_alpha;
            best_result.xmin = *px;
            best_result.D = curr_D;
	    best_n = n - m;
        }

        prev_x = *px;
        px++; m++;
    }

    *result = best_result;
    if (options->finite_size_correction)
        plfit_i_perform_finite_size_correction(result, best_n);
    result->p = plfit_ks_test_one_sample_p(result->D, best_n);
    plfit_log_likelihood_discrete(xs_copy+(n-best_n), best_n,
			result->alpha, result->xmin, &result->L);

    free(xs_copy);

    return PLFIT_SUCCESS;
}

