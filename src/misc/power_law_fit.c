/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_nongraph.h"

#include "igraph_random.h"
#include "igraph_types.h"

#include "plfit/plfit_error.h"
#include "plfit/plfit.h"

#include <math.h>

static const char* igraph_i_plfit_error_message = 0;

static void igraph_i_plfit_error_handler_store(const char *reason, const char *file,
        int line, int plfit_errno) {

    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(plfit_errno);

    igraph_i_plfit_error_message = reason;
}

static void igraph_i_plfit_prepare_continuous_options(
    plfit_continuous_options_t* options, igraph_bool_t finite_size_correction
) {
    plfit_continuous_options_init(options);
    options->p_value_method = PLFIT_P_VALUE_SKIP;
    options->xmin_method = PLFIT_STRATIFIED_SAMPLING;
    options->finite_size_correction = (plfit_bool_t) finite_size_correction;
}

static void igraph_i_plfit_prepare_discrete_options(
    plfit_discrete_options_t* options, igraph_bool_t finite_size_correction
) {
    plfit_discrete_options_init(options);
    options->p_value_method = PLFIT_P_VALUE_SKIP;
    options->finite_size_correction = (plfit_bool_t) finite_size_correction;
}

/* Decides whether to use finite size correction for the given input data */
static igraph_bool_t igraph_i_plfit_should_use_finite_size_correction(const igraph_vector_t* data) {
    return igraph_vector_size(data) < 50;
}

static igraph_error_t igraph_i_handle_plfit_error(int code) {
    switch (code) {
    case PLFIT_SUCCESS:
        return IGRAPH_SUCCESS;

    case PLFIT_FAILURE:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_FAILURE);
        break;

    case PLFIT_EINVAL:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EINVAL);
        break;

    case PLFIT_UNDRFLOW:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EUNDERFLOW); /* LCOV_EXCL_LINE */
        break;

    case PLFIT_OVERFLOW:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EOVERFLOW); /* LCOV_EXCL_LINE */
        break;

    case PLFIT_ENOMEM:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        break;

    case PLFIT_EMAXITER:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_DIVERGED); /* LCOV_EXCL_LINE */
        break;

    default:
        IGRAPH_ERRORF("Unknown error code returned from plfit (%d)", IGRAPH_FAILURE, code);
        break;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup nongraph
 * \function igraph_power_law_fit
 * \brief Fits a power-law distribution to a vector of numbers.
 *
 * This function fits a power-law distribution to a vector containing samples
 * from a distribution (that is assumed to follow a power-law of course). In
 * a power-law distribution, it is generally assumed that P(X=x) is
 * proportional to x<superscript>-alpha</superscript>, where x is a positive number and alpha
 * is greater than 1. In many real-world cases, the power-law behaviour kicks
 * in only above a threshold value \em xmin. The goal of this functions is to
 * determine \em alpha if \em xmin is given, or to determine \em xmin and the
 * corresponding value of \em alpha.
 *
 * </para><para>
 * The function uses the maximum likelihood principle to determine \em alpha
 * for a given \em xmin; in other words, the function will return the \em alpha
 * value for which the probability of drawing the given sample is the highest.
 * When \em xmin is not given in advance, the algorithm will attempt to find
 * the optimal \em xmin value for which the p-value of a Kolmogorov-Smirnov
 * test between the fitted distribution and the original sample is the largest.
 * The function uses the method of Clauset, Shalizi and Newman to calculate the
 * parameters of the fitted distribution. See the following reference for
 * details:
 *
 * </para><para>
 * Aaron Clauset, Cosma R. Shalizi and Mark E.J. Newman: Power-law
 * distributions in empirical data. SIAM Review 51(4):661-703, 2009.
 * https://doi.org/10.1137/070710111
 *
 * \param data vector containing the samples for which a power-law distribution
 *             is to be fitted. Note that you have to provide the \em samples,
 *             not the probability density function or the cumulative
 *             distribution function. For example, if you wish to fit
 *             a power-law to the degrees of a graph, you can use the output of
 *             \ref igraph_degree directly as an input argument to
 *             \ref igraph_power_law_fit
 * \param result the result of the fitting algorithm. See \ref igraph_plfit_result_t
 *             for more details. Note that the p-value of the fit is \em not
 *             calculated by default as it is time-consuming; you need to call
 *             \ref igraph_plfit_result_calculate_p_value() to calculate the
 *             p-value itself
 * \param xmin the minimum value in the sample vector where the power-law
 *             behaviour is expected to kick in. Samples smaller than \c xmin
 *             will be ignored by the algorithm. Pass zero here if you want to
 *             include all the samples. If \c xmin is negative, the algorithm
 *             will attempt to determine its best value automatically.
 * \param force_continuous assume that the samples in the \c data argument come
 *             from a continuous distribution even if the sample vector
 *             contains integer values only (by chance). If this argument is
 *             false, igraph will assume a continuous distribution if at least
 *             one sample is non-integer and assume a discrete distribution
 *             otherwise.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory
 *         \c IGRAPH_EINVAL: one of the arguments is invalid
 *         \c IGRAPH_EOVERFLOW: overflow during the fitting process
 *         \c IGRAPH_EUNDERFLOW: underflow during the fitting process
 *         \c IGRAPH_FAILURE: the underlying algorithm signaled a failure
 *         without returning a more specific error code
 *
 * Time complexity: in the continuous case, O(n log(n)) if \c xmin is given.
 * In the discrete case, the time complexity is dominated by the complexity of
 * the underlying L-BFGS algorithm that is used to optimize alpha. If \c xmin
 * is not given, the time complexity is multiplied by the number of unique
 * samples in the input vector (although it should be faster in practice).
 *
 * \example examples/simple/igraph_power_law_fit.c
 */
igraph_error_t igraph_power_law_fit(
    const igraph_vector_t* data, igraph_plfit_result_t* result,
    igraph_real_t xmin, igraph_bool_t force_continuous
) {
    plfit_error_handler_t* plfit_stored_error_handler;
    plfit_result_t plfit_result;
    plfit_continuous_options_t cont_options;
    plfit_discrete_options_t disc_options;
    igraph_bool_t discrete = force_continuous ? false : true;
    igraph_bool_t finite_size_correction;

    int retval;
    size_t i, n;

    finite_size_correction = igraph_i_plfit_should_use_finite_size_correction(data);
    n = (size_t) igraph_vector_size(data);

    if (discrete) {
        /* Does the vector contain discrete values only? */
        for (i = 0; i < n; i++) {
            if (trunc(VECTOR(*data)[i]) != VECTOR(*data)[i]) {
                discrete = false;
                break;
            }
        }
    }

    RNG_BEGIN();

    plfit_stored_error_handler = plfit_set_error_handler(igraph_i_plfit_error_handler_store);
    if (discrete) {
        igraph_i_plfit_prepare_discrete_options(&disc_options, finite_size_correction);
        if (xmin >= 0) {
            retval = plfit_estimate_alpha_discrete(VECTOR(*data), n, xmin,
                                                   &disc_options, &plfit_result);
        } else {
            retval = plfit_discrete(VECTOR(*data), n, &disc_options, &plfit_result);
        }
    } else {
        igraph_i_plfit_prepare_continuous_options(&cont_options, finite_size_correction);
        if (xmin >= 0) {
            retval = plfit_estimate_alpha_continuous(VECTOR(*data), n, xmin,
                     &cont_options, &plfit_result);
        } else {
            retval = plfit_continuous(VECTOR(*data), n, &cont_options, &plfit_result);
        }
    }
    plfit_set_error_handler(plfit_stored_error_handler);

    RNG_END();

    IGRAPH_CHECK(igraph_i_handle_plfit_error(retval));

    if (result) {
        result->data = data;
        result->continuous = !discrete;
        result->alpha = plfit_result.alpha;
        result->xmin = plfit_result.xmin;
        result->L = plfit_result.L;
        result->D = plfit_result.D;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup nongraph
 * \function igraph_plfit_result_calculate_p_value
 * \brief Calculates the p-value of a fitted power-law model.
 *
 * </para><para>
 * The p-value is calculated by resampling the input data many times in a way
 * that the part below the fitted \c x_min threshold is resampled from the
 * input data itself, while the part above the fitted \c x_min threshold is
 * drawn from the fitted power-law function. A Kolmogorov-Smirnov test is then
 * performed for each resampled dataset and its test statistic is compared with the
 * observed test statistic from the original dataset. The fraction of resampled
 * datasets that have a \em higher test statistic is the returned p-value.
 *
 * </para><para>
 * Note that the precision of the returned p-value depends on the number of
 * resampling attempts. The number of resampling trials is determined by
 * 0.25 divided by the square of the required precision. For instance, a required
 * precision of 0.01 means that 2500 samples will be drawn.
 *
 * </para><para>
 * If igraph is compiled with OpenMP support, this function will use parallel
 * OpenMP threads for the resampling. Each OpenMP thread gets its own instance
 * of a random number generator. However, since the scheduling of OpenMP threads
 * is outside our control, we cannot guarantee how many resampling instances the
 * threads are asked to execute, thus it may happen that the random number
 * generators are used differently between runs. If you want to obtain
 * reproducible results, seed igraph's master RNG appropriately, and force the
 * number of OpenMP threads to 1 early in your program, either by calling
 * <code>omp_set_num_threads(1)</code> or by setting the value of the \c OMP_NUM_THREADS
 * environment variable to 1.
 *
 * \param model The fitted power-law model from the \ref igraph_power_law_fit()
 *        function
 * \param result The calculated p-value is returned here
 * \param precision The desired precision of the p-value. Higher values correspond
 *        to longer calculation time.
 * @return igraph_error_t
 */
igraph_error_t igraph_plfit_result_calculate_p_value(
    const igraph_plfit_result_t* model, igraph_real_t* result, igraph_real_t precision
) {
    int retval;
    plfit_continuous_options_t cont_options;
    plfit_discrete_options_t disc_options;
    plfit_result_t plfit_result;
    plfit_error_handler_t* plfit_stored_error_handler;
    igraph_bool_t finite_size_correction;

    IGRAPH_ASSERT(model != NULL);

    plfit_result.alpha = model->alpha;
    plfit_result.xmin = model->xmin;
    plfit_result.L = model->L;
    plfit_result.D = model->D;

    finite_size_correction = igraph_i_plfit_should_use_finite_size_correction(model->data);

    RNG_BEGIN();

    plfit_stored_error_handler = plfit_set_error_handler(igraph_i_plfit_error_handler_store);
    if (model->continuous) {
        igraph_i_plfit_prepare_continuous_options(&cont_options, finite_size_correction);
        cont_options.p_value_method = PLFIT_P_VALUE_EXACT;
        cont_options.p_value_precision = precision;
        retval = plfit_calculate_p_value_continuous(
            VECTOR(*model->data), (size_t) igraph_vector_size(model->data),
            &cont_options, /* xmin_fixed = */ 0, &plfit_result
        );
    } else {
        igraph_i_plfit_prepare_discrete_options(&disc_options, finite_size_correction);
        disc_options.p_value_method = PLFIT_P_VALUE_EXACT;
        disc_options.p_value_precision = precision;
        retval = plfit_calculate_p_value_discrete(
            VECTOR(*model->data), (size_t) igraph_vector_size(model->data),
            &disc_options, /* xmin_fixed = */ 0, &plfit_result
        );
    }
    plfit_set_error_handler(plfit_stored_error_handler);

    RNG_END();

    IGRAPH_CHECK(igraph_i_handle_plfit_error(retval));

    if (result) {
        *result = plfit_result.p;
    }

    return IGRAPH_SUCCESS;
}
