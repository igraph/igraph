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

#include "core/interruption.h"
#include "plfit/error.h"
#include "plfit/plfit.h"

#include <math.h>

/**
 * \ingroup nongraph
 * \function igraph_running_mean
 * \brief Calculates the running mean of a vector.
 *
 * </para><para>
 * The running mean is defined by the mean of the
 * previous \p binwidth values.
 * \param data The vector containing the data.
 * \param res The vector containing the result. This should be
 *        initialized before calling this function and will be
 *        resized.
 * \param binwidth Integer giving the width of the bin for the running
 *        mean calculation.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the length of
 * the data vector.
 */

int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res,
                        igraph_integer_t binwidth) {

    double sum = 0;
    long int i;

    /* Check */
    if (igraph_vector_size(data) < binwidth) {
        IGRAPH_ERRORF("Data vector length (%ld) smaller than bin width (%" IGRAPH_PRId ").", IGRAPH_EINVAL, igraph_vector_size(data), binwidth);
    }
    if (binwidth < 1) {
        IGRAPH_ERRORF("Bin width for running mean should be at least 1, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, igraph_vector_size(data), binwidth);
    }

    /* Memory for result */

    IGRAPH_CHECK(igraph_vector_resize(res, (long int)(igraph_vector_size(data) - binwidth + 1)));

    /* Initial bin */
    for (i = 0; i < binwidth; i++) {
        sum += VECTOR(*data)[i];
    }

    VECTOR(*res)[0] = sum / binwidth;

    for (i = 1; i < igraph_vector_size(data) - binwidth + 1; i++) {
        IGRAPH_ALLOW_INTERRUPTION();
        sum -= VECTOR(*data)[i - 1];
        sum += VECTOR(*data)[ (long int)(i + binwidth - 1)];
        VECTOR(*res)[i] = sum / binwidth;
    }

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup nongraph
 * \function igraph_convex_hull
 * \brief Determines the convex hull of a given set of points in the 2D plane
 *
 * </para><para>
 * The convex hull is determined by the Graham scan algorithm.
 * See the following reference for details:
 *
 * </para><para>
 * Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford
 * Stein. Introduction to Algorithms, Second Edition. MIT Press and
 * McGraw-Hill, 2001. ISBN 0262032937. Pages 949-955 of section 33.3:
 * Finding the convex hull.
 *
 * \param data vector containing the coordinates. The length of the
 *        vector must be even, since it contains X-Y coordinate pairs.
 * \param resverts the vector containing the result, e.g. the vector of
 *        vertex indices used as the corners of the convex hull. Supply
 *        \c NULL here if you are only interested in the coordinates of
 *        the convex hull corners.
 * \param rescoords the matrix containing the coordinates of the selected
 *        corner vertices. Supply \c NULL here if you are only interested in
 *        the vertex indices.
 * \return Error code:
 *         \c IGRAPH_ENOMEM: not enough memory
 *
 * Time complexity: O(n log(n)) where n is the number of vertices
 *
 * \example examples/simple/igraph_convex_hull.c
 */
int igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_t *resverts,
                       igraph_matrix_t *rescoords) {
    igraph_integer_t no_of_nodes;
    long int i, pivot_idx = 0, last_idx, before_last_idx, next_idx, j;
    igraph_vector_t angles, stack, order;
    igraph_real_t px, py, cp;

    no_of_nodes = (igraph_integer_t) igraph_matrix_nrow(data);
    if (igraph_matrix_ncol(data) != 2) {
        IGRAPH_ERROR("matrix must have 2 columns", IGRAPH_EINVAL);
    }
    if (no_of_nodes == 0) {
        if (resverts != 0) {
            IGRAPH_CHECK(igraph_vector_resize(resverts, 0));
        }
        if (rescoords != 0) {
            IGRAPH_CHECK(igraph_matrix_resize(rescoords, 0, 2));
        }
        /**************************** this is an exit here *********/
        return 0;
    }

    IGRAPH_VECTOR_INIT_FINALLY(&angles, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&stack, 0);

    /* Search for the pivot vertex */
    for (i = 1; i < no_of_nodes; i++) {
        if (MATRIX(*data, i, 1) < MATRIX(*data, pivot_idx, 1)) {
            pivot_idx = i;
        } else if (MATRIX(*data, i, 1) == MATRIX(*data, pivot_idx, 1) &&
                   MATRIX(*data, i, 0) < MATRIX(*data, pivot_idx, 0)) {
            pivot_idx = i;
        }
    }
    px = MATRIX(*data, pivot_idx, 0);
    py = MATRIX(*data, pivot_idx, 1);

    /* Create angle array */
    for (i = 0; i < no_of_nodes; i++) {
        if (i == pivot_idx) {
            /* We can't calculate the angle of the pivot point with itself,
             * so we use 10 here. This way, after sorting the angle vector,
             * the pivot point will always be the first one, since the range
             * of atan2 is -3.14..3.14 */
            VECTOR(angles)[i] = 10;
        } else {
            VECTOR(angles)[i] = atan2(MATRIX(*data, i, 1) - py, MATRIX(*data, i, 0) - px);
        }
    }

    /* Sort points by angles */
    IGRAPH_VECTOR_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_qsort_ind(&angles, &order, 0));

    /* Check if two points have the same angle. If so, keep only the point that
     * is farthest from the pivot */
    j = 0;
    last_idx = (long int) VECTOR(order)[0];
    pivot_idx = (long int) VECTOR(order)[no_of_nodes - 1];
    for (i = 1; i < no_of_nodes; i++) {
        next_idx = (long int) VECTOR(order)[i];
        if (VECTOR(angles)[last_idx] == VECTOR(angles)[next_idx]) {
            /* Keep the vertex that is farther from the pivot, drop the one that is
             * closer */
            px = pow(MATRIX(*data, last_idx, 0) - MATRIX(*data, pivot_idx, 0), 2) +
                 pow(MATRIX(*data, last_idx, 1) - MATRIX(*data, pivot_idx, 1), 2);
            py = pow(MATRIX(*data, next_idx, 0) - MATRIX(*data, pivot_idx, 0), 2) +
                 pow(MATRIX(*data, next_idx, 1) - MATRIX(*data, pivot_idx, 1), 2);
            if (px > py) {
                VECTOR(order)[i] = -1;
            } else {
                VECTOR(order)[j] = -1;
                last_idx = next_idx;
                j = i;
            }
        } else {
            last_idx = next_idx;
            j = i;
        }
    }

    j = 0;
    last_idx = -1;
    before_last_idx = -1;
    while (!igraph_vector_empty(&order)) {
        next_idx = (long int)VECTOR(order)[igraph_vector_size(&order) - 1];
        if (next_idx < 0) {
            /* This vertex should be skipped; was excluded in an earlier step */
            igraph_vector_pop_back(&order);
            continue;
        }
        /* Determine whether we are at a left or right turn */
        if (j < 2) {
            /* Pretend that we are turning into the right direction if we have less
             * than two items in the stack */
            cp = -1;
        } else {
            cp = (MATRIX(*data, last_idx, 0) - MATRIX(*data, before_last_idx, 0)) *
                 (MATRIX(*data, next_idx, 1) - MATRIX(*data, before_last_idx, 1)) -
                 (MATRIX(*data, next_idx, 0) - MATRIX(*data, before_last_idx, 0)) *
                 (MATRIX(*data, last_idx, 1) - MATRIX(*data, before_last_idx, 1));
        }
        /*
        printf("B L N cp: %ld, %ld, %ld, %f [", before_last_idx, last_idx, next_idx, (float)cp);
        for (int k=0; k<j; k++) printf("%ld ", (long)VECTOR(stack)[k]);
        printf("]\n");
        */
        if (cp < 0) {
            /* We are turning into the right direction */
            igraph_vector_pop_back(&order);
            IGRAPH_CHECK(igraph_vector_push_back(&stack, next_idx));
            before_last_idx = last_idx;
            last_idx = next_idx;
            j++;
        } else {
            /* No, skip back and try again in the next iteration */
            igraph_vector_pop_back(&stack);
            j--;
            last_idx = before_last_idx;
            before_last_idx = (j >= 2) ? (long int) VECTOR(stack)[j - 2] : -1;
        }
    }

    /* Create result vector */
    if (resverts != 0) {
        igraph_vector_clear(resverts);
        IGRAPH_CHECK(igraph_vector_append(resverts, &stack));
    }
    if (rescoords != 0) {
        igraph_matrix_select_rows(data, rescoords, &stack);
    }

    /* Free everything */
    igraph_vector_destroy(&order);
    igraph_vector_destroy(&stack);
    igraph_vector_destroy(&angles);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}


static const char* igraph_i_plfit_error_message = 0;

static void igraph_i_plfit_error_handler_store(const char *reason, const char *file,
        int line, int plfit_errno) {

    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(plfit_errno);

    igraph_i_plfit_error_message = reason;
}

/**
 * \ingroup nongraph
 * \function igraph_power_law_fit
 * \brief Fits a power-law distribution to a vector of numbers
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
 * Aaron Clauset, Cosma R .Shalizi and Mark E.J. Newman: Power-law
 * distributions in empirical data. SIAM Review 51(4):661-703, 2009.
 *
 * \param data vector containing the samples for which a power-law distribution
 *             is to be fitted. Note that you have to provide the \em samples,
 *             not the probability density function or the cumulative
 *             distribution function. For example, if you wish to fit
 *             a power-law to the degrees of a graph, you can use the output of
 *             \ref igraph_degree directly as an input argument to
 *             \ref igraph_power_law_fit
 * \param result the result of the fitting algorithm. See \ref igraph_plfit_result_t
 *             for more details.
 * \param xmin the minimum value in the sample vector where the power-law
 *             behaviour is expected to kick in. Samples smaller than \c xmin
 *             will be ignored by the algoritm. Pass zero here if you want to
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
int igraph_power_law_fit(const igraph_vector_t* data, igraph_plfit_result_t* result,
                         igraph_real_t xmin, igraph_bool_t force_continuous) {
    plfit_error_handler_t* plfit_stored_error_handler;
    plfit_result_t plfit_result;
    plfit_continuous_options_t cont_options;
    plfit_discrete_options_t disc_options;
    igraph_bool_t discrete = force_continuous ? 0 : 1;
    igraph_bool_t finite_size_correction;
    int retval;
    size_t i, n;

    n = (size_t) igraph_vector_size(data);
    finite_size_correction = (n < 50);

    if (discrete) {
        /* Does the vector contain discrete values only? */
        for (i = 0; i < n; i++) {
            if ((long int)(VECTOR(*data)[i]) != VECTOR(*data)[i]) {
                discrete = 0;
                break;
            }
        }
    }

    RNG_BEGIN();

    plfit_stored_error_handler = plfit_set_error_handler(igraph_i_plfit_error_handler_store);
    if (discrete) {
        plfit_discrete_options_init(&disc_options);
        /* TODO: approximation method should be switched to PLFIT_P_VALUE_EXACT in igraph 0.9 */
        disc_options.p_value_method = PLFIT_P_VALUE_APPROXIMATE;
        disc_options.finite_size_correction = (plfit_bool_t) finite_size_correction;

        if (xmin >= 0) {
            retval = plfit_estimate_alpha_discrete(VECTOR(*data), n, xmin,
                                                   &disc_options, &plfit_result);
        } else {
            retval = plfit_discrete(VECTOR(*data), n, &disc_options, &plfit_result);
        }
    } else {
        plfit_continuous_options_init(&cont_options);
        /* TODO: approximation method should be switched to PLFIT_P_VALUE_EXACT in igraph 0.9 */
        cont_options.p_value_method = PLFIT_P_VALUE_APPROXIMATE;
        /* TODO: xmin method should be switched to PLFIT_STRATIFIED_SAMPLING in igraph 0.9 */
        cont_options.xmin_method = PLFIT_GSS_OR_LINEAR;
        cont_options.finite_size_correction = (plfit_bool_t) finite_size_correction;

        if (xmin >= 0) {
            retval = plfit_estimate_alpha_continuous(VECTOR(*data), n, xmin,
                     &cont_options, &plfit_result);
        } else {
            retval = plfit_continuous(VECTOR(*data), n, &cont_options, &plfit_result);
        }
    }
    plfit_set_error_handler(plfit_stored_error_handler);

    RNG_END();

    switch (retval) {
    case PLFIT_FAILURE:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_FAILURE);
        break;

    case PLFIT_EINVAL:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EINVAL);
        break;

    case PLFIT_UNDRFLOW:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EUNDERFLOW);
        break;

    case PLFIT_OVERFLOW:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_EOVERFLOW);
        break;

    case PLFIT_ENOMEM:
        IGRAPH_ERROR(igraph_i_plfit_error_message, IGRAPH_ENOMEM);
        break;

    default:
        break;
    }

    if (result) {
        result->continuous = !discrete;
        result->alpha = plfit_result.alpha;
        result->xmin = plfit_result.xmin;
        result->L = plfit_result.L;
        result->D = plfit_result.D;
        result->p = plfit_result.p;
    }

    return 0;
}

