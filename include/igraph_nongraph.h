/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_NONGRAPH_H
#define IGRAPH_NONGRAPH_H

#include "igraph_decls.h"
#include "igraph_error.h"
#include "igraph_matrix.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/**
 * \def IGRAPH_SHORTEST_PATH_EPSILON
 *
 * Relative error threshold used in weighted shortest path calculations
 * to decide whether two shortest paths are of equal length.
 */
#define IGRAPH_SHORTEST_PATH_EPSILON 1e-10

typedef igraph_real_t  igraph_scalar_function_t(const igraph_vector_t *var,
        const igraph_vector_t *par,
        void* extra);
typedef void igraph_vector_function_t(const igraph_vector_t *var,
                                      const igraph_vector_t *par,
                                      igraph_vector_t* res, void* extra);

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

/**
 * \struct igraph_plfit_result_t
 * \brief Result of fitting a power-law distribution to a vector.
 *
 * This data structure contains the result of \ref igraph_power_law_fit(),
 * which tries to fit a power-law distribution to a vector of numbers. The
 * structure contains the following members:
 *
 * \member continuous Whether the fitted power-law distribution was continuous
 *                    or discrete.
 * \member alpha The exponent of the fitted power-law distribution.
 * \member xmin  The minimum value from which the power-law distribution was
 *               fitted. In other words, only the values larger than \c xmin
 *               were used from the input vector.
 * \member L     The log-likelihood of the fitted parameters; in other words,
 *               the probability of observing the input vector given the
 *               parameters.
 * \member D     The test statistic of a Kolmogorov-Smirnov test that compares
 *               the fitted distribution with the input vector. Smaller scores
 *               denote better fit.
 * \member p     The p-value of the Kolmogorov-Smirnov test; \c NaN if it has
 *               not been calculated yet. Small p-values (less than 0.05)
 *               indicate that the test rejected the hypothesis that the
 *               original data could have been drawn from the fitted power-law
 *               distribution.
 * \member data  The vector containing the original input data. May not be valid
 *               any more if the caller already destroyed the vector.
 */
typedef struct igraph_plfit_result_t {
    igraph_bool_t continuous;
    igraph_real_t alpha;
    igraph_real_t xmin;
    igraph_real_t L;
    igraph_real_t D;
    const igraph_vector_t* data;
} igraph_plfit_result_t;

IGRAPH_EXPORT igraph_error_t igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res,
                                      igraph_integer_t binwidth);
IGRAPH_EXPORT igraph_error_t igraph_random_sample(igraph_vector_int_t *res, igraph_integer_t l, igraph_integer_t h,
                                       igraph_integer_t length);
IGRAPH_EXPORT igraph_error_t igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_int_t *resverts,
                                     igraph_matrix_t *rescoords);
IGRAPH_EXPORT igraph_bool_t igraph_almost_equals(double a, double b, double eps);
IGRAPH_EXPORT int igraph_cmp_epsilon(double a, double b, double eps);

IGRAPH_EXPORT igraph_error_t igraph_power_law_fit(
    const igraph_vector_t* vector, igraph_plfit_result_t* result,
    igraph_real_t xmin, igraph_bool_t force_continuous
);
IGRAPH_EXPORT igraph_error_t igraph_plfit_result_calculate_p_value(
    const igraph_plfit_result_t* model, igraph_real_t* result, igraph_real_t precision
);

IGRAPH_EXPORT IGRAPH_DEPRECATED igraph_error_t igraph_zeroin(
    igraph_real_t *ax, igraph_real_t *bx, igraph_real_t (*f)(igraph_real_t x, void *info),
    void *info, igraph_real_t *Tol, int *Maxit, igraph_real_t *res
);

__END_DECLS

#endif
