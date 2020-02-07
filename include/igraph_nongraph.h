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
#include "igraph_constants.h"
#include "igraph_matrix.h"
#include "igraph_types.h"
#include "igraph_vector.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

/**
 * \struct igraph_plfit_result_t
 * \brief Result of fitting a power-law distribution to a vector
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
 * \member p     The p-value of the Kolmogorov-Smirnov test. Small p-values
 *               (less than 0.05) indicate that the test rejected the hypothesis
 *               that the original data could have been drawn from the fitted
 *               power-law distribution.
 */
typedef struct igraph_plfit_result_t {
    igraph_bool_t continuous;
    double alpha;
    double xmin;
    double L;
    double D;
    double p;
} igraph_plfit_result_t;

DECLDIR int igraph_running_mean(const igraph_vector_t *data, igraph_vector_t *res,
                                igraph_integer_t binwidth);
DECLDIR int igraph_fisher_yates_shuffle(igraph_vector_t *seq);
DECLDIR int igraph_random_sample(igraph_vector_t *res, igraph_real_t l, igraph_real_t h,
                                 igraph_integer_t length);
DECLDIR int igraph_convex_hull(const igraph_matrix_t *data, igraph_vector_t *resverts,
                               igraph_matrix_t *rescoords);
DECLDIR int igraph_zeroin(igraph_real_t *ax, igraph_real_t *bx,
                          igraph_real_t (*f)(igraph_real_t x, void *info),
                          void *info, igraph_real_t *Tol, int *Maxit, igraph_real_t *res);
DECLDIR int igraph_bfgs(igraph_vector_t *b, igraph_real_t *Fmin,
                        igraph_scalar_function_t fminfn, igraph_vector_function_t fmingr,
                        int maxit, int trace,
                        igraph_real_t abstol, igraph_real_t reltol, int nREPORT, void *ex,
                        igraph_integer_t *fncount, igraph_integer_t *grcount);
DECLDIR int igraph_power_law_fit(const igraph_vector_t* vector, igraph_plfit_result_t* result,
                                 igraph_real_t xmin, igraph_bool_t force_continuous);

__END_DECLS

#endif
