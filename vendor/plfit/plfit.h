/* plfit.h
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

#ifndef PLFIT_H
#define PLFIT_H

#include <stdlib.h>
#include "plfit_decls.h"
#include "plfit_error.h"
#include "plfit_mt.h"
#include "plfit_sampling.h"
#include "plfit_version.h"

PLFIT_BEGIN_C_DECLS

typedef unsigned short int plfit_bool_t;

typedef enum {
    PLFIT_LINEAR_ONLY,
    PLFIT_STRATIFIED_SAMPLING,
    PLFIT_GSS_OR_LINEAR,
    PLFIT_DEFAULT_CONTINUOUS_METHOD = PLFIT_STRATIFIED_SAMPLING
} plfit_continuous_method_t;

typedef enum {
    PLFIT_LBFGS,
    PLFIT_LINEAR_SCAN,
    PLFIT_PRETEND_CONTINUOUS,
    PLFIT_DEFAULT_DISCRETE_METHOD = PLFIT_LBFGS
} plfit_discrete_method_t;

typedef enum {
    PLFIT_P_VALUE_SKIP,
    PLFIT_P_VALUE_APPROXIMATE,
    PLFIT_P_VALUE_EXACT,
    PLFIT_DEFAULT_P_VALUE_METHOD = PLFIT_P_VALUE_EXACT
} plfit_p_value_method_t;

typedef struct _plfit_result_t {
    double alpha;     /* fitted power-law exponent */
    double xmin;      /* cutoff where the power-law behaviour kicks in */
    double L;         /* log-likelihood of the sample */
    double D;         /* test statistic for the KS test */
    double p;         /* p-value of the KS test */
} plfit_result_t;

/********** structure that holds the options of plfit **********/

typedef struct _plfit_continuous_options_t {
    plfit_bool_t finite_size_correction;
    plfit_continuous_method_t xmin_method;
    plfit_p_value_method_t p_value_method;
    double p_value_precision;
    plfit_mt_rng_t* rng;
} plfit_continuous_options_t;

typedef struct _plfit_discrete_options_t {
    plfit_bool_t finite_size_correction;
    plfit_discrete_method_t alpha_method;
    struct {
        double min;
        double max;
        double step;
    } alpha;
    plfit_p_value_method_t p_value_method;
    double p_value_precision;
    plfit_mt_rng_t* rng;
} plfit_discrete_options_t;

PLFIT_EXPORT int plfit_continuous_options_init(plfit_continuous_options_t* options);
PLFIT_EXPORT int plfit_discrete_options_init(plfit_discrete_options_t* options);

PLFIT_EXPORT extern const plfit_continuous_options_t plfit_continuous_default_options;
PLFIT_EXPORT extern const plfit_discrete_options_t plfit_discrete_default_options;

/********** continuous power law distribution fitting **********/

PLFIT_EXPORT int plfit_log_likelihood_continuous(const double* xs, size_t n, double alpha,
        double xmin, double* l);
PLFIT_EXPORT int plfit_estimate_alpha_continuous(const double* xs, size_t n, double xmin,
        const plfit_continuous_options_t* options, plfit_result_t* result);
PLFIT_EXPORT int plfit_continuous(const double* xs, size_t n,
        const plfit_continuous_options_t* options, plfit_result_t* result);

/*********** discrete power law distribution fitting ***********/

PLFIT_EXPORT int plfit_estimate_alpha_discrete(const double* xs, size_t n, double xmin,
        const plfit_discrete_options_t* options, plfit_result_t *result);
PLFIT_EXPORT int plfit_log_likelihood_discrete(const double* xs, size_t n, double alpha, double xmin, double* l);
PLFIT_EXPORT int plfit_discrete(const double* xs, size_t n, const plfit_discrete_options_t* options,
        plfit_result_t* result);

/***** resampling routines to generate synthetic replicates ****/

PLFIT_EXPORT int plfit_resample_continuous(const double* xs, size_t n, double alpha, double xmin,
        size_t num_samples, plfit_mt_rng_t* rng, double* result);
PLFIT_EXPORT int plfit_resample_discrete(const double* xs, size_t n, double alpha, double xmin,
        size_t num_samples, plfit_mt_rng_t* rng, double* result);

/******** calculating the p-value of a fitted model only *******/

PLFIT_EXPORT int plfit_calculate_p_value_continuous(const double* xs, size_t n,
        const plfit_continuous_options_t* options, plfit_bool_t xmin_fixed,
        plfit_result_t *result);
PLFIT_EXPORT int plfit_calculate_p_value_discrete(const double* xs, size_t n,
        const plfit_discrete_options_t* options, plfit_bool_t xmin_fixed,
        plfit_result_t *result);

/************* calculating descriptive statistics **************/

PLFIT_EXPORT int plfit_moments(const double* data, size_t n, double* mean, double* variance,
        double* skewness, double* kurtosis);

PLFIT_END_C_DECLS

#endif /* PLFIT_H */
