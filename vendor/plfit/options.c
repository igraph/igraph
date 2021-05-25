/* options.c
 *
 * Copyright (C) 2012 Tamas Nepusz
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

#include "error.h"
#include "plfit.h"

const plfit_continuous_options_t plfit_continuous_default_options = {
    /* .finite_size_correction = */ 0,
    /* .xmin_method = */ PLFIT_DEFAULT_CONTINUOUS_METHOD,
    /* .p_value_method = */ PLFIT_DEFAULT_P_VALUE_METHOD,
    /* .p_value_precision = */ 0.01,
    /* .rng = */ 0
};

const plfit_discrete_options_t plfit_discrete_default_options = {
    /* .finite_size_correction = */ 0,
    /* .alpha_method = */ PLFIT_DEFAULT_DISCRETE_METHOD,
    /* .alpha = */ {
        /* .min = */ 1.01,
        /* .max = */ 5,
        /* .step = */ 0.01
    },
    /* .p_value_method = */ PLFIT_DEFAULT_P_VALUE_METHOD,
    /* .p_value_precision = */ 0.01,
    /* .rng = */ 0
};

int plfit_continuous_options_init(plfit_continuous_options_t* options) {
	*options = plfit_continuous_default_options;
	return PLFIT_SUCCESS;
}

int plfit_discrete_options_init(plfit_discrete_options_t* options) {
	*options = plfit_discrete_default_options;
	return PLFIT_SUCCESS;
}
