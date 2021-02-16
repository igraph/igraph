/* gss.c
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

#include <float.h>
#include <math.h>
#include <string.h>
#include "error.h"
#include "gss.h"
#include "platform.h"

/**
 * \def PHI
 *
 * The golden ratio, i.e. 1+sqrt(5)/2
 */
#define PHI 1.618033988749895

/**
 * \def RESPHI
 *
 * Constant defined as 2 - \c PHI
 */
#define RESPHI 0.3819660112501051

/**
 * \const _defparam
 *
 * Default parameters for the GSS algorithm.
 */
static const gss_parameter_t _defparam = {
    /* .epsilon = */  DBL_MIN,
	/* .on_error = */ GSS_ERROR_STOP
};

/**
 * Stores whether the last optimization run triggered a warning or not.
 */
static unsigned short int gss_i_warning_flag = 0;

void gss_parameter_init(gss_parameter_t *param) {
    memcpy(param, &_defparam, sizeof(*param));
}

unsigned short int gss_get_warning_flag() {
	return gss_i_warning_flag;
}

#define TERMINATE {        \
    if (_min) {            \
        *(_min) = min;     \
    }                      \
    if (_fmin) {           \
        *(_fmin) = fmin;   \
    }                      \
}

#define EVALUATE(x, fx) { \
    fx = proc_evaluate(instance, x); \
    if (fmin > fx) { \
        min = x;     \
        fmin = fx;   \
    } \
    if (proc_progress) { \
        retval = proc_progress(instance, x, fx, min, fmin, \
                (a < b) ? a : b, (a < b) ? b : a, k); \
        if (retval) { \
			TERMINATE;            \
            return PLFIT_SUCCESS; \
        } \
    } \
}

int gss(double a, double b, double *_min, double *_fmin,
        gss_evaluate_t proc_evaluate, gss_progress_t proc_progress,
        void* instance, const gss_parameter_t *_param) {
    double c, d, min;
    double fa, fb, fc, fd, fmin;
    int k = 0;
    int retval;
    unsigned short int successful = 1;

    gss_parameter_t param = _param ? (*_param) : _defparam;

	gss_i_warning_flag = 0;

    if (a > b) {
        c = a; a = b; b = c;
    }

    min = a;
    fmin = proc_evaluate(instance, a);

    c = a + RESPHI*(b-a);

    EVALUATE(a, fa);
    EVALUATE(b, fb);
    EVALUATE(c, fc);

    if (fc >= fa || fc >= fb) {
		if (param.on_error == GSS_ERROR_STOP) {
			return PLFIT_FAILURE;
		} else {
			gss_i_warning_flag = 1;
		}
	}

    while (fabs(a-b) > param.epsilon) {
        k++;

        d = c + RESPHI*(b-c);
        EVALUATE(d, fd);

        if (fd >= fa || fd >= fb) {
			if (param.on_error == GSS_ERROR_STOP) {
				successful = 0;
				break;
			} else {
				gss_i_warning_flag = 1;
			}
        }

        if (fc <= fd) {
            b = a; a = d;
        } else {
            a = c; c = d; fc = fd;
        }
    }

    if (successful) {
        c = (a+b) / 2.0;
        k++;
        EVALUATE(c, fc);
		TERMINATE;
    }

    return successful ? PLFIT_SUCCESS : PLFIT_FAILURE;
}
