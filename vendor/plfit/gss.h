/* gss.h
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

#ifndef __GSS_H__
#define __GSS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/**
 * Enum specifying what the search should do when the function is not U-shaped.
 */
typedef enum {
	GSS_ERROR_STOP,              /**< Stop and return an error code */
	GSS_ERROR_WARN               /**< Continue and set the warning flag */
} gss_error_handling_t;

/**
 * Parameter settings for a golden section search.
 */
typedef struct {
    double epsilon;
	gss_error_handling_t on_error;
} gss_parameter_t;

/**
 * Callback interface to provide objective function evaluations for the golden
 * section search.
 *
 * The gss() function calls this function to obtain the values of the objective
 * function when needed. A client program must implement this function to evaluate
 * the value of the objective function, given the location.
 *
 * @param  instance    The user data sent for the gss() function by the client.
 * @param  x           The current value of the variable.
 * @retval double      The value of the objective function for the current
 *                      variable.
 */
typedef double (*gss_evaluate_t)(void *instance, double x);

/**
 * Callback interface to receive the progress of the optimization process for
 * the golden section search.
 *
 * The gss() function calls this function for each iteration. Implementing
 * this function, a client program can store or display the current progress
 * of the optimization process.
 *
 * @param  instance    The user data sent for the gss() function by the client.
 * @param  x           The current value of the variable.
 * @param  fx          The value of the objective function at x.
 * @param  min         The location of the minimum value of the objective
 *                     function found so far.
 * @param  fmin        The minimum value of the objective function found so far.
 * @param  left        The left side of the current bracket.
 * @param  right       The right side of the current bracket.
 * @param  k           The index of the current iteration.
 * @retval int         Zero to continue the optimization process. Returning a
 *                     non-zero value will cancel the optimization process.
 */
typedef int (*gss_progress_t)(void *instance, double x, double fx, double min,
        double fmin, double left, double right, int k);

/**
 * Start a golden section search optimization.
 *
 * @param  a    The left side of the bracket to start from
 * @param  b    The right side of the bracket to start from
 * @param  min  The pointer to the variable that receives the location of the
 *              final value of the objective function. This argument can be set to
 *              \c NULL if the location of the final value of the objective
 *              function is unnecessary.
 * @param  fmin The pointer to the variable that receives the final value of
 *              the objective function. This argument can be st to \c NULL if the
 *              final value of the objective function is unnecessary.
 * @param  proc_evaluate  The callback function to evaluate the objective
 *                        function at a given location.
 * @param  proc_progress  The callback function to receive the progress (the
 *                        last evaluated location, the value of the objective
 *                        function at that location, the width of the current
 *                        bracket, the minimum found so far and the step
 *                        count). This argument can be set to \c NULL if
 *                        a progress report is unnecessary.
 * @param  instance    A user data for the client program. The callback
 *                     functions will receive the value of this argument.
 * @param  param       The pointer to a structure representing parameters for
 *                     GSS algorithm. A client program can set this parameter
 *                     to \c NULL to use the default parameters.
 *                     Call the \ref gss_parameter_init() function to fill a
 *                     structure with the default values.
 * @retval int         The status code. This function returns zero if the
 *                     minimization process terminates without an error. A
 *                     non-zero value indicates an error; in particular,
 *                     \c PLFIT_FAILURE means that the function is not
 *                     U-shaped.
 */
int gss(double a, double b, double *min, double *fmin,
        gss_evaluate_t proc_evaluate, gss_progress_t proc_progress,
        void* instance, const gss_parameter_t *_param);

/**
 * Return the state of the warning flag.
 *
 * The warning flag is 1 if the last optimization was run on a function that
 * was not U-shaped.
 */
unsigned short int gss_get_warning_flag();

/**
 * Initialize GSS parameters to the default values.
 *
 * Call this function to fill a parameter structure with the default values
 * and overwrite parameter values if necessary.
 *
 * @param  param       The pointer to the parameter structure.
 */
void gss_parameter_init(gss_parameter_t *param);

__END_DECLS

#endif /* __GSS_H__ */
