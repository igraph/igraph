/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_progress.h"
#include "config.h"

static IGRAPH_THREAD_LOCAL igraph_progress_handler_t *igraph_i_progress_handler = 0;
static IGRAPH_THREAD_LOCAL char igraph_i_progressmsg_buffer[1000];

/**
 * \function igraph_progress
 * Report progress
 *
 * Note that the usual way to report progress is the \ref IGRAPH_PROGRESS
 * macro, as that takes care of the return value of the progress
 * handler.
 * \param message A string describing the function or algorithm
 *     that is reporting the progress. Current igraph functions
 *     always use the name \p message argument if reporting from the
 *     same function.
 * \param percent Numeric, the percentage that was completed by the
 *     algorithm or function.
 * \param data User-defined data. Current igraph functions that
 *     report progress pass a null pointer here. Users can
 *     write their own progress handlers and functions with progress
 *     reporting, and then pass some meaningfull context here.
 * \return If there is a progress handler installed and
 *     it does not return \c IGRAPH_SUCCESS, then \c IGRAPH_INTERRUPTED
 *     is returned.
 *
 * Time complexity: O(1).
 */

int igraph_progress(const char *message, igraph_real_t percent, void *data) {
    if (igraph_i_progress_handler) {
        if (igraph_i_progress_handler(message, percent, data) != IGRAPH_SUCCESS) {
            return IGRAPH_INTERRUPTED;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_progressf
 * Report progress, printf-like version
 *
 * This is a more flexible version of \ref igraph_progress(), with
 * a printf-like template string. First the template string
 * is filled with the additional arguments and then \ref
 * igraph_progress() is called.
 *
 * </para><para>Note that there is an upper limit for the length of
 * the \p message string, currently 1000 characters.
 * \param message A string describing the function or algorithm
 *     that is reporting the progress. For this function this is a
 *     template string, using the same syntax as the standard
 *     \c libc \c printf function.
 * \param percent Numeric, the percentage that was completed by the
 *     algorithm or function.
 * \param data User-defined data. Current igraph functions that
 *     report progress pass a null pointer here. Users can
 *     write their own progress handlers and functions with progress
 *     reporting, and then pass some meaningfull context here.
 * \param ... Additional argument that were specified in the
 *     \p message argument.
 * \return If there is a progress handler installed and
 *     it does not return \c IGRAPH_SUCCESS, then \c IGRAPH_INTERRUPTED
 *     is returned.
 * \return
 */

int igraph_progressf(const char *message, igraph_real_t percent, void *data,
                     ...) {
    va_list ap;
    va_start(ap, data);
    vsnprintf(igraph_i_progressmsg_buffer,
              sizeof(igraph_i_progressmsg_buffer) / sizeof(char), message, ap);
    return igraph_progress(igraph_i_progressmsg_buffer, percent, data);
}

#ifndef USING_R

/**
 * \function igraph_progress_handler_stderr
 * A simple predefined progress handler
 *
 * This simple progress handler first prints \p message, and then
 * the percentage complete value in a short message to standard error.
 * \param message A string describing the function or algorithm
 *     that is reporting the progress. Current igraph functions
 *     always use the name \p message argument if reporting from the
 *     same function.
 * \param percent Numeric, the percentage that was completed by the
 *     algorithm or function.
 * \param data User-defined data. Current igraph functions that
 *     report progress pass a null pointer here. Users can
 *     write their own progress handlers and functions with progress
 *     reporting, and then pass some meaningfull context here.
 * \return This function always returns with \c IGRAPH_SUCCESS.
 *
 * Time complexity: O(1).
 */

int igraph_progress_handler_stderr(const char *message, igraph_real_t percent,
                                   void* data) {
    IGRAPH_UNUSED(data);
    fputs(message, stderr);
    fprintf(stderr, "%.1f percent ready\n", (double)percent);
    return 0;
}
#endif

/**
 * \function igraph_set_progress_handler
 * Install a progress handler, or remove the current handler
 *
 * There is a single simple predefined progress handler:
 * \ref igraph_progress_handler_stderr().
 * \param new_handler Pointer to a function of type
 *     \ref igraph_progress_handler_t, the progress handler function to
 *     install. To uninstall the current progress handler, this argument
 *     can be a null pointer.
 * \return Pointer to the previously installed progress handler function.
 *
 * Time complexity: O(1).
 */

igraph_progress_handler_t *
igraph_set_progress_handler(igraph_progress_handler_t new_handler) {
    igraph_progress_handler_t *previous_handler = igraph_i_progress_handler;
    igraph_i_progress_handler = new_handler;
    return previous_handler;
}
