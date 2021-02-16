/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_statusbar.h"
#include "igraph_error.h"

#include "config.h"

#include <stdio.h>
#include <stdarg.h>

static IGRAPH_THREAD_LOCAL igraph_status_handler_t *igraph_i_status_handler = 0;

/**
 * \function igraph_status
 * Report status from an igraph function.
 *
 * It calls the installed status handler function, if there is
 * one. Otherwise it does nothing. Note that the standard way to
 * report the status from an igraph function is the
 * \ref IGRAPH_STATUS or \ref IGRAPH_STATUSF macro, as these
 * take care of the termination of the calling function if the
 * status handler returns with \c IGRAPH_INTERRUPTED.
 * \param message The status message.
 * \param data Additional context, with user-defined semantics.
 *        Existing igraph functions pass a null pointer here.
 * \return Error code. If a status handler function was called
 *        and it did not return with \c IGRAPH_SUCCESS, then
 *        \c IGRAPH_INTERRUPTED is returned by \c igraph_status().
 *
 * Time complexity: O(1).
 */

int igraph_status(const char *message, void *data) {
    if (igraph_i_status_handler) {
        if (igraph_i_status_handler(message, data) != IGRAPH_SUCCESS) {
            return IGRAPH_INTERRUPTED;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_statusf
 * Report status, more flexible printf-like version.
 *
 * This is the more flexible version of \ref igraph_status(),
 * that has a syntax similar to the \c printf standard C library function.
 * It substitutes the values of the additional arguments into the
 * \p message template string and calls \ref igraph_status().
 * \param message Status message template string, the syntax is the same
 *        as for the \c printf function.
 * \param data Additional context, with user-defined semantics.
 *        Existing igraph functions pass a null pointer here.
 * \param ... The additional arguments to fill the template given in the
 *        \p message argument.
 * \return Error code. If a status handler function was called
 *        and it did not return with \c IGRAPH_SUCCESS, then
 *        \c IGRAPH_INTERRUPTED is returned by \c igraph_status().
 */

int igraph_statusf(const char *message, void *data, ...) {
    char buffer[300];
    va_list ap;
    va_start(ap, data);
    vsnprintf(buffer, sizeof(buffer) - 1, message, ap);
    return igraph_status(buffer, data);
}

#ifndef USING_R

/**
 * \function igraph_status_handler_stderr
 * A simple predefined status handler function.
 *
 * A simple status handler function, that writes the status
 * message to the standard errror.
 * \param message The status message.
 * \param data Additional context, with user-defined semantics.
 *        Existing igraph functions pass a null pointer here.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_status_handler_stderr(const char *message, void *data) {
    IGRAPH_UNUSED(data);
    fputs(message, stderr);
    return 0;
}
#endif

/**
 * \function igraph_set_status_handler
 * Install of uninstall a status handler function.
 *
 * To uninstall the currently installed status handler, call
 * this function with a null pointer.
 * \param new_handler The status handler function to install.
 * \return The previously installed status handler function.
 *
 * Time complexity: O(1).
 */

igraph_status_handler_t *
igraph_set_status_handler(igraph_status_handler_t new_handler) {
    igraph_status_handler_t *previous_handler = igraph_i_status_handler;
    igraph_i_status_handler = new_handler;
    return previous_handler;
}

