/*
   igraph library.
   Copyright (C) 2010-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef IGRAPH_STATUSBAR_H
#define IGRAPH_STATUSBAR_H

#include "igraph_decls.h"
#include "igraph_error.h"

IGRAPH_BEGIN_C_DECLS

/**
 * \section about_status_handlers Status reporting
 *
 * <para>
 * In addition to the possibility of reporting the progress of an
 * igraph computation via \ref igraph_progress(), it is also possible
 * to report simple status messages from within igraph functions,
 * without having to judge how much of the computation was performed
 * already. For this one needs to install a status handler function.
 * </para>
 *
 * <para>
 * Status handler functions must be of type \ref igraph_status_handler_t
 * and they can be installed by a call to \ref igraph_set_status_handler().
 * Currently there is a simple predefined status handler function,
 * called \ref igraph_status_handler_stderr(), but the user can define
 * new ones.
 * </para>
 *
 * <para>
 * igraph functions report their status via a call to the
 * \ref IGRAPH_STATUS() or the \ref IGRAPH_STATUSF() macro.
 * </para>
 */

/**
 * \typedef igraph_status_handler_t
 *
 * The type of the igraph status handler functions
 * \param message The status message.
 * \param data Additional context, with user-defined semantics.
 *        Existing igraph functions pass a null pointer here.
 * \return Error code. The current calculation will abort if you return anything
 *         else than \c IGRAPH_SUCCESS here.
 */

typedef igraph_error_t igraph_status_handler_t(const char *message, void *data);

IGRAPH_EXPORT extern igraph_status_handler_t igraph_status_handler_stderr;

IGRAPH_EXPORT igraph_status_handler_t *igraph_set_status_handler(igraph_status_handler_t new_handler);

IGRAPH_EXPORT igraph_error_t igraph_status(const char *message, void *data);

/**
 * \define IGRAPH_STATUS
 * Report the status of an igraph function.
 *
 * Typically this function is called only a handful of times from
 * an igraph function. E.g. if an algorithm has three major
 * steps, then it is logical to call it three times, to
 * signal the three major steps.
 * \param message The status message.
 * \param data Additional context, with user-defined semantics.
 *        Existing igraph functions pass a null pointer here.
 * \return If the status handler returns with a value other than
 *        \c IGRAPH_SUCCESS, then the function that called this
 *        macro returns as well, with the same error code, after
 *        cleaning up all allocated memory as needed.
 */

#define IGRAPH_STATUS(message, data) \
    do { \
        IGRAPH_CHECK(igraph_status((message), (data))); \
    } while (0)

IGRAPH_EXPORT igraph_error_t igraph_statusf(const char *message, void *data, ...);

/**
 * \define IGRAPH_STATUSF
 * Report the status from an igraph function
 *
 * This is the more flexible version of \ref IGRAPH_STATUS(),
 * having a printf-like syntax. As this macro takes variable
 * number of arguments, they must be all supplied as a single
 * argument, enclosed in parentheses. \ref igraph_statusf() is then
 * called with the given arguments.
 *
 * \param args The arguments to pass to \ref igraph_statusf().
 * \return If the status handler returns with a value other than
 *        \c IGRAPH_SUCCESS, then the function that called this
 *        macro returns as well, with the same error code, after
 *        cleaning up all allocated memory as needed.
 */

#define IGRAPH_STATUSF(args) \
    do { \
        IGRAPH_CHECK(igraph_statusf args); \
    } while (0)

IGRAPH_END_C_DECLS

#endif
