/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_INTERRUPT_H
#define IGRAPH_INTERRUPT_H

#include "igraph_decls.h"
#include "igraph_error.h"

__BEGIN_DECLS

/* This file contains the igraph interruption handling. */

/**
 * \section interrupthandlers Interruption handlers
 *
 * <para>
 * \a igraph is designed to be embeddable into several higher level
 * languages (R and Python interfaces are included in the original
 * package). Since most higher level languages consider internal \a igraph
 * calls as atomic, interruption requests (like Ctrl-C in Python) must
 * be handled differently depending on the environment \a igraph embeds
 * into.</para>
 * <para>
 * An \emb interruption handler \eme is a function which is called regularly
 * by \a igraph during long calculations. A typical usage of the interruption
 * handler is to check whether the user tried to interrupt the calculation
 * and return an appropriate value to signal this condition. For example,
 * in R, one must call an internal R function regularly to check for
 * interruption requests, and the \a igraph interruption handler is the
 * perfect place to do that.</para>
 * <para>
 * If you are using the plain C interface of \a igraph or if you are
 * allowed to replace the operating system's interruption handler (like
 * SIGINT in Un*x systems), these calls are not of much use to you.</para>
 * <para>
 * The default interruption handler is empty.
 * The \ref igraph_set_interruption_handler() function can be used to set a
 * new interruption handler function of type
 * \ref igraph_interruption_handler_t, see the
 * documentation of this type for details.
 * </para>
 */

/**
 * \section writing_interruption_handlers Writing interruption handlers
 *
 * <para>
 * You can write and install interruption handlers simply by defining a
 * function of type \ref igraph_interruption_handler_t and calling
 * \ref igraph_set_interruption_handler(). This feature is useful for
 * interface writers, because usually this is the only way to allow handling
 * of Ctrl-C and similar keypresses properly.
 * </para>
 * <para>
 * Your interruption handler will be called regularly during long operations
 * (so it is not guaranteed to be called during operations which tend to be
 * short, like adding single edges). An interruption handler accepts no
 * parameters and must return \c IGRAPH_SUCCESS if the calculation should go on. All
 * other return values are considered to be a request for interruption,
 * and the caller function would return a special error code, \c IGRAPH_INTERRUPTED.
 * It is up to your error handler function to handle this error properly.
 * </para>
 */

/**
 * \section writing_functions_interruption_handling Writing \a igraph functions with
 * proper interruption handling
 *
 * <para>
 * There is practically a simple rule that should be obeyed when writing
 * \a igraph functions. If the calculation is expected to take a long time
 * in large graphs (a simple rule of thumb is to assume this for every
 * function with a time complexity of at least O(n^2)), call
 * \ref IGRAPH_ALLOW_INTERRUPTION in regular intervals like every 10th
 * iteration or so.
 * </para>
 */

/**
 * \typedef igraph_interruption_handler_t
 *
 * This is the type of the interruption handler functions.
 *
 * \param data reserved for possible future use
 * \return \c IGRAPH_SUCCESS if the calculation should go on, anything else otherwise.
 */

typedef int igraph_interruption_handler_t (void* data);

/**
 * \function igraph_allow_interruption
 *
 * This is the function which is called (usually via the
 * \ref IGRAPH_ALLOW_INTERRUPTION macro) if \a igraph is checking for interruption
 * requests.
 *
 * \param data reserved for possible future use, now it is always \c NULL
 * \return \c IGRAPH_SUCCESS if the calculation should go on, anything else otherwise.
 */

IGRAPH_EXPORT int igraph_allow_interruption(void* data);

IGRAPH_EXPORT igraph_interruption_handler_t * igraph_set_interruption_handler (igraph_interruption_handler_t * new_handler);

__END_DECLS

#endif
