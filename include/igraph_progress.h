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

#ifndef IGRAPH_PROGRESS_H
#define IGRAPH_PROGRESS_H

#include "igraph_decls.h"
#include "igraph_types.h"

__BEGIN_DECLS

/**
 * \section about_progress_handlers About progress handlers
 *
 * <para>It is often useful to report the progress of some long
 * calculation, to allow the user to follow the computation and
 * guess the total running time. A couple of igraph functions
 * support this at the time of writing, hopefully more will support it
 * in the future.
 * </para>
 *
 * <para>
 * To see the progress of a computation, the user has to install a
 * progress handler, as there is none installed by default.
 * If an igraph function supports progress reporting, then it
 * calls the installed progress handler periodically, and passes a
 * percentage value to it, the percentage of computation already
 * performed. To install a progress handler, you need to call
 * \ref igraph_set_progress_handler(). Currently there is a single
 * pre-defined progress handler, called \ref
 * igraph_progress_handler_stderr().
 * </para>
 */

/**
 * \section writing_progress_handlers Writing progress handlers
 *
 * <para>
 * To write a new progress handler, one needs to create a function of
 * type \ref igraph_progress_handler_t. The new progress handler
 * can then be installed with the \ref igraph_set_progress_handler()
 * function.
 * </para>
 *
 * <para>
 * One can assume that the first progress handler call from a
 * calculation will be call with zero as the \p percentage argument,
 * and the last call from a function will have 100 as the \p
 * percentage argument. Note, however, that if an error happens in the
 * middle of a computation, then the 100 percent call might be
 * omitted.
 * </para>
 */

/**
 * \section igraph_functions_with_progress Writing igraph functions with progress reporting
 *
 * <para>
 * If you want to write a function that uses igraph and supports
 * progress reporting, you need to include \ref igraph_progress()
 * calls in your function, usually via the \ref IGRAPH_PROGRESS()
 * macro.
 * </para>
 *
 * <para>
 * It is good practice to always include a call to \ref
 * igraph_progress() with a zero \p percentage argument, before the
 * computation; and another call with 100 \p percentage value
 * after the computation is completed.
 * </para>
 *
 * <para>
 * It is also good practice \em not to call \ref igraph_progress() too
 * often, as this would slow down the computation. It might not be
 * worth to support progress reporting in functions with linear or
 * log-linear time complexity, as these are fast, even with a large
 * amount of data. For functions with quadratic or higher time
 * complexity make sure that the time complexity of the progress
 * reporting is constant or at least linear. In practice this means
 * having at most O(n) progress checks and at most 100
 * \ref igraph_progress() calls.
 * </para>
 */

/**
 * \section progress_and_threads Multi-threaded programs
 *
 * <para>
 * In multi-threaded programs, each thread has its own progress
 * handler, if thread-local storage is supported and igraph is
 * thread-safe. See the \ref IGRAPH_THREAD_SAFE macro for checking
 * whether an igraph build is thread-safe.
 * </para>
 */

/* -------------------------------------------------- */
/* Progress handlers                                  */
/* -------------------------------------------------- */

/**
 * \typedef igraph_progress_handler_t
 * \brief Type of progress handler functions
 *
 * This is the type of the igraph progress handler functions.
 * There is currently one such predefined function,
 * \ref igraph_progress_handler_stderr(), but the user can
 * write and set up more sophisticated ones.
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
 * \return If the return value of the progress handler is not
 *     \c IGRAPH_SUCCESS, then \ref igraph_progress() returns the
 *     error code \c IGRAPH_INTERRUPTED. The \ref IGRAPH_PROGRESS()
 *     macro frees all memory and finishes the igraph function with
 *     error code \c IGRAPH_INTERRUPTED in this case.
 */

typedef int igraph_progress_handler_t(const char *message, igraph_real_t percent,
                                      void *data);

IGRAPH_EXPORT extern igraph_progress_handler_t igraph_progress_handler_stderr;

IGRAPH_EXPORT igraph_progress_handler_t * igraph_set_progress_handler(igraph_progress_handler_t new_handler);

IGRAPH_EXPORT int igraph_progress(const char *message, igraph_real_t percent, void *data);

IGRAPH_EXPORT int igraph_progressf(const char *message, igraph_real_t percent, void *data,
                                   ...);

/**
 * \define IGRAPH_PROGRESS
 * \brief Report progress.
 *
 * The standard way to report progress from an igraph function
 * \param message A string, a textual message that references the
 *    calculation under progress.
 * \param percent Numeric scalar, the percentage that is complete.
 * \param data User-defined data, this can be used in user-defined
 *    progress handler functions, from user-written igraph functions.
 * \return If the progress handler returns with \c IGRAPH_INTERRUPTED,
 *    then this macro frees up the igraph allocated memory for
 *    temporary data and returns to the caller with \c
 *    IGRAPH_INTERRUPTED.
 */

#define IGRAPH_PROGRESS(message, percent, data) \
    do { \
        if (igraph_progress((message), (percent), (data)) != IGRAPH_SUCCESS) { \
            IGRAPH_FINALLY_FREE(); \
            return IGRAPH_INTERRUPTED; \
        } \
    } while (0)

__END_DECLS

#endif
