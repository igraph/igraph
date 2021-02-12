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

#ifndef IGRAPH_ERROR_H
#define IGRAPH_ERROR_H

#include "igraph_decls.h"

#include <stdarg.h>

__BEGIN_DECLS

/* This file contains the igraph error handling.
 * Most bits are taken literally from the GSL library (with the GSL_
 * prefix renamed to IGRAPH_), as I couldn't find a better way to do
 * them. */

/* IGRAPH_NORETURN indicates to the compiler that a function does not return.
 * There are standard facilities for this, namely _Noreturn in C11 and [[noreturn]] in C++11.
 * However, since igraph is currently compiled with older standards, and since
 * the standard 'noreturn' specification would need to be diferent between C and C++,
 * we do not use these facilities.
 */
#if defined(__GNUC__)
/* Compilers that support the GNU C syntax. Use __noreturn__ instead of 'noreturn' as the latter is a macro in C11. */
#define IGRAPH_NORETURN __attribute__((__noreturn__))
#elif defined(_MSC_VER)
/* Compilers that support the MSVC syntax. */
#define IGRAPH_NORETURN __declspec(noreturn)
#else
#define IGRAPH_NORETURN
#endif

/**
 * \section error_handling_basics Error handling basics
 *
 * <para>\a igraph functions can run into various problems preventing them
 * from normal operation. The user might have supplied invalid arguments,
 * e.g. a non-square matrix when a square-matrix was expected, or the program
 * has run out of memory while some more memory allocation is required, etc.
 * </para>
 *
 * <para>By default \a igraph aborts the program when it runs into an
 * error. While this behavior might be good enough for smaller programs,
 * it is without doubt avoidable in larger projects. Please read further
 * if your project requires more sophisticated error handling. You can
 * safely skip the rest of this chapter otherwise.
 * </para>
 */

/**
 * \section error_handlers Error handlers
 *
 * <para>
 * If \a igraph runs into an error - an invalid argument was supplied
 * to a function, or we've ran out of memory - the control is
 * transferred to the \emb error handler \eme function.
 * </para><para>
 * The default error handler is \ref igraph_error_handler_abort which
 * prints an error message and aborts the program.
 * </para>
 * <para>
 * The \ref igraph_set_error_handler() function can be used to set a new
 * error handler function of type \ref igraph_error_handler_t; see the
 * documentation of this type for details.
 * </para>
 * <para>
 * There are two other predefined error handler functions,
 * \ref igraph_error_handler_ignore and \ref igraph_error_handler_printignore.
 * These deallocate the temporarily allocated memory (more about this
 * later) and return with the error code. The latter also prints an
 * error message. If you use these error handlers you need to take
 * care about possible errors yourself by checking the return value of
 * (almost) every non-void \a igraph function.
 * </para><para>
 * Independently of the error handler installed, all functions in the
 * library do their best to leave their arguments
 * \em semantically unchanged if an error
 * happens. By semantically we mean that the implementation of an
 * object supplied as an argument might change, but its
 * \quote meaning \endquote in most cases does not. The rare occasions
 * when this rule is violated are documented in this manual.
 * </para>
 */

/**
 * \section error_codes Error codes
 *
 * <para>Every \a igraph function which can fail return a
 * single integer error code. Some functions are very simple and
 * cannot run into any error, these may return other types, or
 * \type void as well. The error codes are defined by the
 * \ref igraph_error_type_t enumeration.
 * </para>
 */

/**
 * \section writing_error_handlers Writing error handlers
 *
 * <para>
 * The contents of the rest of this chapter might be useful only
 * for those who want to create an interface to \a igraph from another
 * language. Most readers can safely skip to the next chapter.
 * </para>
 *
 * <para>
 * You can write and install error handlers simply by defining a
 * function of type \ref igraph_error_handler_t and calling
 * \ref igraph_set_error_handler(). This feature is useful for interface
 * writers, as \a igraph will have the chance to
 * signal errors the appropriate way, e.g. the R interface defines an
 * error handler which calls the <function>error()</function>
 * function, as required by R, while the Python interface has an error
 * handler which raises an exception according to the Python way.
 * </para>
 * <para>
 * If you want to write an error handler, your error handler should
 * call \ref IGRAPH_FINALLY_FREE() to deallocate all temporary memory to
 * prevent memory leaks.
 * </para>
 */

/**
 * \section error_handling_internals Error handling internals
 *
 * <para>
 * If an error happens, the functions in the library call the
 * \ref IGRAPH_ERROR() macro with a textual description of the error and an
 * \a igraph error code. This macro calls (through the \ref
 * igraph_error() function) the installed error handler. Another useful
 * macro is \ref IGRAPH_CHECK(). This checks the return value of its
 * argument, which is normally a function call, and calls \ref
 * IGRAPH_ERROR() if it is not \c IGRAPH_SUCCESS.
 * </para>
 */

/**
 * \section deallocating_memory Deallocating memory
 *
 * <para>
 * If a function runs into an error (and the program is not aborted)
 * the error handler should deallocate all temporary memory. This is
 * done by storing the address and the destroy function of all temporary
 * objects in a stack. The \ref IGRAPH_FINALLY function declares an object as
 * temporary by placing its address in the stack. If an \a igraph function returns
 * with success it calls \ref IGRAPH_FINALLY_CLEAN() with the
 * number of objects to remove from the stack. If an error happens
 * however, the error handler should call \ref IGRAPH_FINALLY_FREE() to
 * deallocate each object added to the stack. This means that the
 * temporary objects allocated in the calling function (and etc.) will
 * be freed as well.
 * </para>
 */

/**
 * \section writing_functions_error_handling Writing \a igraph functions with
 * proper error handling
 *
 * <para>
 * There are some simple rules to keep in order to have functions
 * behaving well in erroneous situations. First, check the arguments
 * of the functions and call \ref IGRAPH_ERROR() if they are invalid. Second,
 * call \ref IGRAPH_FINALLY on each dynamically allocated object and call
 * \ref IGRAPH_FINALLY_CLEAN() with the proper argument before returning. Third, use
 * \ref IGRAPH_CHECK on all \a igraph function calls which can generate errors.
 * </para>
 * <para>
 * The size of the stack used for this bookkeeping is fixed, and
 * small. If you want to allocate several objects, write a destroy
 * function which can deallocate all of these. See the
 * <filename>adjlist.c</filename> file in the
 * \a igraph source for an example.
 * </para>
 * <para>
 * For some functions these mechanisms are simply not flexible
 * enough. These functions should define their own error handlers and
 * restore the error handler before they return.
 * </para>
 */

/**
 * \section error_handling_threads Error handling and threads
 *
 * <para>
 * It is likely that the \a igraph error handling
 * method is \em not thread-safe, mainly because of
 * the static global stack which is used to store the address of the
 * temporarily allocated objects. This issue might be addressed in a
 * later version of \a igraph.
 * </para>
 */

/**
 * \typedef igraph_error_handler_t
 * \brief The type of error handler functions.
 *
 * This is the type of the error handler functions.
 * \param reason Textual description of the error.
 * \param file The source file in which the error is noticed.
 * \param line The number of the line in the source file which triggered
 *   the error
 * \param igraph_errno The \a igraph error code.
 */

typedef void igraph_error_handler_t (const char * reason, const char * file,
                                     int line, int igraph_errno);

/**
 * \var igraph_error_handler_abort
 * \brief Abort program in case of error.
 *
 * The default error handler, prints an error message and aborts the
 * program.
 */

IGRAPH_EXPORT igraph_error_handler_t igraph_error_handler_abort;

/**
 * \var igraph_error_handler_ignore
 * \brief Ignore errors.
 *
 * This error handler frees the temporarily allocated memory and returns
 * with the error code.
 */

IGRAPH_EXPORT igraph_error_handler_t igraph_error_handler_ignore;

/**
 * \var igraph_error_handler_printignore
 * \brief Print and ignore errors.
 *
 * Frees temporarily allocated memory, prints an error message to the
 * standard error and returns with the error code.
 */

IGRAPH_EXPORT igraph_error_handler_t igraph_error_handler_printignore;

/**
 * \function igraph_set_error_handler
 * \brief Sets a new error handler.
 *
 * Installs a new error handler. If called with 0, it installs the
 * default error handler (which is currently
 * \ref igraph_error_handler_abort).
 * \param new_handler The error handler function to install.
 * \return The old error handler function. This should be saved and
 *   restored if \p new_handler is not needed any
 *   more.
 */

IGRAPH_EXPORT igraph_error_handler_t* igraph_set_error_handler(igraph_error_handler_t* new_handler);

/**
 * \typedef igraph_error_type_t
 * \brief Error code type.
 * These are the possible values returned by \a igraph functions.
 * Note that these are interesting only if you defined an error handler
 * with \ref igraph_set_error_handler(). Otherwise the program is aborted
 * and the function causing the error never returns.
 *
 * \enumval IGRAPH_SUCCESS The function successfully completed its task.
 * \enumval IGRAPH_FAILURE Something went wrong. You'll almost never
 *    meet this error as normally more specific error codes are used.
 * \enumval IGRAPH_ENOMEM There wasn't enough memory to allocate
 *    on the heap.
 * \enumval IGRAPH_PARSEERROR A parse error was found in a file.
 * \enumval IGRAPH_EINVAL A parameter's value is invalid. E.g. negative
 *    number was specified as the number of vertices.
 * \enumval IGRAPH_EXISTS A graph/vertex/edge attribute is already
 *    installed with the given name.
 * \enumval IGRAPH_EINVEVECTOR Invalid vector of vertex ids. A vertex id
 *    is either negative or bigger than the number of vertices minus one.
 * \enumval IGRAPH_EINVVID Invalid vertex id, negative or too big.
 * \enumval IGRAPH_NONSQUARE A non-square matrix was received while a
 *    square matrix was expected.
 * \enumval IGRAPH_EINVMODE Invalid mode parameter.
 * \enumval IGRAPH_EFILE A file operation failed. E.g. a file doesn't exist,
 *   or the user has no rights to open it.
 * \enumval IGRAPH_UNIMPLEMENTED Attempted to call an unimplemented or
 *   disabled (at compile-time) function.
 * \enumval IGRAPH_DIVERGED A numeric algorithm failed to converge.
 * \enumval IGRAPH_ARPACK_PROD Matrix-vector product failed.
 * \enumval IGRAPH_ARPACK_NPOS N must be positive.
 * \enumval IGRAPH_ARPACK_NEVNPOS NEV must be positive.
 * \enumval IGRAPH_ARPACK_NCVSMALL NCV must be bigger.
 * \enumval IGRAPH_ARPACK_NONPOSI Maximum number of iterations should be positive.
 * \enumval IGRAPH_ARPACK_WHICHINV Invalid WHICH parameter.
 * \enumval IGRAPH_ARPACK_BMATINV Invalid BMAT parameter.
 * \enumval IGRAPH_ARPACK_WORKLSMALL WORKL is too small.
 * \enumval IGRAPH_ARPACK_TRIDERR LAPACK error in tridiagonal eigenvalue calculation.
 * \enumval IGRAPH_ARPACK_ZEROSTART Starting vector is zero.
 * \enumval IGRAPH_ARPACK_MODEINV MODE is invalid.
 * \enumval IGRAPH_ARPACK_MODEBMAT MODE and BMAT are not compatible.
 * \enumval IGRAPH_ARPACK_ISHIFT ISHIFT must be 0 or 1.
 * \enumval IGRAPH_ARPACK_NEVBE NEV and WHICH='BE' are incompatible.
 * \enumval IGRAPH_ARPACK_NOFACT Could not build an Arnoldi factorization.
 * \enumval IGRAPH_ARPACK_FAILED No eigenvalues to sufficient accuracy.
 * \enumval IGRAPH_ARPACK_HOWMNY HOWMNY is invalid.
 * \enumval IGRAPH_ARPACK_HOWMNYS HOWMNY='S' is not implemented.
 * \enumval IGRAPH_ARPACK_EVDIFF Different number of converged Ritz values.
 * \enumval IGRAPH_ARPACK_SHUR Error from calculation of a real Schur form.
 * \enumval IGRAPH_ARPACK_LAPACK LAPACK (dtrevc) error for calculating eigenvectors.
 * \enumval IGRAPH_ARPACK_UNKNOWN Unknown ARPACK error.
 * \enumval IGRAPH_ENEGLOOP Negative loop detected while calculating shortest paths.
 * \enumval IGRAPH_EINTERNAL Internal error, likely a bug in igraph.
 * \enumval IGRAPH_EDIVZERO Big integer division by zero.
 * \enumval IGARPH_GLP_EBOUND GLPK error (GLP_EBOUND).
 * \enumval IGARPH_GLP_EROOT GLPK error (GLP_EROOT).
 * \enumval IGARPH_GLP_ENOPFS GLPK error (GLP_ENOPFS).
 * \enumval IGARPH_GLP_ENODFS GLPK error (GLP_ENODFS).
 * \enumval IGARPH_GLP_EFAIL GLPK error (GLP_EFAIL).
 * \enumval IGARPH_GLP_EMIPGAP GLPK error (GLP_EMIPGAP).
 * \enumval IGARPH_GLP_ETMLIM GLPK error (GLP_ETMLIM).
 * \enumval IGARPH_GLP_ESTOP GLPK error (GLP_ESTOP).
 * \enumval IGRAPH_EATTRIBUTES Attribute handler error. The user is not
 *   expected to find this; it is signalled if some igraph function is
 *   not using the attribute handler interface properly.
 * \enumval IGRAPH_EATTRCOMBINE Unimplemented attribute combination
 *   method for the given attribute type.
 * \enumval IGRAPH_ELAPACK A LAPACK call resulted an error.
 * \enumval IGRAPH_EDRL Internal error in the DrL layout generator.
 * \enumval IGRAPH_EOVERFLOW Integer or double overflow.
 * \enumval IGRAPH_EGLP Internal GLPK error.
 * \enumval IGRAPH_CPUTIME CPU time exceeded.
 * \enumval IGRAPH_EUNDERFLOW Integer or double underflow.
 * \enumval IGRAPH_ERWSTUCK Random walk got stuck.
 */

typedef enum {
    IGRAPH_SUCCESS           = 0,
    IGRAPH_FAILURE           = 1,
    IGRAPH_ENOMEM            = 2,
    IGRAPH_PARSEERROR        = 3,
    IGRAPH_EINVAL            = 4,
    IGRAPH_EXISTS            = 5,
    IGRAPH_EINVEVECTOR       = 6,
    IGRAPH_EINVVID           = 7,
    IGRAPH_NONSQUARE         = 8,
    IGRAPH_EINVMODE          = 9,
    IGRAPH_EFILE             = 10,
    IGRAPH_UNIMPLEMENTED     = 12,
    IGRAPH_INTERRUPTED       = 13,
    IGRAPH_DIVERGED          = 14,
    IGRAPH_ARPACK_PROD       = 15,
    IGRAPH_ARPACK_NPOS       = 16,
    IGRAPH_ARPACK_NEVNPOS    = 17,
    IGRAPH_ARPACK_NCVSMALL   = 18,
    IGRAPH_ARPACK_NONPOSI    = 19,
    IGRAPH_ARPACK_WHICHINV   = 20,
    IGRAPH_ARPACK_BMATINV    = 21,
    IGRAPH_ARPACK_WORKLSMALL = 22,
    IGRAPH_ARPACK_TRIDERR    = 23,
    IGRAPH_ARPACK_ZEROSTART  = 24,
    IGRAPH_ARPACK_MODEINV    = 25,
    IGRAPH_ARPACK_MODEBMAT   = 26,
    IGRAPH_ARPACK_ISHIFT     = 27,
    IGRAPH_ARPACK_NEVBE      = 28,
    IGRAPH_ARPACK_NOFACT     = 29,
    IGRAPH_ARPACK_FAILED     = 30,
    IGRAPH_ARPACK_HOWMNY     = 31,
    IGRAPH_ARPACK_HOWMNYS    = 32,
    IGRAPH_ARPACK_EVDIFF     = 33,
    IGRAPH_ARPACK_SHUR       = 34,
    IGRAPH_ARPACK_LAPACK     = 35,
    IGRAPH_ARPACK_UNKNOWN    = 36,
    IGRAPH_ENEGLOOP          = 37,
    IGRAPH_EINTERNAL         = 38,
    IGRAPH_ARPACK_MAXIT      = 39,
    IGRAPH_ARPACK_NOSHIFT    = 40,
    IGRAPH_ARPACK_REORDER    = 41,
    IGRAPH_EDIVZERO          = 42,
    IGRAPH_GLP_EBOUND        = 43,
    IGRAPH_GLP_EROOT         = 44,
    IGRAPH_GLP_ENOPFS        = 45,
    IGRAPH_GLP_ENODFS        = 46,
    IGRAPH_GLP_EFAIL         = 47,
    IGRAPH_GLP_EMIPGAP       = 48,
    IGRAPH_GLP_ETMLIM        = 49,
    IGRAPH_GLP_ESTOP         = 50,
    IGRAPH_EATTRIBUTES       = 51,
    IGRAPH_EATTRCOMBINE      = 52,
    IGRAPH_ELAPACK           = 53,
    IGRAPH_EDRL              = 54,
    IGRAPH_EOVERFLOW         = 55,
    IGRAPH_EGLP              = 56,
    IGRAPH_CPUTIME           = 57,
    IGRAPH_EUNDERFLOW        = 58,
    IGRAPH_ERWSTUCK          = 59,
    IGRAPH_STOP              = 60  /* undocumented, used internally; signals a request to stop in functions like igraph_i_maximal_cliques_bk */
} igraph_error_type_t;
/* Each enum value above must have a corresponding error string in
 * igraph_i_error_strings[] in igraph_error.c */

/* We use IGRAPH_FILE_BASENAME instead of __FILE__ to ensure that full
 * paths don't leak into the library code. IGRAPH_FILE_BASENAME is set up
 * by the build system when compiling the individual files. However, when
 * including igraph_error.h in user code, this macro is not defined so we
 * fall back to __FILE__ here
 */
#ifndef IGRAPH_FILE_BASENAME
#  define IGRAPH_FILE_BASENAME __FILE__
#endif

/**
 * \define IGRAPH_ERROR
 * \brief Trigger an error.
 *
 * \a igraph functions usually use this macro when they notice an error.
 * It calls
 * \ref igraph_error() with the proper parameters and if that returns
 * the macro returns the "calling" function as well, with the error
 * code. If for some (suspicious) reason you want to call the error
 * handler without returning from the current function, call
 * \ref igraph_error() directly.
 * \param reason Textual description of the error. This should be
 *   something more descriptive than the text associated with the error
 *   code. E.g. if the error code is \c IGRAPH_EINVAL,
 *   its associated text (see  \ref igraph_strerror()) is "Invalid
 *   value" and this string should explain which parameter was invalid
 *   and maybe why.
 * \param igraph_errno The \a igraph error code.
 */

#define IGRAPH_ERROR(reason, igraph_errno) \
    do { \
        igraph_error (reason, IGRAPH_FILE_BASENAME, __LINE__, igraph_errno) ; \
        return igraph_errno ; \
    } while (0)

/**
 * \function igraph_error
 * \brief Triggers an error.
 *
 * \a igraph functions usually call this function (most often via the
 * \ref IGRAPH_ERROR macro) if they notice an error.
 * It calls the currently installed error handler function with the
 * supplied arguments.
 *
 * \param reason Textual description of the error.
 * \param file The source file in which the error was noticed.
 * \param line The number of line in the source file which triggered the
 *   error.
 * \param igraph_errno The \a igraph error code.
 * \return the error code (if it returns)
 *
 * \sa igraph_errorf().
 */

IGRAPH_EXPORT int igraph_error(const char *reason, const char *file, int line,
                               int igraph_errno);

/**
 * \define IGRAPH_ERRORF
 * \brief Triggers an error, with printf-like syntax.
 *
 * \a igraph functions can use this macro when they notice an error and
 * want to pass on extra information to the user about what went wrong.
 * It calls \ref igraph_errorf() with the proper parameters and if that
 * returns the macro returns the "calling" function as well, with the
 * error code. If for some (suspicious) reason you want to call the
 * error handler without returning from the current function, call
 * \ref igraph_errorf() directly.
 * \param reason Textual description of the error, a template string
 *   with the same syntax as the standard printf C library function.
 *   This should be something more descriptive than the text associated
 *   with the error code. E.g. if the error code is \c IGRAPH_EINVAL,
 *   its associated text (see  \ref igraph_strerror()) is "Invalid
 *   value" and this string should explain which parameter was invalid
 *   and maybe what was expected and what was recieved.
 * \param igraph_errno The \a igraph error code.
 * \param ... The additional arguments to be substituted into the
 *   template string.
 */

#define IGRAPH_ERRORF(reason, igraph_errno, ...) \
    do { \
        igraph_errorf(reason, IGRAPH_FILE_BASENAME, __LINE__, \
                      igraph_errno, __VA_ARGS__) ; \
        return igraph_errno; \
    } while (0)

/**
 * \function igraph_errorf
 * \brief Triggers an error, printf-like version.
 *
 * \param reason Textual description of the error, interpreted as
 *               a \c printf format string.
 * \param file The source file in which the error was noticed.
 * \param line The line in the source file which triggered the error.
 * \param igraph_errno The \a igraph error code.
 * \param ... Additional parameters, the values to substitute into the
 *            format string.
 *
 * \sa igraph_error().
 */

IGRAPH_EXPORT int igraph_errorf(const char *reason, const char *file, int line,
                                int igraph_errno, ...);

IGRAPH_EXPORT int igraph_errorvf(const char *reason, const char *file, int line,
                                 int igraph_errno, va_list ap);

/**
 * \function igraph_strerror
 * \brief Textual description of an error.
 *
 * This is a simple utility function, it gives a short general textual
 * description for an \a igraph error code.
 *
 * \param igraph_errno The \a igraph error code.
 * \return pointer to the textual description of the error code.
 */

IGRAPH_EXPORT const char* igraph_strerror(const int igraph_errno);

#define IGRAPH_ERROR_SELECT_2(a,b)       ((a) != IGRAPH_SUCCESS ? (a) : ((b) != IGRAPH_SUCCESS ? (b) : IGRAPH_SUCCESS))
#define IGRAPH_ERROR_SELECT_3(a,b,c)     ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_2(b,c))
#define IGRAPH_ERROR_SELECT_4(a,b,c,d)   ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_3(b,c,d))
#define IGRAPH_ERROR_SELECT_5(a,b,c,d,e) ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_4(b,c,d,e))

/* Now comes the more convenient error handling macro arsenal.
 * Ideas taken from exception.{h,c} by Laurent Deniau see
 * http://cern.ch/Laurent.Deniau/html/oopc/oopc.html#Exceptions for more
 * information. We don't use the exception handling code though.  */

struct igraph_i_protectedPtr {
    int all;
    void *ptr;
    void (*func)(void*);
};

typedef void igraph_finally_func_t (void*);

IGRAPH_EXPORT void IGRAPH_FINALLY_REAL(void (*func)(void*), void* ptr);

/**
 * \function IGRAPH_FINALLY_CLEAN
 * \brief Signals clean deallocation of objects.
 *
 * Removes the specified number of objects from the stack of
 * temporarily allocated objects. Most often this is called just
 * before returning from a function.
 * \param num The number of objects to remove from the bookkeeping
 *   stack.
 */

IGRAPH_EXPORT void IGRAPH_FINALLY_CLEAN(int num);

/**
 * \function IGRAPH_FINALLY_FREE
 * \brief Deallocates all registered objects.
 *
 * Calls the destroy function for all objects in the stack of
 * temporarily allocated objects. This is usually called only from an
 * error handler. It is \em not appropriate to use it
 * instead of destroying each unneeded object of a function, as it
 * destroys the temporary objects of the caller function (and so on)
 * as well.
 */

IGRAPH_EXPORT void IGRAPH_FINALLY_FREE(void);

/**
 * \function IGRAPH_FINALLY_STACK_SIZE
 * \brief The number of registered objects.
 *
 * Returns the number of objects in the stack of temporarily allocated
 * objects. This function is handy if you write an own igraph routine and
 * you want to make sure it handles errors properly. A properly written
 * igraph routine should not leave pointers to temporarily allocated objects
 * in the finally stack, because otherwise an \ref IGRAPH_FINALLY_FREE call
 * in another igraph function would result in freeing these objects as well
 * (and this is really hard to debug, since the error will be not in that
 * function that shows erroneous behaviour). Therefore, it is advised to
 * write your own test cases and examine \ref IGRAPH_FINALLY_STACK_SIZE
 * before and after your test cases - the numbers should be equal.
 */
IGRAPH_EXPORT int IGRAPH_FINALLY_STACK_SIZE(void);

/**
 * \define IGRAPH_FINALLY_STACK_EMPTY
 * \brief Returns true if there are no registered objects, false otherwise.
 *
 * This is just a shorthand notation for checking that
 * \ref IGRAPH_FINALLY_STACK_SIZE() is zero.
 */
#define IGRAPH_FINALLY_STACK_EMPTY (IGRAPH_FINALLY_STACK_SIZE() == 0)

/**
 * \define IGRAPH_FINALLY
 * \brief Registers an object for deallocation.
 * \param func The address of the function which is normally called to
 *   destroy the object.
 * \param ptr Pointer to the object itself.
 *
 * This macro places the address of an object, together with the
 * address of its destructor in a stack. This stack is used if an
 * error happens to deallocate temporarily allocated objects to
 * prevent memory leaks.
 */

#define IGRAPH_FINALLY(func, ptr) \
    do { \
        /* the following branch makes the compiler check the compatibility of \
         * func and ptr to detect cases when we are accidentally invoking an \
         * incorrect destructor function with the pointer */ \
        if (0) { func(ptr); } \
        IGRAPH_FINALLY_REAL((igraph_finally_func_t*)(func), (ptr)); \
    } while (0)

#if !defined(GCC_VERSION_MAJOR) && defined(__GNUC__)
    #define GCC_VERSION_MAJOR  __GNUC__
#endif

#if defined(GCC_VERSION_MAJOR) && (GCC_VERSION_MAJOR >= 3)
    #define IGRAPH_UNLIKELY(a) __builtin_expect((a), 0)
    #define IGRAPH_LIKELY(a)   __builtin_expect((a), 1)
#else
    #define IGRAPH_UNLIKELY(a) a
    #define IGRAPH_LIKELY(a)   a
#endif

#if IGRAPH_VERIFY_FINALLY_STACK == 1
#define IGRAPH_CHECK(a) \
        do { \
            int enter_stack_size = IGRAPH_FINALLY_STACK_SIZE(); \
            int igraph_i_ret=(a); \
            if (IGRAPH_UNLIKELY(igraph_i_ret != 0)) {\
                IGRAPH_ERROR("", igraph_i_ret); \
            } \
            if (IGRAPH_UNLIKELY(enter_stack_size != IGRAPH_FINALLY_STACK_SIZE())) { \
                IGRAPH_ERROR("Non-matching number of IGRAPH_FINALLY and IGRAPH_FINALLY_CLEAN", IGRAPH_FAILURE); \
            } \
        } while (0)
#else
/**
 * \define IGRAPH_CHECK
 * \brief Checks the return value of a function call.
 *
 * \param a An expression, usually a function call.
 *
 * Executes the expression and checks its value. If this is not
 * \c IGRAPH_SUCCESS, it calls \ref IGRAPH_ERROR with
 * the value as the error code. Here is an example usage:
 * \verbatim IGRAPH_CHECK(vector_push_back(&amp;v, 100)); \endverbatim
 *
 * </para><para>There is only one reason to use this macro when writing
 * \a igraph functions. If the user installs an error handler which
 * returns to the auxiliary calling code (like \ref
 * igraph_error_handler_ignore and \ref
 * igraph_error_handler_printignore), and the \a igraph function
 * signalling the error is called from another \a igraph function
 * then we need to make sure that the error is propagated back to
 * the auxiliary (i.e. non-igraph) calling function. This is achieved
 * by using <function>IGRAPH_CHECK</function> on every \a igraph
 * call which can return an error code.
 */
#define IGRAPH_CHECK(a) do { \
        int igraph_i_ret=(a); \
        if (IGRAPH_UNLIKELY(igraph_i_ret != 0)) {\
            IGRAPH_ERROR("", igraph_i_ret); \
        } } while (0)
#endif



/**
 * \section about_igraph_warnings Warning messages
 *
 * <para>
 * \a igraph also supports warning messages in addition to error
 * messages. Warning messages typically do not terminate the
 * program, but they are usually crucial to the user.
 * </para>
 *
 * <para>
 * \a igraph warnings are handled similarly to errors. There is a
 * separate warning handler function that is called whenever
 * an \a igraph function triggers a warning. This handler can be
 * set by the \ref igraph_set_warning_handler() function. There are
 * two predefined simple warning handlers,
 * \ref igraph_warning_handler_ignore() and
 * \ref igraph_warning_handler_print(), the latter being the default.
 * </para>
 *
 * <para>
 * To trigger a warning, \a igraph functions typically use the
 * \ref IGRAPH_WARNING() macro, the \ref igraph_warning() function,
 * or if more flexibility is needed, \ref igraph_warningf().
 * </para>
 */

/**
 * \typedef igraph_warning_handler_t
 * \brief The type of igraph warning handler functions.
 *
 * Currently it is defined to have the same type as
 * \ref igraph_error_handler_t, although the last (error code)
 * argument is not used.
 */

typedef igraph_error_handler_t igraph_warning_handler_t;

/**
 * \function igraph_set_warning_handler
 * \brief Installs a warning handler.
 *
 * Install the supplied warning handler function.
 * \param new_handler The new warning handler function to install.
 *        Supply a null pointer here to uninstall the current
 *        warning handler, without installing a new one.
 * \return The current warning handler function.
 */

IGRAPH_EXPORT igraph_warning_handler_t* igraph_set_warning_handler(igraph_warning_handler_t* new_handler);

IGRAPH_EXPORT extern igraph_warning_handler_t igraph_warning_handler_ignore;
IGRAPH_EXPORT extern igraph_warning_handler_t igraph_warning_handler_print;

/**
 * \function igraph_warning
 * \brief Triggers a warning.
 *
 * Call this function if you want to trigger a warning from within
 * a function that uses \a igraph.
 * \param reason Textual description of the warning.
 * \param file The source file in which the warning was noticed.
 * \param line The number of line in the source file which triggered the
 *         warning.
 * \param igraph_errno Warnings could have potentially error codes as well,
 *        but this is currently not used in igraph.
 * \return The supplied error code.
 */

IGRAPH_EXPORT int igraph_warning(const char *reason, const char *file, int line,
                                 int igraph_errno);

/**
 * \define IGRAPH_WARNINGF
 * \brief Triggers a warning, with printf-like syntax.
 *
 * \a igraph functions can use this macro when they notice a warning and
 * want to pass on extra information to the user about what went wrong.
 * It calls \ref igraph_warningf() with the proper parameters and no
 * error code.
 * \param reason Textual description of the warning, a template string
 *        with the same syntax as the standard printf C library function.
 * \param ... The additional arguments to be substituted into the
 *        template string.
 */

#define IGRAPH_WARNINGF(reason, ...) \
    do { \
        igraph_warningf(reason, IGRAPH_FILE_BASENAME, __LINE__, \
                        -1, __VA_ARGS__); \
    } while (0)



/**
 * \function igraph_warningf
 * \brief Triggers a warning, printf-like version.
 *
 * This function is similar to \ref igraph_warning(), but
 * uses a printf-like syntax. It substitutes the additional arguments
 * into the \p reason template string and calls \ref igraph_warning().
 * \param reason Textual description of the warning, a template string
 *        with the same syntax as the standard printf C library function.
 * \param file The source file in which the warning was noticed.
 * \param line The number of line in the source file which triggered the
 *         warning.
 * \param igraph_errno Warnings could have potentially error codes as well,
 *        but this is currently not used in igraph.
 * \param ... The additional arguments to be substituted into the
 *        template string.
 * \return The supplied error code.
 */

IGRAPH_EXPORT int igraph_warningf(const char *reason, const char *file, int line,
                                  int igraph_errno, ...);

/**
 * \define IGRAPH_WARNING
 * \brief Triggers a warning.
 *
 * This is the usual way of triggering a warning from an igraph
 * function. It calls \ref igraph_warning().
 * \param reason The warning message.
 */

#define IGRAPH_WARNING(reason) \
    do { \
        igraph_warning(reason, IGRAPH_FILE_BASENAME, __LINE__, -1); \
    } while (0)


/**
 * \section fatal_error_handlers Fatal errors
 *
 * <para>
 * In some rare situations, \a igraph may encounter an internal error
 * that cannot be fully handled. In this case, it will call the
 * current fatal error handler. The default fatal error handler
 * simply prints the error and aborts the program.
 * </para>
 *
 * <para>
 * Fatal error handlers do not return. Typically, they might abort the
 * the program immediately, or in the case of the high-level \a igraph
 * interfaces, they might return to the top level using a
 * <code>longjmp()</code>. The fatal error handler is only called when
 * a serious error has occurred, and as a result igraph may be in an
 * inconsistent state. The purpose of returning to the top level is to
 * give the user a chance to save their work instead of aborting immediately.
 * However, the program session should be restarted as soon as possible.
 * </para>
 *
 * <para>
 * Most projects that use \a igraph will use the default fatal error
 * handler.
 * </para>
 */

/**
 * \typedef igraph_fatal_handler_t
 * \brief The type of igraph fatal error handler functions.
 *
 * Functions of this type \em must not return. Typically they
 * call <code>abort()</code> or do a <code>longjmp()</code>.
 *
 * \param reason Textual description of the error.
 * \param file The source file in which the error is noticed.
 * \param line The number of the line in the source file which triggered the error
 */

typedef void igraph_fatal_handler_t (const char *reason, const char *file, int line);

/**
 * \function igraph_set_fatal_handler
 * \brief Installs a fatal error handler.
 *
 * Installs the supplied fatal error handler function.
 *
 * </para><para>
 * Fatal error handler functions \em must not return. Typically, the fatal
 * error handler would either call <code>abort()</code> or <code>longjmp()</code>.
 *
 * \param new_handler The new fatal error handler function to install.
 *        Supply a null pointer here to uninstall the current
 *        fatal error handler, without installing a new one.
 * \return The current fatal error handler function.
 */

IGRAPH_EXPORT igraph_fatal_handler_t* igraph_set_fatal_handler(igraph_fatal_handler_t* new_handler);

/**
 * \var igraph_fatal_handler_abort
 * \brief Abort program in case of fatal error.
 *
 * The default fatal error handler, prints an error message and aborts the program.
 */

IGRAPH_EXPORT igraph_fatal_handler_t igraph_fatal_handler_abort;

/**
 * \function igraph_fatal
 * \brief Triggers a fatal error.
 *
 * This function triggers a fatal error. Typically it is called indirectly through
 * \ref IGRAPH_FATAL() or \ref IGRAPH_ASSERT().
 *
 * \param reason Textual description of the error.
 * \param file The source file in which the error was noticed.
 * \param line The number of line in the source file which triggered the error.
 */

IGRAPH_EXPORT IGRAPH_NORETURN void igraph_fatal(const char *reason, const char *file, int line);

/**
 * \function igraph_fatalf
 * \brief Triggers a fatal error, printf-like syntax.
 *
 * This function is similar to \ref igraph_fatal(), but
 * uses a printf-like syntax. It substitutes the additional arguments
 * into the \p reason template string and calls \ref igraph_fatal().
 *
 * \param reason Textual description of the error.
 * \param file The source file in which the error was noticed.
 * \param line The number of line in the source file which triggered the error.
 * \param ... The additional arguments to be substituted into the template string.
 */

IGRAPH_EXPORT IGRAPH_NORETURN void igraph_fatalf(const char *reason, const char *file, int line, ...);

/**
 * \define IGRAPH_FATALF
 * \brief Triggers a fatal error, with printf-like syntax.
 *
 * \a igraph functions can use this macro when a fatal error occurs and
 * want to pass on extra information to the user about what went wrong.
 * It calls \ref igraph_fatalf() with the proper parameters.
 * \param reason Textual description of the error, a template string
 *        with the same syntax as the standard printf C library function.
 * \param ... The additional arguments to be substituted into the
 *        template string.
 */

#define IGRAPH_FATALF(reason, ...) \
    do { \
        igraph_fatalf(reason, IGRAPH_FILE_BASENAME, __LINE__, \
                      __VA_ARGS__); \
    } while (0)

/**
 * \define IGRAPH_FATAL
 * \brief Triggers a fatal error.
 *
 * This is the usual way of triggering a fatal error from an igraph
 * function. It calls \ref igraph_fatal().
 *
 * </para><para>
 * Use this macro only in situations where the error cannot be handled.
 * The normal way to handle errors is \ref IGRAPH_ERROR().
 *
 * \param reason The error message.
 */

#define IGRAPH_FATAL(reason) \
    do { \
        igraph_fatal(reason, IGRAPH_FILE_BASENAME, __LINE__); \
    } while (0)

/**
 * \define IGRAPH_ASSERT
 * \brief igraph-specific replacement for <code>assert()</code>.
 *
 * This macro is like the standard <code>assert()</code>, but instead of
 * calling <code>abort()</code>, it calls \ref igraph_fatal(). This allows for returning
 * the control to the calling program, e.g. returning to the top level in a high-level
 * \a igraph interface.
 *
 * </para><para>
 * Unlike <code>assert()</code>, <code>IGRAPH_ASSERT()</code> is not disabled
 * when the \c NDEBUG macro is defined.
 *
 * </para><para>
 * This macro is meant for internal use by \a igraph.
 *
 * </para><para>
 * Since a typial fatal error handler does a <code>longjmp()</code>, avoid using this
 * macro in C++ code. With most compilers, destructor will not be called when
 * <code>longjmp()</code> leaves the current scope.
 *
 * \param condition The condition to be checked.
 */

#define IGRAPH_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            igraph_fatal("Assertion failed: " #condition, IGRAPH_FILE_BASENAME, __LINE__); \
        } \
    } while (0)

__END_DECLS

#endif
