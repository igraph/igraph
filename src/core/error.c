/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "config.h"
#include "igraph_error.h"
#include "igraph_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/* Detecting ASan with GCC:
 *   https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
 * Detecting ASan with Clang:
 *   https://clang.llvm.org/docs/AddressSanitizer.html#conditional-compilation-with-has-feature-address-sanitizer
 */
#if defined(__SANITIZE_ADDRESS__)
#  define IGRAPH_SANITIZER_AVAILABLE 1
#elif defined(__has_feature)
#  if __has_feature(address_sanitizer)
#    define IGRAPH_SANITIZER_AVAILABLE 1
#  endif
#endif

#ifdef IGRAPH_SANITIZER_AVAILABLE
#include <sanitizer/asan_interface.h>
#endif

#ifdef USING_R
#include <R.h>
#endif

/***** Helper functions *****/

/* All calls to abort() in this compilation unit must go through igraph_abort(),
 * in order to make it easy for igraph's R interface to not have any reference to abort(),
 * which is disallowed by CRAN.
 *
 * Since the R interface sets its own error / fatal error handlers, this function
 * is never actually called by it.
 *
 * Note that some of the other #ifndef USING_R's in this file are still needed
 * to avoid references to fprintf and stderr.
 */
static IGRAPH_NORETURN void igraph_abort() {
#ifndef USING_R
#ifdef IGRAPH_SANITIZER_AVAILABLE
    fprintf(stderr, "\nStack trace:\n");
    __sanitizer_print_stack_trace();
#endif
    abort();
#else
    /* R's error() function is declared 'noreturn'. We use it here to satisfy the compiler that igraph_abort() does indeed not return. */
    error("igraph_abort() was called. This should never happen. Please report this as an igraph bug, along with steps to reproduce it.");
#endif
}


/***** Handling errors *****/

static IGRAPH_THREAD_LOCAL igraph_error_handler_t *igraph_i_error_handler = 0;
static IGRAPH_THREAD_LOCAL char igraph_i_errormsg_buffer[500];
static IGRAPH_THREAD_LOCAL char igraph_i_warningmsg_buffer[500];
static IGRAPH_THREAD_LOCAL char igraph_i_fatalmsg_buffer[500];

/* Error strings corresponding to each igraph_error_type_t enum value. */
static const char *igraph_i_error_strings[] = {
    /*  0 */ "No error",
    /*  1 */ "Failed",
    /*  2 */ "Out of memory",
    /*  3 */ "Parse error",
    /*  4 */ "Invalid value",
    /*  5 */ "Already exists",
    /*  6 */ "Invalid edge vector",
    /*  7 */ "Invalid vertex id",
    /*  8 */ "Non-square matrix",
    /*  9 */ "Invalid mode",
    /* 10 */ "File operation error",
    /* 11 */ "Unfold infinite iterator",
    /* 12 */ "Unimplemented function call",
    /* 13 */ "Interrupted",
    /* 14 */ "Numeric procedure did not converge",
    /* 15 */ "Matrix-vector product failed",
    /* 16 */ "N must be positive",
    /* 17 */ "NEV must be positive",
    /* 18 */ "NCV must be greater than NEV and less than or equal to N "
    "(and for the non-symmetric solver NCV-NEV >=2 must also hold)",
    /* 19 */ "Maximum number of iterations should be positive",
    /* 20 */ "Invalid WHICH parameter",
    /* 21 */ "Invalid BMAT parameter",
    /* 22 */ "WORKL is too small",
    /* 23 */ "LAPACK error in tridiagonal eigenvalue calculation",
    /* 24 */ "Starting vector is zero",
    /* 25 */ "MODE is invalid",
    /* 26 */ "MODE and BMAT are not compatible",
    /* 27 */ "ISHIFT must be 0 or 1",
    /* 28 */ "NEV and WHICH='BE' are incompatible",
    /* 29 */ "Could not build an Arnoldi factorization",
    /* 30 */ "No eigenvalues to sufficient accuracy",
    /* 31 */ "HOWMNY is invalid",
    /* 32 */ "HOWMNY='S' is not implemented",
    /* 33 */ "Different number of converged Ritz values",
    /* 34 */ "Error from calculation of a real Schur form",
    /* 35 */ "LAPACK (dtrevc) error for calculating eigenvectors",
    /* 36 */ "Unknown ARPACK error",
    /* 37 */ "Negative loop detected while calculating shortest paths",
    /* 38 */ "Internal error, likely a bug in igraph",
    /* 39 */ "Maximum number of iterations reached",
    /* 40 */ "No shifts could be applied during a cycle of the "
    "Implicitly restarted Arnoldi iteration. One possibility "
    "is to increase the size of NCV relative to NEV",
    /* 41 */ "The Schur form computed by LAPACK routine dlahqr "
    "could not be reordered by LAPACK routine dtrsen.",
    /* 42 */ "Big integer division by zero",
    /* 43 */ "GLPK Error, GLP_EBOUND",
    /* 44 */ "GLPK Error, GLP_EROOT",
    /* 45 */ "GLPK Error, GLP_ENOPFS",
    /* 46 */ "GLPK Error, GLP_ENODFS",
    /* 47 */ "GLPK Error, GLP_EFAIL",
    /* 48 */ "GLPK Error, GLP_EMIPGAP",
    /* 49 */ "GLPK Error, GLP_ETMLIM",
    /* 50 */ "GLPK Error, GLP_STOP",
    /* 51 */ "Internal attribute handler error",
    /* 52 */ "Unimplemented attribute combination for this type",
    /* 53 */ "LAPACK call resulted in an error",
    /* 54 */ "Internal DrL error",
    /* 55 */ "Integer or double overflow",
    /* 56 */ "Internal GPLK error",
    /* 57 */ "CPU time exceeded",
    /* 58 */ "Integer or double underflow",
    /* 59 */ "Random walk got stuck",
    /* 60 */ "Search stopped; this error should never be visible to the user, "
    "please report this error along with the steps to reproduce it."
};

const char* igraph_strerror(const int igraph_errno) {
    if (igraph_errno < 0 ||
        ((unsigned long)igraph_errno) >= sizeof(igraph_i_error_strings) / sizeof(char *)) {
        return "Invalid error code; no error string available.";
    }
    return igraph_i_error_strings[igraph_errno];
}

int igraph_error(const char *reason, const char *file, int line,
                 int igraph_errno) {

    if (igraph_i_error_handler) {
        igraph_i_error_handler(reason, file, line, igraph_errno);
#ifndef USING_R
    }  else {
        igraph_error_handler_abort(reason, file, line, igraph_errno);
#endif
    }
    return igraph_errno;
}

int igraph_errorf(const char *reason, const char *file, int line,
                  int igraph_errno, ...) {
    va_list ap;
    va_start(ap, igraph_errno);
    vsnprintf(igraph_i_errormsg_buffer,
              sizeof(igraph_i_errormsg_buffer) / sizeof(char), reason, ap);
    return igraph_error(igraph_i_errormsg_buffer, file, line, igraph_errno);
}

int igraph_errorvf(const char *reason, const char *file, int line,
                   int igraph_errno, va_list ap) {
    vsnprintf(igraph_i_errormsg_buffer,
              sizeof(igraph_i_errormsg_buffer) / sizeof(char), reason, ap);
    return igraph_error(igraph_i_errormsg_buffer, file, line, igraph_errno);
}

#ifndef USING_R
void igraph_error_handler_abort(const char *reason, const char *file,
                                int line, int igraph_errno) {
    fprintf(stderr, "Error at %s:%i : %s - %s.\n",
            file, line, reason, igraph_strerror(igraph_errno));
    igraph_abort();
}
#endif

void igraph_error_handler_ignore(const char *reason, const char *file,
                                 int line, int igraph_errno) {
    IGRAPH_UNUSED(reason);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(igraph_errno);

    IGRAPH_FINALLY_FREE();
}

#ifndef USING_R
void igraph_error_handler_printignore(const char *reason, const char *file,
                                      int line, int igraph_errno) {
    IGRAPH_FINALLY_FREE();
    fprintf(stderr, "Error at %s:%i : %s - %s.\n",
            file, line, reason, igraph_strerror(igraph_errno));
}
#endif

igraph_error_handler_t *igraph_set_error_handler(igraph_error_handler_t *new_handler) {
    igraph_error_handler_t *previous_handler = igraph_i_error_handler;
    igraph_i_error_handler = new_handler;
    return previous_handler;
}


/***** "Finally" stack *****/

IGRAPH_THREAD_LOCAL struct igraph_i_protectedPtr igraph_i_finally_stack[100];

/*
 * Adds another element to the free list
 */

void IGRAPH_FINALLY_REAL(void (*func)(void*), void* ptr) {
    int no = igraph_i_finally_stack[0].all;
    IGRAPH_ASSERT(no < 100);
    IGRAPH_ASSERT(no >= 0);
    igraph_i_finally_stack[no].ptr = ptr;
    igraph_i_finally_stack[no].func = func;
    igraph_i_finally_stack[0].all ++;
    /* printf("--> Finally stack contains now %d elements\n", igraph_i_finally_stack[0].all); */
}

void IGRAPH_FINALLY_CLEAN(int minus) {
    igraph_i_finally_stack[0].all -= minus;
    if (igraph_i_finally_stack[0].all < 0) {
        int left = igraph_i_finally_stack[0].all + minus;
        /* Set to zero in case fatal error handler does a longjmp instead of terminating the process: */
        igraph_i_finally_stack[0].all = 0;
        IGRAPH_FATALF("Corrupt finally stack: trying to pop %d element(s) when only %d left.", minus, left);
    }
    /* printf("<-- Finally stack contains now %d elements\n", igraph_i_finally_stack[0].all); */
}

void IGRAPH_FINALLY_FREE(void) {
    int p;
    /*   printf("[X] Finally stack will be cleaned (contained %d elements)\n", igraph_i_finally_stack[0].all);  */
    for (p = igraph_i_finally_stack[0].all - 1; p >= 0; p--) {
        igraph_i_finally_stack[p].func(igraph_i_finally_stack[p].ptr);
    }
    igraph_i_finally_stack[0].all = 0;
}

int IGRAPH_FINALLY_STACK_SIZE(void) {
    return igraph_i_finally_stack[0].all;
}


/***** Handling warnings *****/

static IGRAPH_THREAD_LOCAL igraph_warning_handler_t *igraph_i_warning_handler = 0;

/**
 * \function igraph_warning_handler_ignore
 * \brief Ignores all warnings.
 *
 * This warning handler function simply ignores all warnings.
 * \param reason Textual description of the warning.
 * \param file The source file in which the warning was noticed.
 * \param line The number of line in the source file which triggered the
 *         warning..
 * \param igraph_errno Warnings could have potentially error codes as well,
 *        but this is currently not used in igraph.
 */

void igraph_warning_handler_ignore(const char *reason, const char *file,
                                   int line, int igraph_errno) {
    IGRAPH_UNUSED(reason);
    IGRAPH_UNUSED(file);
    IGRAPH_UNUSED(line);
    IGRAPH_UNUSED(igraph_errno);
}

#ifndef USING_R

/**
 * \function igraph_warning_handler_print
 * \brief Prints all warnings to the standard error.
 *
 * This warning handler function simply prints all warnings to the
 * standard error.
 * \param reason Textual description of the warning.
 * \param file The source file in which the warning was noticed.
 * \param line The number of line in the source file which triggered the
 *         warning..
 * \param igraph_errno Warnings could have potentially error codes as well,
 *        but this is currently not used in igraph.
 */

void igraph_warning_handler_print(const char *reason, const char *file,
                                  int line, int igraph_errno) {
    IGRAPH_UNUSED(igraph_errno);
    fprintf(stderr, "Warning at %s:%i : %s\n", file, line, reason);
}
#endif

int igraph_warning(const char *reason, const char *file, int line,
                   int igraph_errno) {

    if (igraph_i_warning_handler) {
        igraph_i_warning_handler(reason, file, line, igraph_errno);
#ifndef USING_R
    }  else {
        igraph_warning_handler_print(reason, file, line, igraph_errno);
#endif
    }
    return igraph_errno;
}

int igraph_warningf(const char *reason, const char *file, int line,
                    int igraph_errno, ...) {
    va_list ap;
    va_start(ap, igraph_errno);
    vsnprintf(igraph_i_warningmsg_buffer,
              sizeof(igraph_i_warningmsg_buffer) / sizeof(char), reason, ap);
    return igraph_warning(igraph_i_warningmsg_buffer, file, line,
                          igraph_errno);
}

igraph_warning_handler_t *igraph_set_warning_handler(igraph_warning_handler_t *new_handler) {
    igraph_warning_handler_t *previous_handler = igraph_i_warning_handler;
    igraph_i_warning_handler = new_handler;
    return previous_handler;
}


/***** Handling fatal errors *****/

static IGRAPH_THREAD_LOCAL igraph_fatal_handler_t *igraph_i_fatal_handler = NULL;

igraph_fatal_handler_t *igraph_set_fatal_handler(igraph_fatal_handler_t *new_handler) {
    igraph_fatal_handler_t *previous_handler = igraph_i_fatal_handler;
    igraph_i_fatal_handler = new_handler;
    return previous_handler;
}

#ifndef USING_R
void igraph_fatal_handler_abort(const char *reason, const char *file, int line) {
    fprintf(stderr, "Fatal error at %s:%i : %s\n", file, line, reason);
    igraph_abort();
}
#endif

void igraph_fatal(const char *reason, const char *file, int line) {
    if (igraph_i_fatal_handler) {
        igraph_i_fatal_handler(reason, file, line);
#ifndef USING_R
    }  else {
        igraph_fatal_handler_abort(reason, file, line);
#endif
    }
    /* The following line should never be reached, as fatal error handlers are not supposed to return.
       It is here to satisfy the compiler that this function indeed does not return. */
    igraph_abort();
}

void igraph_fatalf(const char *reason, const char *file, int line, ...) {
    va_list ap;
    va_start(ap, line);
    vsnprintf(igraph_i_fatalmsg_buffer, sizeof(igraph_i_fatalmsg_buffer) / sizeof(char), reason, ap);
    va_end(ap);
    igraph_fatal(igraph_i_fatalmsg_buffer, file, line);
}
