/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "internal/glpk_support.h"

#ifdef HAVE_GLPK

#include "igraph_error.h"

#include "core/interruption.h"

#include <stdio.h>

IGRAPH_THREAD_LOCAL igraph_i_glpk_error_info_t igraph_i_glpk_error_info;

/* glp_at_error() was added in GLPK 4.57. Due to the R interface, we need to
 * support ancient GLPK versions like GLPK 4.38 so we need to guard the
 * invocation of glp_at_error(). Note that this is a temporary workaround only
 * for sake of supporting R 4.1, so it is enabled only if USING_R is defined */
#ifdef USING_R
#  define HAS_GLP_AT_ERROR (GLP_MAJOR_VERSION > 4 || (GLP_MAJOR_VERSION == 4 && GLP_MINOR_VERSION >= 57))
#else
#  define HAS_GLP_AT_ERROR 1
#endif

int igraph_i_glpk_terminal_hook(void *info, const char *s) {
    IGRAPH_UNUSED(info);

    if (igraph_i_interruption_handler &&
        !igraph_i_glpk_error_info.is_interrupted &&
        igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) {
        /* If an interruption has already occurred, do not set another error,
           to avoid an infinite loop between the term_hook (this function)
           and the error_hook. */
        igraph_i_glpk_error_info.is_interrupted = true;
        glp_error("GLPK was interrupted."); /* This dummy message is never printed */
#if HAS_GLP_AT_ERROR
    } else if (glp_at_error()) {
        /* Copy the error messages into a buffer for later reporting */
        /* We must use glp_at_error() instead of igraph_i_glpk_error_info.is_error
         * to determine if a message is an error message, as the reporting function is
         * called before the error function. */
        const size_t n = sizeof(igraph_i_glpk_error_info.msg) / sizeof(char) - 1;
        while (*s != '\0' && igraph_i_glpk_error_info.msg_ptr < igraph_i_glpk_error_info.msg + n) {
            *(igraph_i_glpk_error_info.msg_ptr++) = *(s++);
        }
        *igraph_i_glpk_error_info.msg_ptr = '\0';
#endif
    }

    return 1; /* Non-zero return value signals to GLPK not to print to the terminal */
}

void igraph_i_glpk_error_hook(void *info) {
    IGRAPH_UNUSED(info);
    igraph_i_glpk_error_info.is_error = true;
    glp_free_env();
    longjmp(igraph_i_glpk_error_info.jmp, 1);
}

void igraph_i_glpk_interruption_hook(glp_tree *tree, void *info) {
    IGRAPH_UNUSED(info);

    /* This is a callback function meant to be used with glp_intopt(),
       in order to support interruption. It is essentially a GLPK-compatible
       replacement for IGRAPH_ALLOW_INTERRUPTION().
       Calling glp_ios_terminate() from glp_intopt()'s callback function
       signals to GLPK that it should terminate the optimization and return
       with the code GLP_ESTOP.
    */
    if (igraph_i_interruption_handler) {
        if (igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) {
            glp_ios_terminate(tree);
        }
    }
}

/**
 * \ingroup internal
 * \function igraph_i_glp_delete_prob
 * \brief Safe replacement for glp_delete_prob().
 *
 * This function is meant to be used with IGRAPH_FINALLY()
 * in conjunction with glp_create_prob().
 *
 * When using GLPK, normally glp_delete_prob() is used to free
 * problems created with glp_create_prob(). However, when GLPK
 * encounters an error, the error handler installed by igraph
 * will call glp_free_env() which invalidates all problems.
 * Calling glp_delete_prob() would then lead to a crash.
 * This replacement function avoids this situation by first
 * checking if GLPK is at an error state.
 */
void igraph_i_glp_delete_prob(glp_prob *p) {
    if (! igraph_i_glpk_error_info.is_error) {
        glp_delete_prob(p);
    }
}

igraph_error_t igraph_i_glpk_check(int retval, const char* message) {
    const char *code = "none";
    char message_and_code[4096];
    igraph_error_t ret;

    if (retval == 0) {
        return IGRAPH_SUCCESS;
    }

    /* handle errors */
#define HANDLE_CODE(c)  case c: code = #c; ret = IGRAPH_##c; break;
#define HANDLE_CODE2(c) case c: code = #c; ret = IGRAPH_FAILURE; break;
#define HANDLE_CODE3(c) case c: code = #c; ret = IGRAPH_INTERRUPTED; break;
    switch (retval) {
        HANDLE_CODE(GLP_EBOUND);
        HANDLE_CODE(GLP_EROOT);
        HANDLE_CODE(GLP_ENOPFS);
        HANDLE_CODE(GLP_ENODFS);
        HANDLE_CODE(GLP_EFAIL);
        HANDLE_CODE(GLP_EMIPGAP);
        HANDLE_CODE(GLP_ETMLIM);

        HANDLE_CODE3(GLP_ESTOP);

        HANDLE_CODE2(GLP_EBADB);
        HANDLE_CODE2(GLP_ESING);
        HANDLE_CODE2(GLP_ECOND);
        HANDLE_CODE2(GLP_EOBJLL);
        HANDLE_CODE2(GLP_EOBJUL);
        HANDLE_CODE2(GLP_EITLIM);

    default:
        IGRAPH_ERROR("Unknown GLPK error.", IGRAPH_FAILURE);
    }
#undef HANDLE_CODE
#undef HANDLE_CODE2
#undef HANDLE_CODE3

    snprintf(message_and_code, sizeof(message_and_code) / sizeof(message_and_code[0]),
            "%s (%s)", message, code);
    IGRAPH_ERROR(message_and_code, ret);
}

#else

int igraph_glpk_dummy(void) {
    /* get rid of "ISO C requires a translation unit to contain at least one
     * declaration" warning */
    return 'd' + 'u' + 'm' + 'm' + 'y';
}

#endif
