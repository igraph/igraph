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

int igraph_i_glpk_terminal_hook(void *info, const char *s) {
    IGRAPH_UNUSED(info);

    if (igraph_i_interruption_handler &&
        !igraph_i_glpk_error_info.is_interrupted &&
        igraph_allow_interruption(NULL) != IGRAPH_SUCCESS) {
        /* If an interruption has already occurred, do not set another error,
           to avoid an infinite loop between the term_hook (this function)
           and the error_hook. */
        igraph_i_glpk_error_info.is_interrupted = 1;
        glp_error("GLPK was interrupted."); /* This dummy message is never printed */
    } else if (glp_at_error()) {
        /* Copy the error messages into a buffer for later reporting */
        const size_t n = sizeof(igraph_i_glpk_error_info.msg) / sizeof(char) - 1;
        while (*s != '\0' && igraph_i_glpk_error_info.msg_ptr < igraph_i_glpk_error_info.msg + n) {
            *(igraph_i_glpk_error_info.msg_ptr++) = *(s++);
        }
        *igraph_i_glpk_error_info.msg_ptr = '\0';
    }

    return 1; /* Non-zero return value signals to GLPK not to print to the terminal */
}

void igraph_i_glpk_error_hook(void *info) {
    IGRAPH_UNUSED(info);
    igraph_i_glpk_error_info.is_error = 1;
    glp_free_env();
    longjmp(igraph_i_glpk_error_info.jmp, 1);
}

void igraph_i_glpk_interruption_hook(glp_tree *tree, void *info) {
    IGRAPH_UNUSED(info);

    /* This is a special version of IGRAPH_ALLOW_INTERRUPTION().
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

/* Only delete problem if GLPK is not at an error state.
 * If there is an error, glp_free_env() will be called instead.
 * We do not use glp_at_error() as igraph's old embedded GLPK
 * does not yet have this function. */
void igraph_i_glp_delete_prob(glp_prob *p) {
    if (! igraph_i_glpk_error_info.is_error) {
        glp_delete_prob(p);
    }
}

int igraph_i_glpk_check(int retval, const char* message) {
    char* code = "none";
    char message_and_code[4096];

    if (retval == IGRAPH_SUCCESS) {
        return IGRAPH_SUCCESS;
    }

    /* handle errors */
#define HANDLE_CODE(c) case c: code = #c; retval = IGRAPH_##c; break;
#define HANDLE_CODE2(c) case c: code = #c; retval = IGRAPH_FAILURE; break;
#define HANDLE_CODE3(c) case c: code = #c; retval = IGRAPH_INTERRUPTED; break;
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
        IGRAPH_ERROR("Unknown GLPK error", IGRAPH_FAILURE);
    }
#undef HANDLE_CODE
#undef HANDLE_CODE2
#undef HANDLE_CODE3

    sprintf(message_and_code, "%s (%s)", message, code);
    IGRAPH_ERROR(message_and_code, retval);
}

#else

int igraph_glpk_dummy() {
    /* get rid of "ISO C requires a translation unit to contain at least one
     * declaration" warning */
    return 'd' + 'u' + 'm' + 'm' + 'y';
}

#endif
