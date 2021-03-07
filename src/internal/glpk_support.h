/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_GLPK_SUPPORT_H
#define IGRAPH_GLPK_SUPPORT_H

#include "config.h"

/* Note: only files calling the GLPK routines directly need to
   include this header.
*/

#ifdef HAVE_GLPK

#include <glpk.h>
#include <setjmp.h>

typedef struct igraph_i_glpk_error_info_s {
    jmp_buf jmp;            /* used for bailing when there is a GLPK error */
    int     is_interrupted; /* Boolean; true if there was an interruption */
    int     is_error;       /* Boolean; true if the error hook was called */
    char    msg[4096];      /* GLPK error messages are collected here */
    char   *msg_ptr;        /* Points to the end (null terminator) of msg */
} igraph_i_glpk_error_info_t;

extern IGRAPH_THREAD_LOCAL igraph_i_glpk_error_info_t igraph_i_glpk_error_info;

int igraph_i_glpk_check(int retval, const char* message);
void igraph_i_glpk_interruption_hook(glp_tree *tree, void *info);
void igraph_i_glpk_error_hook(void *info);
int igraph_i_glpk_terminal_hook(void *info, const char *s);
void igraph_i_glp_delete_prob(glp_prob *p);

#define IGRAPH_GLPK_CHECK(func, message) do { \
        int igraph_i_ret = igraph_i_glpk_check(func, message); \
        if (IGRAPH_UNLIKELY(igraph_i_ret != 0)) { \
            return igraph_i_ret; \
        } } while (0)

/* Notes:
 *  - IGRAPH_GLPK_SETUP() must be called in all top-level functions that
 *    use GLPK, before beginning to use any GLPK functins.
 *  - Do not call glp_term_out(OFF) as interruption support relies on
 *    the terminal hook being called.
 *  - This must be a macro and not a function, as jumping into a function
 *    that has already returned with longjmp() is not possible.
 */
#define IGRAPH_GLPK_SETUP() \
    do { \
        glp_error_hook(igraph_i_glpk_error_hook, NULL); \
        glp_term_hook(igraph_i_glpk_terminal_hook, NULL); \
        igraph_i_glpk_error_info.is_interrupted = 0; \
        igraph_i_glpk_error_info.is_error = 0; \
        igraph_i_glpk_error_info.msg_ptr = igraph_i_glpk_error_info.msg; \
        if (setjmp(igraph_i_glpk_error_info.jmp)) { \
            if (igraph_i_glpk_error_info.is_interrupted) { \
                return IGRAPH_INTERRUPTED; \
            } else { \
                if (igraph_i_glpk_error_info.msg_ptr != igraph_i_glpk_error_info.msg) { \
                    while ( *(igraph_i_glpk_error_info.msg_ptr - 1) == '\n' && \
                            igraph_i_glpk_error_info.msg_ptr > igraph_i_glpk_error_info.msg ) { \
                        igraph_i_glpk_error_info.msg_ptr--; \
                    } \
                    *igraph_i_glpk_error_info.msg_ptr = '\0'; \
                    igraph_error(igraph_i_glpk_error_info.msg, IGRAPH_FILE_BASENAME, __LINE__, IGRAPH_EGLP); \
                } \
                return IGRAPH_EGLP; \
            } \
        } \
    } while (0)

#endif

#endif
