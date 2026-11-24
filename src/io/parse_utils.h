/*
   igraph library.
   Copyright (C) 2022-2026  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_PARSE_UTILS_H
#define IGRAPH_PARSE_UTILS_H

#include "igraph_error.h"
#include "igraph_types.h"

/* This macro must be used only in Bison actions, in place of IGRAPH_CHECK(). */
#define IGRAPH_YY_CHECK(expr) \
    do { \
        igraph_error_t igraph_i_ret = (expr); \
        if (IGRAPH_UNLIKELY(igraph_i_ret != IGRAPH_SUCCESS)) { \
            context->igraph_errno = igraph_i_ret; \
            yyerror(&yylloc, context, "failed"); \
            YYABORT; \
        } \
    } while (0)

/* This macro must be used only in Bison actions, in place of IGRAPH_CHECK(). */
/* Note:
 * Don't name macro argument 'igraph_errno' due to use of context->igraph_errno,
 * or 'errno' due to use of #include <errno.h> in parse_utils.c. */
#define IGRAPH_YY_ERRORF(reason, error_code, ...) \
    do { \
        igraph_errorf(reason, IGRAPH_FILE_BASENAME, __LINE__, \
                      error_code, __VA_ARGS__) ; \
        context->igraph_errno = error_code; \
        YYABORT; \
    } while (0)

void igraph_i_trim_whitespace(const char *str, size_t str_len, const char **res, size_t *res_len);

igraph_error_t igraph_i_fskip_whitespace(FILE *file);

igraph_error_t igraph_i_parse_integer(const char *str, size_t length, igraph_int_t *value);
igraph_error_t igraph_i_parse_real(const char *str, size_t length, igraph_real_t *value);

igraph_error_t igraph_i_fget_integer(FILE *file, igraph_int_t *value);
igraph_error_t igraph_i_fget_real(FILE *file, igraph_real_t *value);

#endif /* IGRAPH_PARSE_UTILS_H */
