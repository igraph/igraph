/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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

%{

/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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

#include <stdio.h>
#include <string.h>

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"
#include "config.h"

#include "core/math.h"
#include "io/ncol-header.h"
#include "io/parsers/ncol-parser.h"
#include "io/parsers/ncol-lexer.h"
#include "internal/hacks.h"

/* This macro must be used only in Bison actions, in place of IGRAPH_CHECK(). */
#define IGRAPH_YY_CHECK(expr) \
    do { \
        igraph_error_t igraph_i_ret = (expr); \
        if (IGRAPH_UNLIKELY(igraph_i_ret != IGRAPH_SUCCESS)) { \
            context->igraph_errno = igraph_i_ret; \
            YYABORT; \
        } \
    } while (0)

/* This macro must be used only in Bison actions, in place of IGRAPH_CHECK(). */
#define IGRAPH_YY_ERRORF(reason, errno, ...) \
    do { \
        igraph_errorf(reason, IGRAPH_FILE_BASENAME, __LINE__, \
                      errno, __VA_ARGS__) ; \
        context->igraph_errno = errno; \
        YYABORT; \
    } while (0)

int igraph_ncol_yyerror(YYLTYPE* locp,
                        igraph_i_ncol_parsedata_t *context,
                        const char *s);
igraph_real_t igraph_ncol_get_number(const char *str, yy_size_t len);

#define scanner context->scanner
%}

%pure-parser
/* bison: do not remove the equals sign; macOS XCode ships with bison 2.3, which
 * needs the equals sign */
%name-prefix="igraph_ncol_yy"
%defines
%locations
%error-verbose
%parse-param { igraph_i_ncol_parsedata_t* context }
%lex-param { void *scanner }

%union {
  igraph_integer_t edgenum;
  igraph_real_t weightnum;
}

%type <edgenum>   edgeid
%type <weightnum> weight

%token ALNUM
%token NEWLINE
%token ERROR

%%

input :    /* empty */
         | input NEWLINE
         | input edge
;

edge :   edgeid edgeid NEWLINE        {
           IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1));
           IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $2));
           IGRAPH_YY_CHECK(igraph_vector_push_back(context->weights, 0.0));
       }
       | edgeid edgeid weight NEWLINE {
           IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1));
           IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $2));
           IGRAPH_YY_CHECK(igraph_vector_push_back(context->weights, $3));
           context->has_weights = 1;
       }
;

edgeid : ALNUM  {
  igraph_integer_t trie_id;
  igraph_trie_get2(context->trie,
    igraph_ncol_yyget_text(scanner),
    igraph_ncol_yyget_leng(scanner),
    &trie_id
  );
  $$ = trie_id;
};

weight : ALNUM  { $$=igraph_ncol_get_number(igraph_ncol_yyget_text(scanner),
                        igraph_ncol_yyget_leng(scanner)); } ;

%%

int igraph_ncol_yyerror(YYLTYPE* locp,
            igraph_i_ncol_parsedata_t *context,
            const char *s) {
    snprintf(context->errmsg, sizeof(context->errmsg)/sizeof(char)-1,
            "Parse error in NCOL file, line %i (%s)",
            locp->first_line, s);
    return 0;
}

igraph_real_t igraph_ncol_get_number(const char *str, yy_size_t length) {
    igraph_real_t num;
    char *tmp=IGRAPH_CALLOC(length+1, char);

    strncpy(tmp, str, length);
    tmp[length]='\0';
    sscanf(tmp, "%lf", &num);
    IGRAPH_FREE(tmp);
    return num;
}
