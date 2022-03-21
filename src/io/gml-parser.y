/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_error.h"
#include "igraph_memory.h"
#include "config.h"

#include "io/gml-header.h"
#include "io/gml-tree.h"
#include "io/parsers/gml-parser.h"
#include "io/parsers/gml-lexer.h"
#include "io/parse_utils.h"
#include "internal/hacks.h" /* strcasecmp */

#include <stdio.h>
#include <math.h>
#include <string.h>

int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context,
                       const char *s);
static igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, igraph_gml_string_t *res);
static igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, igraph_gml_string_t *res);
static igraph_error_t igraph_i_gml_make_numeric(const igraph_gml_string_t name, igraph_real_t value, igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_numeric2(igraph_gml_string_t name,
                                          igraph_gml_string_t value,
                                          igraph_real_t sign,
                                          igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_string(const igraph_gml_string_t name,
                                        igraph_gml_string_t value,
                                        igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_list(const igraph_gml_string_t name,
                                      igraph_gml_tree_t *list,
                                      igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2);

#define scanner context->scanner
#define USE(x) /*(x)*/

%}

%pure-parser
/* bison: do not remove the equals sign; macOS XCode ships with bison 2.3, which
 * needs the equals sign */
%name-prefix="igraph_gml_yy"
%defines
%locations
%error-verbose
%parse-param { igraph_i_gml_parsedata_t* context }
%lex-param { void *scanner }

%union {
   igraph_gml_string_t str;
   igraph_gml_tree_t *tree;
   igraph_real_t real;
}

%type <tree>    list;
%type <tree>    keyvalue;
%type <str>     key;
%type <real>    num;
%type <str>     string;

%token STRING
%token NUM
%token <str>    KEYWORD
%token MINUS
%token LISTOPEN
%token LISTCLOSE
%token EOFF
%token ERROR

%destructor { igraph_gml_string_destroy(&$$); } string key KEYWORD;
%destructor { igraph_gml_tree_destroy($$); } list keyvalue;

%%

input:   list      { context->tree=$1; }
       | list EOFF { context->tree=$1; }
;

list:   keyvalue      { $$=$1; }
      | list keyvalue { IGRAPH_YY_CHECK(igraph_i_gml_merge($1, $2)); $$ = $1; };

keyvalue:   key num
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric($1, $2, &$$)); }
          | key string
            { IGRAPH_YY_CHECK(igraph_i_gml_make_string($1, $2, &$$)); }
          | key LISTOPEN list LISTCLOSE
            { IGRAPH_YY_CHECK(igraph_i_gml_make_list($1, $3, &$$)); }
          | key key
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric2($1, $2, 1.0, &$$)); }
          | key MINUS key
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric2($1, $3, -1.0, &$$)); }
;

key: KEYWORD { IGRAPH_YY_CHECK(igraph_i_gml_get_keyword(igraph_gml_yyget_text(scanner),
                               igraph_gml_yyget_leng(scanner),
                               &$$)); USE($1); };
num : NUM {
    igraph_real_t val;
    IGRAPH_YY_CHECK(igraph_i_parse_real(igraph_gml_yyget_text(scanner),
                                        igraph_gml_yyget_leng(scanner),
                                        &val));
    $$=val;
};

string: STRING { IGRAPH_YY_CHECK(igraph_i_gml_get_string(igraph_gml_yyget_text(scanner),
                                         igraph_gml_yyget_leng(scanner),
                                         &$$)); };

%%

int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context,
                       const char *s) {
  snprintf(context->errmsg, sizeof(context->errmsg)/sizeof(char)-1,
           "Parse error in GML file, line %i (%s)",
           locp->first_line, s);
  return 0;
}

static igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, igraph_gml_string_t *res) {
  res->str = strndup(s, len);
  if (!res->str) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  res->len=len;
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, igraph_gml_string_t *res) {
  res->str = strndup(s, len-2);
  if (!res->str) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  res->len=len-2;
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_numeric(const igraph_gml_string_t name, igraph_real_t value, igraph_gml_tree_t **tree) {
  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  if (floor(value)==value) {
    IGRAPH_CHECK(igraph_gml_tree_init_integer(t, name, value));
  } else {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, name, value));
  }

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_numeric2(igraph_gml_string_t name,
                                          igraph_gml_string_t value,
                                          igraph_real_t sign,
                                          igraph_gml_tree_t **tree) {

  char tmp = value.str[value.len];

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  value.str[value.len]='\0';

  /* if v == "inf" or v == "nan", the newly created tree node will take ownership
   * of s. If the creation fails, we need to free s and v as well in order not
   * to leak memory */
  IGRAPH_FINALLY(igraph_gml_string_destroy, &name);
  IGRAPH_FINALLY(igraph_gml_string_destroy, &value);
  if (value.len < 3) {
    /* Both "inf" and "nan" are three characters long.
     * Do not read past their end with strcasecmp(). */
    IGRAPH_ERROR("Error while parsing GML.", IGRAPH_PARSEERROR);
  }
  if (strcasecmp(value.str, "inf") == 0) {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, name, sign * IGRAPH_INFINITY));
  } else if (strcasecmp(value.str, "nan") == 0) {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, name, IGRAPH_NAN));
  } else {
    IGRAPH_ERROR("Error while parsing GML.", IGRAPH_PARSEERROR);
  }

  value.str[value.len]=tmp;
  igraph_gml_string_destroy(&value);

  *tree = t;

  IGRAPH_FINALLY_CLEAN(3); /* +name, +t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_string(const igraph_gml_string_t name,
                                        igraph_gml_string_t value,
                                        igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  /* if igraph_gml_tree_init_string succeeds, the newly created tree node takes
   * ownership of 'value'. If it fails, we need to free 'value' ourselves in order
   * not to leak memory */
  IGRAPH_FINALLY(igraph_gml_string_destroy, &value);
  IGRAPH_CHECK(igraph_gml_tree_init_string(t, name, value));

  IGRAPH_FINALLY_CLEAN(1); /* value */

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_list(const igraph_gml_string_t name,
                                      igraph_gml_tree_t *list,
                                      igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  IGRAPH_CHECK(igraph_gml_tree_init_tree(t, name, list));

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2) {

  IGRAPH_CHECK(igraph_gml_tree_mergedest(t1, t2));
  IGRAPH_FREE(t2);

  return IGRAPH_SUCCESS;
}
