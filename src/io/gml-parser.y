/*
   IGraph library.
   Copyright (C) 2009-2022  The igraph development team

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

#include "io/gml-header.h"
#include "io/gml-tree.h"
#include "io/parsers/gml-parser.h"
#include "io/parsers/gml-lexer.h"
#include "io/parse_utils.h"
#include "internal/hacks.h" /* strcasecmp & strndup */
#include "math/safe_intop.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context,
                       const char *s);
static igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, char **res);
static igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, char **res);
static igraph_error_t igraph_i_gml_make_numeric(const char *name,
                                                int line,
                                                igraph_real_t value,
                                                igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_string(const char *name,
                                               int line,
                                               char *value,
                                               igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_list(const char *name,
                                             int line,
                                             igraph_gml_tree_t *list,
                                             igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_make_empty(igraph_gml_tree_t **tree);
static igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2);

#define scanner context->scanner

%}

%pure-parser
/* bison: do not remove the equals sign; macOS XCode ships with bison 2.3, which
 * needs the equals sign */
%name-prefix="igraph_gml_yy"
%defines
%locations
%error-verbose
%expect 2 /* from list rule */
%parse-param { igraph_i_gml_parsedata_t* context }
%lex-param { void *scanner }

%union {
   char *str;
   igraph_gml_tree_t *tree;
   igraph_real_t real;
}

%type <tree>    list;
%type <tree>    keyvalue;
%type <str>     key;
%type <real>    num;
%type <str>     string;

%token STRING           "string"
%token NUM              "number"
%token <str>    KEYWORD "keyword"
%token LISTOPEN         "["
%token LISTCLOSE        "]"
/* The following ensures that the special $end token is shown with a friendly name
 * even in older Bison versions.
 * See https://www.gnu.org/software/bison/manual/bison.html#Token-I18n for more details. */
%token END 0            "end of file" /* friendly name for $end */
%token ERROR

%destructor { free($$); } string key;
%destructor { igraph_gml_tree_destroy($$); } list keyvalue;

%%

input:   list      { context->tree=$1; };

list:   /* empty */   { IGRAPH_YY_CHECK(igraph_i_gml_make_empty(&$$)); }
      | keyvalue      { $$=$1; } /* redundant and causes shift/reduce conflict, but increases performance */
      | list keyvalue { IGRAPH_YY_CHECK(igraph_i_gml_merge($1, $2)); $$ = $1; };

keyvalue:   key num
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric($1, @1.first_line, $2, &$$)); }
          | key string
            { IGRAPH_YY_CHECK(igraph_i_gml_make_string($1, @1.first_line, $2, &$$)); }
          | key LISTOPEN list LISTCLOSE
            { IGRAPH_YY_CHECK(igraph_i_gml_make_list($1, @1.first_line, $3, &$$)); }
;

key: KEYWORD { IGRAPH_YY_CHECK(igraph_i_gml_get_keyword(igraph_gml_yyget_text(scanner),
                               igraph_gml_yyget_leng(scanner),
                               &$$)); };
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

static igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, char **res) {
  *res = strndup(s, len);
  if (! *res) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, char **res) {
  *res = strndup(s+1, len-2);
  if (! *res) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_numeric(const char *name,
                                                int line,
                                                igraph_real_t value,
                                                igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  /* The GML spec only requires support for 32-bit signed integers,
   * but igraph tries to support the same range as igraph_integer_t,
   * so that it can read/write all graphs it can represent.
   * We treat anything out of that range as real. These values end
   * up as igraph_real_t anyway, as igraph does not currently support
   * integer-typed attributes. */
  igraph_real_t trunc_value = trunc(value);
  if (value == trunc_value && igraph_i_is_real_representable_as_integer(trunc_value)) {
    IGRAPH_CHECK(igraph_gml_tree_init_integer(t, name, line, value));
  } else {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, name, line, value));
  }

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_string(const char *name,
                                               int line,
                                               char *value,
                                               igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  /* if igraph_gml_tree_init_string succeeds, the newly created tree node takes
   * ownership of 'value'. If it fails, we need to free 'value' ourselves in order
   * not to leak memory */
  IGRAPH_FINALLY(igraph_free, value);
  IGRAPH_CHECK(igraph_gml_tree_init_string(t, name, line, value));

  IGRAPH_FINALLY_CLEAN(1); /* value */

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_list(const char *name,
                                             int line,
                                             igraph_gml_tree_t *list,
                                             igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
  }
  IGRAPH_FINALLY(igraph_free, t);

  IGRAPH_CHECK(igraph_gml_tree_init_tree(t, name, line, list));

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_make_empty(igraph_gml_tree_t **tree) {
    igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
    if (!t) {
      IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
    }
    IGRAPH_FINALLY(igraph_free, t);

    IGRAPH_CHECK(igraph_gml_tree_init_empty(t));

    *tree = t;
    IGRAPH_FINALLY_CLEAN(1); /* t */

    return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2) {

  IGRAPH_CHECK(igraph_gml_tree_mergedest(t1, t2));
  IGRAPH_FREE(t2);

  return IGRAPH_SUCCESS;
}
