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

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "igraph_error.h"
#include "igraph_memory.h"
#include "config.h"

#include "core/math.h"
#include "io/gml-header.h"
#include "io/gml-tree.h"
#include "io/parsers/gml-parser.h"
#include "io/parsers/gml-lexer.h"
#include "io/parse_utils.h"
#include "internal/hacks.h" /* strcasecmp */

int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context,
                       const char *s);
igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, void *res);
igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, void *res);
igraph_error_t igraph_i_gml_make_numeric(const char* s, size_t len, double value, igraph_gml_tree_t **tree);
igraph_error_t igraph_i_gml_make_numeric2(char* s, size_t len,
                                          char *v, size_t vlen,
                                          igraph_gml_tree_t **tree);
igraph_error_t igraph_i_gml_make_string(const char* s, size_t len,
                                        char *value, size_t valuelen,
                                        igraph_gml_tree_t **tree);
igraph_error_t igraph_i_gml_make_list(const char* s, size_t len,
                                      igraph_gml_tree_t *list,
                                      igraph_gml_tree_t **tree);
igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2);

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
   struct {
      char *s;
      size_t len;
   } str;
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
%token LISTOPEN
%token LISTCLOSE
%token EOFF
%token ERROR

%destructor { IGRAPH_FREE($$.s); } string key KEYWORD;
%destructor { igraph_gml_tree_destroy($$); } list keyvalue;

%%

input:   list      { context->tree=$1; }
       | list EOFF { context->tree=$1; }
;

list:   keyvalue      { $$=$1; }
      | list keyvalue { IGRAPH_YY_CHECK(igraph_i_gml_merge($1, $2)); $$ = $1; };

keyvalue:   key num
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric($1.s, $1.len, $2, &$$)); }
          | key string
            { IGRAPH_YY_CHECK(igraph_i_gml_make_string($1.s, $1.len, $2.s, $2.len, &$$)); }
          | key LISTOPEN list LISTCLOSE
            { IGRAPH_YY_CHECK(igraph_i_gml_make_list($1.s, $1.len, $3, &$$)); }
          | key key
            { IGRAPH_YY_CHECK(igraph_i_gml_make_numeric2($1.s, $1.len, $2.s, $2.len, &$$)); }
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

igraph_error_t igraph_i_gml_get_keyword(const char *s, size_t len, void *res) {
  struct { char *s; size_t len; } *p=res;
  p->s=IGRAPH_CALLOC(len+1, char);
  if (!p->s) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM);
  }
  memcpy(p->s, s, sizeof(char)*len);
  p->s[len]='\0';
  p->len=len;
  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_get_string(const char *s, size_t len, void *res) {
  struct { char *s; size_t len; } *p=res;
  p->s=IGRAPH_CALLOC(len-1, char);
  if (!p->s) {
    IGRAPH_ERROR("Cannot read GML file.", IGRAPH_ENOMEM);
  }
  memcpy(p->s, s+1, sizeof(char)*(len-2));
  p->s[len-2]='\0';
  p->len=len-2;
  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_make_numeric(const char *s, size_t len, igraph_real_t value, igraph_gml_tree_t **tree) {
  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, t);

  if (floor(value)==value) {
    IGRAPH_CHECK(igraph_gml_tree_init_integer(t, s, len, value));
  } else {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, s, len, value));
  }

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_make_numeric2(char *s, size_t len, 
                                          char *v, size_t vlen, 
                                          igraph_gml_tree_t **tree) {

  char tmp = v[vlen];

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, t);

  v[vlen]='\0';

  /* if v == "inf" or v == "nan", the newly created tree node will take ownership
   * of s. If the creation fails, we need to free s and v as well in order not
   * to leak memory */
  IGRAPH_FINALLY(igraph_free, s);
  IGRAPH_FINALLY(igraph_free, v);
  if (strcasecmp(v, "inf")) {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, s, len, IGRAPH_INFINITY));
  } else if (strcasecmp(v, "nan")) {
    IGRAPH_CHECK(igraph_gml_tree_init_real(t, s, len, IGRAPH_NAN));
  } else {
    IGRAPH_ERROR("Error while parsing GML.", IGRAPH_PARSEERROR);
  }

  v[vlen]=tmp;
  free(v);
  IGRAPH_FINALLY_CLEAN(2); /* s no longer needs to be freed */

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_make_string(const char *s, size_t len,
                                        char *value, size_t valuelen,
                                        igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, t);

  /* if igraph_gml_tree_init_string succeeds, the newly created tree node takes
   * ownership of 'value'. If it fails, we need to free 'value' ourselves in order
   * not to leak memory */
  IGRAPH_FINALLY(igraph_free, value);
  IGRAPH_CHECK(igraph_gml_tree_init_string(t, s, len, value, valuelen));

  IGRAPH_FINALLY_CLEAN(1); /* value */

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_make_list(const char* s, size_t len,
                                      igraph_gml_tree_t *list,
                                      igraph_gml_tree_t **tree) {

  igraph_gml_tree_t *t = IGRAPH_CALLOC(1, igraph_gml_tree_t);
  if (!t) {
    IGRAPH_ERROR("Cannot build GML tree.", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, t);

  IGRAPH_CHECK(igraph_gml_tree_init_tree(t, s, len, list));

  *tree = t;
  IGRAPH_FINALLY_CLEAN(1); /* t */

  return IGRAPH_SUCCESS;
}

igraph_error_t igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2) {

  IGRAPH_CHECK(igraph_gml_tree_mergedest(t1, t2));
  IGRAPH_FREE(t2);

  return IGRAPH_SUCCESS;
}
