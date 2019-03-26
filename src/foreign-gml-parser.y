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
#include "igraph_hacks_internal.h"
#include "igraph_math.h"
#include "igraph_gml_tree.h"
#include "foreign-gml-header.h"
#include "foreign-gml-parser.h"

#define yyscan_t void*

int igraph_gml_yylex(YYSTYPE* lvalp, YYLTYPE* llocp, void *scanner);
int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context, 
		       const char *s);
char *igraph_gml_yyget_text (yyscan_t yyscanner );
int igraph_gml_yyget_leng (yyscan_t yyscanner );
void igraph_i_gml_get_keyword(char *s, int len, void *res);
void igraph_i_gml_get_string(char *s, int len, void *res);
double igraph_i_gml_get_real(char *s, int len);
igraph_gml_tree_t *igraph_i_gml_make_numeric(char* s, int len, double value);
igraph_gml_tree_t *igraph_i_gml_make_numeric2(char* s, int len, 
					      char *v, int vlen);
igraph_gml_tree_t *igraph_i_gml_make_string(char* s, int len, 
					    char *value, int valuelen);
igraph_gml_tree_t *igraph_i_gml_make_list(char* s, int len, 
					  igraph_gml_tree_t *list);
igraph_gml_tree_t *igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2);

#define scanner context->scanner
#define USE(x) /*(x)*/

%}

%pure-parser
%output="y.tab.c"
%name-prefix="igraph_gml_yy"
%defines
%locations
%error-verbose
%parse-param { igraph_i_gml_parsedata_t* context }
%lex-param { void *scanner }

%union {
   struct {
      char *s;
      int len;
   } str;
   void *tree;
   double real;
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

%destructor { igraph_Free($$.s); } string key KEYWORD;
%destructor { igraph_gml_tree_destroy($$); } list keyvalue;

%%

input:   list      { context->tree=$1; } 
       | list EOFF { context->tree=$1; }
;

list:   keyvalue      { $$=$1; } 
      | list keyvalue { $$=igraph_i_gml_merge($1, $2); };

keyvalue:   key num
            { $$=igraph_i_gml_make_numeric($1.s, $1.len, $2); }
	  | key string
            { $$=igraph_i_gml_make_string($1.s, $1.len, $2.s, $2.len); }
          | key LISTOPEN list LISTCLOSE
            { $$=igraph_i_gml_make_list($1.s, $1.len, $3); }
          | key key
            { $$=igraph_i_gml_make_numeric2($1.s, $1.len, $2.s, $2.len); }
;

key: KEYWORD { igraph_i_gml_get_keyword(igraph_gml_yyget_text(scanner), 
					igraph_gml_yyget_leng(scanner), 
					&$$); USE($1) };
num : NUM { $$=igraph_i_gml_get_real(igraph_gml_yyget_text(scanner), 
				     igraph_gml_yyget_leng(scanner)); };

string: STRING { igraph_i_gml_get_string(igraph_gml_yyget_text(scanner), 
					 igraph_gml_yyget_leng(scanner), 
					 &$$); };

%%

int igraph_gml_yyerror(YYLTYPE* locp, igraph_i_gml_parsedata_t *context, 
		       const char *s) {
  snprintf(context->errmsg, sizeof(context->errmsg)/sizeof(char)-1, 
	   "Parse error in GML file, line %i (%s)", 
	   locp->first_line, s);
  return 0;
}

void igraph_i_gml_get_keyword(char *s, int len, void *res) {
  struct { char *s; int len; } *p=res;
  p->s=igraph_Calloc(len+1, char);
  if (!p->s) { 
    igraph_error("Cannot read GML file", __FILE__, __LINE__, IGRAPH_PARSEERROR);
  }
  memcpy(p->s, s, sizeof(char)*len);
  p->s[len]='\0';
  p->len=len;
}

void igraph_i_gml_get_string(char *s, int len, void *res) {
  struct { char *s; int len; } *p=res;
  p->s=igraph_Calloc(len-1, char);
  if (!p->s) { 
    igraph_error("Cannot read GML file", __FILE__, __LINE__, IGRAPH_PARSEERROR);
  }
  memcpy(p->s, s+1, sizeof(char)*(len-2));
  p->s[len-2]='\0';
  p->len=len-2;
}

double igraph_i_gml_get_real(char *s, int len) {
  igraph_real_t num;
  char tmp=s[len];
  s[len]='\0';
  sscanf(s, "%lf", &num);
  s[len]=tmp;
  return num;
} 

igraph_gml_tree_t *igraph_i_gml_make_numeric(char* s, int len, double value) {
  igraph_gml_tree_t *t=igraph_Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  if (floor(value)==value) {
    igraph_gml_tree_init_integer(t, s, len, value);
  } else {
    igraph_gml_tree_init_real(t, s, len, value);
  }
  
  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_numeric2(char* s, int len, 
					      char *v, int vlen) {
  igraph_gml_tree_t *t=igraph_Calloc(1, igraph_gml_tree_t);
  char tmp=v[vlen];
  igraph_real_t value=0;
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  v[vlen]='\0';
  if (strcasecmp(v, "inf")) {
    value=IGRAPH_INFINITY;
  } else if (strcasecmp(v, "nan")) {
    value=IGRAPH_NAN;
  } else {
    igraph_error("Parse error", __FILE__, __LINE__, IGRAPH_PARSEERROR);
  }
  v[vlen]=tmp;
  igraph_gml_tree_init_real(t, s, len, value);  

  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_string(char* s, int len, 
					    char *value, int valuelen) {
  igraph_gml_tree_t *t=igraph_Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_string(t, s, len, value, valuelen);

  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_list(char* s, int len, 
					  igraph_gml_tree_t *list) {
  
  igraph_gml_tree_t *t=igraph_Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_tree(t, s, len, list);

  return t;
}

igraph_gml_tree_t *igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2) {

  igraph_gml_tree_mergedest(t1, t2);
  igraph_Free(t2);

  return t1;
}
