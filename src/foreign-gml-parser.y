/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#include "error.h"
#include "memory.h"
#include "gml_tree.h"

extern int igraph_gml_yylex(void);
extern long int igraph_gml_mylineno;
extern char *igraph_gml_yytext;
extern int igraph_gml_yyleng;
igraph_gml_tree_t *igraph_i_gml_parsed_tree;
int igraph_gml_yyerror(char *s);
void igraph_i_gml_reset_scanner(void);
long int igraph_i_gml_get_integer(char *s, int len);
double igraph_i_gml_get_real(char*s1, int len1, char* s2, int len2, int mant);
igraph_gml_tree_t *igraph_i_gml_make_integer(char* s, int len, long int value);
igraph_gml_tree_t *igraph_i_gml_make_real(char* s, int len, double value);
igraph_gml_tree_t *igraph_i_gml_make_string(char* s, int len, 
					    char *value, int valuelen);
igraph_gml_tree_t *igraph_i_gml_make_list(char* s, int len, 
					  igraph_gml_tree_t *list);
igraph_gml_tree_t *igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2);
%}

%output="y.tab.c"
%name-prefix="igraph_gml_yy"
%defines

%union {
   struct {
      char *s;
      int len;
   } str;
   void *tree;
   long int integer;
   double real;
}

%type <tree>    list;
%type <tree>    keyvalue;
%type <str>     string;
%type <str>     key;
%type <str>     alnumstring;
%type <integer> integer;
%type <integer> sign
%type <str>     digits;
%type <str>     digit;
%type <str>     optdigits;
%type <real>    real;
%type <integer> mantissa;
%type <str>     instrings;
%type <str>     instring;
%type <str>     letter;
%type <str>     allowedchar;

%token WS
%token LETTER
%token NUM
%token OTHER

%%

input: list { igraph_i_gml_parsed_tree=$1; };

list:   keyvalue optwhitespace   { $$=$1; } 
      | keyvalue whitespace list { $$=igraph_i_gml_merge($3, $1); };

keyvalue:   key whitespace integer
              { $$=igraph_i_gml_make_integer($1.s, $1.len, $3); }
          | key whitespace real     
              { $$=igraph_i_gml_make_real($1.s, $1.len, $3); }
	  | key whitespace string   
              { $$=igraph_i_gml_make_string($1.s, $1.len, $3.s, $3.len); }
          | key whitespace '[' optwhitespace list optwhitespace ']'
              { $$=igraph_i_gml_make_list($1.s, $1.len, $5); }
;

key:   letter             { $$=$1; }
     | letter alnumstring { $$.s=$1.s; $$.len=$2.len+1; } ;

alnumstring:   letter               { $$=$1; }
             | NUM                  { $$.s=igraph_gml_yytext; $$.len=1; }
             | alnumstring letter   { $$.s=$1.s; $$.len=$1.len+1; }
             | alnumstring NUM      { $$.s=$1.s; $$.len=$1.len+1; }

integer: sign digits { $$=$1*igraph_i_gml_get_integer($2.s, $2.len); };

digits:   digit         { $$=$1; }
        | digit digits  { $$.s=$1.s; $$.len=$2.len+1; } ;

digit: NUM { $$.s=igraph_gml_yytext; $$.len=1; } ;

real: sign optdigits '.' digits mantissa { 
		 $$=$1*igraph_i_gml_get_real($2.s, $2.len, $4.s, $4.len, $5); } ;

optdigits: /* empty */ { $$.len=0; } | digits { $$=$1; } ;

string: '"' instrings '"' { $$=$2; } ;

sign:   /* empty */ { $$= 1; } 
      | '+'         { $$= 1; } 
      | '-'         { $$=-1; };

mantissa:   /* empty */    { $$ = 0; }
          | e sign digits  { $$ = $2 * igraph_i_gml_get_integer($3.s,$3.len); } ;

e: 'e' | 'E' ;

instrings:   /* empty */          { $$.len=0; } 
           | instring instrings   { $$.s=$1.s; $$.len=$2.len+$1.len; }
;

instring:   allowedchar         { $$.s=igraph_gml_yytext; $$.len=1; }
          | '&' alnumstring ';' { $$.s=$2.s-1; $$.len=$2.len+2; } ;

optwhitespace: /* empty */ | whitespace;

whitespace: WS | WS whitespace;

letter: letter2 { $$.len=1; $$.s=igraph_gml_yytext; } ;

letter2: LETTER | 'e' | 'E';

allowedchar: allowedchar2 { $$.len=1; $$.s=igraph_gml_yytext; } ;

allowedchar2: letter | NUM | OTHER | '[' | ']' | ';' | WS | '\n';

%%

int igraph_gml_yyerror(char *s)
{
  char str[200];  
  igraph_i_gml_reset_scanner();
  snprintf(str, sizeof(str), "Parse error in GML file, line %li", 
	   (long)igraph_gml_mylineno);
  IGRAPH_ERROR(str, IGRAPH_PARSEERROR);
  return IGRAPH_PARSEERROR;
}

long int igraph_i_gml_get_integer(char *s, int len) {
  long int num;
  int tmp=s[len];
  s[len]='\0';
  num=strtol(s, 0, 10);
  s[len]=tmp;
  return num;
}

double igraph_i_gml_get_real(char*s1, int len1, char* s2, int len2, int mant) {
  double fracpart, num=0;
  int tmp;
  
  if (len1!=0) {
    tmp=s1[len1];  
    s1[len1]='\0';
    num=strtod(s1, 0);
    s1[len1]=tmp;
  }
  
  if (len2 != 0) {
    tmp=s2[len2];
    s2[len2]='\0';
    fracpart=strtod(s2, 0);
    s2[len2]=tmp;
    num += fracpart / pow(10, len2);
  }
  
  num *= pow(10, mant);

  return 0.0;
}

igraph_gml_tree_t *igraph_i_gml_make_integer(char* s, int len, long int value) {
  igraph_gml_tree_t *t=Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_integer(t, s, len, value);
  
  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_real(char* s, int len, double value) {
  igraph_gml_tree_t *t=Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_real(t, s, len, value);
  
  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_string(char* s, int len, 
					    char *value, int valuelen) {
  igraph_gml_tree_t *t=Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_string(t, s, len, value, valuelen);

  return t;
}

igraph_gml_tree_t *igraph_i_gml_make_list(char* s, int len, 
					  igraph_gml_tree_t *list) {
  
  igraph_gml_tree_t *t=Calloc(1, igraph_gml_tree_t);
  if (!t) { 
    igraph_error("Cannot build GML tree", __FILE__, __LINE__, IGRAPH_ENOMEM);
    return 0;
  }
  igraph_gml_tree_init_tree(t, s, len, list);

  return t;
}

igraph_gml_tree_t *igraph_i_gml_merge(igraph_gml_tree_t *t1, igraph_gml_tree_t* t2) {

  igraph_gml_tree_mergedest(t1, t2);
  Free(t2);

  return t1;
}
