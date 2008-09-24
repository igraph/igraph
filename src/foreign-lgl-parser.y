/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <string.h>
extern int igraph_lgl_yylex(void);
extern long int igraph_lgl_mylineno;
extern char *igraph_lgl_yytext;
extern int igraph_lgl_yyleng;
char *igraph_i_lgl_errmsg=0;
int igraph_lgl_yyerror(char *s);
#include "types.h" 
#include "memory.h"
#include "error.h"
#include "config.h"
extern igraph_vector_t *igraph_lgl_vector;
extern igraph_vector_t *igraph_lgl_weights;
extern igraph_trie_t *igraph_lgl_trie;
igraph_real_t igraph_lgl_get_number(const char *str, long int len);
void igraph_i_lgl_reset_scanner(void);
long int igraph_lgl_actvertex;
%}

%output="y.tab.c"
%name-prefix="igraph_lgl_yy"
%defines

%union {
  long int edgenum;
  double weightnum;
}

%type <edgenum>   edgeid
%type <weightnum> weight

%token ALNUM
%token NEWLINE
%token HASH

%%

input :    /* empty */
         | input NEWLINE
         | input vertex
;

vertex : vertexdef edges ;

vertexdef : HASH edgeid NEWLINE       { igraph_lgl_actvertex=$2; } ;

edges :   /* empty */ | edge edges ;

edge :   edgeid NEWLINE             { 
             igraph_vector_push_back(igraph_lgl_vector, igraph_lgl_actvertex);
	     igraph_vector_push_back(igraph_lgl_vector, $1);
	     igraph_vector_push_back(igraph_lgl_weights, 0);
           }
       | edgeid weight NEWLINE      { 
             igraph_vector_push_back(igraph_lgl_vector, igraph_lgl_actvertex);
	     igraph_vector_push_back(igraph_lgl_vector, $1);
	     igraph_vector_push_back(igraph_lgl_weights, $2);
           } 
;


edgeid : ALNUM  { igraph_trie_get2(igraph_lgl_trie, 
				   igraph_lgl_yytext, 
				   igraph_lgl_yyleng, &$$); };

weight : ALNUM  { $$=igraph_lgl_get_number(igraph_lgl_yytext, 
					   igraph_lgl_yyleng); } ;

%%

int igraph_lgl_yyerror (char *s)
{
  static char str[300];  
  igraph_i_lgl_reset_scanner();

  snprintf(str, sizeof(str), "Parse error in LGL file, line %li (%s)", 
	   (long)igraph_lgl_mylineno, s);
  igraph_i_lgl_errmsg=str;
  return 0;
}

igraph_real_t igraph_lgl_get_number(const char *str, long int length) {
  igraph_real_t num;
  char *tmp=igraph_Calloc(length+1, char);
  
  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  igraph_Free(tmp);
  return num;
} 
