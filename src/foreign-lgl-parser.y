%{
#include <stdio.h>
#include <string.h>
extern int igraph_lgl_yylex();
extern int mylineno;
extern char *igraph_lgl_yytext;
extern int igraph_lgl_yyleng;
int igraph_lgl_yyerror(char *s);
#include "types.h" 
#include "memory.h"
#include "error.h"
extern vector_t *igraph_lgl_vector;
extern vector_t *igraph_lgl_weights;
extern igraph_trie_t *igraph_lgl_trie;
real_t igraph_lgl_get_number(const char *str, long int len);
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
             vector_push_back(igraph_lgl_vector, igraph_lgl_actvertex);
	     vector_push_back(igraph_lgl_vector, $1);
	     vector_push_back(igraph_lgl_weights, 0);
           }
       | edgeid weight NEWLINE      { 
             vector_push_back(igraph_lgl_vector, igraph_lgl_actvertex);
	     vector_push_back(igraph_lgl_vector, $1);
	     vector_push_back(igraph_lgl_weights, $2);
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
  IGRAPH_ERROR("Parse error", IGRAPH_PARSEERROR);
}

real_t igraph_lgl_get_number(const char *str, long int length) {
  real_t num;
  char *tmp=Calloc(length+1, char);
  
  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  Free(tmp);
  return num;
} 
