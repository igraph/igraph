%{
#include <stdio.h>
#include <string.h>
extern int igraph_ncol_yylex();
extern int igraph_ncol_mylineno;
extern char *igraph_ncol_yytext;
extern int igraph_ncol_yyleng;
int igraph_ncol_yyerror(char *s);
#include "types.h" 
#include "memory.h"
#include "error.h"
extern igraph_vector_t *igraph_ncol_vector;
extern igraph_vector_t *igraph_ncol_weights;
extern igraph_trie_t *igraph_ncol_trie;
real_t igraph_ncol_get_number(const char *str, long int len);
%}

%output="y.tab.c"
%name-prefix="igraph_ncol_yy"
%defines

%union {
  long int edgenum;
  double weightnum;
}

%type <edgenum>   edgeid
%type <weightnum> weight

%token ALNUM
%token NEWLINE

%%

input :    /* empty */
         | input NEWLINE
         | input edge
;

edge :   edgeid edgeid NEWLINE        { 
           igraph_vector_push_back(igraph_ncol_vector, $1);
           igraph_vector_push_back(igraph_ncol_vector, $2);
	   igraph_vector_push_back(igraph_ncol_weights, 0);
       }
       | edgeid edgeid weight NEWLINE { 
           igraph_vector_push_back(igraph_ncol_vector, $1);
           igraph_vector_push_back(igraph_ncol_vector, $2);
	   igraph_vector_push_back(igraph_ncol_weights, $3);
       }
;

edgeid : ALNUM  { igraph_trie_get2(igraph_ncol_trie, 
				   igraph_ncol_yytext, 
				   igraph_ncol_yyleng, &$$); };

weight : ALNUM  { $$=igraph_ncol_get_number(igraph_ncol_yytext, 
					   igraph_ncol_yyleng); } ;

%%

int igraph_ncol_yyerror (char *s)
{
  IGRAPH_ERROR("Parse error", IGRAPH_PARSEERROR);
}

real_t igraph_ncol_get_number(const char *str, long int length) {
  real_t num;
  char *tmp=Calloc(length+1, char);
  
  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  Free(tmp);
  return num;
} 
