%{
#include "foreign-ncol-parser.h"
long int igraph_ncol_mylineno;
%}

%option noyywrap
%option prefix="igraph_ncol_yy"
%option outfile="lex.yy.c"

alnum [^ \t\n]

%%

 /* ------------------------------------------------whitespace------*/
[ \t]*               { }

 /* ---------------------------------------------------newline------*/
\n                   { igraph_ncol_mylineno++; return NEWLINE; }

 /* ----------------------------------------------alphanumeric------*/
{alnum}+             { return ALNUM; }

%%
