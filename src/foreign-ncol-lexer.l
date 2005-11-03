%{
#include "foreign-ncol-parser.h"
long int mylineno;
%}

%option noyywrap
%option prefix="igraph_ncol_yy"
%option outfile="lex.yy.c"

alnum [^ \t\n]

%%

 /* ------------------------------------------------whitespace------*/
[ \t]*               { }

 /* ---------------------------------------------------newline------*/
\n                   { mylineno++; return NEWLINE; }

 /* ----------------------------------------------alphanumeric------*/
{alnum}+             { return ALNUM; }

%%
