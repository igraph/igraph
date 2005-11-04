%{
#include "foreign-lgl-parser.h"
long int mylineno;
%}

%option noyywrap
%option prefix="igraph_lgl_yy"
%option outfile="lex.yy.c"

alnum [^ \t\n#]

%%

 /* --------------------------------------------------hashmark------*/
#                    { return HASH; }

 /* ------------------------------------------------whitespace------*/
[ \t]*               { }

 /* ---------------------------------------------------newline------*/
\n                   { mylineno++; return NEWLINE; }

 /* ----------------------------------------------alphanumeric------*/
{alnum}+             { return ALNUM; }

%%
