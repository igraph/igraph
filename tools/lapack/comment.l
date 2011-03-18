%{
#include <stdio.h>

/* extern FILE *commentFile, *localVarFile, *codeFile; */

%}

whitespace [\n\t ]*
any	   .*

%%
"*/"{whitespace}"/*"         {yytext[0]=yytext[1]=yytext[yyleng-1]=yytext[yyleng-2]=' ';printf("%s",yytext);}
"\n"                         {printf("%s", yytext);}
.                            {printf("%s", yytext);}
