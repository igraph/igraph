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

#include "config.h"
#include <stdlib.h>
#include "foreign-lgl-parser.h"
#include "foreign-lgl-header.h"
#define YY_EXTRA_TYPE igraph_i_lgl_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
%}

%option noyywrap
%option prefix="igraph_lgl_yy"
%option outfile="lex.yy.c"
%option nounput
%option noinput
%option reentrant
%option bison-bridge
%option bison-locations

alnum [^ \t\r\n#]

%%

 /* --------------------------------------------------hashmark------*/
#                    { return HASH; }

 /* ------------------------------------------------whitespace------*/
[ \t]*               { }

 /* ---------------------------------------------------newline------*/
\n\r|\r\n|\n|\r      { return NEWLINE; }

 /* ----------------------------------------------alphanumeric------*/
{alnum}+             { return ALNUM; }

<<EOF>>           { if (yyextra->eof) {
                       yyterminate();
                    } else {
                       yyextra->eof=1;
                       return NEWLINE; 
                    }
                  }

%%
