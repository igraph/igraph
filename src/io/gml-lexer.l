/* 
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "config.h"
#include <stdlib.h>

#include "io/gml-header.h"
#include "io/parsers/gml-parser.h"

#define YY_EXTRA_TYPE igraph_i_gml_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
#define YY_FATAL_ERROR(msg) IGRAPH_FATAL("Error in GML parser: " # msg)
#ifdef USING_R
#define fprintf(file, msg, ...) (1)
#ifdef stdout 
#  undef stdout
#endif
#define stdout 0
#endif
%}

%option noyywrap
%option prefix="igraph_gml_yy"
%option nounput
%option noinput
%option nodefault
%option reentrant
%option bison-bridge
%option bison-locations

digit       [0-9]
whitespace  [ \r\n\t]

%%

^#[^\n\r]*[\n]|[\r]     { /* comments ignored */ }

\"[^\"]*\"                                      { return STRING; }
\-?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? { return NUM; }
[a-zA-Z_][a-zA-Z_0-9]*  { return KEYWORD; }
\[                      { return LISTOPEN; }
\]                      { return LISTCLOSE; }
\n\r|\r\n|\r|\n         { }
{whitespace}            { /* other whitespace ignored */ }

<<EOF>>                 { 
                          if (yyextra->eof) {
			    yyterminate();
			  } else {
			    yyextra->eof=1;
			    return EOFF;
			  }			  
                        }
.                       { return ERROR; }
%%
