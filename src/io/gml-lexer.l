/*
   IGraph library.
   Copyright (C) 2007-2021  The igraph development team

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
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
%option yylineno

digit       [0-9]
whitespace  [ \t\v\f]
newline     [\r\n]

/* Use to parse inf/nan as number only when expecting a value, i.e. after a keyword.
 * Otherwise they are parsed as a keyword. */
%s VALUE

%%

^#[^\0\n\r]*            { /* comments ignored */ }

\"[^\0\"]*\"            { BEGIN(INITIAL); return STRING; }
<VALUE>(\+|\-)?((?i:nan)|(?i:inf)) { BEGIN(INITIAL); return NUM; }
(\+|\-)?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? { BEGIN(INITIAL); return NUM; }
[a-zA-Z_][a-zA-Z_0-9]*  { BEGIN(VALUE); return KEYWORD; }
\[                      {
                          BEGIN(INITIAL);
                          yyextra->depth++;
                          if (yyextra->depth >= 32) {
                            return ERROR;
                          } else {
                            return LISTOPEN;
                          }
                        }
\]                      {
                          yyextra->depth--;
                          return LISTCLOSE;
                        }
{whitespace}+           { /* other whitespace ignored */ }
{newline}+              { yy_set_bol(true); /* set "beginning of line" even after \r */ }

.                       { return ERROR; }
%%
