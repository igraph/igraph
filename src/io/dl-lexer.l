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

#include <stdlib.h>
#include <stdarg.h>

#include "io/dl-header.h"
#include "io/parsers/dl-parser.h"

#define YY_EXTRA_TYPE igraph_i_dl_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
#define YY_FATAL_ERROR(msg) IGRAPH_FATAL("Error in DL parser: " # msg)
#ifdef USING_R
#define fprintf(file, msg, ...) (1)
#ifdef stdout
#  undef stdout
#endif
#define stdout 0
#endif
%}

%option noyywrap
%option prefix="igraph_dl_yy"
%option nounput
%option noinput
%option nodefault
%option reentrant
%option bison-bridge
%option bison-locations
%option yylineno
%option caseless

digit      [0-9]
whitespace [ \t\v\f]

%x LABELM FULLMATRIX EDGELIST NODELIST

%%

<*>\n\r|\r\n|\r|\n             { return NEWLINE; }

dl{whitespace}+                { return DL; }
n{whitespace}*[=]{whitespace}* { return NEQ; }
{digit}+                       { return NUM; }

<INITIAL,LABELM>data:        {
  switch (yyextra->mode) {
  case 0: BEGIN(FULLMATRIX);
    break;
  case 1: BEGIN(EDGELIST);
    break;
  case 2: BEGIN(NODELIST);
    break;
  }
  return DATA; }

labels:  { BEGIN(LABELM); return LABELS; }
<INITIAL,LABELM>labels{whitespace}+embedded:?{whitespace}* {
  return LABELSEMBEDDED; }
format{whitespace}*[=]{whitespace}*fullmatrix{whitespace}* {
  yyextra->mode=0; return FORMATFULLMATRIX; }
format{whitespace}*[=]{whitespace}*edgelist1{whitespace}* {
  yyextra->mode=1; return FORMATEDGELIST1; }
format{whitespace}*[=]{whitespace}*nodelist1{whitespace}* {
  yyextra->mode=2; return FORMATNODELIST1; }

<LABELM>[, ]                               { /* eaten up */ }
<LABELM>[^, \t\n\r\f\v\0]+{whitespace}*    { return LABEL; }

<FULLMATRIX>{digit}{whitespace}*          { return DIGIT; }
<FULLMATRIX>[^ \t\n\r\v\f\0,]+            { return LABEL; }
<FULLMATRIX>{whitespace}                  { }

<EDGELIST>(\+|\-)?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)?  { return NUM; }
<EDGELIST>[^ \t\n\r\v\f\0,]+                                 { return LABEL; }
<EDGELIST>{whitespace}*                                    { }

<NODELIST>{digit}+                      { return NUM; }
<NODELIST>[^ \t\r\n\v\f\0,]+            { return LABEL; }
<NODELIST>{whitespace}*                 { }

{whitespace}+                      { /* eaten up */ }

<<EOF>>                 {
                          if (yyextra->eof) {
                            yyterminate();
                          } else {
                            yyextra->eof=1;
                            BEGIN(INITIAL);
                            return EOFF;
                          }
                        }

<*>. { return 0; }
