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

#include "io/lgl-header.h"
#include "io/parsers/lgl-parser.h"

#define YY_EXTRA_TYPE igraph_i_lgl_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
#define YY_FATAL_ERROR(msg) IGRAPH_FATAL("Error in LGL parser: " # msg)
#ifdef USING_R
#define fprintf(file, msg, ...) (1)
#ifdef stdout
#  undef stdout
#endif
#define stdout 0
#endif
%}

%option noyywrap
%option prefix="igraph_lgl_yy"
%option nounput
%option noinput
%option nodefault
%option reentrant
%option bison-bridge
%option bison-locations
%option yylineno

/* Anything except non-printable (00-1F), space (20), del (7F) and # */
alnum [^\x00-\x20\x7f#]

%s LINE

%%

 /* ------------------------------------------------whitespace------*/
[ \t]+               { /* skip space */ }

 /* --------------------------------------------------hashmark------*/
#                    { BEGIN(LINE); return HASH; }

 /* ----------------------------------------------alphanumeric------*/
{alnum}+             { BEGIN(LINE); return ALNUM; }

 /* ---------------------------------------------------newline------*/
\n\r|\r\n|\n|\r      |
<LINE><<EOF>>        { BEGIN(INITIAL); return NEWLINE; }

 /* ---------------------------------------------anything else------*/
.                 { return ERROR; }

%%
