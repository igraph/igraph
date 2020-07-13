/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
#include "foreign-ncol-header.h"
#include "foreign-ncol-parser.h"
#define YY_EXTRA_TYPE igraph_i_ncol_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
#ifdef USING_R
#define fprintf(file, msg, ...) (1)
#define YY_FATAL_ERROR(msg)                                       \
  igraph_error("Fatal error in NCOL parser: " # msg, __FILE__,    \
               __LINE__, IGRAPH_PARSEERROR);
#ifdef stdout 
#  undef stdout
#endif
#define stdout 0
#endif
%}

%option noyywrap
%option prefix="igraph_ncol_yy"
%option outfile="lex.yy.c"
%option nounput
%option noinput
%option nodefault
%option reentrant
%option bison-bridge
%option bison-locations

alnum [^ \t\n\r]

%%

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

 /* ---------------------------------------------anything else------*/
.                    { return ERROR; }

%%

