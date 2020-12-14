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

#include "io/pajek-header.h"
#include "io/parsers/pajek-parser.h"

#define YY_EXTRA_TYPE igraph_i_pajek_parsedata_t*
#define YY_USER_ACTION yylloc->first_line = yylineno;
#define YY_FATAL_ERROR(msg) IGRAPH_FATAL("Error in Pajek parser: " # msg)
#ifdef USING_R
#define fprintf(file, msg, ...) (1)
#ifdef stdout 
#  undef stdout
#endif
#define stdout 0
#endif
%}

%option noyywrap
%option prefix="igraph_pajek_yy"
%option nounput
%option noinput
%option nodefault
%option reentrant
%option bison-bridge
%option bison-locations

digit [0-9]
word [^ \t\r\n\0]

%%

[ \t]+          { }
%[^\n]*\n[\r]*  { }
%[^\n]*\r[\n]*  { }
\*[Nn][eE][Tt]                         { return NETWORKLINE; }
\*[Nn][Ee][Tt][Ww][Oo][Rr][Kk]         { return NETWORKLINE; }
\*[Vv][Ee][Rr][Tt][Ii][Cc][Ee][Ss]     { return VERTICESLINE; }
\*[Aa][Rr][Cc][Ss]                     { return ARCSLINE; }
\*[Ee][Dd][Gg][Ee][Ss]                 { return EDGESLINE; }
\*[Aa][Rr][Cc][Ss][Ll][Ii][Ss][Tt]     { return ARCSLISTLINE; }
\*[Ee][Dd][Gg][Ee][Ss][Ll][Ii][Ss][Tt] { return EDGESLISTLINE; }
\*[Mm][Aa][Tt][Rr][Ii][Xx]             { return MATRIXLINE; }
\n\r|\r\n|\n|\r   { yyextra->mode=0; return NEWLINE; }
\"[^\"]*\"        { return QSTR; }
\([^\)]*\)        { return PSTR; }
\-?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? { 
                    return NUM; }

[Xx]_[Ff][Aa][Cc][Tt]/[ \t\n\r]  { if (yyextra->mode==1) { return VP_X_FACT; } else { return ALNUM; } }
[Yy]_[Ff][Aa][Cc][Tt]/[ \t\n\r]  { if (yyextra->mode==1) { return VP_Y_FACT; } else { return ALNUM; } }
[Ii][Cc]/[ \t\n\r]               { if (yyextra->mode==1) { return VP_IC; } else { return ALNUM; } }
[Bb][Cc]/[ \t\n\r]               { if (yyextra->mode==1) { return VP_BC; } else { return ALNUM; } }
[Bb][Ww]/[ \t\n\r]               { if (yyextra->mode==1) { return VP_BW; } else { return ALNUM; } }
[Pp][Hh][Ii]/[ \t\n\r]           { if (yyextra->mode==1) { return VP_PHI; } else { return ALNUM; } }
[Rr]/[ \t\n\r]                   { if (yyextra->mode==1) { return VP_R; } else { return ALNUM; } }
[Qq]/[ \t\n\r]                   { if (yyextra->mode==1) { return VP_Q; } else { return ALNUM; } }
[Ff][Oo][Nn][Tt]/[ \t\n\r]       { if (yyextra->mode==1) { return VP_FONT; } else { return ALNUM; } }
[Uu][Rr][Ll]/[ \t\n\r]           { if (yyextra->mode==1) { return VP_URL; } else { return ALNUM; } }

[Cc]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_C; } else { return ALNUM; } }
[Pp]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_P; } else { return ALNUM; } }
[Ss]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_S; } else { return ALNUM; } }
[Aa]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_A; } else { return ALNUM; } }
[Ww]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_W; } else { return ALNUM; } }
[Hh]1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_H1; } else { return ALNUM; } }
[Hh]2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_H2; } else { return ALNUM; } }
[Aa]1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_A1; } else { return ALNUM; } }
[Aa]2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_A2; } else { return ALNUM; } }
[Kk]1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_K1; } else { return ALNUM; } }
[Kk]2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_K2; } else { return ALNUM; } }
[Aa][Pp]/[ \t\n\r]   { if (yyextra->mode==2) { return EP_AP; } else { return ALNUM; } }
[Ll]/[ \t\n\r]       { if (yyextra->mode==2) { return EP_L; } else { return ALNUM; } }
[Ll][Pp]/[ \t\n\r]   { if (yyextra->mode==2) { return EP_LP; } else { return ALNUM; } }

[Ll][Pp][Hh][Ii]/[ \t\n\r] { if (yyextra->mode==1) { return VP_LPHI; } else
                             if (yyextra->mode==2) { return EP_LPHI; } else { return ALNUM; } }
[Ll][Cc]/[ \t\n\r]         { if (yyextra->mode==1) { return VP_LC; } else
                             if (yyextra->mode==2) { return EP_LC; } else { return ALNUM; } }
[Ll][Rr]/[ \t\n\r]         { if (yyextra->mode==1) { return VP_LR; } else
                             if (yyextra->mode==2) { return EP_LR; } else { return ALNUM; } }
[Ll][Aa]/[ \t\n\r]         { if (yyextra->mode==1) { return VP_LA; } else
                             if (yyextra->mode==2) { return EP_LA; } else { return ALNUM; } }
[Ss][Ii][Zz][Ee]/[ \t\n\r] { if (yyextra->mode==1) { return VP_SIZE; } else 
                             if (yyextra->mode==2) { return EP_SIZE; } else { return ALNUM; } }
[Ff][Oo][Ss]/[ \t\n\r]     { if (yyextra->mode==1) { return VP_FOS; } else 
                             if (yyextra->mode==2) { return EP_FOS; } else { return ALNUM; } }
                             
{word}+           { return ALNUM; }

<<EOF>>           { if (yyextra->eof) {
                       yyterminate();
                    } else {
                       yyextra->eof=1;
                       return NEWLINE; 
                    }
                  }

.                 { return ERROR; }

%%
