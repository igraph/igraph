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
%option yylineno
%option caseless

digit [0-9]
word [^ \t\v\f\r\n\0]
whitespace [ \t\v\f]

%x netline

%%

<*>{whitespace}+  { }
%[^\n\r\0]*[\n\r]*  { } /* comments */
\*net             { BEGIN(netline); return NETWORKLINE; }
\*network         { BEGIN(netline); return NETWORKLINE; }
<netline>{whitespace}({word}|{whitespace})* { 
   return NET_TITLE; }
\*vertices        { return VERTICESLINE; }
\*arcs            { return ARCSLINE; }
\*edges           { return EDGESLINE; }
\*arcslist        { return ARCSLISTLINE; }
\*edgeslist       { return EDGESLISTLINE; }
\*matrix          { return MATRIXLINE; }
<*>[\n\r]+        { BEGIN(INITIAL); yyextra->mode=0; return NEWLINE; } /* skip over multiple newlines */
\"[^\"\0]*\"      { return QSTR; }
\([^\)\0]*\)      { return PSTR; }
(\+|\-)?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? {
                    return NUM; }

x_fact/[ \t\n\r]  { if (yyextra->mode==1) { return VP_X_FACT; } else { return ALNUM; } }
y_fact/[ \t\n\r]  { if (yyextra->mode==1) { return VP_Y_FACT; } else { return ALNUM; } }
ic/[ \t\n\r]      { if (yyextra->mode==1) { return VP_IC; } else { return ALNUM; } }
bc/[ \t\n\r]      { if (yyextra->mode==1) { return VP_BC; } else { return ALNUM; } }
bw/[ \t\n\r]      { if (yyextra->mode==1) { return VP_BW; } else { return ALNUM; } }
phi/[ \t\n\r]     { if (yyextra->mode==1) { return VP_PHI; } else { return ALNUM; } }
r/[ \t\n\r]       { if (yyextra->mode==1) { return VP_R; } else { return ALNUM; } }
q/[ \t\n\r]       { if (yyextra->mode==1) { return VP_Q; } else { return ALNUM; } }
font/[ \t\n\r]    { if (yyextra->mode==1) { return VP_FONT; } else { return ALNUM; } }
url/[ \t\n\r]     { if (yyextra->mode==1) { return VP_URL; } else { return ALNUM; } }

c/[ \t\n\r]       { if (yyextra->mode==2) { return EP_C; } else { return ALNUM; } }
p/[ \t\n\r]       { if (yyextra->mode==2) { return EP_P; } else { return ALNUM; } }
s/[ \t\n\r]       { if (yyextra->mode==2) { return EP_S; } else { return ALNUM; } }
a/[ \t\n\r]       { if (yyextra->mode==2) { return EP_A; } else { return ALNUM; } }
w/[ \t\n\r]       { if (yyextra->mode==2) { return EP_W; } else { return ALNUM; } }
h1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_H1; } else { return ALNUM; } }
h2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_H2; } else { return ALNUM; } }
a1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_A1; } else { return ALNUM; } }
a2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_A2; } else { return ALNUM; } }
k1/[ \t\n\r]      { if (yyextra->mode==2) { return EP_K1; } else { return ALNUM; } }
k2/[ \t\n\r]      { if (yyextra->mode==2) { return EP_K2; } else { return ALNUM; } }
ap/[ \t\n\r]      { if (yyextra->mode==2) { return EP_AP; } else { return ALNUM; } }
l/[ \t\n\r]       { if (yyextra->mode==2) { return EP_L; } else { return ALNUM; } }
lp/[ \t\n\r]      { if (yyextra->mode==2) { return EP_LP; } else { return ALNUM; } }

lphi/[ \t\n\r]    { if (yyextra->mode==1) { return VP_LPHI; } else
                             if (yyextra->mode==2) { return EP_LPHI; } else { return ALNUM; } }
lc/[ \t\n\r]      { if (yyextra->mode==1) { return VP_LC; } else
                             if (yyextra->mode==2) { return EP_LC; } else { return ALNUM; } }
lr/[ \t\n\r]      { if (yyextra->mode==1) { return VP_LR; } else
                             if (yyextra->mode==2) { return EP_LR; } else { return ALNUM; } }
la/[ \t\n\r]      { if (yyextra->mode==1) { return VP_LA; } else
                             if (yyextra->mode==2) { return EP_LA; } else { return ALNUM; } }
size/[ \t\n\r]    { if (yyextra->mode==1) { return VP_SIZE; } else
                             if (yyextra->mode==2) { return EP_SIZE; } else { return ALNUM; } }
fos/[ \t\n\r]     { if (yyextra->mode==1) { return VP_FOS; } else
                             if (yyextra->mode==2) { return EP_FOS; } else { return ALNUM; } }

{word}+           { return ALNUM; }

<<EOF>>           { if (yyextra->eof) {
                       yyterminate();
                    } else {
                       yyextra->eof=1;
                       return NEWLINE;
                    }
                  }

<*>.              { return ERROR; }

%%
