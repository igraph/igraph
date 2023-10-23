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
#define YY_USER_INIT BEGIN(unknown) /* we start in the 'unknown' start condition */
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

whitespace [ \t]
digit [0-9]

/* Any use of {newline} below must use yy_set_bol(true) in order to mark the character
   following a single \r as the first on a new line, and allow the ^ pattern to match. 
   This pattern must match single newlines only, in order to follow Pajek's "no newline
   after *Vertices" convention. */
newline (\n|\r|\r\n|\n\r) 

/* Anything except non-printable (00-1F), space (20), del (7F), '"' and '*'. */
word [^\x00-\x20\x7f"*]

any_no_star {word}|\"|{whitespace}
any         {any_no_star}|\*

/* 'unknown' skips text at the beginning of the file, lines below an unknown *Word
 * 'unknown_line' skips the rest of the line after an unknown *Word. */
%x unknown unknown_line
%s vert edge

/* Notes:
  *  - Unquoted '*' characters may only appear at the start of a line-initial word.
  *  - Both LF and CR LF line endings are allowed.
  *  - Pajek files do not allow empty lines after *Vertices (empty lines should signify the end of the file),
  *    therefore we are careful not to skip newlines in the lexer.
  */

%%

 /* Skip all text until the next *Word at the beginning of a line. */
<unknown>^\*                     { BEGIN(INITIAL); yyless(0); yy_set_bol(true); }
<unknown>{any_no_star}{any}*     { } /* match cannot start with a * in order not to take precedence over ^\* above */
<unknown,unknown_line>{newline}  { yy_set_bol(true); }
<unknown_line>{any}*             { BEGIN(unknown); }

^%({any})*{newline}        { yy_set_bol(true); } /* comments */

<*>{whitespace}+   { }

^\*network         { BEGIN(unknown_line); return NETWORKLINE; }

^\*vertices        { BEGIN(vert); return VERTICESLINE; }
^\*arcs            { BEGIN(edge); return ARCSLINE; }
^\*edges           { BEGIN(edge); return EDGESLINE; }
^\*arcslist        { BEGIN(INITIAL); return ARCSLISTLINE; }
^\*edgeslist       { BEGIN(INITIAL);return EDGESLISTLINE; }
^\*matrix          { BEGIN(INITIAL); return MATRIXLINE; }

^\*[a-z]+          { BEGIN(unknown_line); IGRAPH_WARNINGF("Skipping unknown section '%s' on line %d.", yytext, yylineno); }

{newline}          { yy_set_bol(true); return NEWLINE; }

 /* Newlines not allowed in strings. */
\"[^\"\0\n\r]*\"   { return QSTR; }

(\+|\-)?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? { return NUM; }

<vert>{
x_fact  { return VP_X_FACT; }
y_fact  { return VP_Y_FACT; }
ic      { return VP_IC; }
bc      { return VP_BC; }
bw      { return VP_BW; }
phi     { return VP_PHI; }
r       { return VP_R; }
q       { return VP_Q; }
font    { return VP_FONT; }
url     { return VP_URL; }

lphi    { return VP_LPHI; }
lc      { return VP_LC; }
lr      { return VP_LR; }
la      { return VP_LA; }
size    { return VP_SIZE; }
fos     { return VP_FOS; }
}

<edge>{
c       { return EP_C; }
p       { return EP_P; }
s       { return EP_S; }
a       { return EP_A; }
w       { return EP_W; }
h1      { return EP_H1; }
h2      { return EP_H2; }
a1      { return EP_A1; }
a2      { return EP_A2; }
k1      { return EP_K1; }
k2      { return EP_K2; }
ap      { return EP_AP; }
l       { return EP_L; }
lp      { return EP_LP; }

lphi    { return EP_LPHI; }
lc      { return EP_LC; }
lr      { return EP_LR; }
la      { return EP_LA; }
size    { return EP_SIZE; }
fos     { return EP_FOS; }
}

{word}+           { return ALNUM; }

<<EOF>>           { if (yyextra->eof) {
                       yyterminate();
                    } else {
                       yyextra->eof=true;
                       return NEWLINE;
                    }
                  }

<*>.              { return ERROR; }

%%
