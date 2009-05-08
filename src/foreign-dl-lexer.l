/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <stdarg.h>
#include "foreign-dl-header.h"
#include "foreign-dl-parser.h"
extern long int igraph_dl_mylineno;
extern int igraph_dl_mode;
extern int igraph_dl_eof;
extern char *igraph_dl_yytext;
int igraph_i_dl_mode=0;
void igraph_i_dl_reset_scanner() {
  YY_FLUSH_BUFFER;
}
int myprintf(const char *arg, ...) {
/*   va_list argptr; */
/*   va_start(argptr, arg); */
/*   printf(arg, argptr); */
  return 0;
}
%}

%option noyywrap
%option prefix="igraph_dl_yy"
%option outfile="lex.yy.c"
%option nounput
%option noinput

digit      [0-9]
whitespace [ \t\v\f]

%x LABELM FULLMATRIX EDGELIST NODELIST

%% 

<*>\n\r|\r\n|\r|\n         { igraph_dl_mylineno++; myprintf("newline\n"); return NEWLINE; }

[dD][lL]{whitespace}+      { myprintf("dl\n"); return DL; }
[nN]{whitespace}*[=]{whitespace}* {
  myprintf("n=\n"); return NEQ; }
{digit}+                   { myprintf("num\n"); return NUM; }

<INITIAL,LABELM>[dD][aA][tT][aA][:]        { 
  switch (igraph_i_dl_mode) { 
  case 0: BEGIN(FULLMATRIX); 
    myprintf("full matrix mode\n");
    break;
  case 1: BEGIN(EDGELIST);
    break;
  case 2: BEGIN(NODELIST);
    break;
  } 
  myprintf("data\n"); return DATA; }

[lL][aA][bB][eE][lL][sS]:  { BEGIN(LABELM); myprintf("labels\n"); return LABELS; }
<INITIAL,LABELM>[lL][aA][bB][eE][lL][sS]{whitespace}+[eE][mM][bB][eE][dD][dD][eE][dD]:?{whitespace}* {
  myprintf("labels embedded\n"); return LABELSEMBEDDED; }
[fF][oO][rR][mM][aA][tT]{whitespace}*[=]{whitespace}*[fF][uU][lL][lL][mM][aA][tT][rR][iI][xX]{whitespace}* {
  igraph_i_dl_mode=0; myprintf("format fm\n"); return FORMATFULLMATRIX; }
[fF][oO][rR][mM][aA][tT]{whitespace}*[=]{whitespace}*[eE][dD][gG][eE][lL][iI][sS][tT][1]{whitespace}* {
  igraph_i_dl_mode=1; myprintf("format el\n"); return FORMATEDGELIST1; }
[fF][oO][rR][mM][aA][tT]{whitespace}*[=]{whitespace}*[nN][oO][dD][eE][lL][iI][sS][tT][1]{whitespace}* {
  igraph_i_dl_mode=2; myprintf("format nl\n"); return FORMATNODELIST1; }

<LABELM>[, ]                               { myprintf(","); /* eaten up */ }
<LABELM>[^, \t\n\r\f\v]+{whitespace}*      { myprintf("label1\n"); return LABEL; }

<FULLMATRIX>{digit}{whitespace}*          { myprintf("digit"); return DIGIT; }
<FULLMATRIX>[^ \t\n\r\v\f,]+              { myprintf("label"); return LABEL; }
<FULLMATRIX>{whitespace}                  { }

<EDGELIST>\-?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)?  { myprintf("num "); return NUM; }
<EDGELIST>[^ \t\n\r\v\f,]+                                 { myprintf("label "); return LABEL; }
<EDGELIST>{whitespace}*                                    { }

<NODELIST>{digit}+                      { myprintf("num3 "); return NUM; }
<NODELIST>[^ \t\r\n\v\f,]+              { myprintf("label3 "); return LABEL; }
<NODELIST>{whitespace}*                 { }

{whitespace}+                      { /* eaten up */ }

<<EOF>>                 { 
                          if (igraph_dl_eof) {
			    yyterminate();
			  } else {
			    myprintf("EOF");
			    igraph_dl_eof=1;
			    BEGIN(INITIAL);
			    return EOFF;
			  }			  
                        }

<*>. { myprintf("gebasz---%s---", igraph_dl_yytext); return 0; }
