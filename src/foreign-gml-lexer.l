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
#include "foreign-gml-parser.h"
extern long int igraph_gml_mylineno;
extern int igraph_gml_mode;
int igraph_gml_eof;
void igraph_i_gml_reset_scanner() {
  YY_FLUSH_BUFFER;
}
%}

%option noyywrap
%option prefix="igraph_gml_yy"
%option outfile="lex.yy.c"
%option nounput
%option noinput

digit       [0-9]
whitespace  [ \r\n\t]

%%

^#[^\n\r]*[\n]|[\r]     { igraph_gml_mylineno++; /* comments ignored */ }

\"[^\"]*\"                                      { return STRING; }
\-?{digit}+(\.{digit}+)?([eE](\+|\-)?{digit}+)? { return NUM; }
[a-zA-Z_][a-zA-Z_0-9]*  { return KEYWORD; }
\[                      { return LISTOPEN; }
\]                      { return LISTCLOSE; }
\n\r|\r\n|\r|\n         { igraph_gml_mylineno++; }
{whitespace}            { /* other whitespace ignored */ }

<<EOF>>                 { 
                          if (igraph_gml_eof) {
			    yyterminate();
			  } else {
			    igraph_gml_eof=1;
			    return EOFF;
			  }			  
                        }

%%
